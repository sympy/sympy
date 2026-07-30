"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.1 Sine/4.1.1.2 (g cos)^p (a+b sin)^m.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, m, p = symbols('a b c d e m p')

def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_1():
    f = (a*sin(c + d*x) + a)*cos(c + d*x)**7
    F = 8*(a*sin(c + d*x) + a)**5/(5*a**4*d) - 2*(a*sin(c + d*x) + a)**6/(a**5*d) + 6*(a*sin(c + d*x) + a)**7/(7*a**6*d) - (a*sin(c + d*x) + a)**8/(8*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_2():
    f = (a*sin(c + d*x) + a)*cos(c + d*x)**6
    F = 5*a*x/16 + a*sin(c + d*x)*cos(c + d*x)**5/(6*d) + 5*a*sin(c + d*x)*cos(c + d*x)**3/(24*d) + 5*a*sin(c + d*x)*cos(c + d*x)/(16*d) - a*cos(c + d*x)**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_3():
    f = (a*sin(c + d*x) + a)*cos(c + d*x)**5
    F = (a*sin(c + d*x) + a)**4/(a**3*d) - 4*(a*sin(c + d*x) + a)**5/(5*a**4*d) + (a*sin(c + d*x) + a)**6/(6*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_4():
    f = (a*sin(c + d*x) + a)*cos(c + d*x)**4
    F = 3*a*x/8 + a*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 3*a*sin(c + d*x)*cos(c + d*x)/(8*d) - a*cos(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_5():
    f = (a*sin(c + d*x) + a)*cos(c + d*x)**3
    F = 2*(a*sin(c + d*x) + a)**3/(3*a**2*d) - (a*sin(c + d*x) + a)**4/(4*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_6():
    f = (a*sin(c + d*x) + a)*cos(c + d*x)**2
    F = a*x/2 + a*sin(c + d*x)*cos(c + d*x)/(2*d) - a*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_7():
    f = (a*sin(c + d*x) + a)*sec(c + d*x)
    F = -a*log(1 - sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_8():
    f = (a*sin(c + d*x) + a)*sec(c + d*x)**2
    F = a*tan(c + d*x)/d + a*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_9():
    f = (a*sin(c + d*x) + a)*sec(c + d*x)**3
    F = a**2/(2*d*(-a*sin(c + d*x) + a)) + a*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_10():
    f = (a*sin(c + d*x) + a)*sec(c + d*x)**4
    F = a*tan(c + d*x)**3/(3*d) + a*tan(c + d*x)/d + a*sec(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_11():
    f = (a*sin(c + d*x) + a)*sec(c + d*x)**5
    F = a**3/(8*d*(-a*sin(c + d*x) + a)**2) - a**2/(8*d*(a*sin(c + d*x) + a)) + a**2/(4*d*(-a*sin(c + d*x) + a)) + 3*a*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_12():
    f = (a*sin(c + d*x) + a)**2*cos(c + d*x)**6
    F = 45*a**2*x/128 + 3*a**2*sin(c + d*x)*cos(c + d*x)**5/(16*d) + 15*a**2*sin(c + d*x)*cos(c + d*x)**3/(64*d) + 45*a**2*sin(c + d*x)*cos(c + d*x)/(128*d) - 9*a**2*cos(c + d*x)**7/(56*d) - (a**2*sin(c + d*x) + a**2)*cos(c + d*x)**7/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_13():
    f = (a*sin(c + d*x) + a)**2*cos(c + d*x)**5
    F = 4*(a*sin(c + d*x) + a)**5/(5*a**3*d) - 2*(a*sin(c + d*x) + a)**6/(3*a**4*d) + (a*sin(c + d*x) + a)**7/(7*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_14():
    f = (a*sin(c + d*x) + a)**2*cos(c + d*x)**4
    F = 7*a**2*x/16 + 7*a**2*sin(c + d*x)*cos(c + d*x)**3/(24*d) + 7*a**2*sin(c + d*x)*cos(c + d*x)/(16*d) - 7*a**2*cos(c + d*x)**5/(30*d) - (a**2*sin(c + d*x) + a**2)*cos(c + d*x)**5/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_15():
    f = (a*sin(c + d*x) + a)**2*cos(c + d*x)**3
    F = (a*sin(c + d*x) + a)**4/(2*a**2*d) - (a*sin(c + d*x) + a)**5/(5*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_16():
    f = (a*sin(c + d*x) + a)**2*cos(c + d*x)**2
    F = 5*a**2*x/8 + 5*a**2*sin(c + d*x)*cos(c + d*x)/(8*d) - 5*a**2*cos(c + d*x)**3/(12*d) - (a**2*sin(c + d*x) + a**2)*cos(c + d*x)**3/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_17():
    f = (a*sin(c + d*x) + a)**2*cos(c + d*x)
    F = (a*sin(c + d*x) + a)**3/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_18():
    f = (a*sin(c + d*x) + a)**2*sec(c + d*x)
    F = -2*a**2*log(1 - sin(c + d*x))/d - a**2*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_19():
    f = (a*sin(c + d*x) + a)**2*sec(c + d*x)**2
    F = 2*a**4*cos(c + d*x)/(d*(-a**2*sin(c + d*x) + a**2)) - a**2*x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_20():
    f = (a*sin(c + d*x) + a)**2*sec(c + d*x)**3
    F = a**3/(d*(-a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_21():
    f = (a*sin(c + d*x) + a)**2*sec(c + d*x)**4
    F = a**4*cos(c + d*x)/(3*d*(-a**2*sin(c + d*x) + a**2)) + a**4*cos(c + d*x)/(3*d*(-a*sin(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_22():
    f = (a*sin(c + d*x) + a)**2*sec(c + d*x)**5
    F = a**4/(4*d*(-a*sin(c + d*x) + a)**2) + a**3/(4*d*(-a*sin(c + d*x) + a)) + a**2*atanh(sin(c + d*x))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_23():
    f = (a*sin(c + d*x) + a)**2*sec(c + d*x)**6
    F = a**2*tan(c + d*x)**3/(5*d) + 3*a**2*tan(c + d*x)/(5*d) + 2*(a**2*sin(c + d*x) + a**2)*sec(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_24():
    f = (a*sin(c + d*x) + a)**2*sec(c + d*x)**7
    F = a**5/(12*d*(-a*sin(c + d*x) + a)**3) + a**4/(8*d*(-a*sin(c + d*x) + a)**2) - a**3/(16*d*(a*sin(c + d*x) + a)) + 3*a**3/(16*d*(-a*sin(c + d*x) + a)) + a**2*atanh(sin(c + d*x))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_25():
    f = (a*sin(c + d*x) + a)**2*sec(c + d*x)**8
    F = a**2*tan(c + d*x)**5/(7*d) + 10*a**2*tan(c + d*x)**3/(21*d) + 5*a**2*tan(c + d*x)/(7*d) + 2*(a**2*sin(c + d*x) + a**2)*sec(c + d*x)**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_26():
    f = (a*sin(c + d*x) + a)**3*cos(c + d*x)**6
    F = 55*a**3*x/128 + 11*a**3*sin(c + d*x)*cos(c + d*x)**5/(48*d) + 55*a**3*sin(c + d*x)*cos(c + d*x)**3/(192*d) + 55*a**3*sin(c + d*x)*cos(c + d*x)/(128*d) - 11*a**3*cos(c + d*x)**7/(56*d) - a*(a*sin(c + d*x) + a)**2*cos(c + d*x)**7/(9*d) - 11*(a**3*sin(c + d*x) + a**3)*cos(c + d*x)**7/(72*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_27():
    f = (a*sin(c + d*x) + a)**3*cos(c + d*x)**5
    F = 2*(a*sin(c + d*x) + a)**6/(3*a**3*d) - 4*(a*sin(c + d*x) + a)**7/(7*a**4*d) + (a*sin(c + d*x) + a)**8/(8*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_28():
    f = (a*sin(c + d*x) + a)**3*cos(c + d*x)**4
    F = 9*a**3*x/16 + 3*a**3*sin(c + d*x)*cos(c + d*x)**3/(8*d) + 9*a**3*sin(c + d*x)*cos(c + d*x)/(16*d) - 3*a**3*cos(c + d*x)**5/(10*d) - a*(a*sin(c + d*x) + a)**2*cos(c + d*x)**5/(7*d) - 3*(a**3*sin(c + d*x) + a**3)*cos(c + d*x)**5/(14*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_29():
    f = (a*sin(c + d*x) + a)**3*cos(c + d*x)**3
    F = 2*(a*sin(c + d*x) + a)**5/(5*a**2*d) - (a*sin(c + d*x) + a)**6/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_30():
    f = (a*sin(c + d*x) + a)**3*cos(c + d*x)**2
    F = 7*a**3*x/8 + 7*a**3*sin(c + d*x)*cos(c + d*x)/(8*d) - 7*a**3*cos(c + d*x)**3/(12*d) - a*(a*sin(c + d*x) + a)**2*cos(c + d*x)**3/(5*d) - 7*(a**3*sin(c + d*x) + a**3)*cos(c + d*x)**3/(20*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_31():
    f = (a*sin(c + d*x) + a)**3*cos(c + d*x)
    F = (a*sin(c + d*x) + a)**4/(4*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_32():
    f = (a*sin(c + d*x) + a)**3*sec(c + d*x)
    F = -4*a**3*log(1 - sin(c + d*x))/d - a**3*sin(c + d*x)**2/(2*d) - 3*a**3*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_33():
    f = (a*sin(c + d*x) + a)**3*sec(c + d*x)**2
    F = 2*a**5*cos(c + d*x)**3/(d*(-a*sin(c + d*x) + a)**2) - 3*a**3*x + 3*a**3*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_34():
    f = (a*sin(c + d*x) + a)**3*sec(c + d*x)**3
    F = 2*a**4/(d*(-a*sin(c + d*x) + a)) + a**3*log(1 - sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_35():
    f = (a*sin(c + d*x) + a)**3*sec(c + d*x)**4
    F = a**6*cos(c + d*x)**3/(3*d*(-a*sin(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_36():
    f = (a*sin(c + d*x) + a)**3*sec(c + d*x)**5
    F = a**5/(2*d*(-a*sin(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_37():
    f = (a*sin(c + d*x) + a)**3*sec(c + d*x)**6
    F = 2*a**6*cos(c + d*x)/(15*d*(-a**3*sin(c + d*x) + a**3)) + a**6*cos(c + d*x)/(5*d*(-a*sin(c + d*x) + a)**3) + 2*a**5*cos(c + d*x)/(15*d*(-a*sin(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_38():
    f = (a*sin(c + d*x) + a)**3*sec(c + d*x)**7
    F = a**6/(6*d*(-a*sin(c + d*x) + a)**3) + a**5/(8*d*(-a*sin(c + d*x) + a)**2) + a**4/(8*d*(-a*sin(c + d*x) + a)) + a**3*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_39():
    f = (a*sin(c + d*x) + a)**3*sec(c + d*x)**8
    F = 3*a**3*tan(c + d*x)**5/(35*d) + 2*a**3*tan(c + d*x)**3/(7*d) + 3*a**3*tan(c + d*x)/(7*d) + 3*a**3*sec(c + d*x)**5/(35*d) + 2*a*(a*sin(c + d*x) + a)**2*sec(c + d*x)**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_40():
    f = (a*sin(c + d*x) + a)**8*cos(c + d*x)**5
    F = 4*(a*sin(c + d*x) + a)**11/(11*a**3*d) - (a*sin(c + d*x) + a)**12/(3*a**4*d) + (a*sin(c + d*x) + a)**13/(13*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_41():
    f = (a*sin(c + d*x) + a)**8*cos(c + d*x)**4
    F = 4199*a**8*x/1024 + 4199*a**8*sin(c + d*x)*cos(c + d*x)**3/(1536*d) + 4199*a**8*sin(c + d*x)*cos(c + d*x)/(1024*d) - 4199*a**8*cos(c + d*x)**5/(1920*d) - 323*a**3*(a*sin(c + d*x) + a)**5*cos(c + d*x)**5/(1320*d) - 19*a**2*(a*sin(c + d*x) + a)**6*cos(c + d*x)**5/(132*d) - 4199*a**2*(a**2*sin(c + d*x) + a**2)**3*cos(c + d*x)**5/(6336*d) - a*(a*sin(c + d*x) + a)**7*cos(c + d*x)**5/(12*d) - 323*(a**2*sin(c + d*x) + a**2)**4*cos(c + d*x)**5/(792*d) - 4199*(a**4*sin(c + d*x) + a**4)**2*cos(c + d*x)**5/(4032*d) - 4199*(a**8*sin(c + d*x) + a**8)*cos(c + d*x)**5/(2688*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_42():
    f = (a*sin(c + d*x) + a)**8*cos(c + d*x)**3
    F = (a*sin(c + d*x) + a)**10/(5*a**2*d) - (a*sin(c + d*x) + a)**11/(11*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_43():
    f = (a*sin(c + d*x) + a)**8*cos(c + d*x)**2
    F = 2431*a**8*x/256 + 2431*a**8*sin(c + d*x)*cos(c + d*x)/(256*d) - 2431*a**8*cos(c + d*x)**3/(384*d) - 17*a**3*(a*sin(c + d*x) + a)**5*cos(c + d*x)**3/(48*d) - 17*a**2*(a*sin(c + d*x) + a)**6*cos(c + d*x)**3/(90*d) - 2431*a**2*(a**2*sin(c + d*x) + a**2)**3*cos(c + d*x)**3/(2016*d) - a*(a*sin(c + d*x) + a)**7*cos(c + d*x)**3/(10*d) - 221*(a**2*sin(c + d*x) + a**2)**4*cos(c + d*x)**3/(336*d) - 2431*(a**4*sin(c + d*x) + a**4)**2*cos(c + d*x)**3/(1120*d) - 2431*(a**8*sin(c + d*x) + a**8)*cos(c + d*x)**3/(640*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_44():
    f = (a*sin(c + d*x) + a)**8*cos(c + d*x)
    F = (a*sin(c + d*x) + a)**9/(9*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_45():
    f = (a*sin(c + d*x) + a)**8*sec(c + d*x)
    F = -128*a**8*log(1 - sin(c + d*x))/d - 64*a**8*sin(c + d*x)/d - 16*a**5*(a*sin(c + d*x) + a)**3/(3*d) - 4*a**3*(a*sin(c + d*x) + a)**5/(5*d) - a**2*(a*sin(c + d*x) + a)**6/(3*d) - a*(a*sin(c + d*x) + a)**7/(7*d) - 2*(a**2*sin(c + d*x) + a**2)**4/d - 16*(a**4*sin(c + d*x) + a**4)**2/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_46():
    f = (a*sin(c + d*x) + a)**8*sec(c + d*x)**2
    F = 143*a**16*cos(c + d*x)**7/(2*d*(-a**8*sin(c + d*x) + a**8)) + 2*a**15*cos(c + d*x)**13/(d*(-a*sin(c + d*x) + a)**7) + 286*a**14*cos(c + d*x)**9/(3*d*(-a**2*sin(c + d*x) + a**2)**3) + 26*a**13*cos(c + d*x)**11/(d*(-a*sin(c + d*x) + a)**5) - 3003*a**8*x/16 - 1001*a**8*sin(c + d*x)*cos(c + d*x)**3/(8*d) - 3003*a**8*sin(c + d*x)*cos(c + d*x)/(16*d) + 1001*a**8*cos(c + d*x)**5/(10*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_47():
    f = (a*sin(c + d*x) + a)**8*sec(c + d*x)**3
    F = 64*a**9/(d*(-a*sin(c + d*x) + a)) + 192*a**8*log(1 - sin(c + d*x))/d + a**8*sin(c + d*x)**5/(5*d) + 2*a**8*sin(c + d*x)**4/d + 10*a**8*sin(c + d*x)**3/d + 36*a**8*sin(c + d*x)**2/d + 129*a**8*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_48():
    f = (a*sin(c + d*x) + a)**8*sec(c + d*x)**4
    F = -231*a**16*cos(c + d*x)**5/(4*d*(-a**8*sin(c + d*x) + a**8)) + 2*a**15*cos(c + d*x)**11/(3*d*(-a*sin(c + d*x) + a)**7) - 66*a**14*cos(c + d*x)**7/(d*(-a**2*sin(c + d*x) + a**2)**3) - 22*a**13*cos(c + d*x)**9/(3*d*(-a*sin(c + d*x) + a)**5) + 1155*a**8*x/8 + 1155*a**8*sin(c + d*x)*cos(c + d*x)/(8*d) - 385*a**8*cos(c + d*x)**3/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_49():
    f = (a*sin(c + d*x) + a)**8*sec(c + d*x)**5
    F = 16*a**10/(d*(-a*sin(c + d*x) + a)**2) - 80*a**9/(d*(-a*sin(c + d*x) + a)) - 80*a**8*log(1 - sin(c + d*x))/d - a**8*sin(c + d*x)**3/(3*d) - 4*a**8*sin(c + d*x)**2/d - 31*a**8*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_50():
    f = cos(c + d*x)**6/(a*sin(c + d*x) + a)
    F = 3*x/(8*a) + sin(c + d*x)*cos(c + d*x)**3/(4*a*d) + 3*sin(c + d*x)*cos(c + d*x)/(8*a*d) + cos(c + d*x)**5/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_51():
    f = cos(c + d*x)**5/(a*sin(c + d*x) + a)
    F = -2*(-a*sin(c + d*x) + a)**3/(3*a**4*d) + (-a*sin(c + d*x) + a)**4/(4*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_52():
    f = cos(c + d*x)**4/(a*sin(c + d*x) + a)
    F = x/(2*a) + sin(c + d*x)*cos(c + d*x)/(2*a*d) + cos(c + d*x)**3/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_53():
    f = cos(c + d*x)**3/(a*sin(c + d*x) + a)
    F = -sin(c + d*x)**2/(2*a*d) + sin(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_54():
    f = cos(c + d*x)**2/(a*sin(c + d*x) + a)
    F = x/a + cos(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_55():
    f = cos(c + d*x)/(a*sin(c + d*x) + a)
    F = log(sin(c + d*x) + 1)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_56():
    f = sec(c + d*x)/(a*sin(c + d*x) + a)
    F = -1/(2*d*(a*sin(c + d*x) + a)) + atanh(sin(c + d*x))/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_57():
    f = sec(c + d*x)**2/(a*sin(c + d*x) + a)
    F = -sec(c + d*x)/(3*d*(a*sin(c + d*x) + a)) + 2*tan(c + d*x)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_58():
    f = sec(c + d*x)**3/(a*sin(c + d*x) + a)
    F = -a/(8*d*(a*sin(c + d*x) + a)**2) - 1/(4*d*(a*sin(c + d*x) + a)) + 1/(8*d*(-a*sin(c + d*x) + a)) + 3*atanh(sin(c + d*x))/(8*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_59():
    f = sec(c + d*x)**4/(a*sin(c + d*x) + a)
    F = -sec(c + d*x)**3/(5*d*(a*sin(c + d*x) + a)) + 4*tan(c + d*x)**3/(15*a*d) + 4*tan(c + d*x)/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_60():
    f = sec(c + d*x)**5/(a*sin(c + d*x) + a)
    F = -a**2/(24*d*(a*sin(c + d*x) + a)**3) - 3*a/(32*d*(a*sin(c + d*x) + a)**2) + a/(32*d*(-a*sin(c + d*x) + a)**2) - 3/(16*d*(a*sin(c + d*x) + a)) + 1/(8*d*(-a*sin(c + d*x) + a)) + 5*atanh(sin(c + d*x))/(16*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_61():
    f = cos(c + d*x)**8/(a*sin(c + d*x) + a)**2
    F = cos(c + d*x)**7/(6*d*(a**2*sin(c + d*x) + a**2)) + 7*x/(16*a**2) + 7*sin(c + d*x)*cos(c + d*x)**3/(24*a**2*d) + 7*sin(c + d*x)*cos(c + d*x)/(16*a**2*d) + 7*cos(c + d*x)**5/(30*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_62():
    f = cos(c + d*x)**7/(a*sin(c + d*x) + a)**2
    F = -(-a*sin(c + d*x) + a)**4/(2*a**6*d) + (-a*sin(c + d*x) + a)**5/(5*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_63():
    f = cos(c + d*x)**6/(a*sin(c + d*x) + a)**2
    F = cos(c + d*x)**5/(4*d*(a**2*sin(c + d*x) + a**2)) + 5*x/(8*a**2) + 5*sin(c + d*x)*cos(c + d*x)/(8*a**2*d) + 5*cos(c + d*x)**3/(12*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_64():
    f = cos(c + d*x)**5/(a*sin(c + d*x) + a)**2
    F = -(-a*sin(c + d*x) + a)**3/(3*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_65():
    f = cos(c + d*x)**4/(a*sin(c + d*x) + a)**2
    F = cos(c + d*x)**3/(2*d*(a**2*sin(c + d*x) + a**2)) + 3*x/(2*a**2) + 3*cos(c + d*x)/(2*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_66():
    f = cos(c + d*x)**3/(a*sin(c + d*x) + a)**2
    F = 2*log(sin(c + d*x) + 1)/(a**2*d) - sin(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_67():
    f = cos(c + d*x)**2/(a*sin(c + d*x) + a)**2
    F = -2*cos(c + d*x)/(d*(a**2*sin(c + d*x) + a**2)) - x/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_68():
    f = cos(c + d*x)/(a*sin(c + d*x) + a)**2
    F = -1/(d*(a**2*sin(c + d*x) + a**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_69():
    f = sec(c + d*x)/(a*sin(c + d*x) + a)**2
    F = -1/(4*d*(a**2*sin(c + d*x) + a**2)) - 1/(4*d*(a*sin(c + d*x) + a)**2) + atanh(sin(c + d*x))/(4*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_70():
    f = sec(c + d*x)**2/(a*sin(c + d*x) + a)**2
    F = -sec(c + d*x)/(5*d*(a**2*sin(c + d*x) + a**2)) - sec(c + d*x)/(5*d*(a*sin(c + d*x) + a)**2) + 2*tan(c + d*x)/(5*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_71():
    f = sec(c + d*x)**3/(a*sin(c + d*x) + a)**2
    F = -a/(12*d*(a*sin(c + d*x) + a)**3) - 3/(16*d*(a**2*sin(c + d*x) + a**2)) + 1/(16*d*(-a**2*sin(c + d*x) + a**2)) - 1/(8*d*(a*sin(c + d*x) + a)**2) + atanh(sin(c + d*x))/(4*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_72():
    f = sec(c + d*x)**4/(a*sin(c + d*x) + a)**2
    F = -sec(c + d*x)**3/(7*d*(a**2*sin(c + d*x) + a**2)) - sec(c + d*x)**3/(7*d*(a*sin(c + d*x) + a)**2) + 4*tan(c + d*x)**3/(21*a**2*d) + 4*tan(c + d*x)/(7*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_73():
    f = sec(c + d*x)**5/(a*sin(c + d*x) + a)**2
    F = -a**2/(32*d*(a*sin(c + d*x) + a)**4) - a/(16*d*(a*sin(c + d*x) + a)**3) - 5/(32*d*(a**2*sin(c + d*x) + a**2)) + 5/(64*d*(-a**2*sin(c + d*x) + a**2)) - 3/(32*d*(a*sin(c + d*x) + a)**2) + 1/(64*d*(-a*sin(c + d*x) + a)**2) + 15*atanh(sin(c + d*x))/(64*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_74():
    f = cos(c + d*x)**8/(a*sin(c + d*x) + a)**3
    F = 2*cos(c + d*x)**7/(3*a*d*(a*sin(c + d*x) + a)**2) + 7*x/(8*a**3) + 7*sin(c + d*x)*cos(c + d*x)**3/(12*a**3*d) + 7*sin(c + d*x)*cos(c + d*x)/(8*a**3*d) + 7*cos(c + d*x)**5/(15*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_75():
    f = cos(c + d*x)**7/(a*sin(c + d*x) + a)**3
    F = -(-a*sin(c + d*x) + a)**4/(4*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_76():
    f = cos(c + d*x)**6/(a*sin(c + d*x) + a)**3
    F = 2*cos(c + d*x)**5/(a*d*(a*sin(c + d*x) + a)**2) + 5*x/(2*a**3) + 5*sin(c + d*x)*cos(c + d*x)/(2*a**3*d) + 5*cos(c + d*x)**3/(3*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_77():
    f = cos(c + d*x)**5/(a*sin(c + d*x) + a)**3
    F = 4*log(sin(c + d*x) + 1)/(a**3*d) + sin(c + d*x)**2/(2*a**3*d) - 3*sin(c + d*x)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_78():
    f = cos(c + d*x)**4/(a*sin(c + d*x) + a)**3
    F = -2*cos(c + d*x)**3/(a*d*(a*sin(c + d*x) + a)**2) - 3*x/a**3 - 3*cos(c + d*x)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_79():
    f = cos(c + d*x)**3/(a*sin(c + d*x) + a)**3
    F = -2/(d*(a**3*sin(c + d*x) + a**3)) - log(sin(c + d*x) + 1)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_80():
    f = cos(c + d*x)**2/(a*sin(c + d*x) + a)**3
    F = -cos(c + d*x)**3/(3*d*(a*sin(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_81():
    f = cos(c + d*x)/(a*sin(c + d*x) + a)**3
    F = -1/(2*a*d*(a*sin(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_82():
    f = sec(c + d*x)/(a*sin(c + d*x) + a)**3
    F = -1/(8*d*(a**3*sin(c + d*x) + a**3)) - 1/(6*d*(a*sin(c + d*x) + a)**3) - 1/(8*a*d*(a*sin(c + d*x) + a)**2) + atanh(sin(c + d*x))/(8*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_83():
    f = sec(c + d*x)**2/(a*sin(c + d*x) + a)**3
    F = -4*sec(c + d*x)/(35*d*(a**3*sin(c + d*x) + a**3)) - sec(c + d*x)/(7*d*(a*sin(c + d*x) + a)**3) - 4*sec(c + d*x)/(35*a*d*(a*sin(c + d*x) + a)**2) + 8*tan(c + d*x)/(35*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_84():
    f = sec(c + d*x)**3/(a*sin(c + d*x) + a)**3
    F = -a/(16*d*(a*sin(c + d*x) + a)**4) - 1/(8*d*(a**3*sin(c + d*x) + a**3)) + 1/(32*d*(-a**3*sin(c + d*x) + a**3)) - 1/(12*d*(a*sin(c + d*x) + a)**3) - 3/(32*a*d*(a*sin(c + d*x) + a)**2) + 5*atanh(sin(c + d*x))/(32*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_85():
    f = sec(c + d*x)**4/(a*sin(c + d*x) + a)**3
    F = -2*sec(c + d*x)**3/(21*d*(a**3*sin(c + d*x) + a**3)) - sec(c + d*x)**3/(9*d*(a*sin(c + d*x) + a)**3) - 2*sec(c + d*x)**3/(21*a*d*(a*sin(c + d*x) + a)**2) + 8*tan(c + d*x)**3/(63*a**3*d) + 8*tan(c + d*x)/(21*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_86():
    f = sec(c + d*x)**5/(a*sin(c + d*x) + a)**3
    F = -a**2/(40*d*(a*sin(c + d*x) + a)**5) - 3*a/(64*d*(a*sin(c + d*x) + a)**4) - 15/(128*d*(a**3*sin(c + d*x) + a**3)) + 3/(64*d*(-a**3*sin(c + d*x) + a**3)) - 1/(16*d*(a*sin(c + d*x) + a)**3) - 5/(64*a*d*(a*sin(c + d*x) + a)**2) + 1/(128*a*d*(-a*sin(c + d*x) + a)**2) + 21*atanh(sin(c + d*x))/(128*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_87():
    f = cos(c + d*x)**8/(a*sin(c + d*x) + a)**8
    F = 2*cos(c + d*x)/(d*(a**8*sin(c + d*x) + a**8)) - 2*cos(c + d*x)**7/(7*a*d*(a*sin(c + d*x) + a)**7) - 2*cos(c + d*x)**3/(3*a**2*d*(a**2*sin(c + d*x) + a**2)**3) + 2*cos(c + d*x)**5/(5*a**3*d*(a*sin(c + d*x) + a)**5) + x/a**8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_88():
    f = cos(c + d*x)**7/(a*sin(c + d*x) + a)**8
    F = -(-a*sin(c + d*x) + a)**4/(8*d*(a**3*sin(c + d*x) + a**3)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_89():
    f = cos(c + d*x)**6/(a*sin(c + d*x) + a)**8
    F = -cos(c + d*x)**7/(9*d*(a*sin(c + d*x) + a)**8) - cos(c + d*x)**7/(63*a*d*(a*sin(c + d*x) + a)**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_90():
    f = cos(c + d*x)**5/(a*sin(c + d*x) + a)**8
    F = 1/(d*(a**2*sin(c + d*x) + a**2)**4) - 4/(5*a**3*d*(a*sin(c + d*x) + a)**5) - 1/(3*a**5*d*(a*sin(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_91():
    f = cos(c + d*x)**4/(a*sin(c + d*x) + a)**8
    F = -cos(c + d*x)**5/(11*d*(a*sin(c + d*x) + a)**8) - cos(c + d*x)**5/(33*a*d*(a*sin(c + d*x) + a)**7) - 2*cos(c + d*x)**5/(231*a**2*d*(a*sin(c + d*x) + a)**6) - 2*cos(c + d*x)**5/(1155*a**3*d*(a*sin(c + d*x) + a)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_92():
    f = cos(c + d*x)**3/(a*sin(c + d*x) + a)**8
    F = -1/(3*a**2*d*(a*sin(c + d*x) + a)**6) + 1/(5*a**3*d*(a*sin(c + d*x) + a)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_93():
    f = cos(c + d*x)**2/(a*sin(c + d*x) + a)**8
    F = -8*cos(c + d*x)**3/(3003*d*(a**2*sin(c + d*x) + a**2)**4) - cos(c + d*x)**3/(13*d*(a*sin(c + d*x) + a)**8) - 5*cos(c + d*x)**3/(143*a*d*(a*sin(c + d*x) + a)**7) - 8*cos(c + d*x)**3/(9009*a**2*d*(a**2*sin(c + d*x) + a**2)**3) - 20*cos(c + d*x)**3/(1287*a**2*d*(a*sin(c + d*x) + a)**6) - 20*cos(c + d*x)**3/(3003*a**3*d*(a*sin(c + d*x) + a)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_94():
    f = cos(c + d*x)/(a*sin(c + d*x) + a)**8
    F = -1/(7*a*d*(a*sin(c + d*x) + a)**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_95():
    f = sec(c + d*x)/(a*sin(c + d*x) + a)**8
    F = -1/(256*d*(a**8*sin(c + d*x) + a**8)) - 1/(256*d*(a**4*sin(c + d*x) + a**4)**2) - 1/(128*d*(a**2*sin(c + d*x) + a**2)**4) - 1/(16*d*(a*sin(c + d*x) + a)**8) - 1/(28*a*d*(a*sin(c + d*x) + a)**7) - 1/(48*a**2*d*(a*sin(c + d*x) + a)**6) - 1/(80*a**3*d*(a*sin(c + d*x) + a)**5) - 1/(192*a**5*d*(a*sin(c + d*x) + a)**3) + atanh(sin(c + d*x))/(256*a**8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_96():
    f = sec(c + d*x)**2/(a*sin(c + d*x) + a)**8
    F = -64*sec(c + d*x)/(12155*d*(a**8*sin(c + d*x) + a**8)) - 64*sec(c + d*x)/(12155*d*(a**4*sin(c + d*x) + a**4)**2) - 112*sec(c + d*x)/(12155*d*(a**2*sin(c + d*x) + a**2)**4) - sec(c + d*x)/(17*d*(a*sin(c + d*x) + a)**8) - 3*sec(c + d*x)/(85*a*d*(a*sin(c + d*x) + a)**7) - 16*sec(c + d*x)/(2431*a**2*d*(a**2*sin(c + d*x) + a**2)**3) - 24*sec(c + d*x)/(1105*a**2*d*(a*sin(c + d*x) + a)**6) - 168*sec(c + d*x)/(12155*a**3*d*(a*sin(c + d*x) + a)**5) + 128*tan(c + d*x)/(12155*a**8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_97():
    f = sec(c + d*x)**3/(a*sin(c + d*x) + a)**8
    F = -a/(36*d*(a*sin(c + d*x) + a)**9) - 9/(1024*d*(a**8*sin(c + d*x) + a**8)) + 1/(1024*d*(-a**8*sin(c + d*x) + a**8)) - 1/(128*d*(a**4*sin(c + d*x) + a**4)**2) - 3/(256*d*(a**2*sin(c + d*x) + a**2)**4) - 1/(32*d*(a*sin(c + d*x) + a)**8) - 3/(112*a*d*(a*sin(c + d*x) + a)**7) - 1/(48*a**2*d*(a*sin(c + d*x) + a)**6) - 1/(64*a**3*d*(a*sin(c + d*x) + a)**5) - 7/(768*a**5*d*(a*sin(c + d*x) + a)**3) + 5*atanh(sin(c + d*x))/(512*a**8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_98():
    f = sec(c + d*x)**4/(a*sin(c + d*x) + a)**8
    F = -32*sec(c + d*x)**3/(4199*d*(a**8*sin(c + d*x) + a**8)) - 32*sec(c + d*x)**3/(4199*d*(a**4*sin(c + d*x) + a**4)**2) - 48*sec(c + d*x)**3/(4199*d*(a**2*sin(c + d*x) + a**2)**4) - sec(c + d*x)**3/(19*d*(a*sin(c + d*x) + a)**8) - 11*sec(c + d*x)**3/(323*a*d*(a*sin(c + d*x) + a)**7) - 112*sec(c + d*x)**3/(12597*a**2*d*(a**2*sin(c + d*x) + a**2)**3) - 22*sec(c + d*x)**3/(969*a**2*d*(a*sin(c + d*x) + a)**6) - 66*sec(c + d*x)**3/(4199*a**3*d*(a*sin(c + d*x) + a)**5) + 128*tan(c + d*x)**3/(12597*a**8*d) + 128*tan(c + d*x)/(4199*a**8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_99():
    f = sec(c + d*x)**5/(a*sin(c + d*x) + a)**8
    F = -a**2/(80*d*(a*sin(c + d*x) + a)**10) - a/(48*d*(a*sin(c + d*x) + a)**9) - 55/(4096*d*(a**8*sin(c + d*x) + a**8)) + 11/(4096*d*(-a**8*sin(c + d*x) + a**8)) - 45/(4096*d*(a**4*sin(c + d*x) + a**4)**2) + 1/(4096*d*(-a**4*sin(c + d*x) + a**4)**2) - 7/(512*d*(a**2*sin(c + d*x) + a**2)**4) - 3/(128*d*(a*sin(c + d*x) + a)**8) - 5/(224*a*d*(a*sin(c + d*x) + a)**7) - 5/(256*a**2*d*(a*sin(c + d*x) + a)**6) - 21/(1280*a**3*d*(a*sin(c + d*x) + a)**5) - 3/(256*a**5*d*(a*sin(c + d*x) + a)**3) + 33*atanh(sin(c + d*x))/(2048*a**8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_100():
    f = sqrt(a*sin(c + d*x) + a)*cos(c + d*x)**7
    F = 16*(a*sin(c + d*x) + a)**(sympy.S(9)/2)/(9*a**4*d) - 24*(a*sin(c + d*x) + a)**(sympy.S(11)/2)/(11*a**5*d) + 12*(a*sin(c + d*x) + a)**(sympy.S(13)/2)/(13*a**6*d) - 2*(a*sin(c + d*x) + a)**(sympy.S(15)/2)/(15*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_101():
    f = sqrt(a*sin(c + d*x) + a)*cos(c + d*x)**6
    F = -256*a**4*cos(c + d*x)**7/(3003*d*(a*sin(c + d*x) + a)**(sympy.S(7)/2)) - 64*a**3*cos(c + d*x)**7/(429*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) - 24*a**2*cos(c + d*x)**7/(143*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 2*a*cos(c + d*x)**7/(13*d*sqrt(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_102():
    f = sqrt(a*sin(c + d*x) + a)*cos(c + d*x)**5
    F = 8*(a*sin(c + d*x) + a)**(sympy.S(7)/2)/(7*a**3*d) - 8*(a*sin(c + d*x) + a)**(sympy.S(9)/2)/(9*a**4*d) + 2*(a*sin(c + d*x) + a)**(sympy.S(11)/2)/(11*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_103():
    f = sqrt(a*sin(c + d*x) + a)*cos(c + d*x)**4
    F = -64*a**3*cos(c + d*x)**5/(315*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) - 16*a**2*cos(c + d*x)**5/(63*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 2*a*cos(c + d*x)**5/(9*d*sqrt(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_104():
    f = sqrt(a*sin(c + d*x) + a)*cos(c + d*x)**3
    F = 4*(a*sin(c + d*x) + a)**(sympy.S(5)/2)/(5*a**2*d) - 2*(a*sin(c + d*x) + a)**(sympy.S(7)/2)/(7*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_105():
    f = sqrt(a*sin(c + d*x) + a)*cos(c + d*x)**2
    F = -8*a**2*cos(c + d*x)**3/(15*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 2*a*cos(c + d*x)**3/(5*d*sqrt(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_106():
    f = sqrt(a*sin(c + d*x) + a)*cos(c + d*x)
    F = 2*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_107():
    f = sqrt(a*sin(c + d*x) + a)*sec(c + d*x)
    F = sqrt(2)*sqrt(a)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_108():
    f = sqrt(a*sin(c + d*x) + a)*sec(c + d*x)**2
    F = -sqrt(2)*sqrt(a)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(2*d) + sqrt(a*sin(c + d*x) + a)*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_109():
    f = sqrt(a*sin(c + d*x) + a)*sec(c + d*x)**3
    F = 3*sqrt(2)*sqrt(a)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/(8*d) - 3*a/(4*d*sqrt(a*sin(c + d*x) + a)) + sqrt(a*sin(c + d*x) + a)*sec(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_110():
    f = sqrt(a*sin(c + d*x) + a)*sec(c + d*x)**4
    F = -5*sqrt(2)*sqrt(a)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(16*d) - 5*a**2*cos(c + d*x)/(8*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + 5*a*sec(c + d*x)/(6*d*sqrt(a*sin(c + d*x) + a)) + sqrt(a*sin(c + d*x) + a)*sec(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_111():
    f = sqrt(a*sin(c + d*x) + a)*sec(c + d*x)**5
    F = 35*sqrt(2)*sqrt(a)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/(128*d) - 35*a**2/(96*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + 7*a*sec(c + d*x)**2/(16*d*sqrt(a*sin(c + d*x) + a)) - 35*a/(64*d*sqrt(a*sin(c + d*x) + a)) + sqrt(a*sin(c + d*x) + a)*sec(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_112():
    f = sqrt(a*sin(c + d*x) + a)*sec(c + d*x)**6
    F = -63*sqrt(2)*sqrt(a)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(256*d) - 63*a**2*cos(c + d*x)/(128*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 21*a**2*sec(c + d*x)/(80*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + 3*a*sec(c + d*x)**3/(10*d*sqrt(a*sin(c + d*x) + a)) + 21*a*sec(c + d*x)/(32*d*sqrt(a*sin(c + d*x) + a)) + sqrt(a*sin(c + d*x) + a)*sec(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_113():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**7
    F = 16*(a*sin(c + d*x) + a)**(sympy.S(11)/2)/(11*a**4*d) - 24*(a*sin(c + d*x) + a)**(sympy.S(13)/2)/(13*a**5*d) + 4*(a*sin(c + d*x) + a)**(sympy.S(15)/2)/(5*a**6*d) - 2*(a*sin(c + d*x) + a)**(sympy.S(17)/2)/(17*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_114():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**6
    F = -4096*a**5*cos(c + d*x)**7/(45045*d*(a*sin(c + d*x) + a)**(sympy.S(7)/2)) - 1024*a**4*cos(c + d*x)**7/(6435*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) - 128*a**3*cos(c + d*x)**7/(715*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 32*a**2*cos(c + d*x)**7/(195*d*sqrt(a*sin(c + d*x) + a)) - 2*a*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)**7/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_115():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**5
    F = 8*(a*sin(c + d*x) + a)**(sympy.S(9)/2)/(9*a**3*d) - 8*(a*sin(c + d*x) + a)**(sympy.S(11)/2)/(11*a**4*d) + 2*(a*sin(c + d*x) + a)**(sympy.S(13)/2)/(13*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_116():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**4
    F = -256*a**4*cos(c + d*x)**5/(1155*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) - 64*a**3*cos(c + d*x)**5/(231*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 8*a**2*cos(c + d*x)**5/(33*d*sqrt(a*sin(c + d*x) + a)) - 2*a*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)**5/(11*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_117():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**3
    F = 4*(a*sin(c + d*x) + a)**(sympy.S(7)/2)/(7*a**2*d) - 2*(a*sin(c + d*x) + a)**(sympy.S(9)/2)/(9*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_118():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**2
    F = -64*a**3*cos(c + d*x)**3/(105*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 16*a**2*cos(c + d*x)**3/(35*d*sqrt(a*sin(c + d*x) + a)) - 2*a*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)**3/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_119():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)
    F = 2*(a*sin(c + d*x) + a)**(sympy.S(5)/2)/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_120():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)
    F = 2*sqrt(2)*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/d - 2*a*sqrt(a*sin(c + d*x) + a)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_121():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**2
    F = 2*a*sqrt(a*sin(c + d*x) + a)*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_122():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**3
    F = sqrt(2)*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/(4*d) + (a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_123():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**4
    F = -sqrt(2)*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(4*d) + a*sqrt(a*sin(c + d*x) + a)*sec(c + d*x)/(2*d) + (a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_124():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**5
    F = 15*sqrt(2)*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/(64*d) - 15*a**2/(32*d*sqrt(a*sin(c + d*x) + a)) + 5*a*sqrt(a*sin(c + d*x) + a)*sec(c + d*x)**2/(16*d) + (a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_125():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**6
    F = -7*sqrt(2)*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(32*d) - 7*a**3*cos(c + d*x)/(16*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + 7*a**2*sec(c + d*x)/(12*d*sqrt(a*sin(c + d*x) + a)) + 7*a*sqrt(a*sin(c + d*x) + a)*sec(c + d*x)**3/(30*d) + (a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_126():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**5
    F = 8*(a*sin(c + d*x) + a)**(sympy.S(11)/2)/(11*a**3*d) - 8*(a*sin(c + d*x) + a)**(sympy.S(13)/2)/(13*a**4*d) + 2*(a*sin(c + d*x) + a)**(sympy.S(15)/2)/(15*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_127():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**4
    F = -4096*a**5*cos(c + d*x)**5/(15015*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) - 1024*a**4*cos(c + d*x)**5/(3003*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 128*a**3*cos(c + d*x)**5/(429*d*sqrt(a*sin(c + d*x) + a)) - 32*a**2*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)**5/(143*d) - 2*a*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**5/(13*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_128():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**3
    F = 4*(a*sin(c + d*x) + a)**(sympy.S(9)/2)/(9*a**2*d) - 2*(a*sin(c + d*x) + a)**(sympy.S(11)/2)/(11*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_129():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**2
    F = -256*a**4*cos(c + d*x)**3/(315*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 64*a**3*cos(c + d*x)**3/(105*d*sqrt(a*sin(c + d*x) + a)) - 8*a**2*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)**3/(21*d) - 2*a*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**3/(9*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_130():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)
    F = 2*(a*sin(c + d*x) + a)**(sympy.S(7)/2)/(7*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_131():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)
    F = 4*sqrt(2)*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/d - 4*a**2*sqrt(a*sin(c + d*x) + a)/d - 2*a*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_132():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**2
    F = 8*a**2*sqrt(a*sin(c + d*x) + a)*sec(c + d*x)/d - 2*a*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_133():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**3
    F = -sqrt(2)*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/(2*d) + a*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**2/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_134():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**4
    F = 2*a*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_135():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**5
    F = 3*sqrt(2)*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/(32*d) + 3*a*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**2/(16*d) + (a*sin(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_136():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**6
    F = -sqrt(2)*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(8*d) + a**2*sqrt(a*sin(c + d*x) + a)*sec(c + d*x)/(4*d) + a*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**3/(6*d) + (a*sin(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_137():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**7
    F = 35*sqrt(2)*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/(256*d) - 35*a**3/(128*d*sqrt(a*sin(c + d*x) + a)) + 35*a**2*sqrt(a*sin(c + d*x) + a)*sec(c + d*x)**2/(192*d) + 7*a*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**4/(48*d) + (a*sin(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**6/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_138():
    f = (a*sin(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**7
    F = 16*(a*sin(c + d*x) + a)**(sympy.S(15)/2)/(15*a**4*d) - 24*(a*sin(c + d*x) + a)**(sympy.S(17)/2)/(17*a**5*d) + 12*(a*sin(c + d*x) + a)**(sympy.S(19)/2)/(19*a**6*d) - 2*(a*sin(c + d*x) + a)**(sympy.S(21)/2)/(21*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_139():
    f = (a*sin(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**6
    F = -131072*a**7*cos(c + d*x)**7/(969969*d*(a*sin(c + d*x) + a)**(sympy.S(7)/2)) - 32768*a**6*cos(c + d*x)**7/(138567*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) - 12288*a**5*cos(c + d*x)**7/(46189*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 1024*a**4*cos(c + d*x)**7/(4199*d*sqrt(a*sin(c + d*x) + a)) - 64*a**3*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)**7/(323*d) - 48*a**2*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**7/(323*d) - 2*a*(a*sin(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**7/(19*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_140():
    f = (a*sin(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**5
    F = 8*(a*sin(c + d*x) + a)**(sympy.S(13)/2)/(13*a**3*d) - 8*(a*sin(c + d*x) + a)**(sympy.S(15)/2)/(15*a**4*d) + 2*(a*sin(c + d*x) + a)**(sympy.S(17)/2)/(17*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_141():
    f = (a*sin(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**4
    F = -16384*a**6*cos(c + d*x)**5/(45045*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) - 4096*a**5*cos(c + d*x)**5/(9009*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 512*a**4*cos(c + d*x)**5/(1287*d*sqrt(a*sin(c + d*x) + a)) - 128*a**3*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)**5/(429*d) - 8*a**2*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**5/(39*d) - 2*a*(a*sin(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**5/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_142():
    f = (a*sin(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**3
    F = 4*(a*sin(c + d*x) + a)**(sympy.S(11)/2)/(11*a**2*d) - 2*(a*sin(c + d*x) + a)**(sympy.S(13)/2)/(13*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_143():
    f = (a*sin(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**2
    F = -4096*a**5*cos(c + d*x)**3/(3465*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 1024*a**4*cos(c + d*x)**3/(1155*d*sqrt(a*sin(c + d*x) + a)) - 128*a**3*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)**3/(231*d) - 32*a**2*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**3/(99*d) - 2*a*(a*sin(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**3/(11*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_144():
    f = (a*sin(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)
    F = 2*(a*sin(c + d*x) + a)**(sympy.S(9)/2)/(9*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_145():
    f = (a*sin(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)
    F = 8*sqrt(2)*a**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/d - 8*a**3*sqrt(a*sin(c + d*x) + a)/d - 4*a**2*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(3*d) - 2*a*(a*sin(c + d*x) + a)**(sympy.S(5)/2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_146():
    f = (a*sin(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**2
    F = 64*a**3*sqrt(a*sin(c + d*x) + a)*sec(c + d*x)/(3*d) - 16*a**2*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)/(3*d) - 2*a*(a*sin(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_147():
    f = (a*sin(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**3
    F = -3*sqrt(2)*a**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/d + 3*a**3*sqrt(a*sin(c + d*x) + a)/d + a*(a*sin(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**2/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_148():
    f = (a*sin(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**4
    F = -8*a**2*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**3/(3*d) + 2*a*(a*sin(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**3/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_149():
    f = (a*sin(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**5
    F = -sqrt(2)*a**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/(16*d) - a**2*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**2/(8*d) + a*(a*sin(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**4/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_150():
    f = (a*sin(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**6
    F = 2*a*(a*sin(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_151():
    f = (a*sin(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**7
    F = 5*sqrt(2)*a**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/(128*d) + 5*a**2*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**2/(64*d) + 5*a*(a*sin(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**4/(48*d) + (a*sin(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**6/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_152():
    f = (a*sin(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**8
    F = -sqrt(2)*a**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(16*d) + a**3*sqrt(a*sin(c + d*x) + a)*sec(c + d*x)/(8*d) + a**2*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**3/(12*d) + a*(a*sin(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**5/(10*d) + (a*sin(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_153():
    f = (a*sin(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**9
    F = 315*sqrt(2)*a**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/(4096*d) - 315*a**4/(2048*d*sqrt(a*sin(c + d*x) + a)) + 105*a**3*sqrt(a*sin(c + d*x) + a)*sec(c + d*x)**2/(1024*d) + 21*a**2*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**4/(256*d) + 3*a*(a*sin(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**6/(32*d) + (a*sin(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**8/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_154():
    f = (a*sin(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**10
    F = -11*sqrt(2)*a**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(128*d) - 11*a**5*cos(c + d*x)/(64*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + 11*a**4*sec(c + d*x)/(48*d*sqrt(a*sin(c + d*x) + a)) + 11*a**3*sqrt(a*sin(c + d*x) + a)*sec(c + d*x)**3/(120*d) + 11*a**2*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**5/(140*d) + 11*a*(a*sin(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**7/(126*d) + (a*sin(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**9/(9*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_155():
    f = cos(c + d*x)**7/sqrt(a*sin(c + d*x) + a)
    F = 16*(a*sin(c + d*x) + a)**(sympy.S(7)/2)/(7*a**4*d) - 8*(a*sin(c + d*x) + a)**(sympy.S(9)/2)/(3*a**5*d) + 12*(a*sin(c + d*x) + a)**(sympy.S(11)/2)/(11*a**6*d) - 2*(a*sin(c + d*x) + a)**(sympy.S(13)/2)/(13*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_156():
    f = cos(c + d*x)**6/sqrt(a*sin(c + d*x) + a)
    F = -64*a**3*cos(c + d*x)**7/(693*d*(a*sin(c + d*x) + a)**(sympy.S(7)/2)) - 16*a**2*cos(c + d*x)**7/(99*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) - 2*a*cos(c + d*x)**7/(11*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_157():
    f = cos(c + d*x)**5/sqrt(a*sin(c + d*x) + a)
    F = 8*(a*sin(c + d*x) + a)**(sympy.S(5)/2)/(5*a**3*d) - 8*(a*sin(c + d*x) + a)**(sympy.S(7)/2)/(7*a**4*d) + 2*(a*sin(c + d*x) + a)**(sympy.S(9)/2)/(9*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_158():
    f = cos(c + d*x)**4/sqrt(a*sin(c + d*x) + a)
    F = -8*a**2*cos(c + d*x)**5/(35*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) - 2*a*cos(c + d*x)**5/(7*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_159():
    f = cos(c + d*x)**3/sqrt(a*sin(c + d*x) + a)
    F = 4*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(3*a**2*d) - 2*(a*sin(c + d*x) + a)**(sympy.S(5)/2)/(5*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_160():
    f = cos(c + d*x)**2/sqrt(a*sin(c + d*x) + a)
    F = -2*a*cos(c + d*x)**3/(3*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_161():
    f = cos(c + d*x)/sqrt(a*sin(c + d*x) + a)
    F = 2*sqrt(a*sin(c + d*x) + a)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_162():
    f = sec(c + d*x)/sqrt(a*sin(c + d*x) + a)
    F = -1/(d*sqrt(a*sin(c + d*x) + a)) + sqrt(2)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_163():
    f = sec(c + d*x)**2/sqrt(a*sin(c + d*x) + a)
    F = -3*a*cos(c + d*x)/(4*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + sec(c + d*x)/(d*sqrt(a*sin(c + d*x) + a)) - 3*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(8*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_164():
    f = sec(c + d*x)**3/sqrt(a*sin(c + d*x) + a)
    F = -5*a/(12*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + sec(c + d*x)**2/(2*d*sqrt(a*sin(c + d*x) + a)) - 5/(8*d*sqrt(a*sin(c + d*x) + a)) + 5*sqrt(2)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/(16*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_165():
    f = sec(c + d*x)**4/sqrt(a*sin(c + d*x) + a)
    F = -35*a*cos(c + d*x)/(64*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 7*a*sec(c + d*x)/(24*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + sec(c + d*x)**3/(3*d*sqrt(a*sin(c + d*x) + a)) + 35*sec(c + d*x)/(48*d*sqrt(a*sin(c + d*x) + a)) - 35*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(128*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_166():
    f = sec(c + d*x)**5/sqrt(a*sin(c + d*x) + a)
    F = -9*a*sec(c + d*x)**2/(40*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 21*a/(64*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + sec(c + d*x)**4/(4*d*sqrt(a*sin(c + d*x) + a)) + 63*sec(c + d*x)**2/(160*d*sqrt(a*sin(c + d*x) + a)) - 63/(128*d*sqrt(a*sin(c + d*x) + a)) + 63*sqrt(2)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/(256*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_167():
    f = sec(c + d*x)**6/sqrt(a*sin(c + d*x) + a)
    F = -231*a*cos(c + d*x)/(512*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 11*a*sec(c + d*x)**3/(60*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 77*a*sec(c + d*x)/(320*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + sec(c + d*x)**5/(5*d*sqrt(a*sin(c + d*x) + a)) + 11*sec(c + d*x)**3/(40*d*sqrt(a*sin(c + d*x) + a)) + 77*sec(c + d*x)/(128*d*sqrt(a*sin(c + d*x) + a)) - 231*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(1024*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_168():
    f = cos(c + d*x)**7/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = 16*(a*sin(c + d*x) + a)**(sympy.S(5)/2)/(5*a**4*d) - 24*(a*sin(c + d*x) + a)**(sympy.S(7)/2)/(7*a**5*d) + 4*(a*sin(c + d*x) + a)**(sympy.S(9)/2)/(3*a**6*d) - 2*(a*sin(c + d*x) + a)**(sympy.S(11)/2)/(11*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_169():
    f = cos(c + d*x)**6/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = -8*a**2*cos(c + d*x)**7/(63*d*(a*sin(c + d*x) + a)**(sympy.S(7)/2)) - 2*a*cos(c + d*x)**7/(9*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_170():
    f = cos(c + d*x)**5/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = 8*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(3*a**3*d) - 8*(a*sin(c + d*x) + a)**(sympy.S(5)/2)/(5*a**4*d) + 2*(a*sin(c + d*x) + a)**(sympy.S(7)/2)/(7*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_171():
    f = cos(c + d*x)**4/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = -2*a*cos(c + d*x)**5/(5*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_172():
    f = cos(c + d*x)**3/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = 4*sqrt(a*sin(c + d*x) + a)/(a**2*d) - 2*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(3*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_173():
    f = cos(c + d*x)**2/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*cos(c + d*x)/(a*d*sqrt(a*sin(c + d*x) + a)) - 2*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_174():
    f = cos(c + d*x)/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = -2/(a*d*sqrt(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_175():
    f = sec(c + d*x)/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = -1/(3*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 1/(2*a*d*sqrt(a*sin(c + d*x) + a)) + sqrt(2)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_176():
    f = sec(c + d*x)**2/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = -15*cos(c + d*x)/(32*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - sec(c + d*x)/(4*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + 5*sec(c + d*x)/(8*a*d*sqrt(a*sin(c + d*x) + a)) - 15*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(64*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_177():
    f = sec(c + d*x)**3/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = -sec(c + d*x)**2/(5*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 7/(24*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + 7*sec(c + d*x)**2/(20*a*d*sqrt(a*sin(c + d*x) + a)) - 7/(16*a*d*sqrt(a*sin(c + d*x) + a)) + 7*sqrt(2)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/(32*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_178():
    f = sec(c + d*x)**4/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = -105*cos(c + d*x)/(256*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - sec(c + d*x)**3/(6*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 7*sec(c + d*x)/(32*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + sec(c + d*x)**3/(4*a*d*sqrt(a*sin(c + d*x) + a)) + 35*sec(c + d*x)/(64*a*d*sqrt(a*sin(c + d*x) + a)) - 105*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(512*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_179():
    f = sec(c + d*x)**5/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = -sec(c + d*x)**4/(7*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 99*sec(c + d*x)**2/(560*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 33/(128*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + 11*sec(c + d*x)**4/(56*a*d*sqrt(a*sin(c + d*x) + a)) + 99*sec(c + d*x)**2/(320*a*d*sqrt(a*sin(c + d*x) + a)) - 99/(256*a*d*sqrt(a*sin(c + d*x) + a)) + 99*sqrt(2)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/(512*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_180():
    f = sec(c + d*x)**6/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = -3003*cos(c + d*x)/(8192*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - sec(c + d*x)**5/(8*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 143*sec(c + d*x)**3/(960*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 1001*sec(c + d*x)/(5120*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + 13*sec(c + d*x)**5/(80*a*d*sqrt(a*sin(c + d*x) + a)) + 143*sec(c + d*x)**3/(640*a*d*sqrt(a*sin(c + d*x) + a)) + 1001*sec(c + d*x)/(2048*a*d*sqrt(a*sin(c + d*x) + a)) - 3003*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(16384*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_181():
    f = cos(c + d*x)**10/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = -64*a**3*cos(c + d*x)**11/(2145*d*(a*sin(c + d*x) + a)**(sympy.S(11)/2)) - 16*a**2*cos(c + d*x)**11/(195*d*(a*sin(c + d*x) + a)**(sympy.S(9)/2)) - 2*a*cos(c + d*x)**11/(15*d*(a*sin(c + d*x) + a)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_182():
    f = cos(c + d*x)**9/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = 32*(a*sin(c + d*x) + a)**(sympy.S(5)/2)/(5*a**5*d) - 64*(a*sin(c + d*x) + a)**(sympy.S(7)/2)/(7*a**6*d) + 16*(a*sin(c + d*x) + a)**(sympy.S(9)/2)/(3*a**7*d) - 16*(a*sin(c + d*x) + a)**(sympy.S(11)/2)/(11*a**8*d) + 2*(a*sin(c + d*x) + a)**(sympy.S(13)/2)/(13*a**9*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_183():
    f = cos(c + d*x)**8/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = -8*a**2*cos(c + d*x)**9/(99*d*(a*sin(c + d*x) + a)**(sympy.S(9)/2)) - 2*a*cos(c + d*x)**9/(11*d*(a*sin(c + d*x) + a)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_184():
    f = cos(c + d*x)**7/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = 16*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(3*a**4*d) - 24*(a*sin(c + d*x) + a)**(sympy.S(5)/2)/(5*a**5*d) + 12*(a*sin(c + d*x) + a)**(sympy.S(7)/2)/(7*a**6*d) - 2*(a*sin(c + d*x) + a)**(sympy.S(9)/2)/(9*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_185():
    f = cos(c + d*x)**6/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = -2*a*cos(c + d*x)**7/(7*d*(a*sin(c + d*x) + a)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_186():
    f = cos(c + d*x)**5/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = 8*sqrt(a*sin(c + d*x) + a)/(a**3*d) - 8*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(3*a**4*d) + 2*(a*sin(c + d*x) + a)**(sympy.S(5)/2)/(5*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_187():
    f = cos(c + d*x)**4/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = 2*cos(c + d*x)**3/(3*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + 4*cos(c + d*x)/(a**2*d*sqrt(a*sin(c + d*x) + a)) - 4*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_188():
    f = cos(c + d*x)**3/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = -4/(a**2*d*sqrt(a*sin(c + d*x) + a)) - 2*sqrt(a*sin(c + d*x) + a)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_189():
    f = cos(c + d*x)**2/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = -cos(c + d*x)/(a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(2*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_190():
    f = cos(c + d*x)/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = -2/(3*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_191():
    f = sec(c + d*x)/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = -1/(5*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) - 1/(6*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 1/(4*a**2*d*sqrt(a*sin(c + d*x) + a)) + sqrt(2)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/(8*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_192():
    f = sec(c + d*x)**2/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = -sec(c + d*x)/(6*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) - 35*cos(c + d*x)/(128*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 7*sec(c + d*x)/(48*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + 35*sec(c + d*x)/(96*a**2*d*sqrt(a*sin(c + d*x) + a)) - 35*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(256*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_193():
    f = sec(c + d*x)**3/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = -sec(c + d*x)**2/(7*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) - 9*sec(c + d*x)**2/(70*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 3/(16*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + 9*sec(c + d*x)**2/(40*a**2*d*sqrt(a*sin(c + d*x) + a)) - 9/(32*a**2*d*sqrt(a*sin(c + d*x) + a)) + 9*sqrt(2)*atanh(sqrt(2)*sqrt(a*sin(c + d*x) + a)/(2*sqrt(a)))/(64*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_194():
    f = sec(c + d*x)**4/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = -sec(c + d*x)**3/(8*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) - 1155*cos(c + d*x)/(4096*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 11*sec(c + d*x)**3/(96*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 77*sec(c + d*x)/(512*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + 11*sec(c + d*x)**3/(64*a**2*d*sqrt(a*sin(c + d*x) + a)) + 385*sec(c + d*x)/(1024*a**2*d*sqrt(a*sin(c + d*x) + a)) - 1155*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(8192*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_195():
    f = (e*cos(c + d*x))**(sympy.S(7)/2)*(a*sin(c + d*x) + a)
    F = 10*a*e**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(e*cos(c + d*x))) + 10*a*e**3*sqrt(e*cos(c + d*x))*sin(c + d*x)/(21*d) + 2*a*e*(e*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*d) - 2*a*(e*cos(c + d*x))**(sympy.S(9)/2)/(9*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_196():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)*(a*sin(c + d*x) + a)
    F = 6*a*e**2*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) + 2*a*e*(e*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*d) - 2*a*(e*cos(c + d*x))**(sympy.S(7)/2)/(7*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_197():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a)
    F = 2*a*e**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(e*cos(c + d*x))) + 2*a*e*sqrt(e*cos(c + d*x))*sin(c + d*x)/(3*d) - 2*a*(e*cos(c + d*x))**(sympy.S(5)/2)/(5*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_198():
    f = sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)
    F = 2*a*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x))) - 2*a*(e*cos(c + d*x))**(sympy.S(3)/2)/(3*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_199():
    f = (a*sin(c + d*x) + a)/sqrt(e*cos(c + d*x))
    F = 2*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*sqrt(e*cos(c + d*x))) - 2*a*sqrt(e*cos(c + d*x))/(d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_200():
    f = (a*sin(c + d*x) + a)/(e*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*a*sin(c + d*x)/(d*e*sqrt(e*cos(c + d*x))) + 2*a/(d*e*sqrt(e*cos(c + d*x))) - 2*a*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*e**2*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_201():
    f = (a*sin(c + d*x) + a)/(e*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*a*sin(c + d*x)/(3*d*e*(e*cos(c + d*x))**(sympy.S(3)/2)) + 2*a/(3*d*e*(e*cos(c + d*x))**(sympy.S(3)/2)) + 2*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*e**2*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_202():
    f = (a*sin(c + d*x) + a)/(e*cos(c + d*x))**(sympy.S(7)/2)
    F = 2*a*sin(c + d*x)/(5*d*e*(e*cos(c + d*x))**(sympy.S(5)/2)) + 2*a/(5*d*e*(e*cos(c + d*x))**(sympy.S(5)/2)) + 6*a*sin(c + d*x)/(5*d*e**3*sqrt(e*cos(c + d*x))) - 6*a*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*e**4*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_203():
    f = (e*cos(c + d*x))**(sympy.S(7)/2)*(a*sin(c + d*x) + a)**2
    F = 130*a**2*e**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(231*d*sqrt(e*cos(c + d*x))) + 130*a**2*e**3*sqrt(e*cos(c + d*x))*sin(c + d*x)/(231*d) + 26*a**2*e*(e*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(77*d) - 26*a**2*(e*cos(c + d*x))**(sympy.S(9)/2)/(99*d*e) - 2*(e*cos(c + d*x))**(sympy.S(9)/2)*(a**2*sin(c + d*x) + a**2)/(11*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_204():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)*(a*sin(c + d*x) + a)**2
    F = 22*a**2*e**2*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(15*d*sqrt(cos(c + d*x))) + 22*a**2*e*(e*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(45*d) - 22*a**2*(e*cos(c + d*x))**(sympy.S(7)/2)/(63*d*e) - 2*(e*cos(c + d*x))**(sympy.S(7)/2)*(a**2*sin(c + d*x) + a**2)/(9*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_205():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a)**2
    F = 6*a**2*e**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(7*d*sqrt(e*cos(c + d*x))) + 6*a**2*e*sqrt(e*cos(c + d*x))*sin(c + d*x)/(7*d) - 18*a**2*(e*cos(c + d*x))**(sympy.S(5)/2)/(35*d*e) - 2*(e*cos(c + d*x))**(sympy.S(5)/2)*(a**2*sin(c + d*x) + a**2)/(7*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_206():
    f = sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**2
    F = 14*a**2*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) - 14*a**2*(e*cos(c + d*x))**(sympy.S(3)/2)/(15*d*e) - 2*(e*cos(c + d*x))**(sympy.S(3)/2)*(a**2*sin(c + d*x) + a**2)/(5*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_207():
    f = (a*sin(c + d*x) + a)**2/sqrt(e*cos(c + d*x))
    F = 10*a**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(e*cos(c + d*x))) - 10*a**2*sqrt(e*cos(c + d*x))/(3*d*e) - 2*sqrt(e*cos(c + d*x))*(a**2*sin(c + d*x) + a**2)/(3*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_208():
    f = (a*sin(c + d*x) + a)**2/(e*cos(c + d*x))**(sympy.S(3)/2)
    F = 4*a**4*(e*cos(c + d*x))**(sympy.S(3)/2)/(d*e**3*(-a**2*sin(c + d*x) + a**2)) - 6*a**2*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*e**2*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_209():
    f = (a*sin(c + d*x) + a)**2/(e*cos(c + d*x))**(sympy.S(5)/2)
    F = 4*a**4*sqrt(e*cos(c + d*x))/(3*d*e**3*(-a**2*sin(c + d*x) + a**2)) - 2*a**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*e**2*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_210():
    f = (a*sin(c + d*x) + a)**2/(e*cos(c + d*x))**(sympy.S(7)/2)
    F = 2*a**4*(e*cos(c + d*x))**(sympy.S(3)/2)/(5*d*e**5*(-a**2*sin(c + d*x) + a**2)) + 2*a**4*(e*cos(c + d*x))**(sympy.S(3)/2)/(5*d*e**5*(-a*sin(c + d*x) + a)**2) - 2*a**2*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*e**4*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_211():
    f = (a*sin(c + d*x) + a)**2/(e*cos(c + d*x))**(sympy.S(9)/2)
    F = 2*a**2*sin(c + d*x)/(7*d*e**3*(e*cos(c + d*x))**(sympy.S(3)/2)) + 2*a**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(7*d*e**4*sqrt(e*cos(c + d*x))) + (4*a**2*sin(c + d*x) + 4*a**2)/(7*d*e*(e*cos(c + d*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_212():
    f = (a*sin(c + d*x) + a)**2/(e*cos(c + d*x))**(sympy.S(11)/2)
    F = 2*a**2*sin(c + d*x)/(9*d*e**3*(e*cos(c + d*x))**(sympy.S(5)/2)) + 2*a**2*sin(c + d*x)/(3*d*e**5*sqrt(e*cos(c + d*x))) - 2*a**2*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(3*d*e**6*sqrt(cos(c + d*x))) + (4*a**2*sin(c + d*x) + 4*a**2)/(9*d*e*(e*cos(c + d*x))**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_213():
    f = (e*cos(c + d*x))**(sympy.S(7)/2)*(a*sin(c + d*x) + a)**3
    F = 170*a**3*e**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(231*d*sqrt(e*cos(c + d*x))) + 170*a**3*e**3*sqrt(e*cos(c + d*x))*sin(c + d*x)/(231*d) + 34*a**3*e*(e*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(77*d) - 34*a**3*(e*cos(c + d*x))**(sympy.S(9)/2)/(99*d*e) - 2*a*(e*cos(c + d*x))**(sympy.S(9)/2)*(a*sin(c + d*x) + a)**2/(13*d*e) - 34*(e*cos(c + d*x))**(sympy.S(9)/2)*(a**3*sin(c + d*x) + a**3)/(143*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_214():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)*(a*sin(c + d*x) + a)**3
    F = 2*a**3*e**2*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x))) + 2*a**3*e*(e*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(3*d) - 10*a**3*(e*cos(c + d*x))**(sympy.S(7)/2)/(21*d*e) - 2*a*(e*cos(c + d*x))**(sympy.S(7)/2)*(a*sin(c + d*x) + a)**2/(11*d*e) - 10*(e*cos(c + d*x))**(sympy.S(7)/2)*(a**3*sin(c + d*x) + a**3)/(33*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_215():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a)**3
    F = 26*a**3*e**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(e*cos(c + d*x))) + 26*a**3*e*sqrt(e*cos(c + d*x))*sin(c + d*x)/(21*d) - 26*a**3*(e*cos(c + d*x))**(sympy.S(5)/2)/(35*d*e) - 2*a*(e*cos(c + d*x))**(sympy.S(5)/2)*(a*sin(c + d*x) + a)**2/(9*d*e) - 26*(e*cos(c + d*x))**(sympy.S(5)/2)*(a**3*sin(c + d*x) + a**3)/(63*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_216():
    f = sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**3
    F = 22*a**3*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) - 22*a**3*(e*cos(c + d*x))**(sympy.S(3)/2)/(15*d*e) - 2*a*(e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a)**2/(7*d*e) - 22*(e*cos(c + d*x))**(sympy.S(3)/2)*(a**3*sin(c + d*x) + a**3)/(35*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_217():
    f = (a*sin(c + d*x) + a)**3/sqrt(e*cos(c + d*x))
    F = 6*a**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*sqrt(e*cos(c + d*x))) - 6*a**3*sqrt(e*cos(c + d*x))/(d*e) - 2*a*sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**2/(5*d*e) - 6*sqrt(e*cos(c + d*x))*(a**3*sin(c + d*x) + a**3)/(5*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_218():
    f = (a*sin(c + d*x) + a)**3/(e*cos(c + d*x))**(sympy.S(3)/2)
    F = 4*a**5*(e*cos(c + d*x))**(sympy.S(7)/2)/(d*e**5*(-a*sin(c + d*x) + a)**2) - 14*a**3*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*e**2*sqrt(cos(c + d*x))) + 14*a**3*(e*cos(c + d*x))**(sympy.S(3)/2)/(3*d*e**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_219():
    f = (a*sin(c + d*x) + a)**3/(e*cos(c + d*x))**(sympy.S(5)/2)
    F = 4*a**5*(e*cos(c + d*x))**(sympy.S(5)/2)/(3*d*e**5*(-a*sin(c + d*x) + a)**2) - 10*a**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*e**2*sqrt(e*cos(c + d*x))) + 10*a**3*sqrt(e*cos(c + d*x))/(3*d*e**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_220():
    f = (a*sin(c + d*x) + a)**3/(e*cos(c + d*x))**(sympy.S(7)/2)
    F = -6*a**6*(e*cos(c + d*x))**(sympy.S(3)/2)/(5*d*e**5*(-a**3*sin(c + d*x) + a**3)) + 4*a**5*(e*cos(c + d*x))**(sympy.S(3)/2)/(5*d*e**5*(-a*sin(c + d*x) + a)**2) + 6*a**3*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*e**4*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_221():
    f = (a*sin(c + d*x) + a)**3/(e*cos(c + d*x))**(sympy.S(9)/2)
    F = -2*a**6*sqrt(e*cos(c + d*x))/(21*d*e**5*(-a**3*sin(c + d*x) + a**3)) + 4*a**5*sqrt(e*cos(c + d*x))/(7*d*e**5*(-a*sin(c + d*x) + a)**2) - 2*a**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*e**4*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_222():
    f = (a*sin(c + d*x) + a)**3/(e*cos(c + d*x))**(sympy.S(11)/2)
    F = 2*a**6*(e*cos(c + d*x))**(sympy.S(3)/2)/(15*d*e**7*(-a**3*sin(c + d*x) + a**3)) + 2*a**6*(e*cos(c + d*x))**(sympy.S(3)/2)/(9*d*e**7*(-a*sin(c + d*x) + a)**3) + 2*a**5*(e*cos(c + d*x))**(sympy.S(3)/2)/(15*d*e**7*(-a*sin(c + d*x) + a)**2) - 2*a**3*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(15*d*e**6*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_223():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a)**4
    F = 442*a**4*e**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(231*d*sqrt(e*cos(c + d*x))) + 442*a**4*e*sqrt(e*cos(c + d*x))*sin(c + d*x)/(231*d) - 442*a**4*(e*cos(c + d*x))**(sympy.S(5)/2)/(385*d*e) - 2*a*(e*cos(c + d*x))**(sympy.S(5)/2)*(a*sin(c + d*x) + a)**3/(11*d*e) - 34*(e*cos(c + d*x))**(sympy.S(5)/2)*(a**2*sin(c + d*x) + a**2)**2/(99*d*e) - 442*(e*cos(c + d*x))**(sympy.S(5)/2)*(a**4*sin(c + d*x) + a**4)/(693*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_224():
    f = sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**4
    F = 22*a**4*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(3*d*sqrt(cos(c + d*x))) - 22*a**4*(e*cos(c + d*x))**(sympy.S(3)/2)/(9*d*e) - 2*a*(e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a)**3/(9*d*e) - 10*(e*cos(c + d*x))**(sympy.S(3)/2)*(a**2*sin(c + d*x) + a**2)**2/(21*d*e) - 22*(e*cos(c + d*x))**(sympy.S(3)/2)*(a**4*sin(c + d*x) + a**4)/(21*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_225():
    f = (a*sin(c + d*x) + a)**4/sqrt(e*cos(c + d*x))
    F = 78*a**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(7*d*sqrt(e*cos(c + d*x))) - 78*a**4*sqrt(e*cos(c + d*x))/(7*d*e) - 2*a*sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**3/(7*d*e) - 26*sqrt(e*cos(c + d*x))*(a**2*sin(c + d*x) + a**2)**2/(35*d*e) - 78*sqrt(e*cos(c + d*x))*(a**4*sin(c + d*x) + a**4)/(35*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_226():
    f = (a*sin(c + d*x) + a)**4/(e*cos(c + d*x))**(sympy.S(3)/2)
    F = 44*a**8*(e*cos(c + d*x))**(sympy.S(7)/2)/(3*d*e**5*(-a**4*sin(c + d*x) + a**4)) + 4*a**7*(e*cos(c + d*x))**(sympy.S(11)/2)/(d*e**7*(-a*sin(c + d*x) + a)**3) - 154*a**4*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*e**2*sqrt(cos(c + d*x))) - 154*a**4*(e*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(15*d*e**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_227():
    f = (a*sin(c + d*x) + a)**4/(e*cos(c + d*x))**(sympy.S(5)/2)
    F = 12*a**8*(e*cos(c + d*x))**(sympy.S(5)/2)/(d*e**5*(-a**4*sin(c + d*x) + a**4)) + 4*a**7*(e*cos(c + d*x))**(sympy.S(9)/2)/(3*d*e**7*(-a*sin(c + d*x) + a)**3) - 10*a**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*e**2*sqrt(e*cos(c + d*x))) - 10*a**4*sqrt(e*cos(c + d*x))*sin(c + d*x)/(d*e**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_228():
    f = (a*sin(c + d*x) + a)**4/(e*cos(c + d*x))**(sympy.S(7)/2)
    F = -28*a**8*(e*cos(c + d*x))**(sympy.S(3)/2)/(5*d*e**5*(-a**4*sin(c + d*x) + a**4)) + 4*a**7*(e*cos(c + d*x))**(sympy.S(7)/2)/(5*d*e**7*(-a*sin(c + d*x) + a)**3) + 42*a**4*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*e**4*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_229():
    f = (a*sin(c + d*x) + a)**4/(e*cos(c + d*x))**(sympy.S(9)/2)
    F = -20*a**8*sqrt(e*cos(c + d*x))/(21*d*e**5*(-a**4*sin(c + d*x) + a**4)) + 4*a**7*(e*cos(c + d*x))**(sympy.S(5)/2)/(7*d*e**7*(-a*sin(c + d*x) + a)**3) + 10*a**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*e**4*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_230():
    f = (a*sin(c + d*x) + a)**4/(e*cos(c + d*x))**(sympy.S(11)/2)
    F = -2*a**8*(e*cos(c + d*x))**(sympy.S(3)/2)/(15*d*e**7*(-a**4*sin(c + d*x) + a**4)) - 2*a**8*(e*cos(c + d*x))**(sympy.S(3)/2)/(15*d*e**7*(-a**2*sin(c + d*x) + a**2)**2) + 4*a**7*(e*cos(c + d*x))**(sympy.S(3)/2)/(9*d*e**7*(-a*sin(c + d*x) + a)**3) + 2*a**4*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(15*d*e**6*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_231():
    f = (a*sin(c + d*x) + a)**4/(e*cos(c + d*x))**(sympy.S(13)/2)
    F = -2*a**8*sqrt(e*cos(c + d*x))/(77*d*e**7*(-a**4*sin(c + d*x) + a**4)) - 2*a**8*sqrt(e*cos(c + d*x))/(77*d*e**7*(-a**2*sin(c + d*x) + a**2)**2) + 4*a**7*sqrt(e*cos(c + d*x))/(11*d*e**7*(-a*sin(c + d*x) + a)**3) - 2*a**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(77*d*e**6*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_232():
    f = (e*cos(c + d*x))**(sympy.S(11)/2)/(a*sin(c + d*x) + a)
    F = 10*e**6*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*a*d*sqrt(e*cos(c + d*x))) + 10*e**5*sqrt(e*cos(c + d*x))*sin(c + d*x)/(21*a*d) + 2*e**3*(e*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*a*d) + 2*e*(e*cos(c + d*x))**(sympy.S(9)/2)/(9*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_233():
    f = (e*cos(c + d*x))**(sympy.S(9)/2)/(a*sin(c + d*x) + a)
    F = 6*e**4*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*a*d*sqrt(cos(c + d*x))) + 2*e**3*(e*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*a*d) + 2*e*(e*cos(c + d*x))**(sympy.S(7)/2)/(7*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_234():
    f = (e*cos(c + d*x))**(sympy.S(7)/2)/(a*sin(c + d*x) + a)
    F = 2*e**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*a*d*sqrt(e*cos(c + d*x))) + 2*e**3*sqrt(e*cos(c + d*x))*sin(c + d*x)/(3*a*d) + 2*e*(e*cos(c + d*x))**(sympy.S(5)/2)/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_235():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)/(a*sin(c + d*x) + a)
    F = 2*e**2*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(a*d*sqrt(cos(c + d*x))) + 2*e*(e*cos(c + d*x))**(sympy.S(3)/2)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_236():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)/(a*sin(c + d*x) + a)
    F = 2*e**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(a*d*sqrt(e*cos(c + d*x))) + 2*e*sqrt(e*cos(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_237():
    f = sqrt(e*cos(c + d*x))/(a*sin(c + d*x) + a)
    F = -2*(e*cos(c + d*x))**(sympy.S(3)/2)/(d*e*(a*sin(c + d*x) + a)) - 2*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(a*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_238():
    f = 1/(sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a))
    F = -2*sqrt(e*cos(c + d*x))/(3*d*e*(a*sin(c + d*x) + a)) + 2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*a*d*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_239():
    f = 1/((e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a))
    F = -2/(5*d*e*sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)) + 6*sin(c + d*x)/(5*a*d*e*sqrt(e*cos(c + d*x))) - 6*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*a*d*e**2*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_240():
    f = 1/((e*cos(c + d*x))**(sympy.S(5)/2)*(a*sin(c + d*x) + a))
    F = -2/(7*d*e*(e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a)) + 10*sin(c + d*x)/(21*a*d*e*(e*cos(c + d*x))**(sympy.S(3)/2)) + 10*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*a*d*e**2*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_241():
    f = 1/((e*cos(c + d*x))**(sympy.S(7)/2)*(a*sin(c + d*x) + a))
    F = -2/(9*d*e*(e*cos(c + d*x))**(sympy.S(5)/2)*(a*sin(c + d*x) + a)) + 14*sin(c + d*x)/(45*a*d*e*(e*cos(c + d*x))**(sympy.S(5)/2)) + 14*sin(c + d*x)/(15*a*d*e**3*sqrt(e*cos(c + d*x))) - 14*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(15*a*d*e**4*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_242():
    f = (e*cos(c + d*x))**(sympy.S(11)/2)/(a*sin(c + d*x) + a)**2
    F = 4*e*(e*cos(c + d*x))**(sympy.S(9)/2)/(5*d*(a**2*sin(c + d*x) + a**2)) + 6*e**6*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(7*a**2*d*sqrt(e*cos(c + d*x))) + 6*e**5*sqrt(e*cos(c + d*x))*sin(c + d*x)/(7*a**2*d) + 18*e**3*(e*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(35*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_243():
    f = (e*cos(c + d*x))**(sympy.S(9)/2)/(a*sin(c + d*x) + a)**2
    F = 4*e*(e*cos(c + d*x))**(sympy.S(7)/2)/(3*d*(a**2*sin(c + d*x) + a**2)) + 14*e**4*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*a**2*d*sqrt(cos(c + d*x))) + 14*e**3*(e*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(15*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_244():
    f = (e*cos(c + d*x))**(sympy.S(7)/2)/(a*sin(c + d*x) + a)**2
    F = 4*e*(e*cos(c + d*x))**(sympy.S(5)/2)/(d*(a**2*sin(c + d*x) + a**2)) + 10*e**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d*sqrt(e*cos(c + d*x))) + 10*e**3*sqrt(e*cos(c + d*x))*sin(c + d*x)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_245():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)/(a*sin(c + d*x) + a)**2
    F = -4*e*(e*cos(c + d*x))**(sympy.S(3)/2)/(d*(a**2*sin(c + d*x) + a**2)) - 6*e**2*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(a**2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_246():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)/(a*sin(c + d*x) + a)**2
    F = -4*e*sqrt(e*cos(c + d*x))/(3*d*(a**2*sin(c + d*x) + a**2)) - 2*e**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_247():
    f = sqrt(e*cos(c + d*x))/(a*sin(c + d*x) + a)**2
    F = -2*(e*cos(c + d*x))**(sympy.S(3)/2)/(5*d*e*(a**2*sin(c + d*x) + a**2)) - 2*(e*cos(c + d*x))**(sympy.S(3)/2)/(5*d*e*(a*sin(c + d*x) + a)**2) - 2*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*a**2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_248():
    f = 1/(sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**2)
    F = -2*sqrt(e*cos(c + d*x))/(7*d*e*(a**2*sin(c + d*x) + a**2)) - 2*sqrt(e*cos(c + d*x))/(7*d*e*(a*sin(c + d*x) + a)**2) + 2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(7*a**2*d*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_249():
    f = 1/((e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a)**2)
    F = -2/(9*d*e*sqrt(e*cos(c + d*x))*(a**2*sin(c + d*x) + a**2)) - 2/(9*d*e*sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**2) + 2*sin(c + d*x)/(3*a**2*d*e*sqrt(e*cos(c + d*x))) - 2*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(3*a**2*d*e**2*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_250():
    f = 1/((e*cos(c + d*x))**(sympy.S(5)/2)*(a*sin(c + d*x) + a)**2)
    F = -2/(11*d*e*(e*cos(c + d*x))**(sympy.S(3)/2)*(a**2*sin(c + d*x) + a**2)) - 2/(11*d*e*(e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a)**2) + 10*sin(c + d*x)/(33*a**2*d*e*(e*cos(c + d*x))**(sympy.S(3)/2)) + 10*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(33*a**2*d*e**2*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_251():
    f = 1/((e*cos(c + d*x))**(sympy.S(7)/2)*(a*sin(c + d*x) + a)**2)
    F = -2/(13*d*e*(e*cos(c + d*x))**(sympy.S(5)/2)*(a**2*sin(c + d*x) + a**2)) - 2/(13*d*e*(e*cos(c + d*x))**(sympy.S(5)/2)*(a*sin(c + d*x) + a)**2) + 14*sin(c + d*x)/(65*a**2*d*e*(e*cos(c + d*x))**(sympy.S(5)/2)) + 42*sin(c + d*x)/(65*a**2*d*e**3*sqrt(e*cos(c + d*x))) - 42*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(65*a**2*d*e**4*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_252():
    f = (e*cos(c + d*x))**(sympy.S(15)/2)/(a*sin(c + d*x) + a)**3
    F = 4*e*(e*cos(c + d*x))**(sympy.S(13)/2)/(5*a*d*(a*sin(c + d*x) + a)**2) + 26*e**8*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*a**3*d*sqrt(e*cos(c + d*x))) + 26*e**7*sqrt(e*cos(c + d*x))*sin(c + d*x)/(21*a**3*d) + 26*e**5*(e*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(35*a**3*d) + 26*e**3*(e*cos(c + d*x))**(sympy.S(9)/2)/(45*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_253():
    f = (e*cos(c + d*x))**(sympy.S(13)/2)/(a*sin(c + d*x) + a)**3
    F = 4*e*(e*cos(c + d*x))**(sympy.S(11)/2)/(3*a*d*(a*sin(c + d*x) + a)**2) + 22*e**6*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*a**3*d*sqrt(cos(c + d*x))) + 22*e**5*(e*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(15*a**3*d) + 22*e**3*(e*cos(c + d*x))**(sympy.S(7)/2)/(21*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_254():
    f = (e*cos(c + d*x))**(sympy.S(11)/2)/(a*sin(c + d*x) + a)**3
    F = 4*e*(e*cos(c + d*x))**(sympy.S(9)/2)/(a*d*(a*sin(c + d*x) + a)**2) + 6*e**6*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(a**3*d*sqrt(e*cos(c + d*x))) + 6*e**5*sqrt(e*cos(c + d*x))*sin(c + d*x)/(a**3*d) + 18*e**3*(e*cos(c + d*x))**(sympy.S(5)/2)/(5*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_255():
    f = (e*cos(c + d*x))**(sympy.S(9)/2)/(a*sin(c + d*x) + a)**3
    F = -4*e*(e*cos(c + d*x))**(sympy.S(7)/2)/(a*d*(a*sin(c + d*x) + a)**2) - 14*e**4*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(a**3*d*sqrt(cos(c + d*x))) - 14*e**3*(e*cos(c + d*x))**(sympy.S(3)/2)/(3*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_256():
    f = (e*cos(c + d*x))**(sympy.S(7)/2)/(a*sin(c + d*x) + a)**3
    F = -4*e*(e*cos(c + d*x))**(sympy.S(5)/2)/(3*a*d*(a*sin(c + d*x) + a)**2) - 10*e**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*a**3*d*sqrt(e*cos(c + d*x))) - 10*e**3*sqrt(e*cos(c + d*x))/(3*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_257():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)/(a*sin(c + d*x) + a)**3
    F = 6*e*(e*cos(c + d*x))**(sympy.S(3)/2)/(5*d*(a**3*sin(c + d*x) + a**3)) - 4*e*(e*cos(c + d*x))**(sympy.S(3)/2)/(5*a*d*(a*sin(c + d*x) + a)**2) + 6*e**2*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*a**3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_258():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)/(a*sin(c + d*x) + a)**3
    F = 2*e*sqrt(e*cos(c + d*x))/(21*d*(a**3*sin(c + d*x) + a**3)) - 4*e*sqrt(e*cos(c + d*x))/(7*a*d*(a*sin(c + d*x) + a)**2) - 2*e**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*a**3*d*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_259():
    f = sqrt(e*cos(c + d*x))/(a*sin(c + d*x) + a)**3
    F = -2*(e*cos(c + d*x))**(sympy.S(3)/2)/(15*d*e*(a**3*sin(c + d*x) + a**3)) - 2*(e*cos(c + d*x))**(sympy.S(3)/2)/(9*d*e*(a*sin(c + d*x) + a)**3) - 2*(e*cos(c + d*x))**(sympy.S(3)/2)/(15*a*d*e*(a*sin(c + d*x) + a)**2) - 2*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(15*a**3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_260():
    f = 1/(sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**3)
    F = -10*sqrt(e*cos(c + d*x))/(77*d*e*(a**3*sin(c + d*x) + a**3)) - 2*sqrt(e*cos(c + d*x))/(11*d*e*(a*sin(c + d*x) + a)**3) - 10*sqrt(e*cos(c + d*x))/(77*a*d*e*(a*sin(c + d*x) + a)**2) + 10*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(77*a**3*d*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_261():
    f = 1/((e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a)**3)
    F = -14/(117*d*e*sqrt(e*cos(c + d*x))*(a**3*sin(c + d*x) + a**3)) - 2/(13*d*e*sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**3) - 14/(117*a*d*e*sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**2) + 14*sin(c + d*x)/(39*a**3*d*e*sqrt(e*cos(c + d*x))) - 14*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(39*a**3*d*e**2*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_262():
    f = (e*cos(c + d*x))**(sympy.S(15)/2)/(a*sin(c + d*x) + a)**4
    F = 52*e**3*(e*cos(c + d*x))**(sympy.S(9)/2)/(5*d*(a**4*sin(c + d*x) + a**4)) + 4*e*(e*cos(c + d*x))**(sympy.S(13)/2)/(a*d*(a*sin(c + d*x) + a)**3) + 78*e**8*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(7*a**4*d*sqrt(e*cos(c + d*x))) + 78*e**7*sqrt(e*cos(c + d*x))*sin(c + d*x)/(7*a**4*d) + 234*e**5*(e*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(35*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_263():
    f = (e*cos(c + d*x))**(sympy.S(13)/2)/(a*sin(c + d*x) + a)**4
    F = -44*e**3*(e*cos(c + d*x))**(sympy.S(7)/2)/(3*d*(a**4*sin(c + d*x) + a**4)) - 4*e*(e*cos(c + d*x))**(sympy.S(11)/2)/(a*d*(a*sin(c + d*x) + a)**3) - 154*e**6*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*a**4*d*sqrt(cos(c + d*x))) - 154*e**5*(e*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(15*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_264():
    f = (e*cos(c + d*x))**(sympy.S(11)/2)/(a*sin(c + d*x) + a)**4
    F = -12*e**3*(e*cos(c + d*x))**(sympy.S(5)/2)/(d*(a**4*sin(c + d*x) + a**4)) - 4*e*(e*cos(c + d*x))**(sympy.S(9)/2)/(3*a*d*(a*sin(c + d*x) + a)**3) - 10*e**6*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(a**4*d*sqrt(e*cos(c + d*x))) - 10*e**5*sqrt(e*cos(c + d*x))*sin(c + d*x)/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_265():
    f = (e*cos(c + d*x))**(sympy.S(9)/2)/(a*sin(c + d*x) + a)**4
    F = 28*e**3*(e*cos(c + d*x))**(sympy.S(3)/2)/(5*d*(a**4*sin(c + d*x) + a**4)) - 4*e*(e*cos(c + d*x))**(sympy.S(7)/2)/(5*a*d*(a*sin(c + d*x) + a)**3) + 42*e**4*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*a**4*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_266():
    f = (e*cos(c + d*x))**(sympy.S(7)/2)/(a*sin(c + d*x) + a)**4
    F = 20*e**3*sqrt(e*cos(c + d*x))/(21*d*(a**4*sin(c + d*x) + a**4)) - 4*e*(e*cos(c + d*x))**(sympy.S(5)/2)/(7*a*d*(a*sin(c + d*x) + a)**3) + 10*e**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*a**4*d*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_267():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)/(a*sin(c + d*x) + a)**4
    F = 2*e*(e*cos(c + d*x))**(sympy.S(3)/2)/(15*d*(a**4*sin(c + d*x) + a**4)) + 2*e*(e*cos(c + d*x))**(sympy.S(3)/2)/(15*d*(a**2*sin(c + d*x) + a**2)**2) - 4*e*(e*cos(c + d*x))**(sympy.S(3)/2)/(9*a*d*(a*sin(c + d*x) + a)**3) + 2*e**2*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(15*a**4*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_268():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)/(a*sin(c + d*x) + a)**4
    F = 2*e*sqrt(e*cos(c + d*x))/(77*d*(a**4*sin(c + d*x) + a**4)) + 2*e*sqrt(e*cos(c + d*x))/(77*d*(a**2*sin(c + d*x) + a**2)**2) - 4*e*sqrt(e*cos(c + d*x))/(11*a*d*(a*sin(c + d*x) + a)**3) - 2*e**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(77*a**4*d*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_269():
    f = sqrt(e*cos(c + d*x))/(a*sin(c + d*x) + a)**4
    F = -2*(e*cos(c + d*x))**(sympy.S(3)/2)/(39*d*e*(a**4*sin(c + d*x) + a**4)) - 2*(e*cos(c + d*x))**(sympy.S(3)/2)/(39*d*e*(a**2*sin(c + d*x) + a**2)**2) - 2*(e*cos(c + d*x))**(sympy.S(3)/2)/(13*d*e*(a*sin(c + d*x) + a)**4) - 10*(e*cos(c + d*x))**(sympy.S(3)/2)/(117*a*d*e*(a*sin(c + d*x) + a)**3) - 2*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(39*a**4*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_270():
    f = 1/(sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**4)
    F = -2*sqrt(e*cos(c + d*x))/(33*d*e*(a**4*sin(c + d*x) + a**4)) - 2*sqrt(e*cos(c + d*x))/(33*d*e*(a**2*sin(c + d*x) + a**2)**2) - 2*sqrt(e*cos(c + d*x))/(15*d*e*(a*sin(c + d*x) + a)**4) - 14*sqrt(e*cos(c + d*x))/(165*a*d*e*(a*sin(c + d*x) + a)**3) + 2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(33*a**4*d*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_271():
    f = 1/((e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a)**4)
    F = -14/(221*d*e*sqrt(e*cos(c + d*x))*(a**4*sin(c + d*x) + a**4)) - 14/(221*d*e*sqrt(e*cos(c + d*x))*(a**2*sin(c + d*x) + a**2)**2) - 2/(17*d*e*sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**4) - 18/(221*a*d*e*sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**3) + 42*sin(c + d*x)/(221*a**4*d*e*sqrt(e*cos(c + d*x))) - 42*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(221*a**4*d*e**2*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_272():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)*sqrt(a*sin(c + d*x) + a)
    F = -a*(e*cos(c + d*x))**(sympy.S(5)/2)/(2*d*e*sqrt(a*sin(c + d*x) + a)) - 3*e**(sympy.S(3)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(4*d*(sin(c + d*x) + cos(c + d*x) + 1)) + 3*e**(sympy.S(3)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(4*d*(sin(c + d*x) + cos(c + d*x) + 1)) + 3*e*sqrt(e*cos(c + d*x))*sqrt(a*sin(c + d*x) + a)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_273():
    f = sqrt(e*cos(c + d*x))*sqrt(a*sin(c + d*x) + a)
    F = -a*(e*cos(c + d*x))**(sympy.S(3)/2)/(d*e*sqrt(a*sin(c + d*x) + a)) + sqrt(e)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(d*(sin(c + d*x) + cos(c + d*x) + 1)) + sqrt(e)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(d*(sin(c + d*x) + cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_274():
    f = sqrt(a*sin(c + d*x) + a)/sqrt(e*cos(c + d*x))
    F = -2*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(d*sqrt(e)*(sin(c + d*x) + cos(c + d*x) + 1)) + 2*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(d*sqrt(e)*(sin(c + d*x) + cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_275():
    f = sqrt(a*sin(c + d*x) + a)/(e*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*sqrt(a*sin(c + d*x) + a)/(d*e*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_276():
    f = sqrt(a*sin(c + d*x) + a)/(e*cos(c + d*x))**(sympy.S(5)/2)
    F = -2*sqrt(a*sin(c + d*x) + a)/(d*e*(e*cos(c + d*x))**(sympy.S(3)/2)) + 4*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(3*a*d*e*(e*cos(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_277():
    f = sqrt(a*sin(c + d*x) + a)/(e*cos(c + d*x))**(sympy.S(7)/2)
    F = -2*sqrt(a*sin(c + d*x) + a)/(3*d*e*(e*cos(c + d*x))**(sympy.S(5)/2)) + 8*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(3*a*d*e*(e*cos(c + d*x))**(sympy.S(5)/2)) - 16*(a*sin(c + d*x) + a)**(sympy.S(5)/2)/(15*a**2*d*e*(e*cos(c + d*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_278():
    f = sqrt(a*sin(c + d*x) + a)/(e*cos(c + d*x))**(sympy.S(9)/2)
    F = -2*sqrt(a*sin(c + d*x) + a)/(5*d*e*(e*cos(c + d*x))**(sympy.S(7)/2)) - 12*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(5*a*d*e*(e*cos(c + d*x))**(sympy.S(7)/2)) + 16*(a*sin(c + d*x) + a)**(sympy.S(5)/2)/(5*a**2*d*e*(e*cos(c + d*x))**(sympy.S(7)/2)) - 32*(a*sin(c + d*x) + a)**(sympy.S(7)/2)/(35*a**3*d*e*(e*cos(c + d*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_279():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)*(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = -15*a**3*(e*cos(c + d*x))**(sympy.S(7)/2)/(32*d*e*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + 15*a**2*e*(e*cos(c + d*x))**(sympy.S(3)/2)/(64*d*sqrt(a*sin(c + d*x) + a)) - 3*a**2*(e*cos(c + d*x))**(sympy.S(7)/2)/(8*d*e*sqrt(a*sin(c + d*x) + a)) + 45*a*e**(sympy.S(5)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(64*d*(sin(c + d*x) + cos(c + d*x) + 1)) + 45*a*e**(sympy.S(5)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(64*d*(sin(c + d*x) + cos(c + d*x) + 1)) - a*(e*cos(c + d*x))**(sympy.S(7)/2)*sqrt(a*sin(c + d*x) + a)/(4*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_280():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = -7*a**2*(e*cos(c + d*x))**(sympy.S(5)/2)/(12*d*e*sqrt(a*sin(c + d*x) + a)) - 7*a*e**(sympy.S(3)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(8*d*(sin(c + d*x) + cos(c + d*x) + 1)) + 7*a*e**(sympy.S(3)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(8*d*(sin(c + d*x) + cos(c + d*x) + 1)) + 7*a*e*sqrt(e*cos(c + d*x))*sqrt(a*sin(c + d*x) + a)/(8*d) - a*(e*cos(c + d*x))**(sympy.S(5)/2)*sqrt(a*sin(c + d*x) + a)/(3*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_281():
    f = sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = -5*a**2*(e*cos(c + d*x))**(sympy.S(3)/2)/(4*d*e*sqrt(a*sin(c + d*x) + a)) + 5*a*sqrt(e)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(4*d*(sin(c + d*x) + cos(c + d*x) + 1)) + 5*a*sqrt(e)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(4*d*(sin(c + d*x) + cos(c + d*x) + 1)) - a*(e*cos(c + d*x))**(sympy.S(3)/2)*sqrt(a*sin(c + d*x) + a)/(2*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_282():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)/sqrt(e*cos(c + d*x))
    F = -a*sqrt(e*cos(c + d*x))*sqrt(a*sin(c + d*x) + a)/(d*e) - 3*a*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(d*sqrt(e)*(sin(c + d*x) + cos(c + d*x) + 1)) + 3*a*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(d*sqrt(e)*(sin(c + d*x) + cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_283():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)/(e*cos(c + d*x))**(sympy.S(3)/2)
    F = -2*a**2*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(d*e**(sympy.S(3)/2)*(a*sin(c + d*x) + a*cos(c + d*x) + a)) - 2*a**2*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(d*e**(sympy.S(3)/2)*(a*sin(c + d*x) + a*cos(c + d*x) + a)) + 4*a*sqrt(a*sin(c + d*x) + a)/(d*e*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_284():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)/(e*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(3*d*e*(e*cos(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_285():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)/(e*cos(c + d*x))**(sympy.S(7)/2)
    F = 2*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(d*e*(e*cos(c + d*x))**(sympy.S(5)/2)) - 4*(a*sin(c + d*x) + a)**(sympy.S(5)/2)/(5*a*d*e*(e*cos(c + d*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_286():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)/(e*cos(c + d*x))**(sympy.S(9)/2)
    F = -2*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(d*e*(e*cos(c + d*x))**(sympy.S(7)/2)) + 8*(a*sin(c + d*x) + a)**(sympy.S(5)/2)/(3*a*d*e*(e*cos(c + d*x))**(sympy.S(7)/2)) - 16*(a*sin(c + d*x) + a)**(sympy.S(7)/2)/(21*a**2*d*e*(e*cos(c + d*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_287():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)/(e*cos(c + d*x))**(sympy.S(11)/2)
    F = -2*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(3*d*e*(e*cos(c + d*x))**(sympy.S(9)/2)) + 4*(a*sin(c + d*x) + a)**(sympy.S(5)/2)/(a*d*e*(e*cos(c + d*x))**(sympy.S(9)/2)) - 16*(a*sin(c + d*x) + a)**(sympy.S(7)/2)/(5*a**2*d*e*(e*cos(c + d*x))**(sympy.S(9)/2)) + 32*(a*sin(c + d*x) + a)**(sympy.S(9)/2)/(45*a**3*d*e*(e*cos(c + d*x))**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_288():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = -77*a**3*(e*cos(c + d*x))**(sympy.S(5)/2)/(96*d*e*sqrt(a*sin(c + d*x) + a)) - 77*a**2*e**(sympy.S(3)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(64*d*(sin(c + d*x) + cos(c + d*x) + 1)) + 77*a**2*e**(sympy.S(3)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(64*d*(sin(c + d*x) + cos(c + d*x) + 1)) + 77*a**2*e*sqrt(e*cos(c + d*x))*sqrt(a*sin(c + d*x) + a)/(64*d) - 11*a**2*(e*cos(c + d*x))**(sympy.S(5)/2)*sqrt(a*sin(c + d*x) + a)/(24*d*e) - a*(e*cos(c + d*x))**(sympy.S(5)/2)*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(4*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_289():
    f = sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = -15*a**3*(e*cos(c + d*x))**(sympy.S(3)/2)/(8*d*e*sqrt(a*sin(c + d*x) + a)) + 15*a**2*sqrt(e)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(8*d*(sin(c + d*x) + cos(c + d*x) + 1)) + 15*a**2*sqrt(e)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(8*d*(sin(c + d*x) + cos(c + d*x) + 1)) - 3*a**2*(e*cos(c + d*x))**(sympy.S(3)/2)*sqrt(a*sin(c + d*x) + a)/(4*d*e) - a*(e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(3*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_290():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)/sqrt(e*cos(c + d*x))
    F = -7*a**2*sqrt(e*cos(c + d*x))*sqrt(a*sin(c + d*x) + a)/(4*d*e) - 21*a**2*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(4*d*sqrt(e)*(sin(c + d*x) + cos(c + d*x) + 1)) + 21*a**2*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(4*d*sqrt(e)*(sin(c + d*x) + cos(c + d*x) + 1)) - a*sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(2*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_291():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)/(e*cos(c + d*x))**(sympy.S(3)/2)
    F = 5*a**3*(e*cos(c + d*x))**(sympy.S(3)/2)/(d*e**3*sqrt(a*sin(c + d*x) + a)) - 5*a**2*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(d*e**(sympy.S(3)/2)*(sin(c + d*x) + cos(c + d*x) + 1)) - 5*a**2*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(d*e**(sympy.S(3)/2)*(sin(c + d*x) + cos(c + d*x) + 1)) + 4*a*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(d*e*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_292():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)/(e*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*a**2*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(d*e**(sympy.S(5)/2)*(sin(c + d*x) + cos(c + d*x) + 1)) - 2*a**2*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(d*e**(sympy.S(5)/2)*(sin(c + d*x) + cos(c + d*x) + 1)) + 4*a*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(3*d*e*(e*cos(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_293():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)/(e*cos(c + d*x))**(sympy.S(7)/2)
    F = 2*(a*sin(c + d*x) + a)**(sympy.S(5)/2)/(5*d*e*(e*cos(c + d*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_294():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)/(e*cos(c + d*x))**(sympy.S(9)/2)
    F = 2*(a*sin(c + d*x) + a)**(sympy.S(5)/2)/(3*d*e*(e*cos(c + d*x))**(sympy.S(7)/2)) - 4*(a*sin(c + d*x) + a)**(sympy.S(7)/2)/(21*a*d*e*(e*cos(c + d*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_295():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)/(e*cos(c + d*x))**(sympy.S(11)/2)
    F = 2*(a*sin(c + d*x) + a)**(sympy.S(5)/2)/(d*e*(e*cos(c + d*x))**(sympy.S(9)/2)) - 8*(a*sin(c + d*x) + a)**(sympy.S(7)/2)/(5*a*d*e*(e*cos(c + d*x))**(sympy.S(9)/2)) + 16*(a*sin(c + d*x) + a)**(sympy.S(9)/2)/(45*a**2*d*e*(e*cos(c + d*x))**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_296():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)/(e*cos(c + d*x))**(sympy.S(13)/2)
    F = -2*(a*sin(c + d*x) + a)**(sympy.S(5)/2)/(d*e*(e*cos(c + d*x))**(sympy.S(11)/2)) + 4*(a*sin(c + d*x) + a)**(sympy.S(7)/2)/(a*d*e*(e*cos(c + d*x))**(sympy.S(11)/2)) - 16*(a*sin(c + d*x) + a)**(sympy.S(9)/2)/(7*a**2*d*e*(e*cos(c + d*x))**(sympy.S(11)/2)) + 32*(a*sin(c + d*x) + a)**(sympy.S(11)/2)/(77*a**3*d*e*(e*cos(c + d*x))**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_297():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)/sqrt(a*sin(c + d*x) + a)
    F = -a*(e*cos(c + d*x))**(sympy.S(7)/2)/(2*d*e*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + 3*e**(sympy.S(5)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(4*d*(a*sin(c + d*x) + a*cos(c + d*x) + a)) + 3*e**(sympy.S(5)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(4*d*(a*sin(c + d*x) + a*cos(c + d*x) + a)) + e*(e*cos(c + d*x))**(sympy.S(3)/2)/(4*d*sqrt(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_298():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)/sqrt(a*sin(c + d*x) + a)
    F = -e**(sympy.S(3)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(a*d*(sin(c + d*x) + cos(c + d*x) + 1)) + e**(sympy.S(3)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(a*d*(sin(c + d*x) + cos(c + d*x) + 1)) + e*sqrt(e*cos(c + d*x))*sqrt(a*sin(c + d*x) + a)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_299():
    f = sqrt(e*cos(c + d*x))/sqrt(a*sin(c + d*x) + a)
    F = 2*sqrt(e)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(d*(a*sin(c + d*x) + a*cos(c + d*x) + a)) + 2*sqrt(e)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(d*(a*sin(c + d*x) + a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_300():
    f = 1/(sqrt(e*cos(c + d*x))*sqrt(a*sin(c + d*x) + a))
    F = -2*sqrt(e*cos(c + d*x))/(d*e*sqrt(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_301():
    f = 1/((e*cos(c + d*x))**(sympy.S(3)/2)*sqrt(a*sin(c + d*x) + a))
    F = -2/(3*d*e*sqrt(e*cos(c + d*x))*sqrt(a*sin(c + d*x) + a)) + 4*sqrt(a*sin(c + d*x) + a)/(3*a*d*e*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_302():
    f = 1/((e*cos(c + d*x))**(sympy.S(5)/2)*sqrt(a*sin(c + d*x) + a))
    F = -2/(5*d*e*(e*cos(c + d*x))**(sympy.S(3)/2)*sqrt(a*sin(c + d*x) + a)) - 8*sqrt(a*sin(c + d*x) + a)/(5*a*d*e*(e*cos(c + d*x))**(sympy.S(3)/2)) + 16*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(15*a**2*d*e*(e*cos(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_303():
    f = 1/((e*cos(c + d*x))**(sympy.S(7)/2)*sqrt(a*sin(c + d*x) + a))
    F = -2/(7*d*e*(e*cos(c + d*x))**(sympy.S(5)/2)*sqrt(a*sin(c + d*x) + a)) - 4*sqrt(a*sin(c + d*x) + a)/(7*a*d*e*(e*cos(c + d*x))**(sympy.S(5)/2)) + 16*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(7*a**2*d*e*(e*cos(c + d*x))**(sympy.S(5)/2)) - 32*(a*sin(c + d*x) + a)**(sympy.S(5)/2)/(35*a**3*d*e*(e*cos(c + d*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_304():
    f = (e*cos(c + d*x))**(sympy.S(7)/2)/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = e*(e*cos(c + d*x))**(sympy.S(5)/2)/(2*a*d*sqrt(a*sin(c + d*x) + a)) - 5*e**(sympy.S(7)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(4*a**2*d*(sin(c + d*x) + cos(c + d*x) + 1)) + 5*e**(sympy.S(7)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(4*a**2*d*(sin(c + d*x) + cos(c + d*x) + 1)) + 5*e**3*sqrt(e*cos(c + d*x))*sqrt(a*sin(c + d*x) + a)/(4*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_305():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = 3*e**(sympy.S(5)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(d*(a**2*sin(c + d*x) + a**2*cos(c + d*x) + a**2)) + 3*e**(sympy.S(5)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(d*(a**2*sin(c + d*x) + a**2*cos(c + d*x) + a**2)) + e*(e*cos(c + d*x))**(sympy.S(3)/2)/(a*d*sqrt(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_306():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = -2*(e*cos(c + d*x))**(sympy.S(5)/2)/(d*e*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + 2*e**(sympy.S(3)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(a**2*d*(sin(c + d*x) + cos(c + d*x) + 1)) - 2*e**(sympy.S(3)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(a**2*d*(sin(c + d*x) + cos(c + d*x) + 1)) - 2*e*sqrt(e*cos(c + d*x))*sqrt(a*sin(c + d*x) + a)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_307():
    f = sqrt(e*cos(c + d*x))/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = -2*(e*cos(c + d*x))**(sympy.S(3)/2)/(3*d*e*(a*sin(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_308():
    f = 1/(sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**(sympy.S(3)/2))
    F = -2*sqrt(e*cos(c + d*x))/(5*d*e*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 4*sqrt(e*cos(c + d*x))/(5*a*d*e*sqrt(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_309():
    f = 1/((e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a)**(sympy.S(3)/2))
    F = -2/(7*d*e*sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 8/(21*a*d*e*sqrt(e*cos(c + d*x))*sqrt(a*sin(c + d*x) + a)) + 16*sqrt(a*sin(c + d*x) + a)/(21*a**2*d*e*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_310():
    f = 1/((e*cos(c + d*x))**(sympy.S(5)/2)*(a*sin(c + d*x) + a)**(sympy.S(3)/2))
    F = -2/(9*d*e*(e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 4/(15*a*d*e*(e*cos(c + d*x))**(sympy.S(3)/2)*sqrt(a*sin(c + d*x) + a)) - 16*sqrt(a*sin(c + d*x) + a)/(15*a**2*d*e*(e*cos(c + d*x))**(sympy.S(3)/2)) + 32*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(45*a**3*d*e*(e*cos(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_311():
    f = 1/((e*cos(c + d*x))**(sympy.S(7)/2)*(a*sin(c + d*x) + a)**(sympy.S(3)/2))
    F = -2/(11*d*e*(e*cos(c + d*x))**(sympy.S(5)/2)*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 16/(77*a*d*e*(e*cos(c + d*x))**(sympy.S(5)/2)*sqrt(a*sin(c + d*x) + a)) - 32*sqrt(a*sin(c + d*x) + a)/(77*a**2*d*e*(e*cos(c + d*x))**(sympy.S(5)/2)) + 128*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(77*a**3*d*e*(e*cos(c + d*x))**(sympy.S(5)/2)) - 256*(a*sin(c + d*x) + a)**(sympy.S(5)/2)/(385*a**4*d*e*(e*cos(c + d*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_312():
    f = (e*cos(c + d*x))**(sympy.S(9)/2)/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = 21*e**(sympy.S(9)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(4*d*(a**3*sin(c + d*x) + a**3*cos(c + d*x) + a**3)) + 21*e**(sympy.S(9)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(4*d*(a**3*sin(c + d*x) + a**3*cos(c + d*x) + a**3)) + e*(e*cos(c + d*x))**(sympy.S(7)/2)/(2*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + 7*e**3*(e*cos(c + d*x))**(sympy.S(3)/2)/(4*a**2*d*sqrt(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_313():
    f = (e*cos(c + d*x))**(sympy.S(7)/2)/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = -4*e*(e*cos(c + d*x))**(sympy.S(5)/2)/(a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + 5*e**(sympy.S(7)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(a**3*d*(sin(c + d*x) + cos(c + d*x) + 1)) - 5*e**(sympy.S(7)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(a**3*d*(sin(c + d*x) + cos(c + d*x) + 1)) - 5*e**3*sqrt(e*cos(c + d*x))*sqrt(a*sin(c + d*x) + a)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_314():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = -2*e**(sympy.S(5)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*asinh(sqrt(e*cos(c + d*x))/sqrt(e))/(d*(a**3*sin(c + d*x) + a**3*cos(c + d*x) + a**3)) - 2*e**(sympy.S(5)/2)*sqrt(a*sin(c + d*x) + a)*sqrt(cos(c + d*x) + 1)*atan(sqrt(e)*sin(c + d*x)/(sqrt(e*cos(c + d*x))*sqrt(cos(c + d*x) + 1)))/(d*(a**3*sin(c + d*x) + a**3*cos(c + d*x) + a**3)) - 4*e*(e*cos(c + d*x))**(sympy.S(3)/2)/(3*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_315():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = -2*(e*cos(c + d*x))**(sympy.S(5)/2)/(5*d*e*(a*sin(c + d*x) + a)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_316():
    f = sqrt(e*cos(c + d*x))/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = -2*(e*cos(c + d*x))**(sympy.S(3)/2)/(7*d*e*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) - 4*(e*cos(c + d*x))**(sympy.S(3)/2)/(21*a*d*e*(a*sin(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_317():
    f = 1/(sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**(sympy.S(5)/2))
    F = -2*sqrt(e*cos(c + d*x))/(9*d*e*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) - 8*sqrt(e*cos(c + d*x))/(45*a*d*e*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 16*sqrt(e*cos(c + d*x))/(45*a**2*d*e*sqrt(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_318():
    f = 1/((e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a)**(sympy.S(5)/2))
    F = -2/(11*d*e*sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) - 12/(77*a*d*e*sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 16/(77*a**2*d*e*sqrt(e*cos(c + d*x))*sqrt(a*sin(c + d*x) + a)) + 32*sqrt(a*sin(c + d*x) + a)/(77*a**3*d*e*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_319():
    f = 1/((e*cos(c + d*x))**(sympy.S(5)/2)*(a*sin(c + d*x) + a)**(sympy.S(5)/2))
    F = -2/(13*d*e*(e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) - 16/(117*a*d*e*(e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 32/(195*a**2*d*e*(e*cos(c + d*x))**(sympy.S(3)/2)*sqrt(a*sin(c + d*x) + a)) - 128*sqrt(a*sin(c + d*x) + a)/(195*a**3*d*e*(e*cos(c + d*x))**(sympy.S(3)/2)) + 256*(a*sin(c + d*x) + a)**(sympy.S(3)/2)/(585*a**4*d*e*(e*cos(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_320():
    f = (e*cos(c + d*x))**(sympy.S(7)/3)/sqrt(a*sin(c + d*x) + a)
    F = -3*2**(sympy.S(1)/6)*a*(e*cos(c + d*x))**(sympy.S(10)/3)*hyper((sympy.S(-1)/6, sympy.S(5)/3), (sympy.S(8)/3,), sympy.S.Half - sin(c + d*x)/2)/(5*d*e*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*(sin(c + d*x) + 1)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_321():
    f = (e*cos(c + d*x))**(sympy.S(5)/3)/sqrt(a*sin(c + d*x) + a)
    F = -3*2**(sympy.S(5)/6)*a*(e*cos(c + d*x))**(sympy.S(8)/3)*(sin(c + d*x) + 1)**(sympy.S(1)/6)*hyper((sympy.S(1)/6, sympy.S(4)/3), (sympy.S(7)/3,), sympy.S.Half - sin(c + d*x)/2)/(8*d*e*(a*sin(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_322():
    f = (e*cos(c + d*x))**(sympy.S(2)/3)/sqrt(a*sin(c + d*x) + a)
    F = -3*2**(sympy.S(1)/3)*a*(e*cos(c + d*x))**(sympy.S(5)/3)*(sin(c + d*x) + 1)**(sympy.S(2)/3)*hyper((sympy.S(2)/3, sympy.S(5)/6), (sympy.S(11)/6,), sympy.S.Half - sin(c + d*x)/2)/(5*d*e*(a*sin(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_323():
    f = (e*cos(c + d*x))**(sympy.S(1)/3)/sqrt(a*sin(c + d*x) + a)
    F = -3*2**(sympy.S(1)/6)*a*(e*cos(c + d*x))**(sympy.S(4)/3)*(sin(c + d*x) + 1)**(sympy.S(5)/6)*hyper((sympy.S(2)/3, sympy.S(5)/6), (sympy.S(5)/3,), sympy.S.Half - sin(c + d*x)/2)/(4*d*e*(a*sin(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_324():
    f = 1/((e*cos(c + d*x))**(sympy.S(1)/3)*sqrt(a*sin(c + d*x) + a))
    F = -3*2**(sympy.S(5)/6)*(e*cos(c + d*x))**(sympy.S(2)/3)*(sin(c + d*x) + 1)**(sympy.S(1)/6)*hyper((sympy.S(1)/3, sympy.S(7)/6), (sympy.S(4)/3,), sympy.S.Half - sin(c + d*x)/2)/(4*d*e*sqrt(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_325():
    f = 1/((e*cos(c + d*x))**(sympy.S(4)/3)*sqrt(a*sin(c + d*x) + a))
    F = 3*2**(sympy.S(1)/3)*(sin(c + d*x) + 1)**(sympy.S(2)/3)*hyper((sympy.S(-1)/6, sympy.S(5)/3), (sympy.S(5)/6,), sympy.S.Half - sin(c + d*x)/2)/(2*d*e*(e*cos(c + d*x))**(sympy.S(1)/3)*sqrt(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_326():
    f = (e*cos(c + d*x))**p*(a*sin(c + d*x) + a)**8
    F = -2**(p/2 + sympy.S(17)/2)*a**8*(e*cos(c + d*x))**(p + 1)*(sin(c + d*x) + 1)**(-p/2 + sympy.S(-1)/2)*hyper((-p/2 + sympy.S(-15)/2, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*e*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_327():
    f = (e*cos(c + d*x))**p*(a*sin(c + d*x) + a)**3
    F = -2**(p/2 + sympy.S(7)/2)*a**3*(e*cos(c + d*x))**(p + 1)*(sin(c + d*x) + 1)**(-p/2 + sympy.S(-1)/2)*hyper((-p/2 + sympy.S(-5)/2, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*e*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_328():
    f = (e*cos(c + d*x))**p*(a*sin(c + d*x) + a)**2
    F = -2**(p/2 + sympy.S(5)/2)*a**2*(e*cos(c + d*x))**(p + 1)*(sin(c + d*x) + 1)**(-p/2 + sympy.S(-1)/2)*hyper((-p/2 + sympy.S(-3)/2, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*e*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_329():
    f = (e*cos(c + d*x))**p*(a*sin(c + d*x) + a)
    F = -2**(p/2 + sympy.S(3)/2)*a*(e*cos(c + d*x))**(p + 1)*(sin(c + d*x) + 1)**(-p/2 + sympy.S(-1)/2)*hyper((-p/2 + sympy.S(-1)/2, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*e*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_330():
    f = (e*cos(c + d*x))**p/(a*sin(c + d*x) + a)
    F = -2**(p/2 + sympy.S(-1)/2)*(e*cos(c + d*x))**(p + 1)*(sin(c + d*x) + 1)**(-p/2 + sympy.S(-1)/2)*hyper((sympy.S(3)/2 - p/2, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(a*d*e*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_331():
    f = (e*cos(c + d*x))**p/(a*sin(c + d*x) + a)**2
    F = -2**(p/2 + sympy.S(-3)/2)*(e*cos(c + d*x))**(p + 1)*(sin(c + d*x) + 1)**(-p/2 + sympy.S(-1)/2)*hyper((sympy.S(5)/2 - p/2, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(a**2*d*e*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_332():
    f = (e*cos(c + d*x))**p/(a*sin(c + d*x) + a)**3
    F = -2**(p/2 + sympy.S(-5)/2)*(e*cos(c + d*x))**(p + 1)*(sin(c + d*x) + 1)**(-p/2 + sympy.S(-1)/2)*hyper((sympy.S(7)/2 - p/2, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(a**3*d*e*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_333():
    f = (e*cos(c + d*x))**p/(a*sin(c + d*x) + a)**8
    F = -2**(p/2 + sympy.S(-15)/2)*(e*cos(c + d*x))**(p + 1)*(sin(c + d*x) + 1)**(-p/2 + sympy.S(-1)/2)*hyper((sympy.S(17)/2 - p/2, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(a**8*d*e*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_334():
    f = (e*cos(c + d*x))**p*(a*sin(c + d*x) + a)**(sympy.S(7)/2)
    F = -2**(p/2 + 4)*a**4*(e*cos(c + d*x))**(p + 1)*hyper((-p/2 - 3, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*e*(p + 1)*sqrt(a*sin(c + d*x) + a)*(sin(c + d*x) + 1)**(p/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_335():
    f = (e*cos(c + d*x))**p*(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = -2**(p/2 + 3)*a**3*(e*cos(c + d*x))**(p + 1)*hyper((-p/2 - 2, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*e*(p + 1)*sqrt(a*sin(c + d*x) + a)*(sin(c + d*x) + 1)**(p/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_336():
    f = (e*cos(c + d*x))**p*(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = -2**(p/2 + 2)*a**2*(e*cos(c + d*x))**(p + 1)*hyper((-p/2 - 1, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*e*(p + 1)*sqrt(a*sin(c + d*x) + a)*(sin(c + d*x) + 1)**(p/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_337():
    f = (e*cos(c + d*x))**p*sqrt(a*sin(c + d*x) + a)
    F = -2**(p/2 + 1)*a*(e*cos(c + d*x))**(p + 1)*hyper((-p/2, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*e*(p + 1)*sqrt(a*sin(c + d*x) + a)*(sin(c + d*x) + 1)**(p/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_338():
    f = (e*cos(c + d*x))**p/sqrt(a*sin(c + d*x) + a)
    F = -2**(p/2)*a*(e*cos(c + d*x))**(p + 1)*(sin(c + d*x) + 1)**(1 - p/2)*hyper((1 - p/2, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*e*(p + 1)*(a*sin(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_339():
    f = (e*cos(c + d*x))**p/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = -2**(p/2 - 1)*(e*cos(c + d*x))**(p + 1)*(sin(c + d*x) + 1)**(1 - p/2)*hyper((2 - p/2, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*e*(p + 1)*(a*sin(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_340():
    f = (e*cos(c + d*x))**p/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = -2**(p/2 - 2)*(e*cos(c + d*x))**(p + 1)*(sin(c + d*x) + 1)**(1 - p/2)*hyper((3 - p/2, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(a*d*e*(p + 1)*(a*sin(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_341():
    f = (e*cos(c + d*x))**p*(a*sin(c + d*x) + a)**m
    F = -2**(m + p/2 + sympy.S.Half)*a*(e*cos(c + d*x))**(p + 1)*(a*sin(c + d*x) + a)**(m - 1)*(sin(c + d*x) + 1)**(-m - p/2 + sympy.S.Half)*hyper((p/2 + sympy.S.Half, -m - p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*e*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_342():
    f = (a*sin(c + d*x) + a)**m*cos(c + d*x)**7
    F = 8*(a*sin(c + d*x) + a)**(m + 4)/(a**4*d*(m + 4)) - 12*(a*sin(c + d*x) + a)**(m + 5)/(a**5*d*(m + 5)) + 6*(a*sin(c + d*x) + a)**(m + 6)/(a**6*d*(m + 6)) - (a*sin(c + d*x) + a)**(m + 7)/(a**7*d*(m + 7))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_343():
    f = (a*sin(c + d*x) + a)**m*cos(c + d*x)**5
    F = 4*(a*sin(c + d*x) + a)**(m + 3)/(a**3*d*(m + 3)) - 4*(a*sin(c + d*x) + a)**(m + 4)/(a**4*d*(m + 4)) + (a*sin(c + d*x) + a)**(m + 5)/(a**5*d*(m + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_344():
    f = (a*sin(c + d*x) + a)**m*cos(c + d*x)**3
    F = 2*(a*sin(c + d*x) + a)**(m + 2)/(a**2*d*(m + 2)) - (a*sin(c + d*x) + a)**(m + 3)/(a**3*d*(m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_345():
    f = (a*sin(c + d*x) + a)**m*cos(c + d*x)
    F = (a*sin(c + d*x) + a)**(m + 1)/(a*d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_346():
    f = (a*sin(c + d*x) + a)**m*sec(c + d*x)
    F = (a*sin(c + d*x) + a)**m*hyper((1, m), (m + 1,), sin(c + d*x)/2 + sympy.S.Half)/(2*d*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_347():
    f = (a*sin(c + d*x) + a)**m*sec(c + d*x)**3
    F = -a*(a*sin(c + d*x) + a)**(m - 1)*hyper((2, m - 1), (m,), sin(c + d*x)/2 + sympy.S.Half)/(4*d*(1 - m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_348():
    f = (a*sin(c + d*x) + a)**m*sec(c + d*x)**5
    F = -a**2*(a*sin(c + d*x) + a)**(m - 2)*hyper((3, m - 2), (m - 1,), sin(c + d*x)/2 + sympy.S.Half)/(8*d*(2 - m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_349():
    f = (a*sin(c + d*x) + a)**m*cos(c + d*x)**4
    F = -2**(m + sympy.S(5)/2)*a**2*(a*sin(c + d*x) + a)**(m - 2)*(sin(c + d*x) + 1)**(-m + sympy.S(-1)/2)*cos(c + d*x)**5*hyper((sympy.S(5)/2, -m + sympy.S(-3)/2), (sympy.S(7)/2,), sympy.S.Half - sin(c + d*x)/2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_350():
    f = (a*sin(c + d*x) + a)**m*cos(c + d*x)**2
    F = -2**(m + sympy.S(3)/2)*a*(a*sin(c + d*x) + a)**(m - 1)*(sin(c + d*x) + 1)**(-m + sympy.S(-1)/2)*cos(c + d*x)**3*hyper((sympy.S(3)/2, -m + sympy.S(-1)/2), (sympy.S(5)/2,), sympy.S.Half - sin(c + d*x)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_351():
    f = (a*sin(c + d*x) + a)**m*sec(c + d*x)**2
    F = 2**(m + sympy.S(-1)/2)*(a*sin(c + d*x) + a)**m*(sin(c + d*x) + 1)**(sympy.S.Half - m)*hyper((sympy.S(-1)/2, sympy.S(3)/2 - m), (sympy.S.Half,), sympy.S.Half - sin(c + d*x)/2)*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_352():
    f = (a*sin(c + d*x) + a)**m*sec(c + d*x)**4
    F = 2**(m + sympy.S(-3)/2)*(a*sin(c + d*x) + a)**(m + 1)*(sin(c + d*x) + 1)**(sympy.S.Half - m)*hyper((sympy.S(-3)/2, sympy.S(5)/2 - m), (sympy.S(-1)/2,), sympy.S.Half - sin(c + d*x)/2)*sec(c + d*x)**3/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_353():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)*(a*sin(c + d*x) + a)**m
    F = -2**(m + sympy.S(11)/4)*a*(e*cos(c + d*x))**(sympy.S(7)/2)*(a*sin(c + d*x) + a)**(m - 1)*(sin(c + d*x) + 1)**(-m + sympy.S(-3)/4)*hyper((sympy.S(7)/4, -m + sympy.S(-3)/4), (sympy.S(11)/4,), sympy.S.Half - sin(c + d*x)/2)/(7*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_354():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a)**m
    F = -2**(m + sympy.S(9)/4)*a*(e*cos(c + d*x))**(sympy.S(5)/2)*(a*sin(c + d*x) + a)**(m - 1)*(sin(c + d*x) + 1)**(-m + sympy.S(-1)/4)*hyper((sympy.S(5)/4, -m + sympy.S(-1)/4), (sympy.S(9)/4,), sympy.S.Half - sin(c + d*x)/2)/(5*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_355():
    f = sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**m
    F = -2**(m + sympy.S(7)/4)*a*(e*cos(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + a)**(m - 1)*(sin(c + d*x) + 1)**(sympy.S(1)/4 - m)*hyper((sympy.S(3)/4, sympy.S(1)/4 - m), (sympy.S(7)/4,), sympy.S.Half - sin(c + d*x)/2)/(3*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_356():
    f = (a*sin(c + d*x) + a)**m/sqrt(e*cos(c + d*x))
    F = -2**(m + sympy.S(5)/4)*a*sqrt(e*cos(c + d*x))*(a*sin(c + d*x) + a)**(m - 1)*(sin(c + d*x) + 1)**(sympy.S(3)/4 - m)*hyper((sympy.S(1)/4, sympy.S(3)/4 - m), (sympy.S(5)/4,), sympy.S.Half - sin(c + d*x)/2)/(d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_357():
    f = (a*sin(c + d*x) + a)**m/(e*cos(c + d*x))**(sympy.S(3)/2)
    F = 2**(m + sympy.S(3)/4)*(a*sin(c + d*x) + a)**m*(sin(c + d*x) + 1)**(sympy.S(1)/4 - m)*hyper((sympy.S(-1)/4, sympy.S(5)/4 - m), (sympy.S(3)/4,), sympy.S.Half - sin(c + d*x)/2)/(d*e*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_358():
    f = (a*sin(c + d*x) + a)**m/(e*cos(c + d*x))**(sympy.S(5)/2)
    F = 2**(m + sympy.S(1)/4)*(a*sin(c + d*x) + a)**m*(sin(c + d*x) + 1)**(sympy.S(3)/4 - m)*hyper((sympy.S(-3)/4, sympy.S(7)/4 - m), (sympy.S(1)/4,), sympy.S.Half - sin(c + d*x)/2)/(3*d*e*(e*cos(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_359():
    f = (e*cos(c + d*x))**(-m - 4)*(a*sin(c + d*x) + a)**m
    F = -(e*cos(c + d*x))**(-m - 3)*(a*sin(c + d*x) + a)**m/(d*e*(3 - m)) - 3*(e*cos(c + d*x))**(-m - 3)*(a*sin(c + d*x) + a)**(m + 1)/(a*d*e*(1 - m)*(3 - m)) + 6*(e*cos(c + d*x))**(-m - 3)*(a*sin(c + d*x) + a)**(m + 2)/(a**2*d*e*(1 - m**2)*(3 - m)) - 6*(e*cos(c + d*x))**(-m - 3)*(a*sin(c + d*x) + a)**(m + 3)/(a**3*d*e*(m**4 - 10*m**2 + 9))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_360():
    f = (e*cos(c + d*x))**(-m - 3)*(a*sin(c + d*x) + a)**m
    F = -(e*cos(c + d*x))**(-m - 2)*(a*sin(c + d*x) + a)**m/(d*e*(2 - m)) + 2*(e*cos(c + d*x))**(-m - 2)*(a*sin(c + d*x) + a)**(m + 1)/(a*d*e*m*(2 - m)) - 2*(e*cos(c + d*x))**(-m - 2)*(a*sin(c + d*x) + a)**(m + 2)/(a**2*d*e*m*(4 - m**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_361():
    f = (e*cos(c + d*x))**(-m - 2)*(a*sin(c + d*x) + a)**m
    F = -(e*cos(c + d*x))**(-m - 1)*(a*sin(c + d*x) + a)**m/(d*e*(1 - m)) + (e*cos(c + d*x))**(-m - 1)*(a*sin(c + d*x) + a)**(m + 1)/(a*d*e*(1 - m**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_362():
    f = (e*cos(c + d*x))**(-m - 1)*(a*sin(c + d*x) + a)**m
    F = (a*sin(c + d*x) + a)**m/(d*e*m*(e*cos(c + d*x))**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_363():
    f = (a*sin(c + d*x) + a)**m/(e*cos(c + d*x))**m
    F = -2**(m/2 + sympy.S.Half)*a*(e*cos(c + d*x))**(1 - m)*(a*sin(c + d*x) + a)**(m - 1)*(sin(c + d*x) + 1)**(sympy.S.Half - m/2)*hyper((sympy.S.Half - m/2, sympy.S.Half - m/2), (sympy.S(3)/2 - m/2,), sympy.S.Half - sin(c + d*x)/2)/(d*e*(1 - m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_364():
    f = (e*cos(c + d*x))**(1 - m)*(a*sin(c + d*x) + a)**m
    F = 2**(1 - m/2)*(e*cos(c + d*x))**(2 - m)*(1 - sin(c + d*x))**(m/2 - 1)*(a*sin(c + d*x) + a)**m*hyper((m/2, m/2 + 1), (m/2 + 2,), sin(c + d*x)/2 + sympy.S.Half)/(d*e*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_365():
    f = (e*cos(c + d*x))**(2 - m)*(a*sin(c + d*x) + a)**m
    F = -2**(m/2 + sympy.S(3)/2)*a*(e*cos(c + d*x))**(3 - m)*(a*sin(c + d*x) + a)**(m - 1)*(sin(c + d*x) + 1)**(-m/2 + sympy.S(-1)/2)*hyper((sympy.S(3)/2 - m/2, -m/2 + sympy.S(-1)/2), (sympy.S(5)/2 - m/2,), sympy.S.Half - sin(c + d*x)/2)/(d*e*(3 - m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_366():
    f = (e*cos(c + d*x))**(5 - 2*m)*(a*sin(c + d*x) + a)**m
    F = -8*a**3*(e*cos(c + d*x))**(6 - 2*m)*(a*sin(c + d*x) + a)**(m - 3)/(d*e*(5 - m)*(m**2 - 7*m + 12)) - 4*a**2*(e*cos(c + d*x))**(6 - 2*m)*(a*sin(c + d*x) + a)**(m - 2)/(d*e*(m**2 - 9*m + 20)) - a*(e*cos(c + d*x))**(6 - 2*m)*(a*sin(c + d*x) + a)**(m - 1)/(d*e*(5 - m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_367():
    f = (e*cos(c + d*x))**(3 - 2*m)*(a*sin(c + d*x) + a)**m
    F = -2*a**2*(e*cos(c + d*x))**(4 - 2*m)*(a*sin(c + d*x) + a)**(m - 2)/(d*e*(m**2 - 5*m + 6)) - a*(e*cos(c + d*x))**(4 - 2*m)*(a*sin(c + d*x) + a)**(m - 1)/(d*e*(3 - m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_368():
    f = (e*cos(c + d*x))**(1 - 2*m)*(a*sin(c + d*x) + a)**m
    F = -a*(e*cos(c + d*x))**(2 - 2*m)*(a*sin(c + d*x) + a)**(m - 1)/(d*e*(1 - m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_369():
    f = (e*cos(c + d*x))**(-2*m - 1)*(a*sin(c + d*x) + a)**m
    F = (a*sin(c + d*x) + a)**m*hyper((1, -m), (1 - m,), sympy.S.Half - sin(c + d*x)/2)/(2*d*e*m*(e*cos(c + d*x))**(2*m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_370():
    f = (e*cos(c + d*x))**(-2*m - 3)*(a*sin(c + d*x) + a)**m
    F = (e*cos(c + d*x))**(-2*m - 2)*(a*sin(c + d*x) + a)**(m + 1)*hyper((2, -m - 1), (-m,), sympy.S.Half - sin(c + d*x)/2)/(4*a*d*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_371():
    f = (e*cos(c + d*x))**(4 - 2*m)*(a*sin(c + d*x) + a)**m
    F = 2**(sympy.S(5)/2 - m)*(e*cos(c + d*x))**(5 - 2*m)*(1 - sin(c + d*x))**(m + sympy.S(-5)/2)*(a*sin(c + d*x) + a)**m*hyper((sympy.S(5)/2, m + sympy.S(-3)/2), (sympy.S(7)/2,), sin(c + d*x)/2 + sympy.S.Half)/(5*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_372():
    f = (e*cos(c + d*x))**(2 - 2*m)*(a*sin(c + d*x) + a)**m
    F = 2**(sympy.S(3)/2 - m)*(e*cos(c + d*x))**(3 - 2*m)*(1 - sin(c + d*x))**(m + sympy.S(-3)/2)*(a*sin(c + d*x) + a)**m*hyper((sympy.S(3)/2, m + sympy.S(-1)/2), (sympy.S(5)/2,), sin(c + d*x)/2 + sympy.S.Half)/(3*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_373():
    f = (a*sin(c + d*x) + a)**m/(e*cos(c + d*x))**(2*m)
    F = 2**(sympy.S.Half - m)*(e*cos(c + d*x))**(1 - 2*m)*(1 - sin(c + d*x))**(m + sympy.S(-1)/2)*(a*sin(c + d*x) + a)**m*hyper((sympy.S.Half, m + sympy.S.Half), (sympy.S(3)/2,), sin(c + d*x)/2 + sympy.S.Half)/(d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_374():
    f = (e*cos(c + d*x))**(-2*m - 2)*(a*sin(c + d*x) + a)**m
    F = -2**(-m + sympy.S(-1)/2)*(e*cos(c + d*x))**(-2*m - 1)*(1 - sin(c + d*x))**(m + sympy.S.Half)*(a*sin(c + d*x) + a)**m*hyper((sympy.S(-1)/2, m + sympy.S(3)/2), (sympy.S.Half,), sin(c + d*x)/2 + sympy.S.Half)/(d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_375():
    f = (a + b*sin(c + d*x))*cos(c + d*x)**5
    F = a*sin(c + d*x)**5/(5*d) - 2*a*sin(c + d*x)**3/(3*d) + a*sin(c + d*x)/d - b*cos(c + d*x)**6/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_376():
    f = (a + b*sin(c + d*x))*cos(c + d*x)**3
    F = -a*sin(c + d*x)**3/(3*d) + a*sin(c + d*x)/d - b*cos(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_377():
    f = (a + b*sin(c + d*x))*sec(c + d*x)
    F = (a - b)*log(sin(c + d*x) + 1)/(2*d) - (a + b)*log(1 - sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_378():
    f = (a + b*sin(c + d*x))*sec(c + d*x)**3
    F = a*atanh(sin(c + d*x))/(2*d) + (a*sin(c + d*x) + b)*sec(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_379():
    f = (a + b*sin(c + d*x))*sec(c + d*x)**5
    F = 3*a*tan(c + d*x)*sec(c + d*x)/(8*d) + 3*a*atanh(sin(c + d*x))/(8*d) + (a*sin(c + d*x) + b)*sec(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_380():
    f = (a + b*sin(c + d*x))*cos(c + d*x)**4
    F = 3*a*x/8 + a*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 3*a*sin(c + d*x)*cos(c + d*x)/(8*d) - b*cos(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_381():
    f = (a + b*sin(c + d*x))*cos(c + d*x)**2
    F = a*x/2 + a*sin(c + d*x)*cos(c + d*x)/(2*d) - b*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_382():
    f = (a + b*sin(c + d*x))*sec(c + d*x)**2
    F = a*tan(c + d*x)/d + b*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_383():
    f = (a + b*sin(c + d*x))*sec(c + d*x)**4
    F = a*tan(c + d*x)**3/(3*d) + a*tan(c + d*x)/d + b*sec(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_384():
    f = (a + b*sin(c + d*x))*sec(c + d*x)**6
    F = a*tan(c + d*x)**5/(5*d) + 2*a*tan(c + d*x)**3/(3*d) + a*tan(c + d*x)/d + b*sec(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_385():
    f = (a + b*sin(c + d*x))**2*cos(c + d*x)**5
    F = a**2*sin(c + d*x)/d - a*b*cos(c + d*x)**6/(3*d) + b**2*sin(c + d*x)**7/(7*d) + (a**2 - 2*b**2)*sin(c + d*x)**5/(5*d) - (2*a**2 - b**2)*sin(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_386():
    f = (a + b*sin(c + d*x))**2*cos(c + d*x)**3
    F = a*(a + b*sin(c + d*x))**4/(2*b**3*d) - (a + b*sin(c + d*x))**5/(5*b**3*d) - (a + b*sin(c + d*x))**3*(a**2 - b**2)/(3*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_387():
    f = (a + b*sin(c + d*x))**2*cos(c + d*x)
    F = (a + b*sin(c + d*x))**3/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_388():
    f = (a + b*sin(c + d*x))**2*sec(c + d*x)
    F = -b**2*sin(c + d*x)/d + (a - b)**2*log(sin(c + d*x) + 1)/(2*d) - (a + b)**2*log(1 - sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_389():
    f = (a + b*sin(c + d*x))**2*sec(c + d*x)**3
    F = (a + b*sin(c + d*x))*(a*sin(c + d*x) + b)*sec(c + d*x)**2/(2*d) + (a**2 - b**2)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_390():
    f = (a + b*sin(c + d*x))**2*sec(c + d*x)**5
    F = (a + b*sin(c + d*x))*(a*sin(c + d*x) + b)*sec(c + d*x)**4/(4*d) + (3*a**2 - b**2)*atanh(sin(c + d*x))/(8*d) + (2*a*b + (3*a**2 - b**2)*sin(c + d*x))*sec(c + d*x)**2/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_391():
    f = (a + b*sin(c + d*x))**2*cos(c + d*x)**6
    F = -9*a*b*cos(c + d*x)**7/(56*d) - b*(a + b*sin(c + d*x))*cos(c + d*x)**7/(8*d) + x*(5*a**2/16 + 5*b**2/128) + (8*a**2 + b**2)*sin(c + d*x)*cos(c + d*x)**5/(48*d) + (40*a**2 + 5*b**2)*sin(c + d*x)*cos(c + d*x)**3/(192*d) + (40*a**2 + 5*b**2)*sin(c + d*x)*cos(c + d*x)/(128*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_392():
    f = (a + b*sin(c + d*x))**2*cos(c + d*x)**4
    F = -7*a*b*cos(c + d*x)**5/(30*d) - b*(a + b*sin(c + d*x))*cos(c + d*x)**5/(6*d) + x*(3*a**2/8 + b**2/16) + (6*a**2 + b**2)*sin(c + d*x)*cos(c + d*x)**3/(24*d) + (6*a**2 + b**2)*sin(c + d*x)*cos(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_393():
    f = (a + b*sin(c + d*x))**2*cos(c + d*x)**2
    F = -5*a*b*cos(c + d*x)**3/(12*d) - b*(a + b*sin(c + d*x))*cos(c + d*x)**3/(4*d) + x*(a**2/2 + b**2/8) + (4*a**2 + b**2)*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_394():
    f = (a + b*sin(c + d*x))**2*sec(c + d*x)**2
    F = a*b*cos(c + d*x)/d - b**2*x + (a + b*sin(c + d*x))*(a*sin(c + d*x) + b)*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_395():
    f = (a + b*sin(c + d*x))**2*sec(c + d*x)**4
    F = a*b*sec(c + d*x)/(3*d) + (a + b*sin(c + d*x))*(a*sin(c + d*x) + b)*sec(c + d*x)**3/(3*d) + (2*a**2 - b**2)*tan(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_396():
    f = (a + b*sin(c + d*x))**2*sec(c + d*x)**6
    F = a*b*sec(c + d*x)**3/(5*d) + (a + b*sin(c + d*x))*(a*sin(c + d*x) + b)*sec(c + d*x)**5/(5*d) + (4*a**2 - b**2)*tan(c + d*x)**3/(15*d) + (4*a**2 - b**2)*tan(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_397():
    f = (a + b*sin(c + d*x))**2*sec(c + d*x)**8
    F = a*b*sec(c + d*x)**5/(7*d) + (a + b*sin(c + d*x))*(a*sin(c + d*x) + b)*sec(c + d*x)**7/(7*d) + (6*a**2 - b**2)*tan(c + d*x)**5/(35*d) + (6*a**2 - b**2)*tan(c + d*x)/(7*d) + (12*a**2 - 2*b**2)*tan(c + d*x)**3/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_398():
    f = (a + b*sin(c + d*x))**3*cos(c + d*x)**5
    F = -4*a*(a + b*sin(c + d*x))**7/(7*b**5*d) - 4*a*(a + b*sin(c + d*x))**5*(a**2 - b**2)/(5*b**5*d) + (a + b*sin(c + d*x))**8/(8*b**5*d) + (a + b*sin(c + d*x))**6*(3*a**2 - b**2)/(3*b**5*d) + (a + b*sin(c + d*x))**4*(a**2 - b**2)**2/(4*b**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_399():
    f = (a + b*sin(c + d*x))**3*cos(c + d*x)**3
    F = 2*a*(a + b*sin(c + d*x))**5/(5*b**3*d) - (a + b*sin(c + d*x))**6/(6*b**3*d) - (a + b*sin(c + d*x))**4*(a**2 - b**2)/(4*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_400():
    f = (a + b*sin(c + d*x))**3*cos(c + d*x)
    F = (a + b*sin(c + d*x))**4/(4*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_401():
    f = (a + b*sin(c + d*x))**3*sec(c + d*x)
    F = -3*a*b**2*sin(c + d*x)/d - b**3*sin(c + d*x)**2/(2*d) + (a - b)**3*log(sin(c + d*x) + 1)/(2*d) - (a + b)**3*log(1 - sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_402():
    f = (a + b*sin(c + d*x))**3*sec(c + d*x)**3
    F = a*b**2*sin(c + d*x)/(2*d) - (a - 2*b)*(a + b)**2*log(1 - sin(c + d*x))/(4*d) + (a - b)**2*(a + 2*b)*log(sin(c + d*x) + 1)/(4*d) + (a + b*sin(c + d*x))**2*(a*sin(c + d*x) + b)*sec(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_403():
    f = (a + b*sin(c + d*x))**3*sec(c + d*x)**5
    F = 3*a*(a + b*sin(c + d*x))*(a*sin(c + d*x) + b)*sec(c + d*x)**2/(8*d) + 3*a*(a**2 - b**2)*atanh(sin(c + d*x))/(8*d) + (a + b*sin(c + d*x))**3*tan(c + d*x)*sec(c + d*x)**3/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_404():
    f = (a + b*sin(c + d*x))**3*cos(c + d*x)**4
    F = -3*a*b*(a + b*sin(c + d*x))*cos(c + d*x)**5/(14*d) + 3*a*x*(2*a**2 + b**2)/16 + a*(2*a**2 + b**2)*sin(c + d*x)*cos(c + d*x)**3/(8*d) + 3*a*(2*a**2 + b**2)*sin(c + d*x)*cos(c + d*x)/(16*d) - b*(a + b*sin(c + d*x))**2*cos(c + d*x)**5/(7*d) - b*(17*a**2 + 4*b**2)*cos(c + d*x)**5/(70*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_405():
    f = (a + b*sin(c + d*x))**3*cos(c + d*x)**2
    F = -7*a*b*(a + b*sin(c + d*x))*cos(c + d*x)**3/(20*d) + a*x*(4*a**2 + 3*b**2)/8 + a*(4*a**2 + 3*b**2)*sin(c + d*x)*cos(c + d*x)/(8*d) - b*(a + b*sin(c + d*x))**2*cos(c + d*x)**3/(5*d) - b*(27*a**2 + 8*b**2)*cos(c + d*x)**3/(60*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_406():
    f = (a + b*sin(c + d*x))**3*sec(c + d*x)**2
    F = -3*a*b**2*x + a*b**2*sin(c + d*x)*cos(c + d*x)/d + 2*b*(a**2 + b**2)*cos(c + d*x)/d + (a + b*sin(c + d*x))**2*(a*sin(c + d*x) + b)*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_407():
    f = (a + b*sin(c + d*x))**3*sec(c + d*x)**4
    F = 2*a*(a**2 - b**2)*tan(c + d*x)/(3*d) + 2*b*(a**2 - b**2)*sec(c + d*x)/(3*d) + (a + b*sin(c + d*x))**2*(a*sin(c + d*x) + b)*sec(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_408():
    f = (a + b*sin(c + d*x))**3*sec(c + d*x)**6
    F = 2*a*(4*a**2 - 3*b**2)*tan(c + d*x)/(15*d) + 2*b*(2*a**2 - b**2)*sec(c + d*x)/(15*d) + (a + b*sin(c + d*x))**2*(a*sin(c + d*x) + b)*sec(c + d*x)**5/(5*d) + 2*(a + b*sin(c + d*x))*(a*b + (2*a**2 - b**2)*sin(c + d*x))*sec(c + d*x)**3/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_409():
    f = (a + b*sin(c + d*x))**3*sec(c + d*x)**8
    F = 4*a*(2*a**2 - b**2)*tan(c + d*x)**3/(35*d) + 12*a*(2*a**2 - b**2)*tan(c + d*x)/(35*d) + 2*b*(3*a**2 - b**2)*sec(c + d*x)**3/(35*d) + (a + b*sin(c + d*x))**2*(a*sin(c + d*x) + b)*sec(c + d*x)**7/(7*d) + 2*(a + b*sin(c + d*x))*(2*a*b + (3*a**2 - b**2)*sin(c + d*x))*sec(c + d*x)**5/(35*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_410():
    f = (a + b*sin(c + d*x))**3*sec(c + d*x)**10
    F = 2*a*(8*a**2 - 3*b**2)*tan(c + d*x)**5/(105*d) + 4*a*(8*a**2 - 3*b**2)*tan(c + d*x)**3/(63*d) + 2*a*(8*a**2 - 3*b**2)*tan(c + d*x)/(21*d) + 2*b*(4*a**2 - b**2)*sec(c + d*x)**5/(63*d) + (a + b*sin(c + d*x))**2*(a*sin(c + d*x) + b)*sec(c + d*x)**9/(9*d) + 2*(a + b*sin(c + d*x))*(3*a*b + (4*a**2 - b**2)*sin(c + d*x))*sec(c + d*x)**7/(63*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_411():
    f = (a + b*sin(c + d*x))**8*cos(c + d*x)**5
    F = -a*(a + b*sin(c + d*x))**12/(3*b**5*d) - 2*a*(a + b*sin(c + d*x))**10*(a**2 - b**2)/(5*b**5*d) + (a + b*sin(c + d*x))**13/(13*b**5*d) + (a + b*sin(c + d*x))**11*(6*a**2 - 2*b**2)/(11*b**5*d) + (a + b*sin(c + d*x))**9*(a**2 - b**2)**2/(9*b**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_412():
    f = (a + b*sin(c + d*x))**8*cos(c + d*x)**3
    F = a*(a + b*sin(c + d*x))**10/(5*b**3*d) - (a + b*sin(c + d*x))**11/(11*b**3*d) - (a + b*sin(c + d*x))**9*(a**2 - b**2)/(9*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_413():
    f = (a + b*sin(c + d*x))**8*cos(c + d*x)
    F = (a + b*sin(c + d*x))**9/(9*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_414():
    f = (a + b*sin(c + d*x))**8*sec(c + d*x)
    F = -4*a*b**7*sin(c + d*x)**6/(3*d) - 2*a*b**5*(7*a**2 + b**2)*sin(c + d*x)**4/d - 4*a*b**3*(7*a**4 + 7*a**2*b**2 + b**4)*sin(c + d*x)**2/d - b**8*sin(c + d*x)**7/(7*d) - b**6*(28*a**2 + b**2)*sin(c + d*x)**5/(5*d) - b**4*(70*a**4 + 28*a**2*b**2 + b**4)*sin(c + d*x)**3/(3*d) - b**2*(28*a**6 + 70*a**4*b**2 + 28*a**2*b**4 + b**6)*sin(c + d*x)/d + (a - b)**8*log(sin(c + d*x) + 1)/(2*d) - (a + b)**8*log(1 - sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_415():
    f = (a + b*sin(c + d*x))**8*sec(c + d*x)**3
    F = a*b**7*sin(c + d*x)**6/(2*d) + 3*a*b**5*(7*a**2 + 4*b**2)*sin(c + d*x)**4/(2*d) + a*b**3*(35*a**4 + 112*a**2*b**2 + 24*b**4)*sin(c + d*x)**2/(2*d) + 7*b**6*(5*a**2 + b**2)*sin(c + d*x)**5/(10*d) + 7*b**4*(15*a**4 + 20*a**2*b**2 + b**4)*sin(c + d*x)**3/(6*d) + 7*b**2*(3*a**6 + 30*a**4*b**2 + 20*a**2*b**4 + b**6)*sin(c + d*x)/(2*d) - (a - 7*b)*(a + b)**7*log(1 - sin(c + d*x))/(4*d) + (a - b)**7*(a + 7*b)*log(sin(c + d*x) + 1)/(4*d) + (a + b*sin(c + d*x))**7*(a*sin(c + d*x) + b)*sec(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_416():
    f = (a + b*sin(c + d*x))**8*sec(c + d*x)**5
    F = -a*b**7*(-3*a**2/b**2 + 13)*sin(c + d*x)**4/(8*d) + a*b**3*(15*a**4 - 77*a**2*b**2 - 48*b**4)*sin(c + d*x)**2/(4*d) + 5*b**4*(9*a**4 - 42*a**2*b**2 - 7*b**4)*sin(c + d*x)**3/(24*d) + 5*b**2*(6*a**6 - 35*a**4*b**2 - 84*a**2*b**4 - 7*b**6)*sin(c + d*x)/(8*d) + (a - b)**6*(3*a**2 + 18*a*b + 35*b**2)*log(sin(c + d*x) + 1)/(16*d) - (a + b)**6*(3*a**2 - 18*a*b + 35*b**2)*log(1 - sin(c + d*x))/(16*d) + (a + b*sin(c + d*x))**7*(a*sin(c + d*x) + b)*sec(c + d*x)**4/(4*d) - (a + b*sin(c + d*x))**5*(-a*(3*a**2 - 11*b**2)*sin(c + d*x) + b*(a**2 + 7*b**2))*sec(c + d*x)**2/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_417():
    f = (a + b*sin(c + d*x))**8*cos(c + d*x)**2
    F = -17*a*b*(a + b*sin(c + d*x))**6*cos(c + d*x)**3/(90*d) - a*b*(a + b*sin(c + d*x))**4*(112*a**2 + 109*b**2)*cos(c + d*x)**3/(336*d) - 13*a*b*(a + b*sin(c + d*x))**2*(112*a**4 + 348*a**2*b**2 + 101*b**4)*cos(c + d*x)**3/(3360*d) - 11*a*b*(1792*a**6 + 10536*a**4*b**2 + 9588*a**2*b**4 + 1289*b**6)*cos(c + d*x)**3/(40320*d) - b*(a + b*sin(c + d*x))**7*cos(c + d*x)**3/(10*d) - b*(a + b*sin(c + d*x))**5*(64*a**2 + 21*b**2)*cos(c + d*x)**3/(240*d) - b*(a + b*sin(c + d*x))**3*(784*a**4 + 1500*a**2*b**2 + 147*b**4)*cos(c + d*x)**3/(2016*d) - b*(a + b*sin(c + d*x))*(6272*a**6 + 28088*a**4*b**2 + 15956*a**2*b**4 + 735*b**6)*cos(c + d*x)**3/(13440*d) + x*(a**8/2 + 7*a**6*b**2/2 + 35*a**4*b**4/8 + 35*a**2*b**6/32 + 7*b**8/256) + (128*a**8 + 896*a**6*b**2 + 1120*a**4*b**4 + 280*a**2*b**6 + 7*b**8)*sin(c + d*x)*cos(c + d*x)/(256*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_418():
    f = (a + b*sin(c + d*x))**8*sec(c + d*x)**2
    F = a*b*(a + b*sin(c + d*x))**6*cos(c + d*x)/d + a*b*(a + b*sin(c + d*x))**4*(30*a**2 + 113*b**2)*cos(c + d*x)/(30*d) + a*b*(a + b*sin(c + d*x))**2*(40*a**4 + 624*a**2*b**2 + 337*b**4)*cos(c + d*x)/(40*d) + a*b*(40*a**6 + 1664*a**4*b**2 + 2789*a**2*b**4 + 512*b**6)*cos(c + d*x)/(20*d) - 7*b**2*x*(64*a**6 + 240*a**4*b**2 + 120*a**2*b**4 + 5*b**6)/16 + b**2*(80*a**6 + 2248*a**4*b**2 + 2502*a**2*b**4 + 175*b**6)*sin(c + d*x)*cos(c + d*x)/(80*d) + b*(a + b*sin(c + d*x))**5*(6*a**2 + 7*b**2)*cos(c + d*x)/(6*d) + b*(a + b*sin(c + d*x))**3*(120*a**4 + 992*a**2*b**2 + 175*b**4)*cos(c + d*x)/(120*d) + (a + b*sin(c + d*x))**7*(a*sin(c + d*x) + b)*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_419():
    f = (a + b*sin(c + d*x))**8*sec(c + d*x)**4
    F = a*b*(a + b*sin(c + d*x))**4*(2*a**2 - 13*b**2)*cos(c + d*x)/(3*d) + a*b*(a + b*sin(c + d*x))**2*(8*a**4 - 88*a**2*b**2 - 151*b**4)*cos(c + d*x)/(12*d) + a*b*(8*a**6 - 104*a**4*b**2 - 803*a**2*b**4 - 256*b**6)*cos(c + d*x)/(6*d) + 35*b**4*x*(16*a**4 + 16*a**2*b**2 + b**4)/8 + b**2*(16*a**6 - 200*a**4*b**2 - 866*a**2*b**4 - 105*b**6)*sin(c + d*x)*cos(c + d*x)/(24*d) + b*(a + b*sin(c + d*x))**5*(2*a**2 - 7*b**2)*cos(c + d*x)/(3*d) + b*(a + b*sin(c + d*x))**3*(8*a**4 - 72*a**2*b**2 - 35*b**4)*cos(c + d*x)/(12*d) + (a + b*sin(c + d*x))**7*(a*sin(c + d*x) + b)*sec(c + d*x)**3/(3*d) - (a + b*sin(c + d*x))**6*(5*a*b - (2*a**2 - 7*b**2)*sin(c + d*x))*sec(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_420():
    f = (a + b*sin(c + d*x))**8*sec(c + d*x)**6
    F = 4*a*b*(a + b*sin(c + d*x))**4*(2*a**2 + b**2)*cos(c + d*x)/(15*d) + a*b*(a + b*sin(c + d*x))**2*(8*a**4 - 32*a**2*b**2 + 87*b**4)*cos(c + d*x)/(15*d) + 2*a*b*(8*a**6 - 48*a**4*b**2 + 163*a**2*b**4 + 192*b**6)*cos(c + d*x)/(15*d) - 7*b**6*x*(8*a**2 + b**2)/2 + b**2*(16*a**6 - 88*a**4*b**2 + 282*a**2*b**4 + 105*b**6)*sin(c + d*x)*cos(c + d*x)/(30*d) + b*(a + b*sin(c + d*x))**3*(8*a**4 - 16*a**2*b**2 + 35*b**4)*cos(c + d*x)/(15*d) + (a + b*sin(c + d*x))**7*(a*sin(c + d*x) + b)*sec(c + d*x)**5/(5*d) - (a + b*sin(c + d*x))**6*(3*a*b - (4*a**2 - 7*b**2)*sin(c + d*x))*sec(c + d*x)**3/(15*d) - 4*(a + b*sin(c + d*x))**5*(-a*(2*a**2 + b**2)*sin(c + d*x) + b*(4*a**2 - 7*b**2))*sec(c + d*x)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_421():
    f = (a + b*sin(c + d*x))**8*sec(c + d*x)**8
    F = 2*a*b*(a + b*sin(c + d*x))**2*(24*a**4 - 40*a**2*b**2 + 9*b**4)*cos(c + d*x)/(105*d) + 4*a*b*(24*a**6 - 88*a**4*b**2 + 125*a**2*b**4 - 96*b**6)*cos(c + d*x)/(105*d) + b**8*x + b**2*(48*a**6 - 152*a**4*b**2 + 174*a**2*b**4 - 105*b**6)*sin(c + d*x)*cos(c + d*x)/(105*d) + 2*b*(a + b*sin(c + d*x))**3*(24*a**4 + 8*a**2*b**2 - 35*b**4)*cos(c + d*x)/(105*d) + (a + b*sin(c + d*x))**7*(a*sin(c + d*x) + b)*sec(c + d*x)**7/(7*d) - (a + b*sin(c + d*x))**6*(a*b - (6*a**2 - 7*b**2)*sin(c + d*x))*sec(c + d*x)**5/(35*d) - 2*(a + b*sin(c + d*x))**5*(-a*(12*a**2 - 11*b**2)*sin(c + d*x) + b*(6*a**2 - 7*b**2))*sec(c + d*x)**3/(105*d) - 2*(a + b*sin(c + d*x))**4*(3*a*b*(12*a**2 - 11*b**2) - (24*a**4 + 8*a**2*b**2 - 35*b**4)*sin(c + d*x))*sec(c + d*x)/(105*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_422():
    f = (a + b*sin(c + d*x))**8*sec(c + d*x)**10
    F = 128*a**2*(a**2 - b**2)**3*tan(c + d*x)/(315*d) + 128*a*b*(a**2 - b**2)**3*sec(c + d*x)/(315*d) + 16*a*(a + b*sin(c + d*x))**4*(a**2 - b**2)*(a*sin(c + d*x) + b)*sec(c + d*x)**5/(105*d) + 64*a*(a + b*sin(c + d*x))**2*(a**2 - b**2)**2*(a*sin(c + d*x) + b)*sec(c + d*x)**3/(315*d) + (a + b*sin(c + d*x))**7*(a*sin(c + d*x) + b)*sec(c + d*x)**9/(9*d) + (a + b*sin(c + d*x))**6*(a*b + (8*a**2 - 7*b**2)*sin(c + d*x))*sec(c + d*x)**7/(63*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_423():
    f = cos(c + d*x)**5/(a + b*sin(c + d*x))
    F = -a*sin(c + d*x)**3/(3*b**2*d) - a*(a**2 - 2*b**2)*sin(c + d*x)/(b**4*d) + sin(c + d*x)**4/(4*b*d) + (a**2 - 2*b**2)*sin(c + d*x)**2/(2*b**3*d) + (a**2 - b**2)**2*log(a + b*sin(c + d*x))/(b**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_424():
    f = cos(c + d*x)**3/(a + b*sin(c + d*x))
    F = a*sin(c + d*x)/(b**2*d) - sin(c + d*x)**2/(2*b*d) - (a**2 - b**2)*log(a + b*sin(c + d*x))/(b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_425():
    f = cos(c + d*x)/(a + b*sin(c + d*x))
    F = log(a + b*sin(c + d*x))/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_426():
    f = sec(c + d*x)/(a + b*sin(c + d*x))
    F = -b*log(a + b*sin(c + d*x))/(d*(a**2 - b**2)) - log(1 - sin(c + d*x))/(d*(2*a + 2*b)) + log(sin(c + d*x) + 1)/(d*(2*a - 2*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_427():
    f = sec(c + d*x)**3/(a + b*sin(c + d*x))
    F = b**3*log(a + b*sin(c + d*x))/(d*(a**2 - b**2)**2) + (a - 2*b)*log(sin(c + d*x) + 1)/(4*d*(a - b)**2) - (-a*sin(c + d*x) + b)*sec(c + d*x)**2/(d*(2*a**2 - 2*b**2)) - (a + 2*b)*log(1 - sin(c + d*x))/(4*d*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_428():
    f = sec(c + d*x)**5/(a + b*sin(c + d*x))
    F = -b**5*log(a + b*sin(c + d*x))/(d*(a**2 - b**2)**3) - (-a*sin(c + d*x) + b)*sec(c + d*x)**4/(d*(4*a**2 - 4*b**2)) + (a*(3*a**2 - 7*b**2)*sin(c + d*x) + 4*b**3)*sec(c + d*x)**2/(8*d*(a**2 - b**2)**2) - (3*a**2 + 9*a*b + 8*b**2)*log(1 - sin(c + d*x))/(16*d*(a + b)**3) + (3*a**2 - 9*a*b + 8*b**2)*log(sin(c + d*x) + 1)/(16*d*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_429():
    f = cos(c + d*x)**6/(a + b*sin(c + d*x))
    F = a*x*(8*a**4 - 20*a**2*b**2 + 15*b**4)/(8*b**6) + cos(c + d*x)**5/(5*b*d) - (4*a**2 - 3*a*b*sin(c + d*x) - 4*b**2)*cos(c + d*x)**3/(12*b**3*d) + (-a*b*(4*a**2 - 7*b**2)*sin(c + d*x) + 8*(a**2 - b**2)**2)*cos(c + d*x)/(8*b**5*d) - 2*(a**2 - b**2)**(sympy.S(5)/2)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(b**6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_430():
    f = cos(c + d*x)**4/(a + b*sin(c + d*x))
    F = -a*x*(2*a**2 - 3*b**2)/(2*b**4) + cos(c + d*x)**3/(3*b*d) - (2*a**2 - a*b*sin(c + d*x) - 2*b**2)*cos(c + d*x)/(2*b**3*d) + 2*(a**2 - b**2)**(sympy.S(3)/2)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_431():
    f = cos(c + d*x)**2/(a + b*sin(c + d*x))
    F = a*x/b**2 + cos(c + d*x)/(b*d) - 2*sqrt(a**2 - b**2)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_432():
    f = sec(c + d*x)**2/(a + b*sin(c + d*x))
    F = -2*b**2*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(d*(a**2 - b**2)**(sympy.S(3)/2)) - (-a*sin(c + d*x) + b)*sec(c + d*x)/(d*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_433():
    f = sec(c + d*x)**4/(a + b*sin(c + d*x))
    F = 2*b**4*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(d*(a**2 - b**2)**(sympy.S(5)/2)) - (-a*sin(c + d*x) + b)*sec(c + d*x)**3/(d*(3*a**2 - 3*b**2)) + (a*(2*a**2 - 5*b**2)*sin(c + d*x) + 3*b**3)*sec(c + d*x)/(3*d*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_434():
    f = sec(c + d*x)**6/(a + b*sin(c + d*x))
    F = -2*b**6*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(d*(a**2 - b**2)**(sympy.S(7)/2)) - (-a*sin(c + d*x) + b)*sec(c + d*x)**5/(d*(5*a**2 - 5*b**2)) + (a*(4*a**2 - 9*b**2)*sin(c + d*x) + 5*b**3)*sec(c + d*x)**3/(15*d*(a**2 - b**2)**2) - (-a*(8*a**4 - 26*a**2*b**2 + 33*b**4)*sin(c + d*x) + 15*b**5)*sec(c + d*x)/(15*d*(a**2 - b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_435():
    f = cos(c + d*x)**7/(a + b*sin(c + d*x))**2
    F = a*sin(c + d*x)**4/(2*b**3*d) + a*(2*a**2 - 3*b**2)*sin(c + d*x)**2/(b**5*d) + 6*a*(a**2 - b**2)**2*log(a + b*sin(c + d*x))/(b**7*d) - sin(c + d*x)**5/(5*b**2*d) - (a**2 - b**2)*sin(c + d*x)**3/(b**4*d) - (5*a**4 - 9*a**2*b**2 + 3*b**4)*sin(c + d*x)/(b**6*d) + (a**2 - b**2)**3/(b**7*d*(a + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_436():
    f = cos(c + d*x)**5/(a + b*sin(c + d*x))**2
    F = -a*sin(c + d*x)**2/(b**3*d) - 4*a*(a**2 - b**2)*log(a + b*sin(c + d*x))/(b**5*d) + sin(c + d*x)**3/(3*b**2*d) + (3*a**2 - 2*b**2)*sin(c + d*x)/(b**4*d) - (a**2 - b**2)**2/(b**5*d*(a + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_437():
    f = cos(c + d*x)**3/(a + b*sin(c + d*x))**2
    F = 2*a*log(a + b*sin(c + d*x))/(b**3*d) - sin(c + d*x)/(b**2*d) + (a**2 - b**2)/(b**3*d*(a + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_438():
    f = cos(c + d*x)/(a + b*sin(c + d*x))**2
    F = -1/(b*d*(a + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_439():
    f = sec(c + d*x)/(a + b*sin(c + d*x))**2
    F = -2*a*b*log(a + b*sin(c + d*x))/(d*(a**2 - b**2)**2) + b/(d*(a + b*sin(c + d*x))*(a**2 - b**2)) - log(1 - sin(c + d*x))/(2*d*(a + b)**2) + log(sin(c + d*x) + 1)/(2*d*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_440():
    f = sec(c + d*x)**3/(a + b*sin(c + d*x))**2
    F = 4*a*b**3*log(a + b*sin(c + d*x))/(d*(a**2 - b**2)**3) - b*(a**2 + 3*b**2)/(2*d*(a + b*sin(c + d*x))*(a**2 - b**2)**2) + (a - 3*b)*log(sin(c + d*x) + 1)/(4*d*(a - b)**3) - (-a*sin(c + d*x) + b)*sec(c + d*x)**2/(d*(a + b*sin(c + d*x))*(2*a**2 - 2*b**2)) - (a + 3*b)*log(1 - sin(c + d*x))/(4*d*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_441():
    f = sec(c + d*x)**5/(a + b*sin(c + d*x))**2
    F = -6*a*b**5*log(a + b*sin(c + d*x))/(d*(a**2 - b**2)**4) - 3*b*(a**4 - 4*a**2*b**2 - 5*b**4)/(8*d*(a + b*sin(c + d*x))*(a**2 - b**2)**3) - (-a*sin(c + d*x) + b)*sec(c + d*x)**4/(d*(a + b*sin(c + d*x))*(4*a**2 - 4*b**2)) + (3*a*(a**2 - 3*b**2)*sin(c + d*x) + b*(a**2 + 5*b**2))*sec(c + d*x)**2/(8*d*(a + b*sin(c + d*x))*(a**2 - b**2)**2) + (-3*a**2 - 12*a*b - 15*b**2)*log(1 - sin(c + d*x))/(16*d*(a + b)**4) + (3*a**2 - 12*a*b + 15*b**2)*log(sin(c + d*x) + 1)/(16*d*(a - b)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_442():
    f = cos(c + d*x)**6/(a + b*sin(c + d*x))**2
    F = 10*a*(a**2 - b**2)**(sympy.S(3)/2)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(b**6*d) - cos(c + d*x)**5/(b*d*(a + b*sin(c + d*x))) + 5*(4*a - 3*b*sin(c + d*x))*cos(c + d*x)**3/(12*b**3*d) - 5*(8*a*(a**2 - b**2) - b*(4*a**2 - 3*b**2)*sin(c + d*x))*cos(c + d*x)/(8*b**5*d) - x*(40*a**4 - 60*a**2*b**2 + 15*b**4)/(8*b**6)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_443():
    f = cos(c + d*x)**4/(a + b*sin(c + d*x))**2
    F = -6*a*sqrt(a**2 - b**2)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(b**4*d) - cos(c + d*x)**3/(b*d*(a + b*sin(c + d*x))) + 3*(2*a - b*sin(c + d*x))*cos(c + d*x)/(2*b**3*d) + x*(6*a**2 - 3*b**2)/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_444():
    f = cos(c + d*x)**2/(a + b*sin(c + d*x))**2
    F = 2*a*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(b**2*d*sqrt(a**2 - b**2)) - cos(c + d*x)/(b*d*(a + b*sin(c + d*x))) - x/b**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_445():
    f = sec(c + d*x)**2/(a + b*sin(c + d*x))**2
    F = -6*a*b**2*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(d*(a**2 - b**2)**(sympy.S(5)/2)) + b*sec(c + d*x)/(d*(a + b*sin(c + d*x))*(a**2 - b**2)) - (3*a*b - (a**2 + 2*b**2)*sin(c + d*x))*sec(c + d*x)/(d*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_446():
    f = sec(c + d*x)**4/(a + b*sin(c + d*x))**2
    F = 10*a*b**4*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(d*(a**2 - b**2)**(sympy.S(7)/2)) + b*sec(c + d*x)**3/(d*(a + b*sin(c + d*x))*(a**2 - b**2)) - (5*a*b - (a**2 + 4*b**2)*sin(c + d*x))*sec(c + d*x)**3/(3*d*(a**2 - b**2)**2) + (15*a*b**3 + (2*a**4 - 9*a**2*b**2 - 8*b**4)*sin(c + d*x))*sec(c + d*x)/(3*d*(a**2 - b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_447():
    f = cos(c + d*x)**7/(a + b*sin(c + d*x))**3
    F = a*sin(c + d*x)**3/(b**4*d) + a*(10*a**2 - 9*b**2)*sin(c + d*x)/(b**6*d) - 6*a*(a**2 - b**2)**2/(b**7*d*(a + b*sin(c + d*x))) - sin(c + d*x)**4/(4*b**3*d) - (6*a**2 - 3*b**2)*sin(c + d*x)**2/(2*b**5*d) - (15*a**4 - 18*a**2*b**2 + 3*b**4)*log(a + b*sin(c + d*x))/(b**7*d) + (a**2 - b**2)**3/(2*b**7*d*(a + b*sin(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_448():
    f = cos(c + d*x)**5/(a + b*sin(c + d*x))**3
    F = -3*a*sin(c + d*x)/(b**4*d) + 4*a*(a**2 - b**2)/(b**5*d*(a + b*sin(c + d*x))) + sin(c + d*x)**2/(2*b**3*d) + (6*a**2 - 2*b**2)*log(a + b*sin(c + d*x))/(b**5*d) - (a**2 - b**2)**2/(2*b**5*d*(a + b*sin(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_449():
    f = cos(c + d*x)**3/(a + b*sin(c + d*x))**3
    F = -2*a/(b**3*d*(a + b*sin(c + d*x))) - log(a + b*sin(c + d*x))/(b**3*d) + (a**2 - b**2)/(2*b**3*d*(a + b*sin(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_450():
    f = cos(c + d*x)/(a + b*sin(c + d*x))**3
    F = -1/(2*b*d*(a + b*sin(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_451():
    f = sec(c + d*x)/(a + b*sin(c + d*x))**3
    F = 2*a*b/(d*(a + b*sin(c + d*x))*(a**2 - b**2)**2) - b*(3*a**2 + b**2)*log(a + b*sin(c + d*x))/(d*(a**2 - b**2)**3) + b/(d*(a + b*sin(c + d*x))**2*(2*a**2 - 2*b**2)) - log(1 - sin(c + d*x))/(2*d*(a + b)**3) + log(sin(c + d*x) + 1)/(2*d*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_452():
    f = sec(c + d*x)**3/(a + b*sin(c + d*x))**3
    F = -a*b*(a**2 + 11*b**2)/(2*d*(a + b*sin(c + d*x))*(a**2 - b**2)**3) + 2*b**3*(5*a**2 + b**2)*log(a + b*sin(c + d*x))/(d*(a**2 - b**2)**4) - b*(a**2 + 2*b**2)/(2*d*(a + b*sin(c + d*x))**2*(a**2 - b**2)**2) + (a - 4*b)*log(sin(c + d*x) + 1)/(4*d*(a - b)**4) - (-a*sin(c + d*x) + b)*sec(c + d*x)**2/(d*(a + b*sin(c + d*x))**2*(2*a**2 - 2*b**2)) - (a + 4*b)*log(1 - sin(c + d*x))/(4*d*(a + b)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_453():
    f = sec(c + d*x)**5/(a + b*sin(c + d*x))**3
    F = -3*a*b*(a**4 - 6*a**2*b**2 - 27*b**4)/(8*d*(a + b*sin(c + d*x))*(a**2 - b**2)**4) - 3*b**5*(7*a**2 + b**2)*log(a + b*sin(c + d*x))/(d*(a**2 - b**2)**5) - 3*b*(a**4 - 5*a**2*b**2 - 4*b**4)/(8*d*(a + b*sin(c + d*x))**2*(a**2 - b**2)**3) - (-a*sin(c + d*x) + b)*sec(c + d*x)**4/(d*(a + b*sin(c + d*x))**2*(4*a**2 - 4*b**2)) + (a*(3*a**2 - 11*b**2)*sin(c + d*x) + 2*b*(a**2 + 3*b**2))*sec(c + d*x)**2/(8*d*(a + b*sin(c + d*x))**2*(a**2 - b**2)**2) + (-3*a**2 - 15*a*b - 24*b**2)*log(1 - sin(c + d*x))/(16*d*(a + b)**5) + (3*a**2 - 15*a*b + 24*b**2)*log(sin(c + d*x) + 1)/(16*d*(a - b)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_454():
    f = cos(c + d*x)**6/(a + b*sin(c + d*x))**3
    F = 5*a*x*(4*a**2 - 3*b**2)/(2*b**6) - cos(c + d*x)**5/(2*b*d*(a + b*sin(c + d*x))**2) - 5*(4*a + b*sin(c + d*x))*cos(c + d*x)**3/(6*b**3*d*(a + b*sin(c + d*x))) + 5*(4*a**2 - 2*a*b*sin(c + d*x) - b**2)*cos(c + d*x)/(2*b**5*d) - (20*a**4 - 25*a**2*b**2 + 5*b**4)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(b**6*d*sqrt(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_455():
    f = cos(c + d*x)**4/(a + b*sin(c + d*x))**3
    F = -3*a*x/b**4 - cos(c + d*x)**3/(2*b*d*(a + b*sin(c + d*x))**2) - 3*(2*a + b*sin(c + d*x))*cos(c + d*x)/(2*b**3*d*(a + b*sin(c + d*x))) + (6*a**2 - 3*b**2)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(b**4*d*sqrt(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_456():
    f = cos(c + d*x)**2/(a + b*sin(c + d*x))**3
    F = a*cos(c + d*x)/(2*b*d*(a + b*sin(c + d*x))*(a**2 - b**2)) + atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(d*(a**2 - b**2)**(sympy.S(3)/2)) - cos(c + d*x)/(2*b*d*(a + b*sin(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_457():
    f = sec(c + d*x)**2/(a + b*sin(c + d*x))**3
    F = 5*a*b*sec(c + d*x)/(2*d*(a + b*sin(c + d*x))*(a**2 - b**2)**2) - 3*b**2*(4*a**2 + b**2)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(d*(a**2 - b**2)**(sympy.S(7)/2)) + b*sec(c + d*x)/(d*(a + b*sin(c + d*x))**2*(2*a**2 - 2*b**2)) - (-a*(2*a**2 + 13*b**2)*sin(c + d*x) + 3*b*(4*a**2 + b**2))*sec(c + d*x)/(2*d*(a**2 - b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_458():
    f = sec(c + d*x)**4/(a + b*sin(c + d*x))**3
    F = 7*a*b*sec(c + d*x)**3/(2*d*(a + b*sin(c + d*x))*(a**2 - b**2)**2) + 5*b**4*(6*a**2 + b**2)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(d*(a**2 - b**2)**(sympy.S(9)/2)) + b*sec(c + d*x)**3/(d*(a + b*sin(c + d*x))**2*(2*a**2 - 2*b**2)) - (-a*(2*a**2 + 33*b**2)*sin(c + d*x) + 5*b*(6*a**2 + b**2))*sec(c + d*x)**3/(6*d*(a**2 - b**2)**3) + (a*(4*a**4 - 28*a**2*b**2 - 81*b**4)*sin(c + d*x) + 15*b**3*(6*a**2 + b**2))*sec(c + d*x)/(6*d*(a**2 - b**2)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_459():
    f = cos(c + d*x)**7/(a + b*sin(c + d*x))**8
    F = -3*a/(b**7*d*(a + b*sin(c + d*x))**2) - a*(5*a**2 - 3*b**2)/(b**7*d*(a + b*sin(c + d*x))**4) - a*(a**2 - b**2)**2/(b**7*d*(a + b*sin(c + d*x))**6) + 1/(b**7*d*(a + b*sin(c + d*x))) + (5*a**2 - b**2)/(b**7*d*(a + b*sin(c + d*x))**3) + (15*a**4 - 18*a**2*b**2 + 3*b**4)/(5*b**7*d*(a + b*sin(c + d*x))**5) + (a**2 - b**2)**3/(7*b**7*d*(a + b*sin(c + d*x))**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_460():
    f = cos(c + d*x)**5/(a + b*sin(c + d*x))**8
    F = a/(b**5*d*(a + b*sin(c + d*x))**4) + 2*a*(a**2 - b**2)/(3*b**5*d*(a + b*sin(c + d*x))**6) - 1/(3*b**5*d*(a + b*sin(c + d*x))**3) - (6*a**2 - 2*b**2)/(5*b**5*d*(a + b*sin(c + d*x))**5) - (a**2 - b**2)**2/(7*b**5*d*(a + b*sin(c + d*x))**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_461():
    f = cos(c + d*x)**3/(a + b*sin(c + d*x))**8
    F = -a/(3*b**3*d*(a + b*sin(c + d*x))**6) + 1/(5*b**3*d*(a + b*sin(c + d*x))**5) + (a**2 - b**2)/(7*b**3*d*(a + b*sin(c + d*x))**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_462():
    f = cos(c + d*x)/(a + b*sin(c + d*x))**8
    F = -1/(7*b*d*(a + b*sin(c + d*x))**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_463():
    f = sec(c + d*x)/(a + b*sin(c + d*x))**8
    F = -8*a*b*(a**2 + b**2)*(a**4 + 6*a**2*b**2 + b**4)*log(a + b*sin(c + d*x))/(d*(a**2 - b**2)**8) + a*b*(a**2 + 3*b**2)*(3*a**2 + b**2)/(d*(a + b*sin(c + d*x))**2*(a**2 - b**2)**6) + a*b*(a**2 + b**2)/(d*(a + b*sin(c + d*x))**4*(a**2 - b**2)**4) + a*b/(3*d*(a + b*sin(c + d*x))**6*(a**2 - b**2)**2) + b*(7*a**6 + 35*a**4*b**2 + 21*a**2*b**4 + b**6)/(d*(a + b*sin(c + d*x))*(a**2 - b**2)**7) + b*(5*a**4 + 10*a**2*b**2 + b**4)/(3*d*(a + b*sin(c + d*x))**3*(a**2 - b**2)**5) + b*(3*a**2 + b**2)/(5*d*(a + b*sin(c + d*x))**5*(a**2 - b**2)**3) + b/(d*(a + b*sin(c + d*x))**7*(7*a**2 - 7*b**2)) - log(1 - sin(c + d*x))/(2*d*(a + b)**8) + log(sin(c + d*x) + 1)/(2*d*(a - b)**8)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_464():
    f = sec(c + d*x)**3/(a + b*sin(c + d*x))**8
    F = 8*a*b**3*(15*a**6 + 63*a**4*b**2 + 45*a**2*b**4 + 5*b**6)*log(a + b*sin(c + d*x))/(d*(a**2 - b**2)**9) - a*b*(a**6 + 77*a**4*b**2 + 147*a**2*b**4 + 31*b**6)/(2*d*(a + b*sin(c + d*x))**2*(a**2 - b**2)**7) - a*b*(a**4 + 20*a**2*b**2 + 11*b**4)/(2*d*(a + b*sin(c + d*x))**4*(a**2 - b**2)**5) - a*b*(3*a**2 + 13*b**2)/(6*d*(a + b*sin(c + d*x))**6*(a**2 - b**2)**3) - b*(a**8 + 196*a**6*b**2 + 574*a**4*b**4 + 244*a**2*b**6 + 9*b**8)/(2*d*(a + b*sin(c + d*x))*(a**2 - b**2)**8) - b*(3*a**6 + 115*a**4*b**2 + 129*a**2*b**4 + 9*b**6)/(6*d*(a + b*sin(c + d*x))**3*(a**2 - b**2)**6) - b*(5*a**4 + 50*a**2*b**2 + 9*b**4)/(10*d*(a + b*sin(c + d*x))**5*(a**2 - b**2)**4) - b*(7*a**2 + 9*b**2)/(14*d*(a + b*sin(c + d*x))**7*(a**2 - b**2)**2) + (a - 9*b)*log(sin(c + d*x) + 1)/(4*d*(a - b)**9) - (-a*sin(c + d*x) + b)*sec(c + d*x)**2/(d*(a + b*sin(c + d*x))**7*(2*a**2 - 2*b**2)) - (a + 9*b)*log(1 - sin(c + d*x))/(4*d*(a + b)**9)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_465():
    f = cos(c + d*x)**8/(a + b*sin(c + d*x))**8
    F = a*cos(c + d*x)**7/(6*b*d*(a + b*sin(c + d*x))**6*(a**2 - b**2)) - a*(6*a**2 - 11*b**2)*cos(c + d*x)**5/(24*b**3*d*(a + b*sin(c + d*x))**4*(a**2 - b**2)**2) + a*(8*a**4 - 22*a**2*b**2 + 19*b**4)*cos(c + d*x)**3/(16*b**5*d*(a + b*sin(c + d*x))**2*(a**2 - b**2)**3) - a*(16*a**6 - 56*a**4*b**2 + 70*a**2*b**4 - 35*b**6)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(8*b**8*d*(a**2 - b**2)**(sympy.S(7)/2)) - cos(c + d*x)**7/(7*b*d*(a + b*sin(c + d*x))**7) + (6*a**2 + 5*a*b*sin(c + d*x) - 6*b**2)*cos(c + d*x)**5/(30*b**3*d*(a + b*sin(c + d*x))**5*(a**2 - b**2)) - (a*b*(6*a**2 - 11*b**2)*sin(c + d*x) + 8*(a**2 - b**2)**2)*cos(c + d*x)**3/(24*b**5*d*(a + b*sin(c + d*x))**3*(a**2 - b**2)**2) + (a*b*(8*a**4 - 22*a**2*b**2 + 19*b**4)*sin(c + d*x) + 16*(a**2 - b**2)**3)*cos(c + d*x)/(16*b**7*d*(a + b*sin(c + d*x))*(a**2 - b**2)**3) + x/b**8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_466():
    f = cos(c + d*x)**6/(a + b*sin(c + d*x))**8
    F = 5*a*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(8*d*(a**2 - b**2)**(sympy.S(9)/2)) + a*(8*a**4 - 30*a**2*b**2 + 57*b**4)*cos(c + d*x)/(336*b**5*d*(a + b*sin(c + d*x))**2*(a**2 - b**2)**3) + a*(4*a**2 - b**2)*cos(c + d*x)/(168*b**5*d*(a + b*sin(c + d*x))**4*(a**2 - b**2)) - cos(c + d*x)**5/(7*b*d*(a + b*sin(c + d*x))**7) + 5*(2*a + 3*b*sin(c + d*x))*cos(c + d*x)**3/(42*b**3*d*(a + b*sin(c + d*x))**6) + (8*a**6 - 38*a**4*b**2 + 87*a**2*b**4 + 48*b**6)*cos(c + d*x)/(336*b**5*d*(a + b*sin(c + d*x))*(a**2 - b**2)**4) + (4*a**4 - 9*a**2*b**2 + 12*b**4)*cos(c + d*x)/(168*b**5*d*(a + b*sin(c + d*x))**3*(a**2 - b**2)**2) - (4*a**2 + 10*a*b*sin(c + d*x) + 9*b**2)*cos(c + d*x)/(42*b**5*d*(a + b*sin(c + d*x))**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_467():
    f = cos(c + d*x)**4/(a + b*sin(c + d*x))**8
    F = 3*a*(2*a**2 + b**2)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(8*d*(a**2 - b**2)**(sympy.S(11)/2)) - a*(4*a**4 - 36*a**2*b**2 - 73*b**4)*cos(c + d*x)/(560*b**3*d*(a + b*sin(c + d*x))**2*(a**2 - b**2)**4) - a*(2*a**2 - 11*b**2)*cos(c + d*x)/(280*b**3*d*(a + b*sin(c + d*x))**4*(a**2 - b**2)**2) - cos(c + d*x)**3/(7*b*d*(a + b*sin(c + d*x))**7) - (4*a**6 - 40*a**4*b**2 - 247*a**2*b**4 - 32*b**6)*cos(c + d*x)/(560*b**3*d*(a + b*sin(c + d*x))*(a**2 - b**2)**5) - (2*a**4 - 15*a**2*b**2 - 8*b**4)*cos(c + d*x)/(280*b**3*d*(a + b*sin(c + d*x))**3*(a**2 - b**2)**3) - (a**2 - 3*b**2)*cos(c + d*x)/(140*b**3*d*(a + b*sin(c + d*x))**5*(a**2 - b**2)) + (a + 3*b*sin(c + d*x))*cos(c + d*x)/(28*b**3*d*(a + b*sin(c + d*x))**6)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_468():
    f = cos(c + d*x)**2/(a + b*sin(c + d*x))**8
    F = a*(8*a**4 + 20*a**2*b**2 + 5*b**4)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(8*d*(a**2 - b**2)**(sympy.S(13)/2)) + a*(40*a**4 + 718*a**2*b**2 + 397*b**4)*cos(c + d*x)/(1680*b*d*(a + b*sin(c + d*x))**2*(a**2 - b**2)**5) + a*(20*a**2 + 79*b**2)*cos(c + d*x)/(840*b*d*(a + b*sin(c + d*x))**4*(a**2 - b**2)**3) + a*cos(c + d*x)/(42*b*d*(a + b*sin(c + d*x))**6*(a**2 - b**2)) + (40*a**6 + 1518*a**4*b**2 + 1779*a**2*b**4 + 128*b**6)*cos(c + d*x)/(1680*b*d*(a + b*sin(c + d*x))*(a**2 - b**2)**6) + (20*a**4 + 179*a**2*b**2 + 32*b**4)*cos(c + d*x)/(840*b*d*(a + b*sin(c + d*x))**3*(a**2 - b**2)**4) + (5*a**2 + 6*b**2)*cos(c + d*x)/(210*b*d*(a + b*sin(c + d*x))**5*(a**2 - b**2)**2) - cos(c + d*x)/(7*b*d*(a + b*sin(c + d*x))**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_469():
    f = sec(c + d*x)**2/(a + b*sin(c + d*x))**8
    F = -9*a*b**2*(64*a**6 + 336*a**4*b**2 + 280*a**2*b**4 + 35*b**6)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(8*d*(a**2 - b**2)**(sympy.S(17)/2)) + 11*a*b*(280*a**4 + 844*a**2*b**2 + 241*b**4)*sec(c + d*x)/(560*d*(a + b*sin(c + d*x))**2*(a**2 - b**2)**6) + 13*a*b*(28*a**2 + 27*b**2)*sec(c + d*x)/(280*d*(a + b*sin(c + d*x))**4*(a**2 - b**2)**4) + 5*a*b*sec(c + d*x)/(14*d*(a + b*sin(c + d*x))**6*(a**2 - b**2)**2) + b*(9800*a**6 + 41484*a**4*b**2 + 22767*a**2*b**4 + 1024*b**6)*sec(c + d*x)/(560*d*(a + b*sin(c + d*x))*(a**2 - b**2)**7) + b*(700*a**4 + 1317*a**2*b**2 + 128*b**4)*sec(c + d*x)/(280*d*(a + b*sin(c + d*x))**3*(a**2 - b**2)**5) + b*(49*a**2 + 16*b**2)*sec(c + d*x)/(70*d*(a + b*sin(c + d*x))**5*(a**2 - b**2)**3) + b*sec(c + d*x)/(d*(a + b*sin(c + d*x))**7*(7*a**2 - 7*b**2)) - (315*a*b*(64*a**6 + 336*a**4*b**2 + 280*a**2*b**4 + 35*b**6) - (560*a**8 + 42472*a**6*b**2 + 125634*a**4*b**4 + 54511*a**2*b**6 + 2048*b**8)*sin(c + d*x))*sec(c + d*x)/(560*d*(a**2 - b**2)**8)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_470():
    f = sec(c + d*x)**4/(a + b*sin(c + d*x))**8
    F = 165*a*b**4*(32*a**6 + 112*a**4*b**2 + 70*a**2*b**4 + 7*b**6)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(8*d*(a**2 - b**2)**(sympy.S(19)/2)) + 13*a*b*(140*a**4 + 336*a**2*b**2 + 85*b**4)*sec(c + d*x)**3/(112*d*(a + b*sin(c + d*x))**2*(a**2 - b**2)**6) + a*b*(118*a**2 + 103*b**2)*sec(c + d*x)**3/(56*d*(a + b*sin(c + d*x))**4*(a**2 - b**2)**4) + 17*a*b*sec(c + d*x)**3/(42*d*(a + b*sin(c + d*x))**6*(a**2 - b**2)**2) + b*(9212*a**6 + 28420*a**4*b**2 + 12907*a**2*b**4 + 512*b**6)*sec(c + d*x)**3/(112*d*(a + b*sin(c + d*x))*(a**2 - b**2)**7) + b*(882*a**4 + 1421*a**2*b**2 + 128*b**4)*sec(c + d*x)**3/(168*d*(a + b*sin(c + d*x))**3*(a**2 - b**2)**5) + b*(13*a**2 + 4*b**2)*sec(c + d*x)**3/(14*d*(a + b*sin(c + d*x))**5*(a**2 - b**2)**3) + b*sec(c + d*x)**3/(d*(a + b*sin(c + d*x))**7*(7*a**2 - 7*b**2)) - (1155*a*b*(32*a**6 + 112*a**4*b**2 + 70*a**2*b**4 + 7*b**6) - (112*a**8 + 52528*a**6*b**2 + 142902*a**4*b**4 + 57665*a**2*b**6 + 2048*b**8)*sin(c + d*x))*sec(c + d*x)**3/(336*d*(a**2 - b**2)**8) + (3465*a*b**3*(32*a**6 + 112*a**4*b**2 + 70*a**2*b**4 + 7*b**6) + (224*a**10 - 6048*a**8*b**2 - 207332*a**6*b**4 - 413024*a**4*b**6 - 135489*a**2*b**8 - 4096*b**10)*sin(c + d*x))*sec(c + d*x)/(336*d*(a**2 - b**2)**9)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_471():
    f = sqrt(a + b*sin(c + d*x))*cos(c + d*x)**5
    F = -8*a*(a + b*sin(c + d*x))**(sympy.S(9)/2)/(9*b**5*d) - 8*a*(a + b*sin(c + d*x))**(sympy.S(5)/2)*(a**2 - b**2)/(5*b**5*d) + 2*(a + b*sin(c + d*x))**(sympy.S(11)/2)/(11*b**5*d) + (a + b*sin(c + d*x))**(sympy.S(7)/2)*(12*a**2 - 4*b**2)/(7*b**5*d) + 2*(a + b*sin(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)**2/(3*b**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_472():
    f = sqrt(a + b*sin(c + d*x))*cos(c + d*x)**3
    F = 4*a*(a + b*sin(c + d*x))**(sympy.S(5)/2)/(5*b**3*d) - 2*(a + b*sin(c + d*x))**(sympy.S(7)/2)/(7*b**3*d) + (a + b*sin(c + d*x))**(sympy.S(3)/2)*(-2*a**2 + 2*b**2)/(3*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_473():
    f = sqrt(a + b*sin(c + d*x))*cos(c + d*x)
    F = 2*(a + b*sin(c + d*x))**(sympy.S(3)/2)/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_474():
    f = sqrt(a + b*sin(c + d*x))*sec(c + d*x)
    F = -sqrt(a - b)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a - b))/d + sqrt(a + b)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a + b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_475():
    f = sqrt(a + b*sin(c + d*x))*sec(c + d*x)**3
    F = sqrt(a + b*sin(c + d*x))*tan(c + d*x)*sec(c + d*x)/(2*d) + (2*a + b)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a + b))/(4*d*sqrt(a + b)) - (2*a - b)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a - b))/(4*d*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_476():
    f = sqrt(a + b*sin(c + d*x))*sec(c + d*x)**5
    F = sqrt(a + b*sin(c + d*x))*tan(c + d*x)*sec(c + d*x)**3/(4*d) - sqrt(a + b*sin(c + d*x))*(a*b - (6*a**2 - 5*b**2)*sin(c + d*x))*sec(c + d*x)**2/(d*(16*a**2 - 16*b**2)) + (12*a**2 + 18*a*b + 5*b**2)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a + b))/(32*d*(a + b)**(sympy.S(3)/2)) - (12*a**2 - 18*a*b + 5*b**2)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a - b))/(32*d*(a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_477():
    f = sqrt(a + b*sin(c + d*x))*cos(c + d*x)**4
    F = -4*a*sqrt(a + b*sin(c + d*x))*cos(c + d*x)**3/(21*b*d) + 32*a*sqrt((a + b*sin(c + d*x))/(a + b))*(a**4 - 4*a**2*b**2 + 3*b**4)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(315*b**4*d*sqrt(a + b*sin(c + d*x))) + 2*(a + b*sin(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**3/(9*b*d) - 4*sqrt(a + b*sin(c + d*x))*(4*a*(a**2 - 3*b**2) - 3*b*(a**2 + 7*b**2)*sin(c + d*x))*cos(c + d*x)/(315*b**3*d) - sqrt(a + b*sin(c + d*x))*(32*a**4 - 120*a**2*b**2 - 168*b**4)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(315*b**4*d*sqrt((a + b*sin(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_478():
    f = sqrt(a + b*sin(c + d*x))*cos(c + d*x)**2
    F = -4*a*sqrt(a + b*sin(c + d*x))*cos(c + d*x)/(15*b*d) - 4*a*sqrt((a + b*sin(c + d*x))/(a + b))*(a**2 - b**2)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(15*b**2*d*sqrt(a + b*sin(c + d*x))) + 2*(a + b*sin(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)/(5*b*d) + sqrt(a + b*sin(c + d*x))*(4*a**2 + 12*b**2)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(15*b**2*d*sqrt((a + b*sin(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_479():
    f = sqrt(a + b*sin(c + d*x))*sec(c + d*x)**2
    F = a*sqrt((a + b*sin(c + d*x))/(a + b))*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(d*sqrt(a + b*sin(c + d*x))) + sqrt(a + b*sin(c + d*x))*tan(c + d*x)/d - sqrt(a + b*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(d*sqrt((a + b*sin(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_480():
    f = sqrt(a + b*sin(c + d*x))*sec(c + d*x)**4
    F = 2*a*sqrt((a + b*sin(c + d*x))/(a + b))*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(3*d*sqrt(a + b*sin(c + d*x))) + sqrt(a + b*sin(c + d*x))*tan(c + d*x)*sec(c + d*x)**2/(3*d) - sqrt(a + b*sin(c + d*x))*(a*b - (4*a**2 - 3*b**2)*sin(c + d*x))*sec(c + d*x)/(d*(6*a**2 - 6*b**2)) - sqrt(a + b*sin(c + d*x))*(4*a**2 - 3*b**2)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(d*sqrt((a + b*sin(c + d*x))/(a + b))*(6*a**2 - 6*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_481():
    f = (a + b*sin(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**5
    F = -8*a*(a + b*sin(c + d*x))**(sympy.S(11)/2)/(11*b**5*d) - 8*a*(a + b*sin(c + d*x))**(sympy.S(7)/2)*(a**2 - b**2)/(7*b**5*d) + 2*(a + b*sin(c + d*x))**(sympy.S(13)/2)/(13*b**5*d) + (a + b*sin(c + d*x))**(sympy.S(9)/2)*(12*a**2 - 4*b**2)/(9*b**5*d) + 2*(a + b*sin(c + d*x))**(sympy.S(5)/2)*(a**2 - b**2)**2/(5*b**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_482():
    f = (a + b*sin(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**3
    F = 4*a*(a + b*sin(c + d*x))**(sympy.S(7)/2)/(7*b**3*d) - 2*(a + b*sin(c + d*x))**(sympy.S(9)/2)/(9*b**3*d) + (a + b*sin(c + d*x))**(sympy.S(5)/2)*(-2*a**2 + 2*b**2)/(5*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_483():
    f = (a + b*sin(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)
    F = 2*(a + b*sin(c + d*x))**(sympy.S(5)/2)/(5*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_484():
    f = (a + b*sin(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)
    F = -2*b*sqrt(a + b*sin(c + d*x))/d - (a - b)**(sympy.S(3)/2)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a - b))/d + (a + b)**(sympy.S(3)/2)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a + b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_485():
    f = (a + b*sin(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**3
    F = -sqrt(a - b)*(2*a + b)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a - b))/(4*d) + sqrt(a + b)*(2*a - b)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a + b))/(4*d) + sqrt(a + b*sin(c + d*x))*(a*sin(c + d*x) + b)*sec(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_486():
    f = (a + b*sin(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**5
    F = -sqrt(a + b*sin(c + d*x))*(-6*a*sin(c + d*x) + b)*sec(c + d*x)**2/(16*d) + sqrt(a + b*sin(c + d*x))*(a*sin(c + d*x) + b)*sec(c + d*x)**4/(4*d) + (12*a**2 + 6*a*b - 3*b**2)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a + b))/(32*d*sqrt(a + b)) + (-12*a**2 + 6*a*b + 3*b**2)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a - b))/(32*d*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_487():
    f = (a + b*sin(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**4
    F = -32*a*sqrt(a + b*sin(c + d*x))*(a**4 - 6*a**2*b**2 - 27*b**4)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(1155*b**4*d*sqrt((a + b*sin(c + d*x))/(a + b))) - 2*b*sqrt(a + b*sin(c + d*x))*cos(c + d*x)**5/(11*d) + 2*sqrt(a + b*sin(c + d*x))*(a**2 + 28*a*b*sin(c + d*x) + 3*b**2)*cos(c + d*x)**3/(231*b*d) - 4*sqrt(a + b*sin(c + d*x))*(4*a**4 - 21*a**2*b**2 - 3*a*b*(a**2 + 31*b**2)*sin(c + d*x) - 15*b**4)*cos(c + d*x)/(1155*b**3*d) + sqrt((a + b*sin(c + d*x))/(a + b))*(32*a**6 - 200*a**4*b**2 + 48*a**2*b**4 + 120*b**6)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(1155*b**4*d*sqrt(a + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_488():
    f = (a + b*sin(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**2
    F = 4*a*sqrt(a + b*sin(c + d*x))*(3*a**2 + 29*b**2)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(105*b**2*d*sqrt((a + b*sin(c + d*x))/(a + b))) - 2*b*sqrt(a + b*sin(c + d*x))*cos(c + d*x)**3/(7*d) + 2*sqrt(a + b*sin(c + d*x))*(3*a**2 + 24*a*b*sin(c + d*x) + 5*b**2)*cos(c + d*x)/(105*b*d) - sqrt((a + b*sin(c + d*x))/(a + b))*(12*a**4 + 8*a**2*b**2 - 20*b**4)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(105*b**2*d*sqrt(a + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_489():
    f = (a + b*sin(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**2
    F = -a*sqrt(a + b*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(d*sqrt((a + b*sin(c + d*x))/(a + b))) + sqrt((a + b*sin(c + d*x))/(a + b))*(a**2 - b**2)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(d*sqrt(a + b*sin(c + d*x))) + sqrt(a + b*sin(c + d*x))*(a*sin(c + d*x) + b)*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_490():
    f = (a + b*sin(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**4
    F = -2*a*sqrt(a + b*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(3*d*sqrt((a + b*sin(c + d*x))/(a + b))) + sqrt((a + b*sin(c + d*x))/(a + b))*(4*a**2 - b**2)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(6*d*sqrt(a + b*sin(c + d*x))) - sqrt(a + b*sin(c + d*x))*(-4*a*sin(c + d*x) + b)*sec(c + d*x)/(6*d) + sqrt(a + b*sin(c + d*x))*(a*sin(c + d*x) + b)*sec(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_491():
    f = (a + b*sin(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**6
    F = -a*sqrt(a + b*sin(c + d*x))*(32*a**2 - 29*b**2)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(d*sqrt((a + b*sin(c + d*x))/(a + b))*(60*a**2 - 60*b**2)) + sqrt((a + b*sin(c + d*x))/(a + b))*(32*a**2 - 5*b**2)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(60*d*sqrt(a + b*sin(c + d*x))) - sqrt(a + b*sin(c + d*x))*(-8*a*sin(c + d*x) + b)*sec(c + d*x)**3/(30*d) + sqrt(a + b*sin(c + d*x))*(a*sin(c + d*x) + b)*sec(c + d*x)**5/(5*d) - sqrt(a + b*sin(c + d*x))*(-a*(32*a**4 - 61*a**2*b**2 + 29*b**4)*sin(c + d*x) + b*(8*a**4 - 13*a**2*b**2 + 5*b**4))*sec(c + d*x)/(60*d*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_492():
    f = (a + b*sin(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**5
    F = -8*a*(a + b*sin(c + d*x))**(sympy.S(13)/2)/(13*b**5*d) - 8*a*(a + b*sin(c + d*x))**(sympy.S(9)/2)*(a**2 - b**2)/(9*b**5*d) + 2*(a + b*sin(c + d*x))**(sympy.S(15)/2)/(15*b**5*d) + (a + b*sin(c + d*x))**(sympy.S(11)/2)*(12*a**2 - 4*b**2)/(11*b**5*d) + 2*(a + b*sin(c + d*x))**(sympy.S(7)/2)*(a**2 - b**2)**2/(7*b**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_493():
    f = (a + b*sin(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**3
    F = 4*a*(a + b*sin(c + d*x))**(sympy.S(9)/2)/(9*b**3*d) - 2*(a + b*sin(c + d*x))**(sympy.S(11)/2)/(11*b**3*d) + (a + b*sin(c + d*x))**(sympy.S(7)/2)*(-2*a**2 + 2*b**2)/(7*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_494():
    f = (a + b*sin(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)
    F = 2*(a + b*sin(c + d*x))**(sympy.S(7)/2)/(7*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_495():
    f = (a + b*sin(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)
    F = -4*a*b*sqrt(a + b*sin(c + d*x))/d - 2*b*(a + b*sin(c + d*x))**(sympy.S(3)/2)/(3*d) - (a - b)**(sympy.S(5)/2)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a - b))/d + (a + b)**(sympy.S(5)/2)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a + b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_496():
    f = (a + b*sin(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**3
    F = a*b*sqrt(a + b*sin(c + d*x))/(2*d) - (a - b)**(sympy.S(3)/2)*(2*a + 3*b)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a - b))/(4*d) + (a + b)**(sympy.S(3)/2)*(2*a - 3*b)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a + b))/(4*d) + (a + b*sin(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + b)*sec(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_497():
    f = (a + b*sin(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**5
    F = -3*sqrt(a - b)*(4*a**2 + 2*a*b - b**2)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a - b))/(32*d) + 3*sqrt(a + b)*(4*a**2 - 2*a*b - b**2)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a + b))/(32*d) + (a + b*sin(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + b)*sec(c + d*x)**4/(4*d) + 3*sqrt(a + b*sin(c + d*x))*(a*b + (2*a**2 - b**2)*sin(c + d*x))*sec(c + d*x)**2/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_498():
    f = (a + b*sin(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**4
    F = -32*a*b*sqrt(a + b*sin(c + d*x))*cos(c + d*x)**5/(143*d) + 32*a*sqrt((a + b*sin(c + d*x))/(a + b))*(5*a**6 - 45*a**4*b**2 - 53*a**2*b**4 + 93*b**6)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(15015*b**4*d*sqrt(a + b*sin(c + d*x))) - 2*b*(a + b*sin(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**5/(13*d) + 2*sqrt(a + b*sin(c + d*x))*(a*(5*a**2 + 59*b**2) + 7*b*(53*a**2 + 11*b**2)*sin(c + d*x))*cos(c + d*x)**3/(3003*b*d) - 4*sqrt(a + b*sin(c + d*x))*(4*a*(5*a**4 - 40*a**2*b**2 - 93*b**4) - 3*b*(5*a**4 + 430*a**2*b**2 + 77*b**4)*sin(c + d*x))*cos(c + d*x)/(15015*b**3*d) - sqrt(a + b*sin(c + d*x))*(160*a**6 - 1400*a**4*b**2 - 13296*a**2*b**4 - 1848*b**6)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(15015*b**4*d*sqrt((a + b*sin(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_499():
    f = (a + b*sin(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**2
    F = -8*a*b*sqrt(a + b*sin(c + d*x))*cos(c + d*x)**3/(21*d) - 4*a*sqrt((a + b*sin(c + d*x))/(a + b))*(5*a**4 + 22*a**2*b**2 - 27*b**4)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(315*b**2*d*sqrt(a + b*sin(c + d*x))) - 2*b*(a + b*sin(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**3/(9*d) + 2*sqrt(a + b*sin(c + d*x))*(a*(5*a**2 + 27*b**2) + 3*b*(25*a**2 + 7*b**2)*sin(c + d*x))*cos(c + d*x)/(315*b*d) + sqrt(a + b*sin(c + d*x))*(20*a**4 + 408*a**2*b**2 + 84*b**4)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(315*b**2*d*sqrt((a + b*sin(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_500():
    f = (a + b*sin(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**2
    F = a*b*sqrt(a + b*sin(c + d*x))*cos(c + d*x)/d + a*sqrt((a + b*sin(c + d*x))/(a + b))*(a**2 - b**2)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(d*sqrt(a + b*sin(c + d*x))) + (a + b*sin(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + b)*sec(c + d*x)/d - sqrt(a + b*sin(c + d*x))*(a**2 + 3*b**2)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(d*sqrt((a + b*sin(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_501():
    f = (a + b*sin(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**4
    F = 2*a*sqrt((a + b*sin(c + d*x))/(a + b))*(a**2 - b**2)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(3*d*sqrt(a + b*sin(c + d*x))) + (a + b*sin(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + b)*sec(c + d*x)**3/(3*d) + sqrt(a + b*sin(c + d*x))*(a*b + (4*a**2 - 3*b**2)*sin(c + d*x))*sec(c + d*x)/(6*d) - sqrt(a + b*sin(c + d*x))*(4*a**2 - 3*b**2)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(6*d*sqrt((a + b*sin(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_502():
    f = (a + b*sin(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**6
    F = a*sqrt((a + b*sin(c + d*x))/(a + b))*(32*a**2 - 17*b**2)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(60*d*sqrt(a + b*sin(c + d*x))) + (a + b*sin(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + b)*sec(c + d*x)**5/(5*d) + sqrt(a + b*sin(c + d*x))*(5*a*b + (8*a**2 - 3*b**2)*sin(c + d*x))*sec(c + d*x)**3/(30*d) - sqrt(a + b*sin(c + d*x))*(8*a*b*(a**2 - b**2) - (32*a**4 - 41*a**2*b**2 + 9*b**4)*sin(c + d*x))*sec(c + d*x)/(d*(60*a**2 - 60*b**2)) - sqrt(a + b*sin(c + d*x))*(32*a**2 - 9*b**2)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(60*d*sqrt((a + b*sin(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_503():
    f = (a + b*sin(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**8
    F = 2*a*sqrt((a + b*sin(c + d*x))/(a + b))*(8*a**2 - 3*b**2)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(35*d*sqrt(a + b*sin(c + d*x))) + (a + b*sin(c + d*x))**(sympy.S(3)/2)*(a*sin(c + d*x) + b)*sec(c + d*x)**7/(7*d) + 3*sqrt(a + b*sin(c + d*x))*(3*a*b + (4*a**2 - b**2)*sin(c + d*x))*sec(c + d*x)**5/(70*d) - sqrt(a + b*sin(c + d*x))*(4*a*b*(a**2 - b**2) - (32*a**4 - 39*a**2*b**2 + 7*b**4)*sin(c + d*x))*sec(c + d*x)**3/(d*(140*a**2 - 140*b**2)) - sqrt(a + b*sin(c + d*x))*(a*b*(32*a**4 - 59*a**2*b**2 + 27*b**4) - (128*a**6 - 272*a**4*b**2 + 165*a**2*b**4 - 21*b**6)*sin(c + d*x))*sec(c + d*x)/(280*d*(a**2 - b**2)**2) - sqrt(a + b*sin(c + d*x))*(128*a**4 - 144*a**2*b**2 + 21*b**4)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(d*sqrt((a + b*sin(c + d*x))/(a + b))*(280*a**2 - 280*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_504():
    f = cos(c + d*x)**5/sqrt(a + b*sin(c + d*x))
    F = -8*a*(a + b*sin(c + d*x))**(sympy.S(7)/2)/(7*b**5*d) - 8*a*(a + b*sin(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)/(3*b**5*d) + 2*(a + b*sin(c + d*x))**(sympy.S(9)/2)/(9*b**5*d) + (a + b*sin(c + d*x))**(sympy.S(5)/2)*(12*a**2 - 4*b**2)/(5*b**5*d) + 2*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)**2/(b**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_505():
    f = cos(c + d*x)**3/sqrt(a + b*sin(c + d*x))
    F = 4*a*(a + b*sin(c + d*x))**(sympy.S(3)/2)/(3*b**3*d) - 2*(a + b*sin(c + d*x))**(sympy.S(5)/2)/(5*b**3*d) + sqrt(a + b*sin(c + d*x))*(-2*a**2 + 2*b**2)/(b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_506():
    f = cos(c + d*x)/sqrt(a + b*sin(c + d*x))
    F = 2*sqrt(a + b*sin(c + d*x))/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_507():
    f = sec(c + d*x)/sqrt(a + b*sin(c + d*x))
    F = atanh(sqrt(a + b*sin(c + d*x))/sqrt(a + b))/(d*sqrt(a + b)) - atanh(sqrt(a + b*sin(c + d*x))/sqrt(a - b))/(d*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_508():
    f = sec(c + d*x)**3/sqrt(a + b*sin(c + d*x))
    F = -sqrt(a + b*sin(c + d*x))*(-a*sin(c + d*x) + b)*sec(c + d*x)**2/(d*(2*a**2 - 2*b**2)) + (2*a + 3*b)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a + b))/(4*d*(a + b)**(sympy.S(3)/2)) - (2*a - 3*b)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a - b))/(4*d*(a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_509():
    f = sec(c + d*x)**5/sqrt(a + b*sin(c + d*x))
    F = -sqrt(a + b*sin(c + d*x))*(-a*sin(c + d*x) + b)*sec(c + d*x)**4/(d*(4*a**2 - 4*b**2)) - sqrt(a + b*sin(c + d*x))*(-6*a*(a**2 - 2*b**2)*sin(c + d*x) + b*(a**2 - 7*b**2))*sec(c + d*x)**2/(16*d*(a**2 - b**2)**2) + (12*a**2 + 30*a*b + 21*b**2)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a + b))/(32*d*(a + b)**(sympy.S(5)/2)) + (-12*a**2 + 30*a*b - 21*b**2)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a - b))/(32*d*(a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_510():
    f = cos(c + d*x)**4/sqrt(a + b*sin(c + d*x))
    F = -32*a*sqrt(a + b*sin(c + d*x))*(a**2 - 2*b**2)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(35*b**4*d*sqrt((a + b*sin(c + d*x))/(a + b))) + 2*sqrt(a + b*sin(c + d*x))*cos(c + d*x)**3/(7*b*d) - 4*sqrt(a + b*sin(c + d*x))*(4*a**2 - 3*a*b*sin(c + d*x) - 5*b**2)*cos(c + d*x)/(35*b**3*d) + sqrt((a + b*sin(c + d*x))/(a + b))*(32*a**4 - 72*a**2*b**2 + 40*b**4)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(35*b**4*d*sqrt(a + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_511():
    f = cos(c + d*x)**2/sqrt(a + b*sin(c + d*x))
    F = 4*a*sqrt(a + b*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(3*b**2*d*sqrt((a + b*sin(c + d*x))/(a + b))) + 2*sqrt(a + b*sin(c + d*x))*cos(c + d*x)/(3*b*d) - sqrt((a + b*sin(c + d*x))/(a + b))*(4*a**2 - 4*b**2)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(3*b**2*d*sqrt(a + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_512():
    f = sec(c + d*x)**2/sqrt(a + b*sin(c + d*x))
    F = -a*sqrt(a + b*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(d*sqrt((a + b*sin(c + d*x))/(a + b))*(a**2 - b**2)) + sqrt((a + b*sin(c + d*x))/(a + b))*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(d*sqrt(a + b*sin(c + d*x))) - sqrt(a + b*sin(c + d*x))*(-a*sin(c + d*x) + b)*sec(c + d*x)/(d*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_513():
    f = sec(c + d*x)**4/sqrt(a + b*sin(c + d*x))
    F = -2*a*sqrt(a + b*sin(c + d*x))*(a**2 - 2*b**2)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(3*d*sqrt((a + b*sin(c + d*x))/(a + b))*(a**2 - b**2)**2) + sqrt((a + b*sin(c + d*x))/(a + b))*(4*a**2 - 5*b**2)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(d*sqrt(a + b*sin(c + d*x))*(6*a**2 - 6*b**2)) - sqrt(a + b*sin(c + d*x))*(-a*sin(c + d*x) + b)*sec(c + d*x)**3/(d*(3*a**2 - 3*b**2)) - sqrt(a + b*sin(c + d*x))*(-4*a*(a**2 - 2*b**2)*sin(c + d*x) + b*(a**2 - 5*b**2))*sec(c + d*x)/(6*d*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_514():
    f = cos(c + d*x)**5/(a + b*sin(c + d*x))**(sympy.S(3)/2)
    F = -8*a*(a + b*sin(c + d*x))**(sympy.S(5)/2)/(5*b**5*d) - 8*a*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)/(b**5*d) + 2*(a + b*sin(c + d*x))**(sympy.S(7)/2)/(7*b**5*d) + (a + b*sin(c + d*x))**(sympy.S(3)/2)*(12*a**2 - 4*b**2)/(3*b**5*d) - 2*(a**2 - b**2)**2/(b**5*d*sqrt(a + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_515():
    f = cos(c + d*x)**3/(a + b*sin(c + d*x))**(sympy.S(3)/2)
    F = 4*a*sqrt(a + b*sin(c + d*x))/(b**3*d) - 2*(a + b*sin(c + d*x))**(sympy.S(3)/2)/(3*b**3*d) + (2*a**2 - 2*b**2)/(b**3*d*sqrt(a + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_516():
    f = cos(c + d*x)/(a + b*sin(c + d*x))**(sympy.S(3)/2)
    F = -2/(b*d*sqrt(a + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_517():
    f = sec(c + d*x)/(a + b*sin(c + d*x))**(sympy.S(3)/2)
    F = 2*b/(d*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)) + atanh(sqrt(a + b*sin(c + d*x))/sqrt(a + b))/(d*(a + b)**(sympy.S(3)/2)) - atanh(sqrt(a + b*sin(c + d*x))/sqrt(a - b))/(d*(a - b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_518():
    f = sec(c + d*x)**3/(a + b*sin(c + d*x))**(sympy.S(3)/2)
    F = -b*(a**2 + 5*b**2)/(2*d*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)**2) - (-a*sin(c + d*x) + b)*sec(c + d*x)**2/(d*sqrt(a + b*sin(c + d*x))*(2*a**2 - 2*b**2)) + (2*a + 5*b)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a + b))/(4*d*(a + b)**(sympy.S(5)/2)) - (2*a - 5*b)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a - b))/(4*d*(a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_519():
    f = sec(c + d*x)**5/(a + b*sin(c + d*x))**(sympy.S(3)/2)
    F = -3*b*(2*a**4 - 7*a**2*b**2 - 15*b**4)/(16*d*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)**3) - (-a*sin(c + d*x) + b)*sec(c + d*x)**4/(d*sqrt(a + b*sin(c + d*x))*(4*a**2 - 4*b**2)) + (2*a*(3*a**2 - 8*b**2)*sin(c + d*x) + b*(a**2 + 9*b**2))*sec(c + d*x)**2/(16*d*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)**2) + (12*a**2 + 42*a*b + 45*b**2)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a + b))/(32*d*(a + b)**(sympy.S(7)/2)) + (-12*a**2 + 42*a*b - 45*b**2)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a - b))/(32*d*(a - b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_520():
    f = cos(c + d*x)**6/(a + b*sin(c + d*x))**(sympy.S(3)/2)
    F = 16*a*sqrt((a + b*sin(c + d*x))/(a + b))*(32*a**4 - 65*a**2*b**2 + 33*b**4)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(63*b**6*d*sqrt(a + b*sin(c + d*x))) - 2*cos(c + d*x)**5/(b*d*sqrt(a + b*sin(c + d*x))) + 20*sqrt(a + b*sin(c + d*x))*(8*a - 7*b*sin(c + d*x))*cos(c + d*x)**3/(63*b**3*d) - 8*sqrt(a + b*sin(c + d*x))*(a*(32*a**2 - 33*b**2) - 3*b*(8*a**2 - 7*b**2)*sin(c + d*x))*cos(c + d*x)/(63*b**5*d) - sqrt(a + b*sin(c + d*x))*(512*a**4 - 912*a**2*b**2 + 336*b**4)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(63*b**6*d*sqrt((a + b*sin(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_521():
    f = cos(c + d*x)**4/(a + b*sin(c + d*x))**(sympy.S(3)/2)
    F = -32*a*sqrt((a + b*sin(c + d*x))/(a + b))*(a**2 - b**2)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(5*b**4*d*sqrt(a + b*sin(c + d*x))) - 2*cos(c + d*x)**3/(b*d*sqrt(a + b*sin(c + d*x))) + 4*sqrt(a + b*sin(c + d*x))*(4*a - 3*b*sin(c + d*x))*cos(c + d*x)/(5*b**3*d) + sqrt(a + b*sin(c + d*x))*(32*a**2 - 24*b**2)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(5*b**4*d*sqrt((a + b*sin(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_522():
    f = cos(c + d*x)**2/(a + b*sin(c + d*x))**(sympy.S(3)/2)
    F = 4*a*sqrt((a + b*sin(c + d*x))/(a + b))*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(b**2*d*sqrt(a + b*sin(c + d*x))) - 2*cos(c + d*x)/(b*d*sqrt(a + b*sin(c + d*x))) - 4*sqrt(a + b*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(b**2*d*sqrt((a + b*sin(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_523():
    f = sec(c + d*x)**2/(a + b*sin(c + d*x))**(sympy.S(3)/2)
    F = a*sqrt((a + b*sin(c + d*x))/(a + b))*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(d*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)) + 2*b*sec(c + d*x)/(d*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)) - sqrt(a + b*sin(c + d*x))*(4*a*b - (a**2 + 3*b**2)*sin(c + d*x))*sec(c + d*x)/(d*(a**2 - b**2)**2) - sqrt(a + b*sin(c + d*x))*(a**2 + 3*b**2)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(d*sqrt((a + b*sin(c + d*x))/(a + b))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_524():
    f = sec(c + d*x)**4/(a + b*sin(c + d*x))**(sympy.S(3)/2)
    F = 2*a*sqrt((a + b*sin(c + d*x))/(a + b))*(a**2 - 3*b**2)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(3*d*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)**2) + 2*b*sec(c + d*x)**3/(d*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)) - sqrt(a + b*sin(c + d*x))*(8*a*b - (a**2 + 7*b**2)*sin(c + d*x))*sec(c + d*x)**3/(3*d*(a**2 - b**2)**2) - sqrt(a + b*sin(c + d*x))*(a*b*(a**2 - 33*b**2) - (4*a**4 - 15*a**2*b**2 - 21*b**4)*sin(c + d*x))*sec(c + d*x)/(6*d*(a**2 - b**2)**3) - sqrt(a + b*sin(c + d*x))*(4*a**4 - 15*a**2*b**2 - 21*b**4)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(6*d*sqrt((a + b*sin(c + d*x))/(a + b))*(a**2 - b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_525():
    f = cos(c + d*x)**5/(a + b*sin(c + d*x))**(sympy.S(5)/2)
    F = -8*a*(a + b*sin(c + d*x))**(sympy.S(3)/2)/(3*b**5*d) + 8*a*(a**2 - b**2)/(b**5*d*sqrt(a + b*sin(c + d*x))) + 2*(a + b*sin(c + d*x))**(sympy.S(5)/2)/(5*b**5*d) + sqrt(a + b*sin(c + d*x))*(12*a**2 - 4*b**2)/(b**5*d) - 2*(a**2 - b**2)**2/(3*b**5*d*(a + b*sin(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_526():
    f = cos(c + d*x)**3/(a + b*sin(c + d*x))**(sympy.S(5)/2)
    F = -4*a/(b**3*d*sqrt(a + b*sin(c + d*x))) - 2*sqrt(a + b*sin(c + d*x))/(b**3*d) + (2*a**2 - 2*b**2)/(3*b**3*d*(a + b*sin(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_527():
    f = cos(c + d*x)/(a + b*sin(c + d*x))**(sympy.S(5)/2)
    F = -2/(3*b*d*(a + b*sin(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_528():
    f = sec(c + d*x)/(a + b*sin(c + d*x))**(sympy.S(5)/2)
    F = 4*a*b/(d*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)**2) + 2*b/(d*(a + b*sin(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)) + atanh(sqrt(a + b*sin(c + d*x))/sqrt(a + b))/(d*(a + b)**(sympy.S(5)/2)) - atanh(sqrt(a + b*sin(c + d*x))/sqrt(a - b))/(d*(a - b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_529():
    f = sec(c + d*x)**3/(a + b*sin(c + d*x))**(sympy.S(5)/2)
    F = -a*b*(a**2 + 19*b**2)/(2*d*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)**3) - b*(3*a**2 + 7*b**2)/(6*d*(a + b*sin(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)**2) - (-a*sin(c + d*x) + b)*sec(c + d*x)**2/(d*(a + b*sin(c + d*x))**(sympy.S(3)/2)*(2*a**2 - 2*b**2)) + (2*a + 7*b)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a + b))/(4*d*(a + b)**(sympy.S(7)/2)) - (2*a - 7*b)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a - b))/(4*d*(a - b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_530():
    f = sec(c + d*x)**5/(a + b*sin(c + d*x))**(sympy.S(5)/2)
    F = -a*b*(3*a**4 - 16*a**2*b**2 - 127*b**4)/(8*d*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)**4) - b*(18*a**4 - 81*a**2*b**2 - 77*b**4)/(48*d*(a + b*sin(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)**3) - (-a*sin(c + d*x) + b)*sec(c + d*x)**4/(d*(a + b*sin(c + d*x))**(sympy.S(3)/2)*(4*a**2 - 4*b**2)) + (2*a*(3*a**2 - 10*b**2)*sin(c + d*x) + b*(3*a**2 + 11*b**2))*sec(c + d*x)**2/(16*d*(a + b*sin(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)**2) + (12*a**2 + 54*a*b + 77*b**2)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a + b))/(32*d*(a + b)**(sympy.S(9)/2)) - (12*a**2 - 54*a*b + 77*b**2)*atanh(sqrt(a + b*sin(c + d*x))/sqrt(a - b))/(32*d*(a - b)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_531():
    f = cos(c + d*x)**8/(a + b*sin(c + d*x))**(sympy.S(5)/2)
    F = -128*a*sqrt(a + b*sin(c + d*x))*(4*a**2 - 3*b**2)*(8*a**2 - 9*b**2)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(99*b**8*d*sqrt((a + b*sin(c + d*x))/(a + b))) - 2*cos(c + d*x)**7/(3*b*d*(a + b*sin(c + d*x))**(sympy.S(3)/2)) - 28*(12*a + b*sin(c + d*x))*cos(c + d*x)**5/(33*b**3*d*sqrt(a + b*sin(c + d*x))) + 40*sqrt(a + b*sin(c + d*x))*(32*a**2 - 28*a*b*sin(c + d*x) - 3*b**2)*cos(c + d*x)**3/(99*b**5*d) - 16*sqrt(a + b*sin(c + d*x))*(128*a**4 - 144*a**2*b**2 - 3*a*b*(32*a**2 - 31*b**2)*sin(c + d*x) + 15*b**4)*cos(c + d*x)/(99*b**7*d) + sqrt((a + b*sin(c + d*x))/(a + b))*(4096*a**6 - 8704*a**4*b**2 + 5088*a**2*b**4 - 480*b**6)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(99*b**8*d*sqrt(a + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_532():
    f = cos(c + d*x)**6/(a + b*sin(c + d*x))**(sympy.S(5)/2)
    F = 16*a*sqrt(a + b*sin(c + d*x))*(32*a**2 - 29*b**2)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(21*b**6*d*sqrt((a + b*sin(c + d*x))/(a + b))) - 2*cos(c + d*x)**5/(3*b*d*(a + b*sin(c + d*x))**(sympy.S(3)/2)) - 20*(8*a + b*sin(c + d*x))*cos(c + d*x)**3/(21*b**3*d*sqrt(a + b*sin(c + d*x))) + 8*sqrt(a + b*sin(c + d*x))*(32*a**2 - 24*a*b*sin(c + d*x) - 5*b**2)*cos(c + d*x)/(21*b**5*d) - sqrt((a + b*sin(c + d*x))/(a + b))*(512*a**4 - 592*a**2*b**2 + 80*b**4)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(21*b**6*d*sqrt(a + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_533():
    f = cos(c + d*x)**4/(a + b*sin(c + d*x))**(sympy.S(5)/2)
    F = -32*a*sqrt(a + b*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(3*b**4*d*sqrt((a + b*sin(c + d*x))/(a + b))) - 2*cos(c + d*x)**3/(3*b*d*(a + b*sin(c + d*x))**(sympy.S(3)/2)) - 4*(4*a + b*sin(c + d*x))*cos(c + d*x)/(3*b**3*d*sqrt(a + b*sin(c + d*x))) + sqrt((a + b*sin(c + d*x))/(a + b))*(32*a**2 - 8*b**2)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(3*b**4*d*sqrt(a + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_534():
    f = cos(c + d*x)**2/(a + b*sin(c + d*x))**(sympy.S(5)/2)
    F = 4*a*cos(c + d*x)/(3*b*d*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)) + 4*a*sqrt(a + b*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(3*b**2*d*sqrt((a + b*sin(c + d*x))/(a + b))*(a**2 - b**2)) - 2*cos(c + d*x)/(3*b*d*(a + b*sin(c + d*x))**(sympy.S(3)/2)) - 4*sqrt((a + b*sin(c + d*x))/(a + b))*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(3*b**2*d*sqrt(a + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_535():
    f = sec(c + d*x)**2/(a + b*sin(c + d*x))**(sympy.S(5)/2)
    F = 16*a*b*sec(c + d*x)/(3*d*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)**2) - a*sqrt(a + b*sin(c + d*x))*(3*a**2 + 29*b**2)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(3*d*sqrt((a + b*sin(c + d*x))/(a + b))*(a**2 - b**2)**3) + 2*b*sec(c + d*x)/(d*(a + b*sin(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)) + sqrt((a + b*sin(c + d*x))/(a + b))*(3*a**2 + 5*b**2)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(3*d*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)**2) - sqrt(a + b*sin(c + d*x))*(-a*(3*a**2 + 29*b**2)*sin(c + d*x) + b*(27*a**2 + 5*b**2))*sec(c + d*x)/(3*d*(a**2 - b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_536():
    f = sec(c + d*x)**4/(a + b*sin(c + d*x))**(sympy.S(5)/2)
    F = 8*a*b*sec(c + d*x)**3/(d*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)**2) - 2*a*sqrt(a + b*sin(c + d*x))*(a**4 - 6*a**2*b**2 - 27*b**4)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(3*d*sqrt((a + b*sin(c + d*x))/(a + b))*(a**2 - b**2)**4) + 2*b*sec(c + d*x)**3/(d*(a + b*sin(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)) + sqrt((a + b*sin(c + d*x))/(a + b))*(4*a**4 - 21*a**2*b**2 - 15*b**4)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(6*d*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)**3) - sqrt(a + b*sin(c + d*x))*(-a*(a**2 + 31*b**2)*sin(c + d*x) + b*(29*a**2 + 3*b**2))*sec(c + d*x)**3/(3*d*(a**2 - b**2)**3) - sqrt(a + b*sin(c + d*x))*(-4*a*(a**4 - 6*a**2*b**2 - 27*b**4)*sin(c + d*x) + b*(a**4 - 114*a**2*b**2 - 15*b**4))*sec(c + d*x)/(6*d*(a**2 - b**2)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_537():
    f = (e*cos(c + d*x))**(sympy.S(7)/2)*(a + b*sin(c + d*x))
    F = 10*a*e**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(e*cos(c + d*x))) + 10*a*e**3*sqrt(e*cos(c + d*x))*sin(c + d*x)/(21*d) + 2*a*e*(e*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*d) - 2*b*(e*cos(c + d*x))**(sympy.S(9)/2)/(9*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_538():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)*(a + b*sin(c + d*x))
    F = 6*a*e**2*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) + 2*a*e*(e*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*d) - 2*b*(e*cos(c + d*x))**(sympy.S(7)/2)/(7*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_539():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)*(a + b*sin(c + d*x))
    F = 2*a*e**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(e*cos(c + d*x))) + 2*a*e*sqrt(e*cos(c + d*x))*sin(c + d*x)/(3*d) - 2*b*(e*cos(c + d*x))**(sympy.S(5)/2)/(5*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_540():
    f = sqrt(e*cos(c + d*x))*(a + b*sin(c + d*x))
    F = 2*a*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x))) - 2*b*(e*cos(c + d*x))**(sympy.S(3)/2)/(3*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_541():
    f = (a + b*sin(c + d*x))/sqrt(e*cos(c + d*x))
    F = 2*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*sqrt(e*cos(c + d*x))) - 2*b*sqrt(e*cos(c + d*x))/(d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_542():
    f = (a + b*sin(c + d*x))/(e*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*a*sin(c + d*x)/(d*e*sqrt(e*cos(c + d*x))) - 2*a*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*e**2*sqrt(cos(c + d*x))) + 2*b/(d*e*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_543():
    f = (a + b*sin(c + d*x))/(e*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*a*sin(c + d*x)/(3*d*e*(e*cos(c + d*x))**(sympy.S(3)/2)) + 2*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*e**2*sqrt(e*cos(c + d*x))) + 2*b/(3*d*e*(e*cos(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_544():
    f = (a + b*sin(c + d*x))/(e*cos(c + d*x))**(sympy.S(7)/2)
    F = 2*a*sin(c + d*x)/(5*d*e*(e*cos(c + d*x))**(sympy.S(5)/2)) + 6*a*sin(c + d*x)/(5*d*e**3*sqrt(e*cos(c + d*x))) - 6*a*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*e**4*sqrt(cos(c + d*x))) + 2*b/(5*d*e*(e*cos(c + d*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_545():
    f = (e*cos(c + d*x))**(sympy.S(7)/2)*(a + b*sin(c + d*x))**2
    F = -26*a*b*(e*cos(c + d*x))**(sympy.S(9)/2)/(99*d*e) - 2*b*(e*cos(c + d*x))**(sympy.S(9)/2)*(a + b*sin(c + d*x))/(11*d*e) + e**4*(110*a**2 + 20*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(231*d*sqrt(e*cos(c + d*x))) + e**3*sqrt(e*cos(c + d*x))*(110*a**2 + 20*b**2)*sin(c + d*x)/(231*d) + e*(e*cos(c + d*x))**(sympy.S(5)/2)*(22*a**2 + 4*b**2)*sin(c + d*x)/(77*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_546():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)*(a + b*sin(c + d*x))**2
    F = -22*a*b*(e*cos(c + d*x))**(sympy.S(7)/2)/(63*d*e) - 2*b*(e*cos(c + d*x))**(sympy.S(7)/2)*(a + b*sin(c + d*x))/(9*d*e) + e**2*sqrt(e*cos(c + d*x))*(18*a**2 + 4*b**2)*elliptic_e(c/2 + d*x/2, 2)/(15*d*sqrt(cos(c + d*x))) + e*(e*cos(c + d*x))**(sympy.S(3)/2)*(18*a**2 + 4*b**2)*sin(c + d*x)/(45*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_547():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)*(a + b*sin(c + d*x))**2
    F = -18*a*b*(e*cos(c + d*x))**(sympy.S(5)/2)/(35*d*e) - 2*b*(e*cos(c + d*x))**(sympy.S(5)/2)*(a + b*sin(c + d*x))/(7*d*e) + e**2*(14*a**2 + 4*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(e*cos(c + d*x))) + e*sqrt(e*cos(c + d*x))*(14*a**2 + 4*b**2)*sin(c + d*x)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_548():
    f = sqrt(e*cos(c + d*x))*(a + b*sin(c + d*x))**2
    F = -14*a*b*(e*cos(c + d*x))**(sympy.S(3)/2)/(15*d*e) - 2*b*(e*cos(c + d*x))**(sympy.S(3)/2)*(a + b*sin(c + d*x))/(5*d*e) + sqrt(e*cos(c + d*x))*(10*a**2 + 4*b**2)*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_549():
    f = (a + b*sin(c + d*x))**2/sqrt(e*cos(c + d*x))
    F = -10*a*b*sqrt(e*cos(c + d*x))/(3*d*e) - 2*b*sqrt(e*cos(c + d*x))*(a + b*sin(c + d*x))/(3*d*e) + (6*a**2 + 4*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_550():
    f = (a + b*sin(c + d*x))**2/(e*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*a*b*(e*cos(c + d*x))**(sympy.S(3)/2)/(d*e**3) + (a + b*sin(c + d*x))*(2*a*sin(c + d*x) + 2*b)/(d*e*sqrt(e*cos(c + d*x))) - sqrt(e*cos(c + d*x))*(2*a**2 + 4*b**2)*elliptic_e(c/2 + d*x/2, 2)/(d*e**2*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_551():
    f = (a + b*sin(c + d*x))**2/(e*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*a*b*sqrt(e*cos(c + d*x))/(3*d*e**3) + (a + b*sin(c + d*x))*(2*a*sin(c + d*x) + 2*b)/(3*d*e*(e*cos(c + d*x))**(sympy.S(3)/2)) + (2*a**2 - 4*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*e**2*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_552():
    f = (a + b*sin(c + d*x))**2/(e*cos(c + d*x))**(sympy.S(7)/2)
    F = 2*a*b/(5*d*e**3*sqrt(e*cos(c + d*x))) + (a + b*sin(c + d*x))*(2*a*sin(c + d*x) + 2*b)/(5*d*e*(e*cos(c + d*x))**(sympy.S(5)/2)) + (6*a**2 - 4*b**2)*sin(c + d*x)/(5*d*e**3*sqrt(e*cos(c + d*x))) - sqrt(e*cos(c + d*x))*(6*a**2 - 4*b**2)*elliptic_e(c/2 + d*x/2, 2)/(5*d*e**4*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_553():
    f = (e*cos(c + d*x))**(sympy.S(7)/2)*(a + b*sin(c + d*x))**3
    F = -34*a*b*(e*cos(c + d*x))**(sympy.S(9)/2)*(a + b*sin(c + d*x))/(143*d*e) + 10*a*e**4*(11*a**2 + 6*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(231*d*sqrt(e*cos(c + d*x))) + 10*a*e**3*sqrt(e*cos(c + d*x))*(11*a**2 + 6*b**2)*sin(c + d*x)/(231*d) + 2*a*e*(e*cos(c + d*x))**(sympy.S(5)/2)*(11*a**2 + 6*b**2)*sin(c + d*x)/(77*d) - 2*b*(e*cos(c + d*x))**(sympy.S(9)/2)*(a + b*sin(c + d*x))**2/(13*d*e) - 2*b*(e*cos(c + d*x))**(sympy.S(9)/2)*(177*a**2 + 44*b**2)/(1287*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_554():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)*(a + b*sin(c + d*x))**3
    F = -10*a*b*(e*cos(c + d*x))**(sympy.S(7)/2)*(a + b*sin(c + d*x))/(33*d*e) + 2*a*e**2*sqrt(e*cos(c + d*x))*(3*a**2 + 2*b**2)*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) + 2*a*e*(e*cos(c + d*x))**(sympy.S(3)/2)*(3*a**2 + 2*b**2)*sin(c + d*x)/(15*d) - 2*b*(e*cos(c + d*x))**(sympy.S(7)/2)*(a + b*sin(c + d*x))**2/(11*d*e) - 2*b*(e*cos(c + d*x))**(sympy.S(7)/2)*(43*a**2 + 12*b**2)/(231*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_555():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)*(a + b*sin(c + d*x))**3
    F = -26*a*b*(e*cos(c + d*x))**(sympy.S(5)/2)*(a + b*sin(c + d*x))/(63*d*e) + 2*a*e**2*(7*a**2 + 6*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(e*cos(c + d*x))) + 2*a*e*sqrt(e*cos(c + d*x))*(7*a**2 + 6*b**2)*sin(c + d*x)/(21*d) - 2*b*(e*cos(c + d*x))**(sympy.S(5)/2)*(a + b*sin(c + d*x))**2/(9*d*e) - 2*b*(e*cos(c + d*x))**(sympy.S(5)/2)*(89*a**2 + 28*b**2)/(315*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_556():
    f = sqrt(e*cos(c + d*x))*(a + b*sin(c + d*x))**3
    F = -22*a*b*(e*cos(c + d*x))**(sympy.S(3)/2)*(a + b*sin(c + d*x))/(35*d*e) + 2*a*sqrt(e*cos(c + d*x))*(5*a**2 + 6*b**2)*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) - 2*b*(e*cos(c + d*x))**(sympy.S(3)/2)*(a + b*sin(c + d*x))**2/(7*d*e) - 2*b*(e*cos(c + d*x))**(sympy.S(3)/2)*(57*a**2 + 20*b**2)/(105*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_557():
    f = (a + b*sin(c + d*x))**3/sqrt(e*cos(c + d*x))
    F = -6*a*b*sqrt(e*cos(c + d*x))*(a + b*sin(c + d*x))/(5*d*e) + 2*a*(a**2 + 2*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*sqrt(e*cos(c + d*x))) - 2*b*sqrt(e*cos(c + d*x))*(a + b*sin(c + d*x))**2/(5*d*e) - 2*b*sqrt(e*cos(c + d*x))*(11*a**2 + 4*b**2)/(5*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_558():
    f = (a + b*sin(c + d*x))**3/(e*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*a*b*(e*cos(c + d*x))**(sympy.S(3)/2)*(a + b*sin(c + d*x))/(d*e**3) - 2*a*sqrt(e*cos(c + d*x))*(a**2 + 6*b**2)*elliptic_e(c/2 + d*x/2, 2)/(d*e**2*sqrt(cos(c + d*x))) + 2*b*(e*cos(c + d*x))**(sympy.S(3)/2)*(3*a**2 + 4*b**2)/(3*d*e**3) + (a + b*sin(c + d*x))**2*(2*a*sin(c + d*x) + 2*b)/(d*e*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_559():
    f = (a + b*sin(c + d*x))**3/(e*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*a*b*sqrt(e*cos(c + d*x))*(a + b*sin(c + d*x))/(3*d*e**3) + 2*a*(a**2 - 6*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*e**2*sqrt(e*cos(c + d*x))) + 2*b*sqrt(e*cos(c + d*x))*(a**2 + 4*b**2)/(3*d*e**3) + (a + b*sin(c + d*x))**2*(2*a*sin(c + d*x) + 2*b)/(3*d*e*(e*cos(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_560():
    f = (a + b*sin(c + d*x))**3/(e*cos(c + d*x))**(sympy.S(7)/2)
    F = -6*a*sqrt(e*cos(c + d*x))*(a**2 - 2*b**2)*elliptic_e(c/2 + d*x/2, 2)/(5*d*e**4*sqrt(cos(c + d*x))) + 2*b*(e*cos(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 4*b**2)/(5*d*e**5) + (a + b*sin(c + d*x))**2*(2*a*sin(c + d*x) + 2*b)/(5*d*e*(e*cos(c + d*x))**(sympy.S(5)/2)) - (2*a + 2*b*sin(c + d*x))*(a*b - (3*a**2 - 4*b**2)*sin(c + d*x))/(5*d*e**3*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_561():
    f = (a + b*sin(c + d*x))**3/(e*cos(c + d*x))**(sympy.S(9)/2)
    F = 2*a*(5*a**2 - 6*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*e**4*sqrt(e*cos(c + d*x))) + 2*b*sqrt(e*cos(c + d*x))*(5*a**2 - 4*b**2)/(21*d*e**5) + (a + b*sin(c + d*x))**2*(2*a*sin(c + d*x) + 2*b)/(7*d*e*(e*cos(c + d*x))**(sympy.S(7)/2)) + (2*a + 2*b*sin(c + d*x))*(a*b + (5*a**2 - 4*b**2)*sin(c + d*x))/(21*d*e**3*(e*cos(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_562():
    f = (e*cos(c + d*x))**(sympy.S(7)/2)*(a + b*sin(c + d*x))**4
    F = -14*a*b*(e*cos(c + d*x))**(sympy.S(9)/2)*(a + b*sin(c + d*x))**2/(65*d*e) - 34*a*b*(e*cos(c + d*x))**(sympy.S(9)/2)*(53*a**2 + 38*b**2)/(6435*d*e) - 2*b*(e*cos(c + d*x))**(sympy.S(9)/2)*(a + b*sin(c + d*x))**3/(15*d*e) - 2*b*(e*cos(c + d*x))**(sympy.S(9)/2)*(a + b*sin(c + d*x))*(93*a**2 + 26*b**2)/(715*d*e) + e**4*(110*a**4 + 120*a**2*b**2 + 8*b**4)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(231*d*sqrt(e*cos(c + d*x))) + e**3*sqrt(e*cos(c + d*x))*(110*a**4 + 120*a**2*b**2 + 8*b**4)*sin(c + d*x)/(231*d) + e*(e*cos(c + d*x))**(sympy.S(5)/2)*(110*a**4 + 120*a**2*b**2 + 8*b**4)*sin(c + d*x)/(385*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_563():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)*(a + b*sin(c + d*x))**4
    F = -38*a*b*(e*cos(c + d*x))**(sympy.S(7)/2)*(a + b*sin(c + d*x))**2/(143*d*e) - 10*a*b*(e*cos(c + d*x))**(sympy.S(7)/2)*(115*a**2 + 94*b**2)/(3003*d*e) - 2*b*(e*cos(c + d*x))**(sympy.S(7)/2)*(a + b*sin(c + d*x))**3/(13*d*e) - 2*b*(e*cos(c + d*x))**(sympy.S(7)/2)*(a + b*sin(c + d*x))*(73*a**2 + 22*b**2)/(429*d*e) + e**2*sqrt(e*cos(c + d*x))*(78*a**4 + 104*a**2*b**2 + 8*b**4)*elliptic_e(c/2 + d*x/2, 2)/(65*d*sqrt(cos(c + d*x))) + e*(e*cos(c + d*x))**(sympy.S(3)/2)*(78*a**4 + 104*a**2*b**2 + 8*b**4)*sin(c + d*x)/(195*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_564():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)*(a + b*sin(c + d*x))**4
    F = -34*a*b*(e*cos(c + d*x))**(sympy.S(5)/2)*(a + b*sin(c + d*x))**2/(99*d*e) - 26*a*b*(e*cos(c + d*x))**(sympy.S(5)/2)*(79*a**2 + 74*b**2)/(3465*d*e) - 2*b*(e*cos(c + d*x))**(sympy.S(5)/2)*(a + b*sin(c + d*x))**3/(11*d*e) - 2*b*(e*cos(c + d*x))**(sympy.S(5)/2)*(a + b*sin(c + d*x))*(167*a**2 + 54*b**2)/(693*d*e) + e**2*(154*a**4 + 264*a**2*b**2 + 24*b**4)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(231*d*sqrt(e*cos(c + d*x))) + e*sqrt(e*cos(c + d*x))*(154*a**4 + 264*a**2*b**2 + 24*b**4)*sin(c + d*x)/(231*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_565():
    f = sqrt(e*cos(c + d*x))*(a + b*sin(c + d*x))**4
    F = -10*a*b*(e*cos(c + d*x))**(sympy.S(3)/2)*(a + b*sin(c + d*x))**2/(21*d*e) - 22*a*b*(e*cos(c + d*x))**(sympy.S(3)/2)*(17*a**2 + 18*b**2)/(315*d*e) - 2*b*(e*cos(c + d*x))**(sympy.S(3)/2)*(a + b*sin(c + d*x))**3/(9*d*e) - 2*b*(e*cos(c + d*x))**(sympy.S(3)/2)*(a + b*sin(c + d*x))*(41*a**2 + 14*b**2)/(105*d*e) + sqrt(e*cos(c + d*x))*(30*a**4 + 72*a**2*b**2 + 8*b**4)*elliptic_e(c/2 + d*x/2, 2)/(15*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_566():
    f = (a + b*sin(c + d*x))**4/sqrt(e*cos(c + d*x))
    F = -26*a*b*sqrt(e*cos(c + d*x))*(a + b*sin(c + d*x))**2/(35*d*e) - 6*a*b*sqrt(e*cos(c + d*x))*(31*a**2 + 34*b**2)/(35*d*e) - 2*b*sqrt(e*cos(c + d*x))*(a + b*sin(c + d*x))**3/(7*d*e) - 2*b*sqrt(e*cos(c + d*x))*(a + b*sin(c + d*x))*(29*a**2 + 10*b**2)/(35*d*e) + (14*a**4 + 56*a**2*b**2 + 8*b**4)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(7*d*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_567():
    f = (a + b*sin(c + d*x))**4/(e*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*a*b*(e*cos(c + d*x))**(sympy.S(3)/2)*(a + b*sin(c + d*x))**2/(d*e**3) + 2*a*b*(e*cos(c + d*x))**(sympy.S(3)/2)*(15*a**2 + 62*b**2)/(15*d*e**3) + 2*b*(e*cos(c + d*x))**(sympy.S(3)/2)*(a + b*sin(c + d*x))*(5*a**2 + 6*b**2)/(5*d*e**3) + (a + b*sin(c + d*x))**3*(2*a*sin(c + d*x) + 2*b)/(d*e*sqrt(e*cos(c + d*x))) - sqrt(e*cos(c + d*x))*(10*a**4 + 120*a**2*b**2 + 24*b**4)*elliptic_e(c/2 + d*x/2, 2)/(5*d*e**2*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_568():
    f = (a + b*sin(c + d*x))**4/(e*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*a*b*sqrt(e*cos(c + d*x))*(a + b*sin(c + d*x))**2/(3*d*e**3) + 2*a*b*sqrt(e*cos(c + d*x))*(a**2 + 14*b**2)/(3*d*e**3) + 2*b*sqrt(e*cos(c + d*x))*(a + b*sin(c + d*x))*(a**2 + 2*b**2)/(3*d*e**3) + (a + b*sin(c + d*x))**3*(2*a*sin(c + d*x) + 2*b)/(3*d*e*(e*cos(c + d*x))**(sympy.S(3)/2)) + (2*a**4 - 24*a**2*b**2 - 8*b**4)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*e**2*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_569():
    f = (a + b*sin(c + d*x))**4/(e*cos(c + d*x))**(sympy.S(7)/2)
    F = 2*a*b*(e*cos(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 10*b**2)/(5*d*e**5) + 6*b*(e*cos(c + d*x))**(sympy.S(3)/2)*(a + b*sin(c + d*x))*(a**2 - 2*b**2)/(5*d*e**5) + (a + b*sin(c + d*x))**3*(2*a*sin(c + d*x) + 2*b)/(5*d*e*(e*cos(c + d*x))**(sympy.S(5)/2)) - 6*(a + b*sin(c + d*x))**2*(a*b - (a**2 - 2*b**2)*sin(c + d*x))/(5*d*e**3*sqrt(e*cos(c + d*x))) - sqrt(e*cos(c + d*x))*(6*a**4 - 24*a**2*b**2 - 24*b**4)*elliptic_e(c/2 + d*x/2, 2)/(5*d*e**4*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_570():
    f = (a + b*sin(c + d*x))**4/(e*cos(c + d*x))**(sympy.S(9)/2)
    F = 10*a*b*sqrt(e*cos(c + d*x))*(a**2 - 2*b**2)/(21*d*e**5) + 2*b*sqrt(e*cos(c + d*x))*(a + b*sin(c + d*x))*(5*a**2 - 6*b**2)/(21*d*e**5) + (a + b*sin(c + d*x))**3*(2*a*sin(c + d*x) + 2*b)/(7*d*e*(e*cos(c + d*x))**(sympy.S(7)/2)) - 2*(a + b*sin(c + d*x))**2*(a*b - (5*a**2 - 6*b**2)*sin(c + d*x))/(21*d*e**3*(e*cos(c + d*x))**(sympy.S(3)/2)) + (10*a**4 - 24*a**2*b**2 + 24*b**4)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*e**4*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_571():
    f = (a + b*sin(c + d*x))**4/(e*cos(c + d*x))**(sympy.S(11)/2)
    F = 2*a*b*(e*cos(c + d*x))**(sympy.S(3)/2)*(21*a**2 - 22*b**2)/(45*d*e**7) + (a + b*sin(c + d*x))**3*(2*a*sin(c + d*x) + 2*b)/(9*d*e*(e*cos(c + d*x))**(sympy.S(9)/2)) + 2*(a + b*sin(c + d*x))**2*(a*b + (7*a**2 - 6*b**2)*sin(c + d*x))/(45*d*e**3*(e*cos(c + d*x))**(sympy.S(5)/2)) - (2*a + 2*b*sin(c + d*x))*(-a*(21*a**2 - 22*b**2)*sin(c + d*x) + b*(7*a**2 - 6*b**2))/(45*d*e**5*sqrt(e*cos(c + d*x))) - sqrt(e*cos(c + d*x))*(14*a**4 - 24*a**2*b**2 + 8*b**4)*elliptic_e(c/2 + d*x/2, 2)/(15*d*e**6*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_572():
    f = (e*cos(c + d*x))**(sympy.S(11)/2)/(a + b*sin(c + d*x))
    F = (Integer(-1) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(9) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('b'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(9) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('b'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(9) * (Integer(2))**(Integer(-1))))) * ((Integer(9) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('a') * ((Integer(21) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(49) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(33) * (Symbol('b'))**(Integer(4)))) * (Symbol('e'))**(Integer(6)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(21) * (Symbol('b'))**(Integer(6)) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * (Symbol('e'))**(Integer(6)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('b'))**(Integer(6)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * (Symbol('e'))**(Integer(6)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('b'))**(Integer(6)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('e'))**(Integer(3)) * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Integer(7) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(35) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (Symbol('e'))**(Integer(5)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Integer(21) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2))) + (Integer(-1) * (Symbol('a') * Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(12) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(21) * (Symbol('b'))**(Integer(5)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_573():
    f = (e*cos(c + d*x))**(sympy.S(9)/2)/(a + b*sin(c + d*x))
    F = (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(7) * (Integer(2))**(Integer(-1))))) * ((Integer(7) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(8) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(4)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(5) * (Symbol('b'))**(Integer(4)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('e'))**(Integer(5)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('b'))**(Integer(5)) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('e'))**(Integer(5)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('b'))**(Integer(5)) * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('e'))**(Integer(3)) * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(5) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(15) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_574():
    f = (e*cos(c + d*x))**(sympy.S(7)/2)/(a + b*sin(c + d*x))
    F = (Integer(-1) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(5) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(4)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(4)) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('e'))**(Integer(4)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('b'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('e'))**(Integer(4)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('b'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) + (Integer(-1) * (Symbol('a') * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_575():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)/(a + b*sin(c + d*x))
    F = (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('a') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(3)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('b'))**(Integer(3)) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(3)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('b'))**(Integer(3)) * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_576():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)/(a + b*sin(c + d*x))
    F = (Integer(-1) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('a') * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_577():
    f = sqrt(e*cos(c + d*x))/(a + b*sin(c + d*x))
    F = ((sympy.sqrt(Symbol('e')) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('b')) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('e')) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('b')) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * Symbol('d')))**(Integer(-1)))) + ((Symbol('a') * Symbol('e') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Symbol('b') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Symbol('a') * Symbol('e') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Symbol('b') * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_578():
    f = 1/(sqrt(e*cos(c + d*x))*(a + b*sin(c + d*x)))
    F = (Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('d') * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('d') * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((Symbol('a') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_579():
    f = 1/((e*cos(c + d*x))**(sympy.S(3)/2)*(a + b*sin(c + d*x)))
    F = (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_580():
    f = 1/((e*cos(c + d*x))**(sympy.S(5)/2)*(a + b*sin(c + d*x)))
    F = (Integer(-1) * (((Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_581():
    f = 1/((e*cos(c + d*x))**(sympy.S(7)/2)*(a + b*sin(c + d*x)))
    F = (((Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(9) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(9) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(8) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(5) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('e'))**(Integer(4)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Symbol('a') * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))) * Symbol('d') * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Symbol('a') * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))) * Symbol('d') * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(5) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * ((Integer(5) * (Symbol('b'))**(Integer(3))) + (Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(8) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(5) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_582():
    f = (e*cos(c + d*x))**(sympy.S(11)/2)/(a + b*sin(c + d*x))**2
    F = ((Integer(-9) * Symbol('a') * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * Symbol('a') * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * ((Integer(21) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(28) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(5) * (Symbol('b'))**(Integer(4)))) * (Symbol('e'))**(Integer(6)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(7) * (Symbol('b'))**(Integer(6)) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(9) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('e'))**(Integer(6)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(6)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(9) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('e'))**(Integer(6)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(6)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(9) * (Symbol('e'))**(Integer(3)) * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Integer(7) * Symbol('a')) + (Integer(-1) * (Integer(5) * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(35) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(9) * (Integer(2))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('e'))**(Integer(5)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Integer(21) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) + (Integer(-1) * (Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(7) * (Symbol('b'))**(Integer(5)) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_583():
    f = (e*cos(c + d*x))**(sympy.S(9)/2)/(a + b*sin(c + d*x))**2
    F = ((Integer(7) * Symbol('a') * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * Symbol('a') * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + ((Integer(7) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(4)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(5) * (Symbol('b'))**(Integer(4)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(5)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(5)) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(7) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(5)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(5)) * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(7) * (Symbol('e'))**(Integer(3)) * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(5) * Symbol('a')) + (Integer(-1) * (Integer(3) * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(15) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(7) * (Integer(2))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_584():
    f = (e*cos(c + d*x))**(sympy.S(7)/2)/(a + b*sin(c + d*x))**2
    F = ((Integer(-5) * Symbol('a') * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('a') * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + ((Integer(5) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(4)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(4)) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(4)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(4)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(5) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * Symbol('a')) + (Integer(-1) * (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_585():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)/(a + b*sin(c + d*x))**2
    F = ((Integer(3) * Symbol('a') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('e'))**(Integer(3)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('e'))**(Integer(3)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_586():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)/(a + b*sin(c + d*x))**2
    F = ((Integer(-1) * (Symbol('a') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('b') * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_587():
    f = sqrt(e*cos(c + d*x))/(a + b*sin(c + d*x))**2
    F = ((Integer(-1) * (Symbol('a') * sympy.sqrt(Symbol('e')) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt(Symbol('e')) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * Symbol('e') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * Symbol('e') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Symbol('b') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_588():
    f = 1/(sqrt(e*cos(c + d*x))*(a + b*sin(c + d*x))**2)
    F = ((Integer(3) * Symbol('a') * sympy.sqrt(Symbol('b')) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * Symbol('d') * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + ((Integer(3) * Symbol('a') * sympy.sqrt(Symbol('b')) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * Symbol('d') * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Symbol('b') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_589():
    f = 1/((e*cos(c + d*x))**(sympy.S(3)/2)*(a + b*sin(c + d*x))**2)
    F = ((Integer(-5) * Symbol('a') * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(9) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(5) * Symbol('a') * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(9) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(3) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Symbol('b') * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * (((Integer(5) * Symbol('a') * Symbol('b')) + (Integer(-1) * (((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(3) * (Symbol('b'))**(Integer(2)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_590():
    f = 1/((e*cos(c + d*x))**(sympy.S(5)/2)*(a + b*sin(c + d*x))**2)
    F = ((Integer(7) * Symbol('a') * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(11) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(7) * Symbol('a') * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(11) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(5) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(3) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(7) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Symbol('b') * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * (((Integer(7) * Symbol('a') * Symbol('b')) + (Integer(-1) * (((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(5) * (Symbol('b'))**(Integer(2)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(3) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_591():
    f = 1/((e*cos(c + d*x))**(sympy.S(7)/2)*(a + b*sin(c + d*x))**2)
    F = ((Integer(-9) * Symbol('a') * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(13) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(9) * Symbol('a') * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(13) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * ((Integer(2) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(10) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(4))))) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(5) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * Symbol('d') * (Symbol('e'))**(Integer(4)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(9) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))) * Symbol('d') * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(9) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))) * Symbol('d') * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Symbol('b') * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * (((Integer(9) * Symbol('a') * Symbol('b')) + (Integer(-1) * (((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(7) * (Symbol('b'))**(Integer(2)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(5) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * ((Integer(15) * Symbol('a') * (Symbol('b'))**(Integer(3))) + (((Integer(2) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(10) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(4))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(5) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * Symbol('d') * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_592():
    f = (e*cos(c + d*x))**(sympy.S(13)/2)/(a + b*sin(c + d*x))**3
    F = ((Integer(-11) * ((Integer(9) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(11) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(2) * (Symbol('b'))**(Integer(4)))) * (Symbol('e'))**((Integer(13) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(13) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * Symbol('d')))**(Integer(-1))) + ((Integer(11) * ((Integer(9) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(11) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(2) * (Symbol('b'))**(Integer(4)))) * (Symbol('e'))**((Integer(13) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(13) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * Symbol('d')))**(Integer(-1))) + ((Integer(11) * Symbol('a') * ((Integer(45) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(37) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(6)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(20) * (Symbol('b'))**(Integer(6)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(11) * Symbol('a') * ((Integer(9) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(11) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(2) * (Symbol('b'))**(Integer(4)))) * (Symbol('e'))**(Integer(7)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(7)) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(11) * Symbol('a') * ((Integer(9) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(11) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(2) * (Symbol('b'))**(Integer(4)))) * (Symbol('e'))**(Integer(7)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(7)) * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(11) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(11) * (Symbol('e'))**(Integer(3)) * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(7) * (Integer(2))**(Integer(-1)))) * ((Integer(9) * Symbol('a')) + (Integer(2) * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(28) * (Symbol('b'))**(Integer(3)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(11) * (Symbol('e'))**(Integer(5)) * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(5) * ((Integer(9) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2)))))) + (Integer(-1) * (Integer(27) * Symbol('a') * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(60) * (Symbol('b'))**(Integer(5)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_593():
    f = (e*cos(c + d*x))**(sympy.S(11)/2)/(a + b*sin(c + d*x))**3
    F = ((Integer(9) * ((Integer(7) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(9) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(2) * (Symbol('b'))**(Integer(4)))) * (Symbol('e'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + ((Integer(9) * ((Integer(7) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(9) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(2) * (Symbol('b'))**(Integer(4)))) * (Symbol('e'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * Symbol('a') * ((Integer(21) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(13) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(6)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(6)) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * Symbol('a') * ((Integer(7) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(9) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(2) * (Symbol('b'))**(Integer(4)))) * (Symbol('e'))**(Integer(6)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(6)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(9) * Symbol('a') * ((Integer(7) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(9) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(2) * (Symbol('b'))**(Integer(4)))) * (Symbol('e'))**(Integer(6)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(6)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(9) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(9) * (Symbol('e'))**(Integer(3)) * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Integer(7) * Symbol('a')) + (Integer(2) * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(20) * (Symbol('b'))**(Integer(3)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('e'))**(Integer(5)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2)))))) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(5)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_594():
    f = (e*cos(c + d*x))**(sympy.S(9)/2)/(a + b*sin(c + d*x))**3
    F = ((Integer(7) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(35) * Symbol('a') * (Symbol('e'))**(Integer(4)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(7) * Symbol('a') * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(5)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(5)) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(7) * Symbol('a') * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(5)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(5)) * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(7) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(7) * (Symbol('e'))**(Integer(3)) * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(5) * Symbol('a')) + (Integer(2) * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(12) * (Symbol('b'))**(Integer(3)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_595():
    f = (e*cos(c + d*x))**(sympy.S(7)/2)/(a + b*sin(c + d*x))**3
    F = ((Integer(-5) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(15) * Symbol('a') * (Symbol('e'))**(Integer(4)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(5) * Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(4)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(5) * Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(4)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * Symbol('a')) + (Integer(2) * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_596():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)/(a + b*sin(c + d*x))**3
    F = ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(3)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(3)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(4) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_597():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)/(a + b*sin(c + d*x))**3
    F = ((((Symbol('a'))**(Integer(2)) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + ((((Symbol('a'))**(Integer(2)) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(2) * Symbol('b') * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + ((Symbol('a') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(4) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_598():
    f = sqrt(e*cos(c + d*x))/(a + b*sin(c + d*x))**3
    F = ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(Symbol('e')) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('b')) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(9) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(Symbol('e')) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('b')) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(9) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + ((Integer(5) * Symbol('a') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * Symbol('e') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * Symbol('e') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Symbol('b') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Integer(5) * Symbol('a') * Symbol('b') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_599():
    f = 1/(sqrt(e*cos(c + d*x))*(a + b*sin(c + d*x))**3)
    F = ((Integer(-3) * sympy.sqrt(Symbol('b')) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(11) * (Integer(4))**(Integer(-1)))) * Symbol('d') * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(Symbol('b')) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(11) * (Integer(4))**(Integer(-1)))) * Symbol('d') * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(7) * Symbol('a') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(3) * Symbol('a') * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Symbol('b') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Integer(7) * Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_600():
    f = 1/((e*cos(c + d*x))**(sympy.S(3)/2)*(a + b*sin(c + d*x))**3)
    F = ((Integer(5) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(13) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(13) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * ((Integer(8) * (Symbol('a'))**(Integer(2))) + (Integer(37) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('a') * Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('a') * Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Symbol('b') * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Integer(9) * Symbol('a') * Symbol('b')) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * (((Integer(5) * Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2))))) + (Integer(-1) * (Symbol('a') * ((Integer(8) * (Symbol('a'))**(Integer(2))) + (Integer(37) * (Symbol('b'))**(Integer(2)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_601():
    f = 1/((e*cos(c + d*x))**(sympy.S(5)/2)*(a + b*sin(c + d*x))**3)
    F = ((Integer(-7) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Integer(9) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(15) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Integer(9) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(15) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('a') * ((Integer(8) * (Symbol('a'))**(Integer(2))) + (Integer(69) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(12) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * ((Integer(9) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * ((Integer(9) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Symbol('b') * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Integer(11) * Symbol('a') * Symbol('b')) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * (((Integer(7) * Symbol('b') * ((Integer(9) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2))))) + (Integer(-1) * (Symbol('a') * ((Integer(8) * (Symbol('a'))**(Integer(2))) + (Integer(69) * (Symbol('b'))**(Integer(2)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(12) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * Symbol('d') * Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_602():
    f = 1/((e*cos(c + d*x))**(sympy.S(7)/2)*(a + b*sin(c + d*x))**3)
    F = ((Integer(9) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * ((Integer(11) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(17) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * ((Integer(11) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(8) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(17) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('a') * ((Integer(8) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(64) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Integer(139) * (Symbol('b'))**(Integer(4))))) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(20) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(4)) * Symbol('d') * (Symbol('e'))**(Integer(4)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(9) * Symbol('a') * (Symbol('b'))**(Integer(3)) * ((Integer(11) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(4)) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))) * Symbol('d') * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(9) * Symbol('a') * (Symbol('b'))**(Integer(3)) * ((Integer(11) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(4)) * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))) * Symbol('d') * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Symbol('b') * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Integer(13) * Symbol('a') * Symbol('b')) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * (((Integer(9) * Symbol('b') * ((Integer(11) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2))))) + (Integer(-1) * (Symbol('a') * ((Integer(8) * (Symbol('a'))**(Integer(2))) + (Integer(109) * (Symbol('b'))**(Integer(2)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(20) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * Symbol('d') * Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * ((Integer(15) * (Symbol('b'))**(Integer(3)) * ((Integer(11) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2))))) + (Symbol('a') * ((Integer(8) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(64) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Integer(139) * (Symbol('b'))**(Integer(4))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(20) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(4)) * Symbol('d') * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_603():
    f = (e*cos(c + d*x))**(sympy.S(15)/2)/(a + b*sin(c + d*x))**4
    F = ((Integer(39) * Symbol('a') * ((Integer(11) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(17) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(6) * (Symbol('b'))**(Integer(4)))) * (Symbol('e'))**((Integer(15) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(15) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + ((Integer(39) * Symbol('a') * ((Integer(11) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(17) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(6) * (Symbol('b'))**(Integer(4)))) * (Symbol('e'))**((Integer(15) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(15) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + ((Integer(13) * ((Integer(231) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(203) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(20) * (Symbol('b'))**(Integer(4)))) * (Symbol('e'))**(Integer(8)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(56) * (Symbol('b'))**(Integer(8)) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(39) * (Symbol('a'))**(Integer(2)) * ((Integer(11) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(17) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(6) * (Symbol('b'))**(Integer(4)))) * (Symbol('e'))**(Integer(8)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(8)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(39) * (Symbol('a'))**(Integer(2)) * ((Integer(11) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(17) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(6) * (Symbol('b'))**(Integer(4)))) * (Symbol('e'))**(Integer(8)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(8)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(13) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * Symbol('b') * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(13) * (Symbol('e'))**(Integer(3)) * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(9) * (Integer(2))**(Integer(-1)))) * ((Integer(11) * Symbol('a')) + (Integer(4) * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(84) * (Symbol('b'))**(Integer(3)) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(39) * (Symbol('e'))**(Integer(5)) * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Integer(77) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(20) * (Symbol('b'))**(Integer(2)))) + (Integer(22) * Symbol('a') * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(280) * (Symbol('b'))**(Integer(5)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(13) * (Symbol('e'))**(Integer(7)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Integer(21) * Symbol('a') * ((Integer(11) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(6) * (Symbol('b'))**(Integer(2)))))) + (Integer(-1) * (Symbol('b') * ((Integer(77) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(20) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(56) * (Symbol('b'))**(Integer(7)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_604():
    f = (e*cos(c + d*x))**(sympy.S(13)/2)/(a + b*sin(c + d*x))**4
    F = ((Integer(77) * Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**((Integer(13) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(13) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(77) * Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**((Integer(13) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(13) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(77) * ((Integer(15) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(6)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(40) * (Symbol('b'))**(Integer(6)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(77) * (Symbol('a'))**(Integer(2)) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(7)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(7)) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(77) * (Symbol('a'))**(Integer(2)) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(7)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(7)) * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(11) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * Symbol('b') * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(11) * (Symbol('e'))**(Integer(3)) * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(7) * (Integer(2))**(Integer(-1)))) * ((Integer(9) * Symbol('a')) + (Integer(4) * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(60) * (Symbol('b'))**(Integer(3)) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(77) * (Symbol('e'))**(Integer(5)) * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(15) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2)))) + (Integer(6) * Symbol('a') * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(120) * (Symbol('b'))**(Integer(5)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_605():
    f = (e*cos(c + d*x))**(sympy.S(11)/2)/(a + b*sin(c + d*x))**4
    F = (Integer(-1) * ((Integer(15) * Symbol('a') * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(6) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(15) * Symbol('a') * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(6) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * ((Integer(21) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(6)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(6)) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(15) * (Symbol('a'))**(Integer(2)) * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(6) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(6)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(6)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(15) * (Symbol('a'))**(Integer(2)) * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(6) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(6)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(6)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(9) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * Symbol('b') * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('e'))**(Integer(3)) * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Integer(7) * Symbol('a')) + (Integer(4) * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Symbol('e'))**(Integer(5)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Integer(21) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2)))) + (Integer(14) * Symbol('a') * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(5)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_606():
    f = (e*cos(c + d*x))**(sympy.S(9)/2)/(a + b*sin(c + d*x))**4
    F = ((Integer(7) * Symbol('a') * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(6) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * Symbol('a') * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(6) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + ((Integer(7) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(4)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * (Symbol('a'))**(Integer(2)) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(6) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(5)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(5)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(7) * (Symbol('a'))**(Integer(2)) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(6) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(5)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(5)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(7) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * Symbol('b') * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))))**(Integer(-1)))) + ((Integer(7) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(3)) * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * (Symbol('e'))**(Integer(3)) * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(5) * Symbol('a')) + (Integer(4) * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(12) * (Symbol('b'))**(Integer(3)) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_607():
    f = (e*cos(c + d*x))**(sympy.S(7)/2)/(a + b*sin(c + d*x))**4
    F = ((Integer(-5) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + ((Integer(5) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(4)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(24) * (Symbol('b'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(4)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(4)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * Symbol('b') * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(24) * (Symbol('b'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(5) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * Symbol('a')) + (Integer(4) * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(12) * (Symbol('b'))**(Integer(3)) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_608():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)/(a + b*sin(c + d*x))**4
    F = ((Integer(-1) * (Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(6) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))))) * ((Integer(16) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(9) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(6) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(9) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + ((((Symbol('a'))**(Integer(2)) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(6) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(3)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(6) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(3)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * Symbol('b') * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))))**(Integer(-1)))) + ((Symbol('a') * Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(4) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((((Symbol('a'))**(Integer(2)) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * Symbol('e') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(8) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_609():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)/(a + b*sin(c + d*x))**4
    F = ((Integer(-1) * (Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(6) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))))) * ((Integer(16) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(11) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(6) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(11) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(24) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(6) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(6) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(3) * Symbol('b') * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))))**(Integer(-1)))) + ((Symbol('a') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(12) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(24) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_610():
    f = sqrt(e*cos(c + d*x))/(a + b*sin(c + d*x))**4
    F = ((Integer(-5) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(Symbol('e')) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('b')) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(13) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + ((Integer(5) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(Symbol('e')) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('b')) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(13) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + ((((Integer(11) * (Symbol('a'))**(Integer(2))) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Integer(5) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * Symbol('e') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(16) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(5) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * Symbol('e') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(16) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Symbol('b') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * Symbol('a') * Symbol('b') * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * Symbol('e') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * ((Integer(11) * (Symbol('a'))**(Integer(2))) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * ((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * Symbol('d') * Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_611():
    f = 1/(sqrt(e*cos(c + d*x))*(a + b*sin(c + d*x))**4)
    F = ((Integer(7) * Symbol('a') * sympy.sqrt(Symbol('b')) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(6) * (Symbol('b'))**(Integer(2)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(16) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(15) * (Integer(4))**(Integer(-1)))) * Symbol('d') * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + ((Integer(7) * Symbol('a') * sympy.sqrt(Symbol('b')) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(6) * (Symbol('b'))**(Integer(2)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(16) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(15) * (Integer(4))**(Integer(-1)))) * Symbol('d') * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((((Integer(57) * (Symbol('a'))**(Integer(2))) + (Integer(20) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(24) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(7) * (Symbol('a'))**(Integer(2)) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(6) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(16) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(7) * (Symbol('a'))**(Integer(2)) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(6) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(16) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b') * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Symbol('b') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))))**(Integer(-1))) + ((Integer(11) * Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(12) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * Symbol('e') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * ((Integer(57) * (Symbol('a'))**(Integer(2))) + (Integer(20) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(24) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * Symbol('d') * Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_612():
    f = 1/((e*cos(c + d*x))**(sympy.S(3)/2)*(a + b*sin(c + d*x))**4)
    F = ((Integer(-15) * Symbol('a') * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(6) * (Symbol('b'))**(Integer(2)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(16) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(17) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(15) * Symbol('a') * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(6) * (Symbol('b'))**(Integer(2)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (((((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(16) * (((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))**((Integer(17) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(16) * (Symbol('a'))**(Integer(4))) + (Integer(151) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2))) + (Integer(28) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(4)) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(15) * (Symbol('a'))**(Integer(2)) * Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(6) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(16) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(4)) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(15) * (Symbol('a'))**(Integer(2)) * Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(6) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(16) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(4)) * (Symbol('b') + sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Symbol('b') * ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))))**(Integer(-1))) + ((Integer(13) * Symbol('a') * Symbol('b')) * ((Integer(12) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * ((Integer(89) * (Symbol('a'))**(Integer(2))) + (Integer(28) * (Symbol('b'))**(Integer(2))))) * ((Integer(24) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * (((Integer(15) * Symbol('a') * Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(6) * (Symbol('b'))**(Integer(2))))) + (Integer(-1) * (((Integer(16) * (Symbol('a'))**(Integer(4))) + (Integer(151) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2))) + (Integer(28) * (Symbol('b'))**(Integer(4)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(4)) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_613():
    f = (e*cos(c + d*x))**p*(a + b*sin(c + d*x))**3
    F = -a*b*(e*cos(c + d*x))**(p + 1)*(a + b*sin(c + d*x))*(p + 5)/(d*e*(p + 2)*(p + 3)) - a*(e*cos(c + d*x))**(p + 1)*(a**2*(p + 2) + 3*b**2)*sin(c + d*x)*hyper((sympy.S.Half, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(d*e*(p + 1)*(p + 2)*sqrt(sin(c + d*x)**2)) - b*(e*cos(c + d*x))**(p + 1)*(a + b*sin(c + d*x))**2/(d*e*(p + 3)) - b*(e*cos(c + d*x))**(p + 1)*(a**2*(p**2 + 6*p + 11) + 2*b**2*(p + 2))/(d*e*(p + 1)*(p + 2)*(p + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_614():
    f = (e*cos(c + d*x))**p*(a + b*sin(c + d*x))**2
    F = -a*b*(e*cos(c + d*x))**(p + 1)*(p + 3)/(d*e*(p + 1)*(p + 2)) - b*(e*cos(c + d*x))**(p + 1)*(a + b*sin(c + d*x))/(d*e*(p + 2)) - (e*cos(c + d*x))**(p + 1)*(a**2*(p + 2) + b**2)*sin(c + d*x)*hyper((sympy.S.Half, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(d*e*(p + 1)*(p + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_615():
    f = (e*cos(c + d*x))**p*(a + b*sin(c + d*x))
    F = -a*(e*cos(c + d*x))**(p + 1)*sin(c + d*x)*hyper((sympy.S.Half, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(d*e*(p + 1)*sqrt(sin(c + d*x)**2)) - b*(e*cos(c + d*x))**(p + 1)/(d*e*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_616():
    f = (e*cos(c + d*x))**p/(a + b*sin(c + d*x))
    F = -e*(e*cos(c + d*x))**(p - 1)*(b*(sin(c + d*x) + 1)/(a + b*sin(c + d*x)))**(sympy.S.Half - p/2)*(-b*(1 - sin(c + d*x))/(a + b*sin(c + d*x)))**(sympy.S.Half - p/2)*appellf1(1 - p, sympy.S.Half - p/2, sympy.S.Half - p/2, 2 - p, (a - b)/(a + b*sin(c + d*x)), (a + b)/(a + b*sin(c + d*x)))/(b*d*(1 - p))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_617():
    f = (e*cos(c + d*x))**p/(a + b*sin(c + d*x))**2
    F = -e*(e*cos(c + d*x))**(p - 1)*(b*(sin(c + d*x) + 1)/(a + b*sin(c + d*x)))**(sympy.S.Half - p/2)*(-b*(1 - sin(c + d*x))/(a + b*sin(c + d*x)))**(sympy.S.Half - p/2)*appellf1(2 - p, sympy.S.Half - p/2, sympy.S.Half - p/2, 3 - p, (a - b)/(a + b*sin(c + d*x)), (a + b)/(a + b*sin(c + d*x)))/(b*d*(2 - p)*(a + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_618():
    f = (e*cos(c + d*x))**p/(a + b*sin(c + d*x))**3
    F = -e*(e*cos(c + d*x))**(p - 1)*(b*(sin(c + d*x) + 1)/(a + b*sin(c + d*x)))**(sympy.S.Half - p/2)*(-b*(1 - sin(c + d*x))/(a + b*sin(c + d*x)))**(sympy.S.Half - p/2)*appellf1(3 - p, sympy.S.Half - p/2, sympy.S.Half - p/2, 4 - p, (a - b)/(a + b*sin(c + d*x)), (a + b)/(a + b*sin(c + d*x)))/(b*d*(3 - p)*(a + b*sin(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_619():
    f = (e*cos(c + d*x))**p/(a + b*sin(c + d*x))**8
    F = -e*(e*cos(c + d*x))**(p - 1)*(b*(sin(c + d*x) + 1)/(a + b*sin(c + d*x)))**(sympy.S.Half - p/2)*(-b*(1 - sin(c + d*x))/(a + b*sin(c + d*x)))**(sympy.S.Half - p/2)*appellf1(8 - p, sympy.S.Half - p/2, sympy.S.Half - p/2, 9 - p, (a - b)/(a + b*sin(c + d*x)), (a + b)/(a + b*sin(c + d*x)))/(b*d*(8 - p)*(a + b*sin(c + d*x))**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_620():
    f = (e*cos(c + d*x))**p*(a + b*sin(c + d*x))**(sympy.S(5)/2)
    F = 2*e*(e*cos(c + d*x))**(p - 1)*(1 - (a + b*sin(c + d*x))/(a - b))**(sympy.S.Half - p/2)*(1 - (a + b*sin(c + d*x))/(a + b))**(sympy.S.Half - p/2)*(a + b*sin(c + d*x))**(sympy.S(7)/2)*appellf1(sympy.S(7)/2, sympy.S.Half - p/2, sympy.S.Half - p/2, sympy.S(9)/2, (a + b*sin(c + d*x))/(a - b), (a + b*sin(c + d*x))/(a + b))/(7*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_621():
    f = (e*cos(c + d*x))**p*(a + b*sin(c + d*x))**(sympy.S(3)/2)
    F = 2*e*(e*cos(c + d*x))**(p - 1)*(1 - (a + b*sin(c + d*x))/(a - b))**(sympy.S.Half - p/2)*(1 - (a + b*sin(c + d*x))/(a + b))**(sympy.S.Half - p/2)*(a + b*sin(c + d*x))**(sympy.S(5)/2)*appellf1(sympy.S(5)/2, sympy.S.Half - p/2, sympy.S.Half - p/2, sympy.S(7)/2, (a + b*sin(c + d*x))/(a - b), (a + b*sin(c + d*x))/(a + b))/(5*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_622():
    f = (e*cos(c + d*x))**p*sqrt(a + b*sin(c + d*x))
    F = 2*e*(e*cos(c + d*x))**(p - 1)*(1 - (a + b*sin(c + d*x))/(a - b))**(sympy.S.Half - p/2)*(1 - (a + b*sin(c + d*x))/(a + b))**(sympy.S.Half - p/2)*(a + b*sin(c + d*x))**(sympy.S(3)/2)*appellf1(sympy.S(3)/2, sympy.S.Half - p/2, sympy.S.Half - p/2, sympy.S(5)/2, (a + b*sin(c + d*x))/(a - b), (a + b*sin(c + d*x))/(a + b))/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_623():
    f = (e*cos(c + d*x))**p/sqrt(a + b*sin(c + d*x))
    F = 2*e*(e*cos(c + d*x))**(p - 1)*(1 - (a + b*sin(c + d*x))/(a - b))**(sympy.S.Half - p/2)*(1 - (a + b*sin(c + d*x))/(a + b))**(sympy.S.Half - p/2)*sqrt(a + b*sin(c + d*x))*appellf1(sympy.S.Half, sympy.S.Half - p/2, sympy.S.Half - p/2, sympy.S(3)/2, (a + b*sin(c + d*x))/(a - b), (a + b*sin(c + d*x))/(a + b))/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_624():
    f = (e*cos(c + d*x))**p/(a + b*sin(c + d*x))**(sympy.S(3)/2)
    F = -2*e*(e*cos(c + d*x))**(p - 1)*(1 - (a + b*sin(c + d*x))/(a - b))**(sympy.S.Half - p/2)*(1 - (a + b*sin(c + d*x))/(a + b))**(sympy.S.Half - p/2)*appellf1(sympy.S(-1)/2, sympy.S.Half - p/2, sympy.S.Half - p/2, sympy.S.Half, (a + b*sin(c + d*x))/(a - b), (a + b*sin(c + d*x))/(a + b))/(b*d*sqrt(a + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_625():
    f = (e*cos(c + d*x))**p/(a + b*sin(c + d*x))**(sympy.S(5)/2)
    F = -2*e*(e*cos(c + d*x))**(p - 1)*(1 - (a + b*sin(c + d*x))/(a - b))**(sympy.S.Half - p/2)*(1 - (a + b*sin(c + d*x))/(a + b))**(sympy.S.Half - p/2)*appellf1(sympy.S(-3)/2, sympy.S.Half - p/2, sympy.S.Half - p/2, sympy.S(-1)/2, (a + b*sin(c + d*x))/(a - b), (a + b*sin(c + d*x))/(a + b))/(3*b*d*(a + b*sin(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_626():
    f = (e*cos(c + d*x))**p*(a + b*sin(c + d*x))**m
    F = e*(e*cos(c + d*x))**(p - 1)*(1 - (a + b*sin(c + d*x))/(a - b))**(sympy.S.Half - p/2)*(1 - (a + b*sin(c + d*x))/(a + b))**(sympy.S.Half - p/2)*(a + b*sin(c + d*x))**(m + 1)*appellf1(m + 1, sympy.S.Half - p/2, sympy.S.Half - p/2, m + 2, (a + b*sin(c + d*x))/(a - b), (a + b*sin(c + d*x))/(a + b))/(b*d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_627():
    f = (a + b*sin(c + d*x))**m*cos(c + d*x)**7
    F = 6*a*(a + b*sin(c + d*x))**(m + 2)*(a**2 - b**2)**2/(b**7*d*(m + 2)) + 4*a*(a + b*sin(c + d*x))**(m + 4)*(5*a**2 - 3*b**2)/(b**7*d*(m + 4)) + 6*a*(a + b*sin(c + d*x))**(m + 6)/(b**7*d*(m + 6)) - (a + b*sin(c + d*x))**(m + 1)*(a**2 - b**2)**3/(b**7*d*(m + 1)) - (a + b*sin(c + d*x))**(m + 3)*(15*a**4 - 18*a**2*b**2 + 3*b**4)/(b**7*d*(m + 3)) - (a + b*sin(c + d*x))**(m + 5)*(15*a**2 - 3*b**2)/(b**7*d*(m + 5)) - (a + b*sin(c + d*x))**(m + 7)/(b**7*d*(m + 7))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_628():
    f = (a + b*sin(c + d*x))**m*cos(c + d*x)**5
    F = -4*a*(a + b*sin(c + d*x))**(m + 2)*(a**2 - b**2)/(b**5*d*(m + 2)) - 4*a*(a + b*sin(c + d*x))**(m + 4)/(b**5*d*(m + 4)) + (a + b*sin(c + d*x))**(m + 1)*(a**2 - b**2)**2/(b**5*d*(m + 1)) + (a + b*sin(c + d*x))**(m + 3)*(6*a**2 - 2*b**2)/(b**5*d*(m + 3)) + (a + b*sin(c + d*x))**(m + 5)/(b**5*d*(m + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_629():
    f = (a + b*sin(c + d*x))**m*cos(c + d*x)**3
    F = 2*a*(a + b*sin(c + d*x))**(m + 2)/(b**3*d*(m + 2)) - (a + b*sin(c + d*x))**(m + 1)*(a**2 - b**2)/(b**3*d*(m + 1)) - (a + b*sin(c + d*x))**(m + 3)/(b**3*d*(m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_630():
    f = (a + b*sin(c + d*x))**m*cos(c + d*x)
    F = (a + b*sin(c + d*x))**(m + 1)/(b*d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_631():
    f = (a + b*sin(c + d*x))**m*sec(c + d*x)
    F = (a + b*sin(c + d*x))**(m + 1)*hyper((1, m + 1), (m + 2,), (a + b*sin(c + d*x))/(a + b))/(d*(2*a + 2*b)*(m + 1)) - (a + b*sin(c + d*x))**(m + 1)*hyper((1, m + 1), (m + 2,), (a + b*sin(c + d*x))/(a - b))/(d*(2*a - 2*b)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_632():
    f = (a + b*sin(c + d*x))**m*sec(c + d*x)**3
    F = -(a + b*sin(c + d*x))**(m + 1)*(-a*sin(c + d*x) + b)*sec(c + d*x)**2/(d*(2*a**2 - 2*b**2)) + (a + b*sin(c + d*x))**(m + 1)*(a - b*m + b)*hyper((1, m + 1), (m + 2,), (a + b*sin(c + d*x))/(a + b))/(4*d*(a + b)**2*(m + 1)) - (a - b*(1 - m))*(a + b*sin(c + d*x))**(m + 1)*hyper((1, m + 1), (m + 2,), (a + b*sin(c + d*x))/(a - b))/(4*d*(a - b)**2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_633():
    f = (a + b*sin(c + d*x))**m*sec(c + d*x)**5
    F = -(a + b*sin(c + d*x))**(m + 1)*(-a*sin(c + d*x) + b)*sec(c + d*x)**4/(d*(4*a**2 - 4*b**2)) + (a + b*sin(c + d*x))**(m + 1)*(a*(3*a**2 - b**2*(5 - 2*m))*sin(c + d*x) + b*(-a**2*(m + 1) + b**2*(3 - m)))*sec(c + d*x)**2/(8*d*(a**2 - b**2)**2) + (a + b*sin(c + d*x))**(m + 1)*(3*a**2 + 3*a*b*(2 - m) + b**2*(m**2 - 4*m + 3))*hyper((1, m + 1), (m + 2,), (a + b*sin(c + d*x))/(a + b))/(16*d*(a + b)**3*(m + 1)) - (a + b*sin(c + d*x))**(m + 1)*(3*a**2 - 3*a*b*(2 - m) + b**2*(m**2 - 4*m + 3))*hyper((1, m + 1), (m + 2,), (a + b*sin(c + d*x))/(a - b))/(16*d*(a - b)**3*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_634():
    f = (a + b*sin(c + d*x))**m*cos(c + d*x)**4
    F = (a + b*sin(c + d*x))**(m + 1)*cos(c + d*x)**3*appellf1(m + 1, sympy.S(-3)/2, sympy.S(-3)/2, m + 2, (a + b*sin(c + d*x))/(a - b), (a + b*sin(c + d*x))/(a + b))/(b*d*(1 - (a + b*sin(c + d*x))/(a - b))**(sympy.S(3)/2)*(1 - (a + b*sin(c + d*x))/(a + b))**(sympy.S(3)/2)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_635():
    f = (a + b*sin(c + d*x))**m*cos(c + d*x)**2
    F = (a + b*sin(c + d*x))**(m + 1)*cos(c + d*x)*appellf1(m + 1, sympy.S(-1)/2, sympy.S(-1)/2, m + 2, (a + b*sin(c + d*x))/(a - b), (a + b*sin(c + d*x))/(a + b))/(b*d*sqrt(1 - (a + b*sin(c + d*x))/(a - b))*sqrt(1 - (a + b*sin(c + d*x))/(a + b))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_636():
    f = (a + b*sin(c + d*x))**m*sec(c + d*x)**2
    F = (1 - (a + b*sin(c + d*x))/(a - b))**(sympy.S(3)/2)*(1 - (a + b*sin(c + d*x))/(a + b))**(sympy.S(3)/2)*(a + b*sin(c + d*x))**(m + 1)*appellf1(m + 1, sympy.S(3)/2, sympy.S(3)/2, m + 2, (a + b*sin(c + d*x))/(a - b), (a + b*sin(c + d*x))/(a + b))*sec(c + d*x)**3/(b*d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_637():
    f = (a + b*sin(c + d*x))**m*sec(c + d*x)**4
    F = (1 - (a + b*sin(c + d*x))/(a - b))**(sympy.S(5)/2)*(1 - (a + b*sin(c + d*x))/(a + b))**(sympy.S(5)/2)*(a + b*sin(c + d*x))**(m + 1)*appellf1(m + 1, sympy.S(5)/2, sympy.S(5)/2, m + 2, (a + b*sin(c + d*x))/(a - b), (a + b*sin(c + d*x))/(a + b))*sec(c + d*x)**5/(b*d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_638():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)*(a + b*sin(c + d*x))**m
    F = e*(e*cos(c + d*x))**(sympy.S(3)/2)*(a + b*sin(c + d*x))**(m + 1)*appellf1(m + 1, sympy.S(-3)/4, sympy.S(-3)/4, m + 2, (a + b*sin(c + d*x))/(a - b), (a + b*sin(c + d*x))/(a + b))/(b*d*(1 - (a + b*sin(c + d*x))/(a - b))**(sympy.S(3)/4)*(1 - (a + b*sin(c + d*x))/(a + b))**(sympy.S(3)/4)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_639():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)*(a + b*sin(c + d*x))**m
    F = e*sqrt(e*cos(c + d*x))*(a + b*sin(c + d*x))**(m + 1)*appellf1(m + 1, sympy.S(-1)/4, sympy.S(-1)/4, m + 2, (a + b*sin(c + d*x))/(a - b), (a + b*sin(c + d*x))/(a + b))/(b*d*(1 - (a + b*sin(c + d*x))/(a - b))**(sympy.S(1)/4)*(1 - (a + b*sin(c + d*x))/(a + b))**(sympy.S(1)/4)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_640():
    f = sqrt(e*cos(c + d*x))*(a + b*sin(c + d*x))**m
    F = e*(1 - (a + b*sin(c + d*x))/(a - b))**(sympy.S(1)/4)*(1 - (a + b*sin(c + d*x))/(a + b))**(sympy.S(1)/4)*(a + b*sin(c + d*x))**(m + 1)*appellf1(m + 1, sympy.S(1)/4, sympy.S(1)/4, m + 2, (a + b*sin(c + d*x))/(a - b), (a + b*sin(c + d*x))/(a + b))/(b*d*sqrt(e*cos(c + d*x))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_641():
    f = (a + b*sin(c + d*x))**m/sqrt(e*cos(c + d*x))
    F = e*(1 - (a + b*sin(c + d*x))/(a - b))**(sympy.S(3)/4)*(1 - (a + b*sin(c + d*x))/(a + b))**(sympy.S(3)/4)*(a + b*sin(c + d*x))**(m + 1)*appellf1(m + 1, sympy.S(3)/4, sympy.S(3)/4, m + 2, (a + b*sin(c + d*x))/(a - b), (a + b*sin(c + d*x))/(a + b))/(b*d*(e*cos(c + d*x))**(sympy.S(3)/2)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_642():
    f = (a + b*sin(c + d*x))**m/(e*cos(c + d*x))**(sympy.S(3)/2)
    F = e*(1 - (a + b*sin(c + d*x))/(a - b))**(sympy.S(5)/4)*(1 - (a + b*sin(c + d*x))/(a + b))**(sympy.S(5)/4)*(a + b*sin(c + d*x))**(m + 1)*appellf1(m + 1, sympy.S(5)/4, sympy.S(5)/4, m + 2, (a + b*sin(c + d*x))/(a - b), (a + b*sin(c + d*x))/(a + b))/(b*d*(e*cos(c + d*x))**(sympy.S(5)/2)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_643():
    f = (a + b*sin(c + d*x))**m/(e*cos(c + d*x))**(sympy.S(5)/2)
    F = e*(1 - (a + b*sin(c + d*x))/(a - b))**(sympy.S(7)/4)*(1 - (a + b*sin(c + d*x))/(a + b))**(sympy.S(7)/4)*(a + b*sin(c + d*x))**(m + 1)*appellf1(m + 1, sympy.S(7)/4, sympy.S(7)/4, m + 2, (a + b*sin(c + d*x))/(a - b), (a + b*sin(c + d*x))/(a + b))/(b*d*(e*cos(c + d*x))**(sympy.S(7)/2)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_644():
    f = (e*cos(c + d*x))**(-m - 4)*(a + b*sin(c + d*x))**m
    F = -2**(sympy.S(3)/2 - m/2)*a*b*(e*cos(c + d*x))**(-m - 1)*((a + b)*(sin(c + d*x) + 1)/(a + b*sin(c + d*x)))**(m/2 + sympy.S.Half)*(a + b*sin(c + d*x))**(m + 1)*hyper((-m/2 + sympy.S(-1)/2, m/2 + sympy.S.Half), (sympy.S.Half - m/2,), (1 - sin(c + d*x))*(a - b)/(2*a + 2*b*sin(c + d*x)))/(d*e**3*(a - b)**2*(a + b)*(m + 1)*(m + 3)) - 2**(-m/2 + sympy.S(-1)/2)*a*(e*cos(c + d*x))**(-m - 3)*((a + b)*(sin(c + d*x) + 1)/(a + b*sin(c + d*x)))**(m/2 + sympy.S(3)/2)*(1 - sin(c + d*x))**2*(a + b*sin(c + d*x))**(m + 1)*(a**2*(m + 2) + 2*a*b - b**2)*hyper((sympy.S.Half - m/2, m/2 + sympy.S(3)/2), (sympy.S(3)/2 - m/2,), (1 - sin(c + d*x))*(a - b)/(2*a + 2*b*sin(c + d*x)))/(d*e*(1 - m)*(a - b)*(a + b)**3*(m + 3)) + a*(e*cos(c + d*x))**(-m - 3)*(1 - sin(c + d*x))*(a + b*sin(c + d*x))**(m + 1)*(a*(m + 2) + 3*b)*(sin(c + d*x) + 1)/(d*e*(a - b)*(a + b)**2*(m + 1)*(m + 3)) + a*(e*cos(c + d*x))**(-m - 3)*(a + b*sin(c + d*x))**(m + 1)*(sin(c + d*x) + 1)/(d*e*(a**2 - b**2)*(m + 3)) + 2*b*(e*cos(c + d*x))**(-m - 1)*(a + b*sin(c + d*x))**(m + 1)/(d*e**3*(a - b)**2*(m + 1)*(m + 3)) - (e*cos(c + d*x))**(-m - 3)*(a + b*sin(c + d*x))**(m + 1)/(d*e*(a - b)*(m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_645():
    f = (e*cos(c + d*x))**(-m - 3)*(a + b*sin(c + d*x))**m
    F = (a + b*sin(c + d*x))**(m + 1)*(sin(c + d*x) - 1)*(sin(c + d*x) + 1)*sec(c + d*x)**4/(d*e**3*(e*cos(c + d*x))**m*(a - b)*(m + 2)) - ((a + b)*(sin(c + d*x) + 1)/((a - b)*(sin(c + d*x) - 1)))**(m/2 - 1)*(a + b*sin(c + d*x))**(m + 1)*(a**2*(m + 1) - b**2)*(sin(c + d*x) + 1)**3*hyper((m/2, m + 1), (m + 2,), -(2*a + 2*b*sin(c + d*x))/((a - b)*(sin(c + d*x) - 1)))*sec(c + d*x)**4/(d*e**3*m*(e*cos(c + d*x))**m*(a - b)**3*(m + 1)) + (a + b*sin(c + d*x))**(m + 1)*(a*(m + 2) - 2*b)*(sin(c + d*x) - 1)*(sin(c + d*x) + 1)**2*sec(c + d*x)**4/(d*e**3*m*(e*cos(c + d*x))**m*(a - b)**2*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_646():
    f = (e*cos(c + d*x))**(-m - 2)*(a + b*sin(c + d*x))**m
    F = 2**(sympy.S.Half - m/2)*a*(e*cos(c + d*x))**(-m - 1)*((a + b)*(sin(c + d*x) + 1)/(a + b*sin(c + d*x)))**(m/2 + sympy.S.Half)*(a + b*sin(c + d*x))**(m + 1)*hyper((-m/2 + sympy.S(-1)/2, m/2 + sympy.S.Half), (sympy.S.Half - m/2,), (1 - sin(c + d*x))*(a - b)/(2*a + 2*b*sin(c + d*x)))/(d*e*(a**2 - b**2)*(m + 1)) - (e*cos(c + d*x))**(-m - 1)*(a + b*sin(c + d*x))**(m + 1)/(d*e*(a - b)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_647():
    f = (e*cos(c + d*x))**(-m - 1)*(a + b*sin(c + d*x))**m
    F = e*(e*cos(c + d*x))**(-m - 2)*(-(1 - sin(c + d*x))*(a - b)/((a + b)*(sin(c + d*x) + 1)))**(m/2)*(1 - sin(c + d*x))*(a + b*sin(c + d*x))**(m + 1)*hyper((m + 1, m/2 + 1), (m + 2,), (2*a + 2*b*sin(c + d*x))/((a + b)*(sin(c + d*x) + 1)))/(d*(a + b)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_648():
    f = (a + b*sin(c + d*x))**m/(e*cos(c + d*x))**m
    F = e*(e*cos(c + d*x))**(-m - 1)*(1 - (a + b*sin(c + d*x))/(a - b))**(m/2 + sympy.S.Half)*(1 - (a + b*sin(c + d*x))/(a + b))**(m/2 + sympy.S.Half)*(a + b*sin(c + d*x))**(m + 1)*appellf1(m + 1, m/2 + sympy.S.Half, m/2 + sympy.S.Half, m + 2, (a + b*sin(c + d*x))/(a - b), (a + b*sin(c + d*x))/(a + b))/(b*d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_649():
    f = (e*cos(c + d*x))**(1 - m)*(a + b*sin(c + d*x))**m
    F = e*(1 - (a + b*sin(c + d*x))/(a - b))**(m/2)*(1 - (a + b*sin(c + d*x))/(a + b))**(m/2)*(a + b*sin(c + d*x))**(m + 1)*appellf1(m + 1, m/2, m/2, m + 2, (a + b*sin(c + d*x))/(a - b), (a + b*sin(c + d*x))/(a + b))/(b*d*(e*cos(c + d*x))**m*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_2_g_cos_pow_p_a_plus_b_sin_pow_m_650():
    f = (e*cos(c + d*x))**(2 - m)*(a + b*sin(c + d*x))**m
    F = e*(e*cos(c + d*x))**(1 - m)*(1 - (a + b*sin(c + d*x))/(a - b))**(m/2 + sympy.S(-1)/2)*(1 - (a + b*sin(c + d*x))/(a + b))**(m/2 + sympy.S(-1)/2)*(a + b*sin(c + d*x))**(m + 1)*appellf1(m + 1, m/2 + sympy.S(-1)/2, m/2 + sympy.S(-1)/2, m + 2, (a + b*sin(c + d*x))/(a - b), (a + b*sin(c + d*x))/(a + b))/(b*d*(m + 1))
    assert integrate(f, x) == F

