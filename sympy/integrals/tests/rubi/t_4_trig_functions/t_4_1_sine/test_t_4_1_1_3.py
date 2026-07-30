"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.1 Sine/4.1.1.3 (g tan)^p (a+b sin)^m.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, f, g, m, p = symbols('a b c d e f g m p')

def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_1():
    f = (a*sin(c + d*x) + a)*tan(c + d*x)**5
    F = a**3/(8*d*(-a*sin(c + d*x) + a)**2) + a**2/(8*d*(a*sin(c + d*x) + a)) - a**2/(d*(-a*sin(c + d*x) + a)) - 23*a*log(1 - sin(c + d*x))/(16*d) + 7*a*log(sin(c + d*x) + 1)/(16*d) - a*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_2():
    f = (a*sin(c + d*x) + a)*tan(c + d*x)**3
    F = a**2/(2*d*(-a*sin(c + d*x) + a)) + 5*a*log(1 - sin(c + d*x))/(4*d) - a*log(sin(c + d*x) + 1)/(4*d) + a*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_3():
    f = (a*sin(c + d*x) + a)*tan(c + d*x)
    F = -a*log(1 - sin(c + d*x))/d - a*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_4():
    f = (a*sin(c + d*x) + a)*cot(c + d*x)
    F = a*log(sin(c + d*x))/d + a*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_5():
    f = (a*sin(c + d*x) + a)*cot(c + d*x)**3
    F = -a*log(sin(c + d*x))/d - a*sin(c + d*x)/d - a*csc(c + d*x)**2/(2*d) - a*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_6():
    f = (a*sin(c + d*x) + a)*cot(c + d*x)**5
    F = a*log(sin(c + d*x))/d + a*sin(c + d*x)/d - a*csc(c + d*x)**4/(4*d) - a*csc(c + d*x)**3/(3*d) + a*csc(c + d*x)**2/d + 2*a*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_7():
    f = (a*sin(c + d*x) + a)*cot(c + d*x)**7
    F = -a*log(sin(c + d*x))/d - a*sin(c + d*x)/d - a*csc(c + d*x)**6/(6*d) - a*csc(c + d*x)**5/(5*d) + 3*a*csc(c + d*x)**4/(4*d) + a*csc(c + d*x)**3/d - 3*a*csc(c + d*x)**2/(2*d) - 3*a*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_8():
    f = (a*sin(c + d*x) + a)*tan(c + d*x)**6
    F = -a*x + a*cos(c + d*x)/d + a*tan(c + d*x)**5/(5*d) - a*tan(c + d*x)**3/(3*d) + a*tan(c + d*x)/d + a*sec(c + d*x)**5/(5*d) - a*sec(c + d*x)**3/d + 3*a*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_9():
    f = (a*sin(c + d*x) + a)*tan(c + d*x)**4
    F = a*x - a*cos(c + d*x)/d + a*tan(c + d*x)**3/(3*d) - a*tan(c + d*x)/d + a*sec(c + d*x)**3/(3*d) - 2*a*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_10():
    f = (a*sin(c + d*x) + a)*tan(c + d*x)**2
    F = -a*x + a*cos(c + d*x)/d + a*cos(c + d*x)/(d*(1 - sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_11():
    f = (a*sin(c + d*x) + a)*cot(c + d*x)**2
    F = -a*x + a*cos(c + d*x)/d - a*cot(c + d*x)/d - a*atanh(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_12():
    f = (a*sin(c + d*x) + a)*cot(c + d*x)**4
    F = a*x - a*cos(c + d*x)*cot(c + d*x)**2/(2*d) - 3*a*cos(c + d*x)/(2*d) - a*cot(c + d*x)**3/(3*d) + a*cot(c + d*x)/d + 3*a*atanh(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_13():
    f = (a*sin(c + d*x) + a)*cot(c + d*x)**6
    F = -a*x - a*cos(c + d*x)*cot(c + d*x)**4/(4*d) + 5*a*cos(c + d*x)*cot(c + d*x)**2/(8*d) + 15*a*cos(c + d*x)/(8*d) - a*cot(c + d*x)**5/(5*d) + a*cot(c + d*x)**3/(3*d) - a*cot(c + d*x)/d - 15*a*atanh(cos(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_14():
    f = (a*sin(c + d*x) + a)**2*tan(c + d*x)**5
    F = a**4/(4*d*(-a*sin(c + d*x) + a)**2) - 9*a**3/(4*d*(-a*sin(c + d*x) + a)) - 31*a**2*log(1 - sin(c + d*x))/(8*d) - a**2*log(sin(c + d*x) + 1)/(8*d) - a**2*sin(c + d*x)**2/(2*d) - 2*a**2*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_15():
    f = (a*sin(c + d*x) + a)**2*tan(c + d*x)**3
    F = a**3/(d*(-a*sin(c + d*x) + a)) + 3*a**2*log(1 - sin(c + d*x))/d + a**2*sin(c + d*x)**2/(2*d) + 2*a**2*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_16():
    f = (a*sin(c + d*x) + a)**2*tan(c + d*x)
    F = -2*a**2*log(1 - sin(c + d*x))/d - a**2*sin(c + d*x)**2/(2*d) - 2*a**2*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_17():
    f = (a*sin(c + d*x) + a)**2*cot(c + d*x)**3
    F = -(a*sin(c + d*x) + a)**4*csc(c + d*x)**2/(2*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_18():
    f = (a*sin(c + d*x) + a)**2*cot(c + d*x)**7
    F = 2*a**2*log(sin(c + d*x))/d - a**2*sin(c + d*x)**2/(2*d) - 2*a**2*sin(c + d*x)/d - a**2*csc(c + d*x)**6/(6*d) - 2*a**2*csc(c + d*x)**5/(5*d) + a**2*csc(c + d*x)**4/(2*d) + 2*a**2*csc(c + d*x)**3/d - 6*a**2*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_19():
    f = (a*sin(c + d*x) + a)**2*tan(c + d*x)**6
    F = -9*a**2*x/2 - a**2*sin(c + d*x)**2*tan(c + d*x)**5/(2*d) + 2*a**2*cos(c + d*x)/d + 9*a**2*tan(c + d*x)**5/(10*d) - 3*a**2*tan(c + d*x)**3/(2*d) + 9*a**2*tan(c + d*x)/(2*d) + 2*a**2*sec(c + d*x)**5/(5*d) - 2*a**2*sec(c + d*x)**3/d + 6*a**2*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_20():
    f = (a*sin(c + d*x) + a)**2*tan(c + d*x)**4
    F = a**4*sin(c + d*x)**3*cos(c + d*x)/(3*d*(-a*sin(c + d*x) + a)**2) + 7*a**2*x/2 - 7*a**2*sin(c + d*x)*cos(c + d*x)/(2*d) - 16*a**2*cos(c + d*x)/(3*d) - 8*a**2*sin(c + d*x)**2*cos(c + d*x)/(3*d*(1 - sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_21():
    f = (a*sin(c + d*x) + a)**2*tan(c + d*x)**2
    F = -5*a**2*x/2 + a**2*sin(c + d*x)*cos(c + d*x)/(2*d) + 2*a**2*cos(c + d*x)/d + 2*a**2*cos(c + d*x)/(d*(1 - sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_22():
    f = (a*sin(c + d*x) + a)**2
    F = 3*a**2*x/2 - a**2*sin(c + d*x)*cos(c + d*x)/(2*d) - 2*a**2*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_23():
    f = (a*sin(c + d*x) + a)**2*cot(c + d*x)**2
    F = -a**2*x/2 + a**2*sin(c + d*x)*cos(c + d*x)/(2*d) + 2*a**2*cos(c + d*x)/d - a**2*cot(c + d*x)/d - 2*a**2*atanh(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_24():
    f = (a*sin(c + d*x) + a)**2*cot(c + d*x)**4
    F = -a**2*x/2 - a**2*sin(c + d*x)*cos(c + d*x)/(2*d) - 2*a**2*cos(c + d*x)/d - a**2*cot(c + d*x)**3/(3*d) - a**2*cot(c + d*x)*csc(c + d*x)/d + 3*a**2*atanh(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_25():
    f = (a*sin(c + d*x) + a)**3*tan(c + d*x)**7
    F = a**6/(6*d*(-a*sin(c + d*x) + a)**3) - 13*a**5/(8*d*(-a*sin(c + d*x) + a)**2) + 71*a**4/(8*d*(-a*sin(c + d*x) + a)) + 209*a**3*log(1 - sin(c + d*x))/(16*d) - a**3*log(sin(c + d*x) + 1)/(16*d) + a**3*sin(c + d*x)**3/(3*d) + 3*a**3*sin(c + d*x)**2/(2*d) + 7*a**3*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_26():
    f = (a*sin(c + d*x) + a)**3*tan(c + d*x)**3
    F = 2*a**4/(d*(-a*sin(c + d*x) + a)) + 7*a**3*log(1 - sin(c + d*x))/d + a**3*sin(c + d*x)**3/(3*d) + 3*a**3*sin(c + d*x)**2/(2*d) + 5*a**3*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_27():
    f = (a*sin(c + d*x) + a)**3*tan(c + d*x)
    F = -4*a**3*log(1 - sin(c + d*x))/d - a**3*sin(c + d*x)**3/(3*d) - 3*a**3*sin(c + d*x)**2/(2*d) - 4*a**3*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_28():
    f = (a*sin(c + d*x) + a)**3*cot(c + d*x)**3
    F = 2*a**3*log(sin(c + d*x))/d - a**3*sin(c + d*x)**3/(3*d) - 3*a**3*sin(c + d*x)**2/(2*d) - 2*a**3*sin(c + d*x)/d - a**3*csc(c + d*x)**2/(2*d) - 3*a**3*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_29():
    f = (a*sin(c + d*x) + a)**3*tan(c + d*x)**6
    F = 23*a**6*sin(c + d*x)**3*cos(c + d*x)/(3*d*(-a**3*sin(c + d*x) + a**3)) + a**6*sin(c + d*x)**5*cos(c + d*x)/(5*d*(-a*sin(c + d*x) + a)**3) - 13*a**5*sin(c + d*x)**4*cos(c + d*x)/(15*d*(-a*sin(c + d*x) + a)**2) - 23*a**3*x/2 + 23*a**3*sin(c + d*x)*cos(c + d*x)/(2*d) - 136*a**3*cos(c + d*x)**3/(15*d) + 136*a**3*cos(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_30():
    f = (a*sin(c + d*x) + a)**3*tan(c + d*x)**4
    F = 17*a**3*x/2 - 3*a**3*sin(c + d*x)*cos(c + d*x)/(2*d) + a**3*cos(c + d*x)**3/(3*d) - 6*a**3*cos(c + d*x)/d - 25*a**3*cos(c + d*x)/(3*d*(1 - sin(c + d*x))) + 2*a**3*cos(c + d*x)/(3*d*(1 - sin(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_31():
    f = (a*sin(c + d*x) + a)**3*tan(c + d*x)**2
    F = -11*a**3*x/2 + 3*a**3*sin(c + d*x)*cos(c + d*x)/(2*d) - a**3*cos(c + d*x)**3/(3*d) + 5*a**3*cos(c + d*x)/d + 4*a**3*cos(c + d*x)/(d*(1 - sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_32():
    f = (a*sin(c + d*x) + a)**3
    F = 5*a**3*x/2 - 3*a**3*sin(c + d*x)*cos(c + d*x)/(2*d) + a**3*cos(c + d*x)**3/(3*d) - 4*a**3*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_33():
    f = (a*sin(c + d*x) + a)**3*cot(c + d*x)**2
    F = a**3*x/2 + 3*a**3*sin(c + d*x)*cos(c + d*x)/(2*d) - a**3*cos(c + d*x)**3/(3*d) + 3*a**3*cos(c + d*x)/d - a**3*cot(c + d*x)/d - 3*a**3*atanh(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_34():
    f = (a*sin(c + d*x) + a)**4*tan(c + d*x)**5
    F = a**6/(d*(-a*sin(c + d*x) + a)**2) - 11*a**5/(d*(-a*sin(c + d*x) + a)) - 25*a**4*log(1 - sin(c + d*x))/d - a**4*sin(c + d*x)**4/(4*d) - 4*a**4*sin(c + d*x)**3/(3*d) - 9*a**4*sin(c + d*x)**2/(2*d) - 16*a**4*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_35():
    f = (a*sin(c + d*x) + a)**4*tan(c + d*x)**3
    F = 4*a**5/(d*(-a*sin(c + d*x) + a)) + 16*a**4*log(1 - sin(c + d*x))/d + a**4*sin(c + d*x)**4/(4*d) + 4*a**4*sin(c + d*x)**3/(3*d) + 4*a**4*sin(c + d*x)**2/d + 12*a**4*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_36():
    f = (a*sin(c + d*x) + a)**4*tan(c + d*x)
    F = -8*a**4*log(1 - sin(c + d*x))/d - a**4*sin(c + d*x)**4/(4*d) - 4*a**4*sin(c + d*x)**3/(3*d) - 7*a**4*sin(c + d*x)**2/(2*d) - 8*a**4*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_37():
    f = (a*sin(c + d*x) + a)**4*cot(c + d*x)**3
    F = 5*a**4*log(sin(c + d*x))/d - a**4*sin(c + d*x)**4/(4*d) - 4*a**4*sin(c + d*x)**3/(3*d) - 5*a**4*sin(c + d*x)**2/(2*d) - a**4*csc(c + d*x)**2/(2*d) - 4*a**4*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_38():
    f = (a*sin(c + d*x) + a)**4*tan(c + d*x)**4
    F = 163*a**4*x/8 - a**4*sin(c + d*x)**3*cos(c + d*x)/(4*d) - 35*a**4*sin(c + d*x)*cos(c + d*x)/(8*d) + 4*a**4*cos(c + d*x)**3/(3*d) - 16*a**4*cos(c + d*x)/d - 56*a**4*cos(c + d*x)/(3*d*(1 - sin(c + d*x))) + 4*a**4*cos(c + d*x)/(3*d*(1 - sin(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_39():
    f = (a*sin(c + d*x) + a)**4*tan(c + d*x)**2
    F = -95*a**4*x/8 + a**4*sin(c + d*x)**3*cos(c + d*x)/(4*d) + 31*a**4*sin(c + d*x)*cos(c + d*x)/(8*d) - 4*a**4*cos(c + d*x)**3/(3*d) + 12*a**4*cos(c + d*x)/d + 8*a**4*cos(c + d*x)/(d*(1 - sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_40():
    f = (a*sin(c + d*x) + a)**4
    F = 35*a**4*x/8 - a**4*sin(c + d*x)**3*cos(c + d*x)/(4*d) - 27*a**4*sin(c + d*x)*cos(c + d*x)/(8*d) + 4*a**4*cos(c + d*x)**3/(3*d) - 8*a**4*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_41():
    f = (a*sin(c + d*x) + a)**4*cot(c + d*x)**2
    F = 17*a**4*x/8 + a**4*sin(c + d*x)**3*cos(c + d*x)/(4*d) + 23*a**4*sin(c + d*x)*cos(c + d*x)/(8*d) - 4*a**4*cos(c + d*x)**3/(3*d) + 4*a**4*cos(c + d*x)/d - a**4*cot(c + d*x)/d - 4*a**4*atanh(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_42():
    f = (a*sin(c + d*x) + a)**4*cot(c + d*x)**4
    F = -61*a**4*x/8 - a**4*sin(c + d*x)**3*cos(c + d*x)/(4*d) - 19*a**4*sin(c + d*x)*cos(c + d*x)/(8*d) + 4*a**4*cos(c + d*x)**3/(3*d) - a**4*cot(c + d*x)**3/(3*d) - 2*a**4*cot(c + d*x)*csc(c + d*x)/d - 5*a**4*cot(c + d*x)/d + 2*a**4*atanh(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_43():
    f = (a*sin(c + d*x) + a)**4*cot(c + d*x)**6
    F = 97*a**4*x/8 + a**4*sin(c + d*x)**3*cos(c + d*x)/(4*d) + 15*a**4*sin(c + d*x)*cos(c + d*x)/(8*d) - 4*a**4*cos(c + d*x)**3/(3*d) - 4*a**4*cos(c + d*x)/d - a**4*cot(c + d*x)**5/(5*d) - 5*a**4*cot(c + d*x)**3/(3*d) - a**4*cot(c + d*x)*csc(c + d*x)**3/d + 5*a**4*cot(c + d*x)*csc(c + d*x)/(2*d) + 10*a**4*cot(c + d*x)/d + 5*a**4*atanh(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_44():
    f = tan(c + d*x)**7/(a*sin(c + d*x) + a)
    F = tan(c + d*x)**8/(8*a*d) - tan(c + d*x)**7*sec(c + d*x)/(8*a*d) + 7*tan(c + d*x)**5*sec(c + d*x)/(48*a*d) - 35*tan(c + d*x)**3*sec(c + d*x)/(192*a*d) + 35*tan(c + d*x)*sec(c + d*x)/(128*a*d) - 35*atanh(sin(c + d*x))/(128*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_45():
    f = tan(c + d*x)**5/(a*sin(c + d*x) + a)
    F = tan(c + d*x)**6/(6*a*d) - tan(c + d*x)**5*sec(c + d*x)/(6*a*d) + 5*tan(c + d*x)**3*sec(c + d*x)/(24*a*d) - 5*tan(c + d*x)*sec(c + d*x)/(16*a*d) + 5*atanh(sin(c + d*x))/(16*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_46():
    f = tan(c + d*x)**3/(a*sin(c + d*x) + a)
    F = tan(c + d*x)**4/(4*a*d) - tan(c + d*x)**3*sec(c + d*x)/(4*a*d) + 3*tan(c + d*x)*sec(c + d*x)/(8*a*d) - 3*atanh(sin(c + d*x))/(8*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_47():
    f = cot(c + d*x)/(a*sin(c + d*x) + a)
    F = -log(sin(c + d*x) + 1)/(a*d) + log(sin(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_48():
    f = cot(c + d*x)**3/(a*sin(c + d*x) + a)
    F = -csc(c + d*x)**2/(2*a*d) + csc(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_49():
    f = cot(c + d*x)**5/(a*sin(c + d*x) + a)
    F = -cot(c + d*x)**4/(4*a*d) + csc(c + d*x)**3/(3*a*d) - csc(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_50():
    f = cot(c + d*x)**7/(a*sin(c + d*x) + a)
    F = -cot(c + d*x)**6/(6*a*d) + csc(c + d*x)**5/(5*a*d) - 2*csc(c + d*x)**3/(3*a*d) + csc(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_51():
    f = cot(c + d*x)**9/(a*sin(c + d*x) + a)
    F = -cot(c + d*x)**8/(8*a*d) + csc(c + d*x)**7/(7*a*d) - 3*csc(c + d*x)**5/(5*a*d) + csc(c + d*x)**3/(a*d) - csc(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_52():
    f = tan(c + d*x)**6/(a*sin(c + d*x) + a)
    F = tan(c + d*x)**7/(7*a*d) - sec(c + d*x)**7/(7*a*d) + 3*sec(c + d*x)**5/(5*a*d) - sec(c + d*x)**3/(a*d) + sec(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_53():
    f = tan(c + d*x)**4/(a*sin(c + d*x) + a)
    F = tan(c + d*x)**5/(5*a*d) - sec(c + d*x)**5/(5*a*d) + 2*sec(c + d*x)**3/(3*a*d) - sec(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_54():
    f = tan(c + d*x)**2/(a*sin(c + d*x) + a)
    F = tan(c + d*x)**3/(3*a*d) - sec(c + d*x)**3/(3*a*d) + sec(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_55():
    f = 1/(a*sin(c + d*x) + a)
    F = -cos(c + d*x)/(d*(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_56():
    f = cot(c + d*x)**2/(a*sin(c + d*x) + a)
    F = -cot(c + d*x)/(a*d) + atanh(cos(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_57():
    f = cot(c + d*x)**4/(a*sin(c + d*x) + a)
    F = -cot(c + d*x)**3/(3*a*d) + cot(c + d*x)*csc(c + d*x)/(2*a*d) - atanh(cos(c + d*x))/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_58():
    f = cot(c + d*x)**6/(a*sin(c + d*x) + a)
    F = -cot(c + d*x)**5/(5*a*d) + cot(c + d*x)**3*csc(c + d*x)/(4*a*d) - 3*cot(c + d*x)*csc(c + d*x)/(8*a*d) + 3*atanh(cos(c + d*x))/(8*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_59():
    f = cot(c + d*x)**8/(a*sin(c + d*x) + a)
    F = -cot(c + d*x)**7/(7*a*d) + cot(c + d*x)**5*csc(c + d*x)/(6*a*d) - 5*cot(c + d*x)**3*csc(c + d*x)/(24*a*d) + 5*cot(c + d*x)*csc(c + d*x)/(16*a*d) - 5*atanh(cos(c + d*x))/(16*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_60():
    f = tan(c + d*x)**7/(a*sin(c + d*x) + a)**2
    F = a**3/(80*d*(a*sin(c + d*x) + a)**5) - 5*a**2/(64*d*(a*sin(c + d*x) + a)**4) + 19*a/(96*d*(a*sin(c + d*x) + a)**3) + a/(192*d*(-a*sin(c + d*x) + a)**3) + 35/(256*d*(a**2*sin(c + d*x) + a**2)) + 21/(256*d*(-a**2*sin(c + d*x) + a**2)) - 1/(4*d*(a*sin(c + d*x) + a)**2) - 1/(32*d*(-a*sin(c + d*x) + a)**2) - 7*atanh(sin(c + d*x))/(128*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_61():
    f = tan(c + d*x)**5/(a*sin(c + d*x) + a)**2
    F = a**2/(32*d*(a*sin(c + d*x) + a)**4) - 7*a/(48*d*(a*sin(c + d*x) + a)**3) - 5/(32*d*(a**2*sin(c + d*x) + a**2)) - 5/(64*d*(-a**2*sin(c + d*x) + a**2)) + 1/(4*d*(a*sin(c + d*x) + a)**2) + 1/(64*d*(-a*sin(c + d*x) + a)**2) + 5*atanh(sin(c + d*x))/(64*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_62():
    f = tan(c + d*x)**3/(a*sin(c + d*x) + a)**2
    F = a/(12*d*(a*sin(c + d*x) + a)**3) + 3/(16*d*(a**2*sin(c + d*x) + a**2)) + 1/(16*d*(-a**2*sin(c + d*x) + a**2)) - 1/(4*d*(a*sin(c + d*x) + a)**2) - atanh(sin(c + d*x))/(8*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_63():
    f = tan(c + d*x)/(a*sin(c + d*x) + a)**2
    F = -1/(4*d*(a**2*sin(c + d*x) + a**2)) + 1/(4*d*(a*sin(c + d*x) + a)**2) + atanh(sin(c + d*x))/(4*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_64():
    f = cot(c + d*x)/(a*sin(c + d*x) + a)**2
    F = 1/(d*(a**2*sin(c + d*x) + a**2)) - log(sin(c + d*x) + 1)/(a**2*d) + log(sin(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_65():
    f = cot(c + d*x)**3/(a*sin(c + d*x) + a)**2
    F = -2*log(sin(c + d*x) + 1)/(a**2*d) + 2*log(sin(c + d*x))/(a**2*d) - csc(c + d*x)**2/(2*a**2*d) + 2*csc(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_66():
    f = cot(c + d*x)**5/(a*sin(c + d*x) + a)**2
    F = -csc(c + d*x)**4/(4*a**2*d) + 2*csc(c + d*x)**3/(3*a**2*d) - csc(c + d*x)**2/(2*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_67():
    f = cot(c + d*x)**7/(a*sin(c + d*x) + a)**2
    F = -csc(c + d*x)**6/(6*a**2*d) + 2*csc(c + d*x)**5/(5*a**2*d) - 2*csc(c + d*x)**3/(3*a**2*d) + csc(c + d*x)**2/(2*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_68():
    f = cot(c + d*x)**9/(a*sin(c + d*x) + a)**2
    F = -csc(c + d*x)**8/(8*a**2*d) + 2*csc(c + d*x)**7/(7*a**2*d) + csc(c + d*x)**6/(6*a**2*d) - 4*csc(c + d*x)**5/(5*a**2*d) + csc(c + d*x)**4/(4*a**2*d) + 2*csc(c + d*x)**3/(3*a**2*d) - csc(c + d*x)**2/(2*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_69():
    f = cot(c + d*x)**11/(a*sin(c + d*x) + a)**2
    F = -csc(c + d*x)**10/(10*a**2*d) + 2*csc(c + d*x)**9/(9*a**2*d) + csc(c + d*x)**8/(4*a**2*d) - 6*csc(c + d*x)**7/(7*a**2*d) + 6*csc(c + d*x)**5/(5*a**2*d) - csc(c + d*x)**4/(2*a**2*d) - 2*csc(c + d*x)**3/(3*a**2*d) + csc(c + d*x)**2/(2*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_70():
    f = cot(c + d*x)**13/(a*sin(c + d*x) + a)**2
    F = -csc(c + d*x)**12/(12*a**2*d) + 2*csc(c + d*x)**11/(11*a**2*d) + 3*csc(c + d*x)**10/(10*a**2*d) - 8*csc(c + d*x)**9/(9*a**2*d) - csc(c + d*x)**8/(4*a**2*d) + 12*csc(c + d*x)**7/(7*a**2*d) - csc(c + d*x)**6/(3*a**2*d) - 8*csc(c + d*x)**5/(5*a**2*d) + 3*csc(c + d*x)**4/(4*a**2*d) + 2*csc(c + d*x)**3/(3*a**2*d) - csc(c + d*x)**2/(2*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_71():
    f = tan(c + d*x)**5/(a*sin(c + d*x) + a)**3
    F = a**2/(40*d*(a*sin(c + d*x) + a)**5) - 7*a/(64*d*(a*sin(c + d*x) + a)**4) - 5/(128*d*(a**3*sin(c + d*x) + a**3)) - 1/(32*d*(-a**3*sin(c + d*x) + a**3)) + 1/(6*d*(a*sin(c + d*x) + a)**3) - 5/(64*a*d*(a*sin(c + d*x) + a)**2) + 1/(128*a*d*(-a*sin(c + d*x) + a)**2) + atanh(sin(c + d*x))/(128*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_72():
    f = tan(c + d*x)**3/(a*sin(c + d*x) + a)**3
    F = a/(16*d*(a*sin(c + d*x) + a)**4) + 1/(16*d*(a**3*sin(c + d*x) + a**3)) + 1/(32*d*(-a**3*sin(c + d*x) + a**3)) - 1/(6*d*(a*sin(c + d*x) + a)**3) + 3/(32*a*d*(a*sin(c + d*x) + a)**2) - atanh(sin(c + d*x))/(32*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_73():
    f = tan(c + d*x)/(a*sin(c + d*x) + a)**3
    F = -1/(8*d*(a**3*sin(c + d*x) + a**3)) + 1/(6*d*(a*sin(c + d*x) + a)**3) - 1/(8*a*d*(a*sin(c + d*x) + a)**2) + atanh(sin(c + d*x))/(8*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_74():
    f = cot(c + d*x)/(a*sin(c + d*x) + a)**3
    F = 1/(d*(a**3*sin(c + d*x) + a**3)) + 1/(2*a*d*(a*sin(c + d*x) + a)**2) - log(sin(c + d*x) + 1)/(a**3*d) + log(sin(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_75():
    f = cot(c + d*x)**3/(a*sin(c + d*x) + a)**3
    F = 2/(d*(a**3*sin(c + d*x) + a**3)) - 5*log(sin(c + d*x) + 1)/(a**3*d) + 5*log(sin(c + d*x))/(a**3*d) - csc(c + d*x)**2/(2*a**3*d) + 3*csc(c + d*x)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_76():
    f = cot(c + d*x)**5/(a*sin(c + d*x) + a)**3
    F = -4*log(sin(c + d*x) + 1)/(a**3*d) + 4*log(sin(c + d*x))/(a**3*d) - csc(c + d*x)**4/(4*a**3*d) + csc(c + d*x)**3/(a**3*d) - 2*csc(c + d*x)**2/(a**3*d) + 4*csc(c + d*x)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_77():
    f = cot(c + d*x)**7/(a*sin(c + d*x) + a)**3
    F = -csc(c + d*x)**6/(6*a**3*d) + 3*csc(c + d*x)**5/(5*a**3*d) - 3*csc(c + d*x)**4/(4*a**3*d) + csc(c + d*x)**3/(3*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_78():
    f = cot(c + d*x)**9/(a*sin(c + d*x) + a)**3
    F = -csc(c + d*x)**8/(8*a**3*d) + 3*csc(c + d*x)**7/(7*a**3*d) - csc(c + d*x)**6/(3*a**3*d) - 2*csc(c + d*x)**5/(5*a**3*d) + 3*csc(c + d*x)**4/(4*a**3*d) - csc(c + d*x)**3/(3*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_79():
    f = cot(c + d*x)**11/(a*sin(c + d*x) + a)**3
    F = -csc(c + d*x)**10/(10*a**3*d) + csc(c + d*x)**9/(3*a**3*d) - csc(c + d*x)**8/(8*a**3*d) - 5*csc(c + d*x)**7/(7*a**3*d) + 5*csc(c + d*x)**6/(6*a**3*d) + csc(c + d*x)**5/(5*a**3*d) - 3*csc(c + d*x)**4/(4*a**3*d) + csc(c + d*x)**3/(3*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_80():
    f = cot(c + d*x)**13/(a*sin(c + d*x) + a)**3
    F = -csc(c + d*x)**12/(12*a**3*d) + 3*csc(c + d*x)**11/(11*a**3*d) - 8*csc(c + d*x)**9/(9*a**3*d) + 3*csc(c + d*x)**8/(4*a**3*d) + 6*csc(c + d*x)**7/(7*a**3*d) - 4*csc(c + d*x)**6/(3*a**3*d) + 3*csc(c + d*x)**4/(4*a**3*d) - csc(c + d*x)**3/(3*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_81():
    f = tan(c + d*x)**5/(a*sin(c + d*x) + a)**4
    F = a**2/(48*d*(a*sin(c + d*x) + a)**6) - 7*a/(80*d*(a*sin(c + d*x) + a)**5) - 1/(256*d*(a**4*sin(c + d*x) + a**4)) - 3/(256*d*(-a**4*sin(c + d*x) + a**4)) - 5/(256*d*(a**2*sin(c + d*x) + a**2)**2) + 1/(256*d*(-a**2*sin(c + d*x) + a**2)**2) + 1/(8*d*(a*sin(c + d*x) + a)**4) - 5/(96*a*d*(a*sin(c + d*x) + a)**3) - atanh(sin(c + d*x))/(128*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_82():
    f = tan(c + d*x)**3/(a*sin(c + d*x) + a)**4
    F = a/(20*d*(a*sin(c + d*x) + a)**5) + 1/(64*d*(a**4*sin(c + d*x) + a**4)) + 1/(64*d*(-a**4*sin(c + d*x) + a**4)) + 1/(32*d*(a**2*sin(c + d*x) + a**2)**2) - 1/(8*d*(a*sin(c + d*x) + a)**4) + 1/(16*a*d*(a*sin(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_83():
    f = tan(c + d*x)/(a*sin(c + d*x) + a)**4
    F = -1/(16*d*(a**4*sin(c + d*x) + a**4)) - 1/(16*d*(a**2*sin(c + d*x) + a**2)**2) + 1/(8*d*(a*sin(c + d*x) + a)**4) - 1/(12*a*d*(a*sin(c + d*x) + a)**3) + atanh(sin(c + d*x))/(16*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_84():
    f = cot(c + d*x)**3/(a*sin(c + d*x) + a)**4
    F = 5/(d*(a**4*sin(c + d*x) + a**4)) + 1/(d*(a**2*sin(c + d*x) + a**2)**2) - 9*log(sin(c + d*x) + 1)/(a**4*d) + 9*log(sin(c + d*x))/(a**4*d) - csc(c + d*x)**2/(2*a**4*d) + 4*csc(c + d*x)/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_85():
    f = cot(c + d*x)**7/(a*sin(c + d*x) + a)**4
    F = -8*log(sin(c + d*x) + 1)/(a**4*d) + 8*log(sin(c + d*x))/(a**4*d) - csc(c + d*x)**6/(6*a**4*d) + 4*csc(c + d*x)**5/(5*a**4*d) - 7*csc(c + d*x)**4/(4*a**4*d) + 8*csc(c + d*x)**3/(3*a**4*d) - 4*csc(c + d*x)**2/(a**4*d) + 8*csc(c + d*x)/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_86():
    f = tan(c + d*x)**2/(a*sin(c + d*x) + a)**4
    F = 8*tan(c + d*x)**9/(9*a**4*d) + 16*tan(c + d*x)**7/(7*a**4*d) + 9*tan(c + d*x)**5/(5*a**4*d) + tan(c + d*x)**3/(3*a**4*d) - 8*sec(c + d*x)**9/(9*a**4*d) + 12*sec(c + d*x)**7/(7*a**4*d) - 4*sec(c + d*x)**5/(5*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_87():
    f = cot(c + d*x)**6/(a*sin(c + d*x) + a)**4
    F = -41*cot(c + d*x)**5/(5*a**4*d) - 27*cot(c + d*x)**3/(a**4*d) + 9*cot(c + d*x)*csc(c + d*x)**3/(a**4*d) + 27*cot(c + d*x)*csc(c + d*x)/(2*a**4*d) - 40*cot(c + d*x)/(a**4*d) + 27*atanh(cos(c + d*x))/(2*a**4*d) + 8*cot(c + d*x)*csc(c + d*x)**4/(a**4*d*(sin(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_88():
    f = sqrt(a*sin(e + f*x) + a)*tan(e + f*x)**2
    F = -sqrt(2)*sqrt(a)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(2*f) + 5*sqrt(a*sin(e + f*x) + a)*sec(e + f*x)/f - 2*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*sec(e + f*x)/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_89():
    f = sqrt(a*sin(e + f*x) + a)*cot(e + f*x)**2
    F = -sqrt(a)*atanh(sqrt(a)*cos(e + f*x)/sqrt(a*sin(e + f*x) + a))/f + 3*a*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)) - sqrt(a*sin(e + f*x) + a)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_90():
    f = sqrt(a*sin(e + f*x) + a)*cot(e + f*x)**4
    F = 11*sqrt(a)*atanh(sqrt(a)*cos(e + f*x)/sqrt(a*sin(e + f*x) + a))/(8*f) - 2*a*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)) - a*cot(e + f*x)*csc(e + f*x)/(12*f*sqrt(a*sin(e + f*x) + a)) + 11*a*cot(e + f*x)/(8*f*sqrt(a*sin(e + f*x) + a)) - sqrt(a*sin(e + f*x) + a)*cot(e + f*x)*csc(e + f*x)**2/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_91():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)*tan(e + f*x)**2
    F = 11*a**2*cos(e + f*x)/(3*f*sqrt(a*sin(e + f*x) + a)) + 7*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*sec(e + f*x)/(3*f) - 2*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*sec(e + f*x)/(3*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_92():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)*cot(e + f*x)**2
    F = -3*a**(sympy.S(3)/2)*atanh(sqrt(a)*cos(e + f*x)/sqrt(a*sin(e + f*x) + a))/f + 11*a**2*cos(e + f*x)/(3*f*sqrt(a*sin(e + f*x) + a)) + 5*a*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(3*f) - (a*sin(e + f*x) + a)**(sympy.S(3)/2)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_93():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)*cot(e + f*x)**4
    F = 37*a**(sympy.S(3)/2)*atanh(sqrt(a)*cos(e + f*x)/sqrt(a*sin(e + f*x) + a))/(8*f) - 8*a**2*cos(e + f*x)/(3*f*sqrt(a*sin(e + f*x) + a)) + 29*a**2*cot(e + f*x)/(24*f*sqrt(a*sin(e + f*x) + a)) - 2*a*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(3*f) - a*sqrt(a*sin(e + f*x) + a)*cot(e + f*x)*csc(e + f*x)/(4*f) - (a*sin(e + f*x) + a)**(sympy.S(3)/2)*cot(e + f*x)*csc(e + f*x)**2/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_94():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)*tan(e + f*x)**2
    F = 124*a**3*cos(e + f*x)/(15*f*sqrt(a*sin(e + f*x) + a)) + 31*a**2*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(15*f) + 9*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*sec(e + f*x)/(5*f) - 2*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*sec(e + f*x)/(5*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_95():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)*cot(e + f*x)**2
    F = -5*a**(sympy.S(5)/2)*atanh(sqrt(a)*cos(e + f*x)/sqrt(a*sin(e + f*x) + a))/f + 49*a**3*cos(e + f*x)/(15*f*sqrt(a*sin(e + f*x) + a)) + 31*a**2*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(15*f) + 7*a*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(5*f) - (a*sin(e + f*x) + a)**(sympy.S(5)/2)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_96():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)*cot(e + f*x)**4
    F = 55*a**(sympy.S(5)/2)*atanh(sqrt(a)*cos(e + f*x)/sqrt(a*sin(e + f*x) + a))/(8*f) - 9*a**3*cos(e + f*x)/(40*f*sqrt(a*sin(e + f*x) + a)) - 16*a**2*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(15*f) + 17*a**2*sqrt(a*sin(e + f*x) + a)*cot(e + f*x)/(24*f) - 2*a*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(5*f) - 5*a*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cot(e + f*x)*csc(e + f*x)/(12*f) - (a*sin(e + f*x) + a)**(sympy.S(5)/2)*cot(e + f*x)*csc(e + f*x)**2/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_97():
    f = tan(e + f*x)**2/sqrt(a*sin(e + f*x) + a)
    F = -sec(e + f*x)/(2*f*sqrt(a*sin(e + f*x) + a)) + 3*sqrt(a*sin(e + f*x) + a)*sec(e + f*x)/(4*a*f) + 5*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(8*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_98():
    f = cot(e + f*x)**2/sqrt(a*sin(e + f*x) + a)
    F = -cot(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)) + atanh(sqrt(a)*cos(e + f*x)/sqrt(a*sin(e + f*x) + a))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_99():
    f = cot(e + f*x)**4/sqrt(a*sin(e + f*x) + a)
    F = -cot(e + f*x)*csc(e + f*x)**2/(3*f*sqrt(a*sin(e + f*x) + a)) + cot(e + f*x)*csc(e + f*x)/(12*f*sqrt(a*sin(e + f*x) + a)) + 9*cot(e + f*x)/(8*f*sqrt(a*sin(e + f*x) + a)) - 7*atanh(sqrt(a)*cos(e + f*x)/sqrt(a*sin(e + f*x) + a))/(8*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_100():
    f = tan(e + f*x)**2/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = cos(e + f*x)/(32*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - sec(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + 5*sec(e + f*x)/(8*a*f*sqrt(a*sin(e + f*x) + a)) + sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(64*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_101():
    f = cot(e + f*x)**2/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -cot(e + f*x)/(a*f*sqrt(a*sin(e + f*x) + a)) + 3*atanh(sqrt(a)*cos(e + f*x)/sqrt(a*sin(e + f*x) + a))/(a**(sympy.S(3)/2)*f) - 2*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_102():
    f = cot(e + f*x)**4/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = 11*cot(e + f*x)*csc(e + f*x)/(12*a*f*sqrt(a*sin(e + f*x) + a)) - cot(e + f*x)/(8*a*f*sqrt(a*sin(e + f*x) + a)) - sqrt(a*sin(e + f*x) + a)*cot(e + f*x)*csc(e + f*x)**2/(3*a**2*f) - atanh(sqrt(a)*cos(e + f*x)/sqrt(a*sin(e + f*x) + a))/(8*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_103():
    f = tan(e + f*x)**2/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -sec(e + f*x)/(6*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) - 11*cos(e + f*x)/(128*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + 17*sec(e + f*x)/(48*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + 11*sec(e + f*x)/(96*a**2*f*sqrt(a*sin(e + f*x) + a)) - 11*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(256*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_104():
    f = cot(e + f*x)**2/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -2*cos(e + f*x)/(a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - cot(e + f*x)/(a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + 5*atanh(sqrt(a)*cos(e + f*x)/sqrt(a*sin(e + f*x) + a))/(a**(sympy.S(5)/2)*f) - 7*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(2*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_105():
    f = cot(e + f*x)**4/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -cot(e + f*x)*csc(e + f*x)**2/(3*a**2*f*sqrt(a*sin(e + f*x) + a)) + 13*cot(e + f*x)*csc(e + f*x)/(12*a**2*f*sqrt(a*sin(e + f*x) + a)) - 19*cot(e + f*x)/(8*a**2*f*sqrt(a*sin(e + f*x) + a)) + 45*atanh(sqrt(a)*cos(e + f*x)/sqrt(a*sin(e + f*x) + a))/(8*a**(sympy.S(5)/2)*f) - 4*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_106():
    f = (a*sin(e + f*x) + a)**(sympy.S(1)/3)*tan(e + f*x)**4
    F = -3*a**2*sin(e + f*x)**2*tan(e + f*x)/(f*(-a*sin(e + f*x) + a)*(a*sin(e + f*x) + a)**(sympy.S(2)/3)) + 3*a**2*sin(e + f*x)*tan(e + f*x)/(2*f*(-a*sin(e + f*x) + a)*(a*sin(e + f*x) + a)**(sympy.S(2)/3)) + 361*(1 - sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(1)/3)*sec(e + f*x)/(63*f) + (1 - sin(e + f*x))*(361 + 361*sqrt(3))*(a*sin(e + f*x) + a)**(sympy.S(2)/3)*sec(e + f*x)/(63*f*(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (1 + sqrt(3))*(a*sin(e + f*x) + a)**(sympy.S(1)/3))) - 361*(a*sin(e + f*x) + a)**(sympy.S(1)/3)*sec(e + f*x)/(126*f) - (-142*a**2*sin(e + f*x) + 65*a**2)*sec(e + f*x)/(42*f*(-a*sin(e + f*x) + a)*(a*sin(e + f*x) + a)**(sympy.S(2)/3)) - 361*2**(sympy.S(1)/3)*3**(sympy.S(1)/4)*sqrt((2**(sympy.S(2)/3)*a**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*a**(sympy.S(1)/3)*(a*sin(e + f*x) + a)**(sympy.S(1)/3) + (a*sin(e + f*x) + a)**(sympy.S(2)/3))/(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (1 + sqrt(3))*(a*sin(e + f*x) + a)**(sympy.S(1)/3))**2)*(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (a*sin(e + f*x) + a)**(sympy.S(1)/3))*(a*sin(e + f*x) + a)**(sympy.S(2)/3)*elliptic_e(acos((2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (1 - sqrt(3))*(a*sin(e + f*x) + a)**(sympy.S(1)/3))/(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (1 + sqrt(3))*(a*sin(e + f*x) + a)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)*sec(e + f*x)/(63*a**(sympy.S(2)/3)*f*sqrt(-(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (a*sin(e + f*x) + a)**(sympy.S(1)/3))*(a*sin(e + f*x) + a)**(sympy.S(1)/3)/(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (1 + sqrt(3))*(a*sin(e + f*x) + a)**(sympy.S(1)/3))**2)) - 2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*sqrt((2**(sympy.S(2)/3)*a**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*a**(sympy.S(1)/3)*(a*sin(e + f*x) + a)**(sympy.S(1)/3) + (a*sin(e + f*x) + a)**(sympy.S(2)/3))/(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (1 + sqrt(3))*(a*sin(e + f*x) + a)**(sympy.S(1)/3))**2)*(361 - 361*sqrt(3))*(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (a*sin(e + f*x) + a)**(sympy.S(1)/3))*(a*sin(e + f*x) + a)**(sympy.S(2)/3)*elliptic_f(acos((2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (1 - sqrt(3))*(a*sin(e + f*x) + a)**(sympy.S(1)/3))/(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (1 + sqrt(3))*(a*sin(e + f*x) + a)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)*sec(e + f*x)/(378*a**(sympy.S(2)/3)*f*sqrt(-(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (a*sin(e + f*x) + a)**(sympy.S(1)/3))*(a*sin(e + f*x) + a)**(sympy.S(1)/3)/(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (1 + sqrt(3))*(a*sin(e + f*x) + a)**(sympy.S(1)/3))**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_107():
    f = (a*sin(e + f*x) + a)**(sympy.S(1)/3)*tan(e + f*x)**2
    F = -5*2**(sympy.S(5)/6)*a*(sin(e + f*x) + 1)**(sympy.S(1)/6)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(3)/2,), sympy.S.Half - sin(e + f*x)/2)/(6*f*(a*sin(e + f*x) + a)**(sympy.S(2)/3)) + 7*(a*sin(e + f*x) + a)**(sympy.S(1)/3)*sec(e + f*x)/f - 3*(a*sin(e + f*x) + a)**(sympy.S(4)/3)*sec(e + f*x)/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_108():
    f = (a*sin(e + f*x) + a)**(sympy.S(1)/3)*cot(e + f*x)**2
    F = 6*sqrt(2)*sqrt(1 - sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(7)/3)*appellf1(sympy.S(11)/6, sympy.S(-1)/2, 2, sympy.S(17)/6, sin(e + f*x)/2 + sympy.S.Half, sin(e + f*x) + 1)*sec(e + f*x)/(11*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_109():
    f = (a*sin(e + f*x) + a)**(sympy.S(1)/3)*cot(e + f*x)**4
    F = 12*sqrt(2)*sqrt(1 - sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(10)/3)*appellf1(sympy.S(17)/6, sympy.S(-3)/2, 4, sympy.S(23)/6, sin(e + f*x)/2 + sympy.S.Half, sin(e + f*x) + 1)*sec(e + f*x)/(17*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_110():
    f = tan(e + f*x)**4/(a*sin(e + f*x) + a)**(sympy.S(1)/3)
    F = 3*a**2*sin(e + f*x)**2*tan(e + f*x)/(f*(-a*sin(e + f*x) + a)*(a*sin(e + f*x) + a)**(sympy.S(4)/3)) + 3*a**2*sin(e + f*x)*tan(e + f*x)/(4*f*(-a*sin(e + f*x) + a)*(a*sin(e + f*x) + a)**(sympy.S(4)/3)) - 973*(1 - sin(e + f*x))*sec(e + f*x)/(495*f*(a*sin(e + f*x) + a)**(sympy.S(1)/3)) + 973*sec(e + f*x)/(396*f*(a*sin(e + f*x) + a)**(sympy.S(1)/3)) - (356*a*sin(e + f*x) + 95*a)*sec(e + f*x)/(132*f*(1 - sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(4)/3)) + 973*2**(sympy.S(2)/3)*3**(sympy.S(3)/4)*sqrt((2**(sympy.S(2)/3)*a**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*a**(sympy.S(1)/3)*(a*sin(e + f*x) + a)**(sympy.S(1)/3) + (a*sin(e + f*x) + a)**(sympy.S(2)/3))/(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (1 + sqrt(3))*(a*sin(e + f*x) + a)**(sympy.S(1)/3))**2)*(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (a*sin(e + f*x) + a)**(sympy.S(1)/3))*(a*sin(e + f*x) + a)**(sympy.S(2)/3)*elliptic_f(acos((2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (1 - sqrt(3))*(a*sin(e + f*x) + a)**(sympy.S(1)/3))/(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (1 + sqrt(3))*(a*sin(e + f*x) + a)**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)*sec(e + f*x)/(2970*a**(sympy.S(4)/3)*f*sqrt(-(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (a*sin(e + f*x) + a)**(sympy.S(1)/3))*(a*sin(e + f*x) + a)**(sympy.S(1)/3)/(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (1 + sqrt(3))*(a*sin(e + f*x) + a)**(sympy.S(1)/3))**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_111():
    f = tan(e + f*x)**2/(a*sin(e + f*x) + a)**(sympy.S(1)/3)
    F = -3*sec(e + f*x)/(5*f*(a*sin(e + f*x) + a)**(sympy.S(1)/3)) + 11*2**(sympy.S(1)/6)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(3)/2,), sympy.S.Half - sin(e + f*x)/2)/(15*f*(a*sin(e + f*x) + a)**(sympy.S(1)/3)*(sin(e + f*x) + 1)**(sympy.S(1)/6)) + 4*(a*sin(e + f*x) + a)**(sympy.S(2)/3)*sec(e + f*x)/(5*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_112():
    f = cot(e + f*x)**2/(a*sin(e + f*x) + a)**(sympy.S(1)/3)
    F = 6*sqrt(2)*sqrt(1 - sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/3)*appellf1(sympy.S(7)/6, sympy.S(-1)/2, 2, sympy.S(13)/6, sin(e + f*x)/2 + sympy.S.Half, sin(e + f*x) + 1)*sec(e + f*x)/(7*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_113():
    f = cot(e + f*x)**4/(a*sin(e + f*x) + a)**(sympy.S(1)/3)
    F = 12*sqrt(2)*sqrt(1 - sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(8)/3)*appellf1(sympy.S(13)/6, sympy.S(-3)/2, 4, sympy.S(19)/6, sin(e + f*x)/2 + sympy.S.Half, sin(e + f*x) + 1)*sec(e + f*x)/(13*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_114():
    f = (g*tan(e + f*x))**p*(a*sin(e + f*x) + a)**3
    F = a**3*(g*tan(e + f*x))**(p + 1)*(cos(e + f*x)**2)**(p/2 + sympy.S.Half)*sin(e + f*x)**3*hyper((p/2 + sympy.S.Half, p/2 + 2), (p/2 + 3,), sin(e + f*x)**2)/(f*g*(p + 4)) + 3*a**3*(g*tan(e + f*x))**(p + 1)*(cos(e + f*x)**2)**(p/2 + sympy.S.Half)*sin(e + f*x)*hyper((p/2 + sympy.S.Half, p/2 + 1), (p/2 + 2,), sin(e + f*x)**2)/(f*g*(p + 2)) + a**3*(g*tan(e + f*x))**(p + 1)*hyper((1, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*g*(p + 1)) + 3*a**3*(g*tan(e + f*x))**(p + 3)*hyper((2, p/2 + sympy.S(3)/2), (p/2 + sympy.S(5)/2,), -tan(e + f*x)**2)/(f*g**3*(p + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_115():
    f = (g*tan(e + f*x))**p*(a*sin(e + f*x) + a)**2
    F = 2*a**2*(g*tan(e + f*x))**(p + 1)*(cos(e + f*x)**2)**(p/2 + sympy.S.Half)*sin(e + f*x)*hyper((p/2 + sympy.S.Half, p/2 + 1), (p/2 + 2,), sin(e + f*x)**2)/(f*g*(p + 2)) + a**2*(g*tan(e + f*x))**(p + 1)*hyper((1, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*g*(p + 1)) + a**2*(g*tan(e + f*x))**(p + 3)*hyper((2, p/2 + sympy.S(3)/2), (p/2 + sympy.S(5)/2,), -tan(e + f*x)**2)/(f*g**3*(p + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_116():
    f = (g*tan(e + f*x))**p*(a*sin(e + f*x) + a)
    F = a*(g*tan(e + f*x))**(p + 1)*(cos(e + f*x)**2)**(p/2 + sympy.S.Half)*sin(e + f*x)*hyper((p/2 + sympy.S.Half, p/2 + 1), (p/2 + 2,), sin(e + f*x)**2)/(f*g*(p + 2)) + a*(g*tan(e + f*x))**(p + 1)*hyper((1, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*g*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_117():
    f = (g*tan(e + f*x))**p/(a*sin(e + f*x) + a)
    F = (g*tan(e + f*x))**(p + 1)/(a*f*g*(p + 1)) - (g*tan(e + f*x))**(p + 2)*(cos(e + f*x)**2)**(p/2 + sympy.S(3)/2)*hyper((p/2 + 1, p/2 + sympy.S(3)/2), (p/2 + 2,), sin(e + f*x)**2)*sec(e + f*x)/(a*f*g**2*(p + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_118():
    f = (g*tan(e + f*x))**p/(a*sin(e + f*x) + a)**2
    F = (g*tan(e + f*x))**(p + 1)/(a**2*f*g*(p + 1)) - 2*(g*tan(e + f*x))**(p + 2)*(cos(e + f*x)**2)**(p/2 + sympy.S(5)/2)*hyper((p/2 + 1, p/2 + sympy.S(5)/2), (p/2 + 2,), sin(e + f*x)**2)*sec(e + f*x)**3/(a**2*f*g**2*(p + 2)) + 2*(g*tan(e + f*x))**(p + 3)/(a**2*f*g**3*(p + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_119():
    f = (g*tan(e + f*x))**p/(a*sin(e + f*x) + a)**3
    F = (g*tan(e + f*x))**(p + 1)/(a**3*f*g*(p + 1)) - 3*(g*tan(e + f*x))**(p + 2)*(cos(e + f*x)**2)**(p/2 + sympy.S(7)/2)*hyper((p/2 + 1, p/2 + sympy.S(7)/2), (p/2 + 2,), sin(e + f*x)**2)*sec(e + f*x)**5/(a**3*f*g**2*(p + 2)) + 5*(g*tan(e + f*x))**(p + 3)/(a**3*f*g**3*(p + 3)) - (g*tan(e + f*x))**(p + 4)*(cos(e + f*x)**2)**(p/2 + sympy.S(7)/2)*hyper((p/2 + 2, p/2 + sympy.S(7)/2), (p/2 + 3,), sin(e + f*x)**2)*sec(e + f*x)**3/(a**3*f*g**4*(p + 4)) + 4*(g*tan(e + f*x))**(p + 5)/(a**3*f*g**5*(p + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_120():
    f = (g*tan(e + f*x))**p*(a*sin(e + f*x) + a)**m
    F = (g*tan(e + f*x))**(p + 1)*(1 - sin(e + f*x))**(p/2 + sympy.S.Half)*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(-m + p/2 + sympy.S.Half)*appellf1(p + 1, p/2 + sympy.S.Half, -m + p/2 + sympy.S.Half, p + 2, sin(e + f*x), -sin(e + f*x))/(f*g*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_121():
    f = (a*sin(e + f*x) + a)**m*tan(e + f*x)**3
    F = -a**2*(a*sin(e + f*x) + a)**(m - 1)*sin(e + f*x)**2/(f*m*(-a*sin(e + f*x) + a)) + a*(m + 4)*(a*sin(e + f*x) + a)**(m - 1)*hyper((1, m - 1), (m,), sin(e + f*x)/2 + sympy.S.Half)/(4*f*(1 - m)) + (a*sin(e + f*x) + a)**(m - 1)*(2*a*m*sin(e + f*x) + a*(-m**2 - 3*m + 2))/(2*f*m*(1 - m)*(1 - sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_122():
    f = (a*sin(e + f*x) + a)**m*tan(e + f*x)
    F = -(a*sin(e + f*x) + a)**m/(2*f*m) + (a*sin(e + f*x) + a)**(m + 1)*hyper((1, m + 1), (m + 2,), sin(e + f*x)/2 + sympy.S.Half)/(4*a*f*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_123():
    f = (a*sin(e + f*x) + a)**m*cot(e + f*x)
    F = -(a*sin(e + f*x) + a)**(m + 1)*hyper((1, m + 1), (m + 2,), sin(e + f*x) + 1)/(a*f*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_124():
    f = (a*sin(e + f*x) + a)**m*cot(e + f*x)**3
    F = -(2 - m)*(a*sin(e + f*x) + a)**(m + 2)*hyper((2, m + 2), (m + 3,), sin(e + f*x) + 1)/(2*a**2*f*(m + 2)) - (a*sin(e + f*x) + a)**(m + 2)*csc(e + f*x)**2/(2*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_125():
    f = (a*sin(e + f*x) + a)**m*cot(e + f*x)**5
    F = (9 - m)*(a*sin(e + f*x) + a)**(m + 3)*csc(e + f*x)**3/(12*a**3*f) - (a*sin(e + f*x) + a)**(m + 3)*csc(e + f*x)**4/(4*a**3*f) - (a*sin(e + f*x) + a)**(m + 3)*(m**2 - 9*m + 12)*hyper((3, m + 3), (m + 4,), sin(e + f*x) + 1)/(12*a**3*f*(m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_126():
    f = (a*sin(e + f*x) + a)**m*tan(e + f*x)**4
    F = 2**(m + sympy.S(-3)/2)*(1 - sin(e + f*x))*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(sympy.S.Half - m)*(m**4 + 6*m**3 - 7*m**2 - 12*m + 9)*hyper((sympy.S.Half, sympy.S(5)/2 - m), (sympy.S(3)/2,), sympy.S.Half - sin(e + f*x)/2)*sec(e + f*x)/(3*f*m*(1 - m)) + a**2*(a*sin(e + f*x) + a)**(m - 1)*sin(e + f*x)*tan(e + f*x)/(f*(1 - m)*(-a*sin(e + f*x) + a)) - a**2*(a*sin(e + f*x) + a)**(m - 1)*sin(e + f*x)**2*tan(e + f*x)/(f*m*(-a*sin(e + f*x) + a)) - (a*sin(e + f*x) + a)**(m - 1)*(-a*(-m**3 - 8*m**2 - 6*m + 9)*sin(e + f*x) + a*(-m**3 - 7*m**2 - m + 6))*sec(e + f*x)/(3*f*m*(1 - m)*(1 - sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_127():
    f = (a*sin(e + f*x) + a)**m*tan(e + f*x)**2
    F = 2**(m + sympy.S(-1)/2)*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(sympy.S.Half - m)*(-m**2 - m + 1)*hyper((sympy.S(-1)/2, sympy.S(3)/2 - m), (sympy.S.Half,), sympy.S.Half - sin(e + f*x)/2)*sec(e + f*x)/(f*m*(1 - m)) + (a*sin(e + f*x) + a)**m*sec(e + f*x)/(f*m*(1 - m)) - (a*sin(e + f*x) + a)**(m + 1)*sec(e + f*x)/(a*f*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_128():
    f = (a*sin(e + f*x) + a)**m
    F = -2**(m + sympy.S.Half)*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(-m + sympy.S(-1)/2)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - sin(e + f*x)/2)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_129():
    f = (a*sin(e + f*x) + a)**m*cot(e + f*x)**2
    F = 2*sqrt(2)*sqrt(1 - sin(e + f*x))*(a*sin(e + f*x) + a)**(m + 2)*appellf1(m + sympy.S(3)/2, sympy.S(-1)/2, 2, m + sympy.S(5)/2, sin(e + f*x)/2 + sympy.S.Half, sin(e + f*x) + 1)*sec(e + f*x)/(a**2*f*(2*m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_130():
    f = (a*sin(e + f*x) + a)**m*cot(e + f*x)**4
    F = 4*sqrt(2)*sqrt(1 - sin(e + f*x))*(a*sin(e + f*x) + a)**(m + 3)*appellf1(m + sympy.S(5)/2, sympy.S(-3)/2, 4, m + sympy.S(7)/2, sin(e + f*x)/2 + sympy.S.Half, sin(e + f*x) + 1)*sec(e + f*x)/(a**3*f*(2*m + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_131():
    f = (a + b*sin(c + d*x))*tan(c + d*x)**3
    F = 3*b*sin(c + d*x)/(2*d) + (a + b*sin(c + d*x))*tan(c + d*x)**2/(2*d) + (2*a - 3*b)*log(sin(c + d*x) + 1)/(4*d) + (2*a + 3*b)*log(1 - sin(c + d*x))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_132():
    f = (a + b*sin(c + d*x))*tan(c + d*x)
    F = -b*sin(c + d*x)/d - (a - b)*log(sin(c + d*x) + 1)/(2*d) - (a + b)*log(1 - sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_133():
    f = (a + b*sin(c + d*x))*cot(c + d*x)
    F = a*log(sin(c + d*x))/d + b*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_134():
    f = (a + b*sin(c + d*x))*cot(c + d*x)**3
    F = -a*log(sin(c + d*x))/d - a*csc(c + d*x)**2/(2*d) - b*sin(c + d*x)/d - b*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_135():
    f = (a + b*sin(c + d*x))*cot(c + d*x)**5
    F = a*log(sin(c + d*x))/d - a*csc(c + d*x)**4/(4*d) + a*csc(c + d*x)**2/d + b*sin(c + d*x)/d - b*csc(c + d*x)**3/(3*d) + 2*b*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_136():
    f = (a + b*sin(c + d*x))*tan(c + d*x)**4
    F = a*x + a*tan(c + d*x)**3/(3*d) - a*tan(c + d*x)/d - b*cos(c + d*x)/d + b*sec(c + d*x)**3/(3*d) - 2*b*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_137():
    f = (a + b*sin(c + d*x))*tan(c + d*x)**2
    F = -a*x + a*tan(c + d*x)/d + b*cos(c + d*x)/d + b*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_138():
    f = (a + b*sin(c + d*x))*cot(c + d*x)**2
    F = -a*x - a*cot(c + d*x)/d + b*cos(c + d*x)/d - b*atanh(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_139():
    f = (a + b*sin(c + d*x))*cot(c + d*x)**4
    F = a*x - a*cot(c + d*x)**3/(3*d) + a*cot(c + d*x)/d - b*cos(c + d*x)*cot(c + d*x)**2/(2*d) - 3*b*cos(c + d*x)/(2*d) + 3*b*atanh(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_140():
    f = (a + b*sin(c + d*x))*cot(c + d*x)**6
    F = -a*x - a*cot(c + d*x)**5/(5*d) + a*cot(c + d*x)**3/(3*d) - a*cot(c + d*x)/d - b*cos(c + d*x)*cot(c + d*x)**4/(4*d) + 5*b*cos(c + d*x)*cot(c + d*x)**2/(8*d) + 15*b*cos(c + d*x)/(8*d) - 15*b*atanh(cos(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_141():
    f = (a + b*sin(c + d*x))**2*tan(c + d*x)**3
    F = 2*a*b*sin(c + d*x)/d + b**2*sin(c + d*x)**2/(2*d) + (a - 2*b)*(a - b)*log(sin(c + d*x) + 1)/(2*d) + (a + b)*(a + 2*b)*log(1 - sin(c + d*x))/(2*d) + (a + b*sin(c + d*x))**2*sec(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_142():
    f = (a + b*sin(c + d*x))**2*tan(c + d*x)
    F = -2*a*b*sin(c + d*x)/d - b**2*sin(c + d*x)**2/(2*d) - (a - b)**2*log(sin(c + d*x) + 1)/(2*d) - (a + b)**2*log(1 - sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_143():
    f = (a + b*sin(c + d*x))**2*cot(c + d*x)
    F = a**2*log(sin(c + d*x))/d + 2*a*b*sin(c + d*x)/d + b**2*sin(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_144():
    f = (a + b*sin(c + d*x))**2*cot(c + d*x)**3
    F = -a**2*csc(c + d*x)**2/(2*d) - 2*a*b*sin(c + d*x)/d - 2*a*b*csc(c + d*x)/d - b**2*sin(c + d*x)**2/(2*d) - (a**2 - b**2)*log(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_145():
    f = (a + b*sin(c + d*x))**2*cot(c + d*x)**5
    F = -a**2*csc(c + d*x)**4/(4*d) + 2*a*b*sin(c + d*x)/d - 2*a*b*csc(c + d*x)**3/(3*d) + 4*a*b*csc(c + d*x)/d + b**2*sin(c + d*x)**2/(2*d) + (a**2 - 2*b**2)*log(sin(c + d*x))/d + (2*a**2 - b**2)*csc(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_146():
    f = (a + b*sin(c + d*x))**2*tan(c + d*x)**4
    F = a**2*x + a**2*tan(c + d*x)**3/(3*d) - a**2*tan(c + d*x)/d - 2*a*b*cos(c + d*x)/d + 2*a*b*sec(c + d*x)**3/(3*d) - 4*a*b*sec(c + d*x)/d + 5*b**2*x/2 - b**2*sin(c + d*x)**2*tan(c + d*x)**3/(2*d) + 5*b**2*tan(c + d*x)**3/(6*d) - 5*b**2*tan(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_147():
    f = (a + b*sin(c + d*x))**2*tan(c + d*x)**2
    F = -a**2*x + a**2*tan(c + d*x)/d + 2*a*b*cos(c + d*x)/d + 2*a*b*sec(c + d*x)/d - 3*b**2*x/2 - b**2*sin(c + d*x)**2*tan(c + d*x)/(2*d) + 3*b**2*tan(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_148():
    f = (a + b*sin(c + d*x))**2*cot(c + d*x)**2
    F = -a**2*x - a**2*cot(c + d*x)/d + 2*a*b*cos(c + d*x)/d - 2*a*b*atanh(cos(c + d*x))/d + b**2*x/2 + b**2*sin(c + d*x)*cos(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_149():
    f = (a + b*sin(c + d*x))**2*cot(c + d*x)**4
    F = a**2*x - a**2*cot(c + d*x)**3/(3*d) + a**2*cot(c + d*x)/d - a*b*cos(c + d*x)*cot(c + d*x)**2/d - 3*a*b*cos(c + d*x)/d + 3*a*b*atanh(cos(c + d*x))/d - 3*b**2*x/2 + b**2*cos(c + d*x)**2*cot(c + d*x)/(2*d) - 3*b**2*cot(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_150():
    f = (a + b*sin(c + d*x))**2*cot(c + d*x)**6
    F = -a**2*x - a**2*cot(c + d*x)**5/(5*d) + a**2*cot(c + d*x)**3/(3*d) - a**2*cot(c + d*x)/d - a*b*cos(c + d*x)*cot(c + d*x)**4/(2*d) + 5*a*b*cos(c + d*x)*cot(c + d*x)**2/(4*d) + 15*a*b*cos(c + d*x)/(4*d) - 15*a*b*atanh(cos(c + d*x))/(4*d) + 5*b**2*x/2 + b**2*cos(c + d*x)**2*cot(c + d*x)**3/(2*d) - 5*b**2*cot(c + d*x)**3/(6*d) + 5*b**2*cot(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_151():
    f = (a + b*sin(c + d*x))**3*tan(c + d*x)**3
    F = 3*a*b**2*sin(c + d*x)**2/(2*d) + b**3*sin(c + d*x)**3/(3*d) + b*(6*a**2 + 5*b**2)*sin(c + d*x)/(2*d) + (a - b)**2*(2*a - 5*b)*log(sin(c + d*x) + 1)/(4*d) + (a + b)**2*(2*a + 5*b)*log(1 - sin(c + d*x))/(4*d) + (a + b*sin(c + d*x))**3*sec(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_152():
    f = (a + b*sin(c + d*x))**3*tan(c + d*x)
    F = -3*a*b**2*sin(c + d*x)**2/(2*d) - b**3*sin(c + d*x)**3/(3*d) - b*(3*a**2 + b**2)*sin(c + d*x)/d - (a - b)**3*log(sin(c + d*x) + 1)/(2*d) - (a + b)**3*log(1 - sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_153():
    f = (a + b*sin(c + d*x))**3*cot(c + d*x)
    F = a**3*log(sin(c + d*x))/d + 3*a**2*b*sin(c + d*x)/d + 3*a*b**2*sin(c + d*x)**2/(2*d) + b**3*sin(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_154():
    f = (a + b*sin(c + d*x))**3*cot(c + d*x)**3
    F = -a**3*csc(c + d*x)**2/(2*d) - 3*a**2*b*csc(c + d*x)/d - 3*a*b**2*sin(c + d*x)**2/(2*d) - a*(a**2 - 3*b**2)*log(sin(c + d*x))/d - b**3*sin(c + d*x)**3/(3*d) - b*(3*a**2 - b**2)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_155():
    f = (a + b*sin(c + d*x))**3*cot(c + d*x)**5
    F = -a**3*csc(c + d*x)**4/(4*d) - a**2*b*csc(c + d*x)**3/d + 3*a*b**2*sin(c + d*x)**2/(2*d) + a*(a**2 - 6*b**2)*log(sin(c + d*x))/d + a*(2*a**2 - 3*b**2)*csc(c + d*x)**2/(2*d) + b**3*sin(c + d*x)**3/(3*d) + b*(3*a**2 - 2*b**2)*sin(c + d*x)/d + b*(6*a**2 - b**2)*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_156():
    f = (a + b*sin(c + d*x))**3*tan(c + d*x)**4
    F = a**3*x + a**3*tan(c + d*x)**3/(3*d) - a**3*tan(c + d*x)/d - 3*a**2*b*cos(c + d*x)/d + a**2*b*sec(c + d*x)**3/d - 6*a**2*b*sec(c + d*x)/d + 15*a*b**2*x/2 - 3*a*b**2*sin(c + d*x)**2*tan(c + d*x)**3/(2*d) + 5*a*b**2*tan(c + d*x)**3/(2*d) - 15*a*b**2*tan(c + d*x)/(2*d) + b**3*cos(c + d*x)**3/(3*d) - 3*b**3*cos(c + d*x)/d + b**3*sec(c + d*x)**3/(3*d) - 3*b**3*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_157():
    f = (a + b*sin(c + d*x))**3*tan(c + d*x)**2
    F = -a**3*x + a**3*tan(c + d*x)/d + 3*a**2*b*cos(c + d*x)/d + 3*a**2*b*sec(c + d*x)/d - 9*a*b**2*x/2 - 3*a*b**2*sin(c + d*x)**2*tan(c + d*x)/(2*d) + 9*a*b**2*tan(c + d*x)/(2*d) - b**3*cos(c + d*x)**3/(3*d) + 2*b**3*cos(c + d*x)/d + b**3*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_158():
    f = (a + b*sin(c + d*x))**3*cot(c + d*x)**2
    F = -a**3*x - a**3*cot(c + d*x)/d + 3*a**2*b*cos(c + d*x)/d - 3*a**2*b*atanh(cos(c + d*x))/d + 3*a*b**2*x/2 + 3*a*b**2*sin(c + d*x)*cos(c + d*x)/(2*d) - b**3*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_159():
    f = (a + b*sin(c + d*x))**3*cot(c + d*x)**4
    F = a**3*x - a**3*cot(c + d*x)**3/(3*d) + a**3*cot(c + d*x)/d - 3*a**2*b*cos(c + d*x)*cot(c + d*x)**2/(2*d) - 9*a**2*b*cos(c + d*x)/(2*d) + 9*a**2*b*atanh(cos(c + d*x))/(2*d) - 9*a*b**2*x/2 + 3*a*b**2*cos(c + d*x)**2*cot(c + d*x)/(2*d) - 9*a*b**2*cot(c + d*x)/(2*d) + b**3*cos(c + d*x)**3/(3*d) + b**3*cos(c + d*x)/d - b**3*atanh(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_160():
    f = (a + b*sin(c + d*x))**3*cot(c + d*x)**6
    F = -a**3*x - a**3*cot(c + d*x)**5/(5*d) + a**3*cot(c + d*x)**3/(3*d) - a**3*cot(c + d*x)/d - 3*a**2*b*cos(c + d*x)*cot(c + d*x)**4/(4*d) + 15*a**2*b*cos(c + d*x)*cot(c + d*x)**2/(8*d) + 45*a**2*b*cos(c + d*x)/(8*d) - 45*a**2*b*atanh(cos(c + d*x))/(8*d) + 15*a*b**2*x/2 + 3*a*b**2*cos(c + d*x)**2*cot(c + d*x)**3/(2*d) - 5*a*b**2*cot(c + d*x)**3/(2*d) + 15*a*b**2*cot(c + d*x)/(2*d) - b**3*cos(c + d*x)**3*cot(c + d*x)**2/(2*d) - 5*b**3*cos(c + d*x)**3/(6*d) - 5*b**3*cos(c + d*x)/(2*d) + 5*b**3*atanh(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_161():
    f = tan(c + d*x)**5/(a + b*sin(c + d*x))
    F = a**5*log(a + b*sin(c + d*x))/(d*(a**2 - b**2)**3) + (a - b*sin(c + d*x))*sec(c + d*x)**4/(d*(4*a**2 - 4*b**2)) - (4*a*(2*a**2 - b**2) - b*(9*a**2 - 5*b**2)*sin(c + d*x))*sec(c + d*x)**2/(8*d*(a**2 - b**2)**2) - (8*a**2 + 9*a*b + 3*b**2)*log(1 - sin(c + d*x))/(16*d*(a + b)**3) - (8*a**2 - 9*a*b + 3*b**2)*log(sin(c + d*x) + 1)/(16*d*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_162():
    f = tan(c + d*x)**3/(a + b*sin(c + d*x))
    F = -a**3*log(a + b*sin(c + d*x))/(d*(a**2 - b**2)**2) + (a - b*sin(c + d*x))*sec(c + d*x)**2/(d*(2*a**2 - 2*b**2)) + (2*a + b)*log(1 - sin(c + d*x))/(4*d*(a + b)**2) + (2*a - b)*log(sin(c + d*x) + 1)/(4*d*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_163():
    f = tan(c + d*x)/(a + b*sin(c + d*x))
    F = a*log(a + b*sin(c + d*x))/(d*(a**2 - b**2)) - log(1 - sin(c + d*x))/(d*(2*a + 2*b)) - log(sin(c + d*x) + 1)/(d*(2*a - 2*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_164():
    f = cot(c + d*x)/(a + b*sin(c + d*x))
    F = -log(a + b*sin(c + d*x))/(a*d) + log(sin(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_165():
    f = cot(c + d*x)**3/(a + b*sin(c + d*x))
    F = -csc(c + d*x)**2/(2*a*d) + b*csc(c + d*x)/(a**2*d) + (a**2 - b**2)*log(a + b*sin(c + d*x))/(a**3*d) - (a**2 - b**2)*log(sin(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_166():
    f = cot(c + d*x)**5/(a + b*sin(c + d*x))
    F = -csc(c + d*x)**4/(4*a*d) + b*csc(c + d*x)**3/(3*a**2*d) + (2*a**2 - b**2)*csc(c + d*x)**2/(2*a**3*d) - b*(2*a**2 - b**2)*csc(c + d*x)/(a**4*d) - (a**2 - b**2)**2*log(a + b*sin(c + d*x))/(a**5*d) + (a**2 - b**2)**2*log(sin(c + d*x))/(a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_167():
    f = tan(c + d*x)**4/(a + b*sin(c + d*x))
    F = 2*a**4*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(d*(a**2 - b**2)**(sympy.S(5)/2)) - a**3*tan(c + d*x)/(d*(a**2 - b**2)**2) + a**2*b*sec(c + d*x)/(d*(a**2 - b**2)**2) + a*tan(c + d*x)**3/(d*(3*a**2 - 3*b**2)) - b*sec(c + d*x)**3/(d*(3*a**2 - 3*b**2)) + b*sec(c + d*x)/(d*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_168():
    f = tan(c + d*x)**2/(a + b*sin(c + d*x))
    F = -2*a**2*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(d*(a**2 - b**2)**(sympy.S(3)/2)) + a*tan(c + d*x)/(d*(a**2 - b**2)) - b*sec(c + d*x)/(d*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_169():
    f = cot(c + d*x)**2/(a + b*sin(c + d*x))
    F = -cot(c + d*x)/(a*d) + b*atanh(cos(c + d*x))/(a**2*d) - 2*sqrt(a**2 - b**2)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_170():
    f = cot(c + d*x)**4/(a + b*sin(c + d*x))
    F = -cot(c + d*x)*csc(c + d*x)**2/(3*a*d) + b*cot(c + d*x)*csc(c + d*x)/(2*a**2*d) + (4*a**2 - 3*b**2)*cot(c + d*x)/(3*a**3*d) - b*(3*a**2 - 2*b**2)*atanh(cos(c + d*x))/(2*a**4*d) + 2*(a**2 - b**2)**(sympy.S(3)/2)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_171():
    f = cot(c + d*x)**6/(a + b*sin(c + d*x))
    F = a*cot(c + d*x)*csc(c + d*x)**2/(2*b**2*d) - cot(c + d*x)*csc(c + d*x)/(b*d) - cot(c + d*x)*csc(c + d*x)**4/(5*a*d) + b*cot(c + d*x)*csc(c + d*x)**3/(4*a**2*d) - (15*a**4 - 22*a**2*b**2 + 10*b**4)*cot(c + d*x)*csc(c + d*x)**2/(30*a**3*b**2*d) + (8*a**4 - 9*a**2*b**2 + 4*b**4)*cot(c + d*x)*csc(c + d*x)/(8*a**4*b*d) - (23*a**4 - 35*a**2*b**2 + 15*b**4)*cot(c + d*x)/(15*a**5*d) + b*(15*a**4 - 20*a**2*b**2 + 8*b**4)*atanh(cos(c + d*x))/(8*a**6*d) - 2*(a**2 - b**2)**(sympy.S(5)/2)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(a**6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_172():
    f = tan(c + d*x)**5/(a + b*sin(c + d*x))**2
    F = -a**5/(d*(a + b*sin(c + d*x))*(a**2 - b**2)**3) + a**4*(a**2 + 5*b**2)*log(a + b*sin(c + d*x))/(d*(a**2 - b**2)**4) - a*(4*a + b)*log(1 - sin(c + d*x))/(8*d*(a + b)**4) - a*(4*a - b)*log(sin(c + d*x) + 1)/(8*d*(a - b)**4) + (a**2 - 2*a*b*sin(c + d*x) + b**2)*sec(c + d*x)**4/(4*d*(a**2 - b**2)**2) - (4*a**4 + 6*a**2*b**2 - a*b*(9*a**2 - b**2)*sin(c + d*x) - 2*b**4)*sec(c + d*x)**2/(4*d*(a**2 - b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_173():
    f = tan(c + d*x)**3/(a + b*sin(c + d*x))**2
    F = a**3/(d*(a + b*sin(c + d*x))*(a**2 - b**2)**2) - a**2*(a**2 + 3*b**2)*log(a + b*sin(c + d*x))/(d*(a**2 - b**2)**3) + a*log(1 - sin(c + d*x))/(2*d*(a + b)**3) + a*log(sin(c + d*x) + 1)/(2*d*(a - b)**3) + (a**2 - 2*a*b*sin(c + d*x) + b**2)*sec(c + d*x)**2/(2*d*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_174():
    f = tan(c + d*x)/(a + b*sin(c + d*x))**2
    F = -a/(d*(a + b*sin(c + d*x))*(a**2 - b**2)) + (a**2 + b**2)*log(a + b*sin(c + d*x))/(d*(a**2 - b**2)**2) - log(1 - sin(c + d*x))/(2*d*(a + b)**2) - log(sin(c + d*x) + 1)/(2*d*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_175():
    f = cot(c + d*x)/(a + b*sin(c + d*x))**2
    F = 1/(a*d*(a + b*sin(c + d*x))) - log(a + b*sin(c + d*x))/(a**2*d) + log(sin(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_176():
    f = cot(c + d*x)**3/(a + b*sin(c + d*x))**2
    F = -csc(c + d*x)**2/(2*a**2*d) + 2*b*csc(c + d*x)/(a**3*d) - (a**2 - b**2)/(a**3*d*(a + b*sin(c + d*x))) + (a**2 - 3*b**2)*log(a + b*sin(c + d*x))/(a**4*d) - (a**2 - 3*b**2)*log(sin(c + d*x))/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_177():
    f = cot(c + d*x)**5/(a + b*sin(c + d*x))**2
    F = -csc(c + d*x)**4/(4*a**2*d) + 2*b*csc(c + d*x)**3/(3*a**3*d) + (2*a**2 - 3*b**2)*csc(c + d*x)**2/(2*a**4*d) - 4*b*(a**2 - b**2)*csc(c + d*x)/(a**5*d) + (a**2 - b**2)**2/(a**5*d*(a + b*sin(c + d*x))) - (a**4 - 6*a**2*b**2 + 5*b**4)*log(a + b*sin(c + d*x))/(a**6*d) + (a**4 - 6*a**2*b**2 + 5*b**4)*log(sin(c + d*x))/(a**6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_178():
    f = tan(c + d*x)**4/(a + b*sin(c + d*x))**2
    F = 2*a**5*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(d*(a**2 - b**2)**(sympy.S(7)/2)) + a**4*b*cos(c + d*x)/(d*(a + b*sin(c + d*x))*(a**2 - b**2)**3) + 8*a**3*b**2*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(d*(a**2 - b**2)**(sympy.S(7)/2)) - cos(c + d*x)/(12*d*(a - b)**2*(sin(c + d*x) + 1)) - cos(c + d*x)/(12*d*(a - b)**2*(sin(c + d*x) + 1)**2) + (3*a - b)*cos(c + d*x)/(4*d*(a - b)**3*(sin(c + d*x) + 1)) + cos(c + d*x)/(12*d*(1 - sin(c + d*x))*(a + b)**2) - (3*a + b)*cos(c + d*x)/(4*d*(1 - sin(c + d*x))*(a + b)**3) + cos(c + d*x)/(12*d*(1 - sin(c + d*x))**2*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_179():
    f = tan(c + d*x)**2/(a + b*sin(c + d*x))**2
    F = -2*a**3*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(d*(a**2 - b**2)**(sympy.S(5)/2)) - a**2*b*cos(c + d*x)/(d*(a + b*sin(c + d*x))*(a**2 - b**2)**2) - 4*a*b**2*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(d*(a**2 - b**2)**(sympy.S(5)/2)) - cos(c + d*x)/(2*d*(a - b)**2*(sin(c + d*x) + 1)) + cos(c + d*x)/(2*d*(1 - sin(c + d*x))*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_180():
    f = cot(c + d*x)**2/(a + b*sin(c + d*x))**2
    F = cot(c + d*x)/(a*d*(a + b*sin(c + d*x))) - 2*cot(c + d*x)/(a**2*d) + 2*b*atanh(cos(c + d*x))/(a**3*d) - (2*a**2 - 4*b**2)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(a**3*d*sqrt(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_181():
    f = cot(c + d*x)**4/(a + b*sin(c + d*x))**2
    F = -cot(c + d*x)*csc(c + d*x)**2/(3*a*d*(a + b*sin(c + d*x))) + (3*a**2 - 4*b**2)*cot(c + d*x)*csc(c + d*x)/(3*a**2*b*d*(a + b*sin(c + d*x))) - (a**2 - 2*b**2)*cot(c + d*x)*csc(c + d*x)/(a**3*b*d) + (7*a**2 - 12*b**2)*cot(c + d*x)/(3*a**4*d) - b*(3*a**2 - 4*b**2)*atanh(cos(c + d*x))/(a**5*d) + (2*a**4 - 10*a**2*b**2 + 8*b**4)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(a**5*d*sqrt(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_182():
    f = cot(c + d*x)**6/(a + b*sin(c + d*x))**2
    F = a*cot(c + d*x)*csc(c + d*x)**2/(6*b**2*d*(a + b*sin(c + d*x))) - cot(c + d*x)*csc(c + d*x)/(2*b*d*(a + b*sin(c + d*x))) - cot(c + d*x)*csc(c + d*x)**4/(5*a*d*(a + b*sin(c + d*x))) + 3*b*cot(c + d*x)*csc(c + d*x)**3/(10*a**2*d*(a + b*sin(c + d*x))) + (2*a**4 - 12*a**2*b**2 + 9*b**4)*cot(c + d*x)*csc(c + d*x)**2/(6*a**3*b**2*d*(a + b*sin(c + d*x))) - (15*a**4 - 82*a**2*b**2 + 60*b**4)*cot(c + d*x)*csc(c + d*x)**2/(30*a**4*b**2*d) + (4*a**4 - 17*a**2*b**2 + 12*b**4)*cot(c + d*x)*csc(c + d*x)/(4*a**5*b*d) - (38*a**4 - 135*a**2*b**2 + 90*b**4)*cot(c + d*x)/(15*a**6*d) + b*(15*a**4 - 40*a**2*b**2 + 24*b**4)*atanh(cos(c + d*x))/(4*a**7*d) - (a**2 - b**2)**(sympy.S(3)/2)*(2*a**2 - 12*b**2)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_183():
    f = tan(c + d*x)**5/(a + b*sin(c + d*x))**3
    F = -a**5/(2*d*(a + b*sin(c + d*x))**2*(a**2 - b**2)**3) - a**4*(a**2 + 5*b**2)/(d*(a + b*sin(c + d*x))*(a**2 - b**2)**4) + a**3*(a**4 + 13*a**2*b**2 + 10*b**4)*log(a + b*sin(c + d*x))/(d*(a**2 - b**2)**5) + (a*(a**2 + 3*b**2) - b*(3*a**2 + b**2)*sin(c + d*x))*sec(c + d*x)**4/(4*d*(a**2 - b**2)**3) - (8*a**3*(a**2 + 5*b**2) - b*(27*a**4 + 22*a**2*b**2 - b**4)*sin(c + d*x))*sec(c + d*x)**2/(8*d*(a**2 - b**2)**4) - (8*a**2 - 5*a*b - b**2)*log(1 - sin(c + d*x))/(16*d*(a + b)**5) - (8*a**2 + 5*a*b - b**2)*log(sin(c + d*x) + 1)/(16*d*(a - b)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_184():
    f = tan(c + d*x)**3/(a + b*sin(c + d*x))**3
    F = a**3/(2*d*(a + b*sin(c + d*x))**2*(a**2 - b**2)**2) + a**2*(a**2 + 3*b**2)/(d*(a + b*sin(c + d*x))*(a**2 - b**2)**3) - a*(a**4 + 8*a**2*b**2 + 3*b**4)*log(a + b*sin(c + d*x))/(d*(a**2 - b**2)**4) + (a*(a**2 + 3*b**2) - b*(3*a**2 + b**2)*sin(c + d*x))*sec(c + d*x)**2/(2*d*(a**2 - b**2)**3) + (2*a - b)*log(1 - sin(c + d*x))/(4*d*(a + b)**4) + (2*a + b)*log(sin(c + d*x) + 1)/(4*d*(a - b)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_185():
    f = tan(c + d*x)/(a + b*sin(c + d*x))**3
    F = a*(a**2 + 3*b**2)*log(a + b*sin(c + d*x))/(d*(a**2 - b**2)**3) - a/(d*(a + b*sin(c + d*x))**2*(2*a**2 - 2*b**2)) - (a**2 + b**2)/(d*(a + b*sin(c + d*x))*(a**2 - b**2)**2) - log(1 - sin(c + d*x))/(2*d*(a + b)**3) - log(sin(c + d*x) + 1)/(2*d*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_186():
    f = cot(c + d*x)/(a + b*sin(c + d*x))**3
    F = 1/(2*a*d*(a + b*sin(c + d*x))**2) + 1/(a**2*d*(a + b*sin(c + d*x))) - log(a + b*sin(c + d*x))/(a**3*d) + log(sin(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_187():
    f = cot(c + d*x)**3/(a + b*sin(c + d*x))**3
    F = -csc(c + d*x)**2/(2*a**3*d) - (a**2 - b**2)/(2*a**3*d*(a + b*sin(c + d*x))**2) + 3*b*csc(c + d*x)/(a**4*d) - (a**2 - 3*b**2)/(a**4*d*(a + b*sin(c + d*x))) + (a**2 - 6*b**2)*log(a + b*sin(c + d*x))/(a**5*d) - (a**2 - 6*b**2)*log(sin(c + d*x))/(a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_188():
    f = cot(c + d*x)**5/(a + b*sin(c + d*x))**3
    F = -csc(c + d*x)**4/(4*a**3*d) + b*csc(c + d*x)**3/(a**4*d) + (a**2 - 3*b**2)*csc(c + d*x)**2/(a**5*d) + (a**2 - b**2)**2/(2*a**5*d*(a + b*sin(c + d*x))**2) - 2*b*(3*a**2 - 5*b**2)*csc(c + d*x)/(a**6*d) + (a**4 - 6*a**2*b**2 + 5*b**4)/(a**6*d*(a + b*sin(c + d*x))) - (a**4 - 12*a**2*b**2 + 15*b**4)*log(a + b*sin(c + d*x))/(a**7*d) + (a**4 - 12*a**2*b**2 + 15*b**4)*log(sin(c + d*x))/(a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_189():
    f = tan(c + d*x)**4/(a + b*sin(c + d*x))**3
    F = 3*a**5*b*cos(c + d*x)/(2*d*(a + b*sin(c + d*x))*(a**2 - b**2)**4) + 8*a**4*b**2*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(d*(a**2 - b**2)**(sympy.S(9)/2)) + a**4*b*cos(c + d*x)/(2*d*(a + b*sin(c + d*x))**2*(a**2 - b**2)**3) + a**4*(2*a**2 + b**2)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(d*(a**2 - b**2)**(sympy.S(9)/2)) + 4*a**3*b**3*cos(c + d*x)/(d*(a + b*sin(c + d*x))*(a**2 - b**2)**4) + 12*a**2*b**2*(a**2 + b**2)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(d*(a**2 - b**2)**(sympy.S(9)/2)) + 3*a*cos(c + d*x)/(4*d*(a - b)**4*(sin(c + d*x) + 1)) - 3*a*cos(c + d*x)/(4*d*(1 - sin(c + d*x))*(a + b)**4) - cos(c + d*x)/(12*d*(a - b)**3*(sin(c + d*x) + 1)) - cos(c + d*x)/(12*d*(a - b)**3*(sin(c + d*x) + 1)**2) + cos(c + d*x)/(12*d*(1 - sin(c + d*x))*(a + b)**3) + cos(c + d*x)/(12*d*(1 - sin(c + d*x))**2*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_190():
    f = tan(c + d*x)**2/(a + b*sin(c + d*x))**3
    F = -3*a**3*b*cos(c + d*x)/(2*d*(a + b*sin(c + d*x))*(a**2 - b**2)**3) - 4*a**2*b**2*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(d*(a**2 - b**2)**(sympy.S(7)/2)) - a**2*b*cos(c + d*x)/(2*d*(a + b*sin(c + d*x))**2*(a**2 - b**2)**2) - a**2*(2*a**2 + b**2)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(d*(a**2 - b**2)**(sympy.S(7)/2)) - 2*a*b**3*cos(c + d*x)/(d*(a + b*sin(c + d*x))*(a**2 - b**2)**3) - 2*b**2*(3*a**2 + b**2)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(d*(a**2 - b**2)**(sympy.S(7)/2)) - cos(c + d*x)/(2*d*(a - b)**3*(sin(c + d*x) + 1)) + cos(c + d*x)/(2*d*(1 - sin(c + d*x))*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_191():
    f = cot(c + d*x)**2/(a + b*sin(c + d*x))**3
    F = cot(c + d*x)/(2*a*d*(a + b*sin(c + d*x))**2) + (2*a**2 - 3*b**2)*cot(c + d*x)/(2*a**2*d*(a + b*sin(c + d*x))*(a**2 - b**2)) - (5*a**2 - 6*b**2)*cot(c + d*x)/(2*a**3*d*(a**2 - b**2)) + 3*b*atanh(cos(c + d*x))/(a**4*d) - (2*a**4 - 9*a**2*b**2 + 6*b**4)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(a**4*d*(a**2 - b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_192():
    f = cot(c + d*x)**4/(a + b*sin(c + d*x))**3
    F = -cot(c + d*x)*csc(c + d*x)**2/(3*a*d*(a + b*sin(c + d*x))**2) + (3*a**2 - 5*b**2)*cot(c + d*x)*csc(c + d*x)/(6*a**2*b*d*(a + b*sin(c + d*x))**2) + (3*a**2 - 20*b**2)*cot(c + d*x)*csc(c + d*x)/(6*a**3*b*d*(a + b*sin(c + d*x))) - (a**2 - 5*b**2)*cot(c + d*x)*csc(c + d*x)/(a**4*b*d) + (17*a**2 - 60*b**2)*cot(c + d*x)/(6*a**5*d) - b*(9*a**2 - 20*b**2)*atanh(cos(c + d*x))/(2*a**6*d) + (2*a**4 - 19*a**2*b**2 + 20*b**4)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(a**6*d*sqrt(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_193():
    f = cot(c + d*x)**6/(a + b*sin(c + d*x))**3
    F = a*cot(c + d*x)*csc(c + d*x)**2/(12*b**2*d*(a + b*sin(c + d*x))**2) - cot(c + d*x)*csc(c + d*x)/(3*b*d*(a + b*sin(c + d*x))**2) - cot(c + d*x)*csc(c + d*x)**4/(5*a*d*(a + b*sin(c + d*x))**2) + 7*b*cot(c + d*x)*csc(c + d*x)**3/(20*a**2*d*(a + b*sin(c + d*x))**2) + (5*a**4 - 60*a**2*b**2 + 63*b**4)*cot(c + d*x)*csc(c + d*x)**2/(60*a**3*b**2*d*(a + b*sin(c + d*x))**2) + (4*a**4 - 54*a**2*b**2 + 63*b**4)*cot(c + d*x)*csc(c + d*x)**2/(12*a**4*b**2*d*(a + b*sin(c + d*x))) - (15*a**4 - 187*a**2*b**2 + 210*b**4)*cot(c + d*x)*csc(c + d*x)**2/(30*a**5*b**2*d) + (8*a**4 - 79*a**2*b**2 + 84*b**4)*cot(c + d*x)*csc(c + d*x)/(8*a**6*b*d) - (91*a**4 - 645*a**2*b**2 + 630*b**4)*cot(c + d*x)/(30*a**7*d) + b*(45*a**4 - 200*a**2*b**2 + 168*b**4)*atanh(cos(c + d*x))/(8*a**8*d) - sqrt(a**2 - b**2)*(2*a**4 - 29*a**2*b**2 + 42*b**4)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(a**8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_194():
    f = (g*tan(e + f*x))**p*(a + b*sin(e + f*x))**3
    F = a**3*(g*tan(e + f*x))**(p + 1)*hyper((1, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*g*(p + 1)) + 3*a**2*b*(g*tan(e + f*x))**(p + 1)*(cos(e + f*x)**2)**(p/2 + sympy.S.Half)*sin(e + f*x)*hyper((p/2 + sympy.S.Half, p/2 + 1), (p/2 + 2,), sin(e + f*x)**2)/(f*g*(p + 2)) + 3*a*b**2*(g*tan(e + f*x))**(p + 3)*hyper((2, p/2 + sympy.S(3)/2), (p/2 + sympy.S(5)/2,), -tan(e + f*x)**2)/(f*g**3*(p + 3)) + b**3*(g*tan(e + f*x))**(p + 1)*(cos(e + f*x)**2)**(p/2 + sympy.S.Half)*sin(e + f*x)**3*hyper((p/2 + sympy.S.Half, p/2 + 2), (p/2 + 3,), sin(e + f*x)**2)/(f*g*(p + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_195():
    f = (g*tan(e + f*x))**p*(a + b*sin(e + f*x))**2
    F = a**2*(g*tan(e + f*x))**(p + 1)*hyper((1, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*g*(p + 1)) + 2*a*b*(g*tan(e + f*x))**(p + 1)*(cos(e + f*x)**2)**(p/2 + sympy.S.Half)*sin(e + f*x)*hyper((p/2 + sympy.S.Half, p/2 + 1), (p/2 + 2,), sin(e + f*x)**2)/(f*g*(p + 2)) + b**2*(g*tan(e + f*x))**(p + 3)*hyper((2, p/2 + sympy.S(3)/2), (p/2 + sympy.S(5)/2,), -tan(e + f*x)**2)/(f*g**3*(p + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_196():
    f = (g*tan(e + f*x))**p*(a + b*sin(e + f*x))
    F = a*(g*tan(e + f*x))**(p + 1)*hyper((1, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*g*(p + 1)) + b*(g*tan(e + f*x))**(p + 1)*(cos(e + f*x)**2)**(p/2 + sympy.S.Half)*sin(e + f*x)*hyper((p/2 + sympy.S.Half, p/2 + 1), (p/2 + 2,), sin(e + f*x)**2)/(f*g*(p + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_3_g_tan_pow_p_a_plus_b_sin_pow_m_197():
    f = (g*tan(e + f*x))**p*(a + b*sin(e + f*x))**m
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m')) * ((Symbol('g') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))**(Symbol('p'))), x)
    assert integrate(f, x) == F

