"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.6 Cosecant/4.6.1.3 (d cos)^n (a+b csc)^m.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b = symbols('a b')

def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_3_d_cos_pow_n_a_plus_b_csc_pow_m_1():
    f = cos(x)**4/(a*csc(x) + a)
    F = -x/(8*a) + sin(x)*cos(x)**3/(4*a) - sin(x)*cos(x)/(8*a) - cos(x)**3/(3*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_3_d_cos_pow_n_a_plus_b_csc_pow_m_2():
    f = cos(x)**3/(a*csc(x) + a)
    F = -sin(x)**3/(3*a) + sin(x)**2/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_3_d_cos_pow_n_a_plus_b_csc_pow_m_3():
    f = cos(x)**2/(a*csc(x) + a)
    F = -x/(2*a) + sin(x)*cos(x)/(2*a) - cos(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_3_d_cos_pow_n_a_plus_b_csc_pow_m_4():
    f = cos(x)/(a*csc(x) + a)
    F = -log(sin(x) + 1)/a + sin(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_3_d_cos_pow_n_a_plus_b_csc_pow_m_5():
    f = sec(x)/(a*csc(x) + a)
    F = -tan(x)*sec(x)/(2*a) + atanh(sin(x))/(2*a) + sec(x)**2/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_3_d_cos_pow_n_a_plus_b_csc_pow_m_6():
    f = sec(x)**2/(a*csc(x) + a)
    F = -tan(x)**3/(3*a) + sec(x)**3/(3*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_3_d_cos_pow_n_a_plus_b_csc_pow_m_7():
    f = sec(x)**3/(a*csc(x) + a)
    F = -tan(x)*sec(x)**3/(4*a) + tan(x)*sec(x)/(8*a) + atanh(sin(x))/(8*a) + sec(x)**4/(4*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_3_d_cos_pow_n_a_plus_b_csc_pow_m_8():
    f = sec(x)**4/(a*csc(x) + a)
    F = -tan(x)**5/(5*a) - tan(x)**3/(3*a) + sec(x)**5/(5*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_3_d_cos_pow_n_a_plus_b_csc_pow_m_9():
    f = cos(x)**4/(a + b*csc(x))
    F = -(-3*a*sin(x) + 4*b)*cos(x)**3/(12*a**2) - (-a*(3*a**2 - 4*b**2)*sin(x) + 8*b*(a**2 - b**2))*cos(x)/(8*a**4) + 2*b*(a**2 - b**2)**(sympy.S(3)/2)*atanh((a + b*tan(x/2))/sqrt(a**2 - b**2))/a**5 + x*(3*a**4 - 12*a**2*b**2 + 8*b**4)/(8*a**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_3_d_cos_pow_n_a_plus_b_csc_pow_m_10():
    f = cos(x)**3/(a + b*csc(x))
    F = -sin(x)**3/(3*a) + b*sin(x)**2/(2*a**2) + (a**2 - b**2)*sin(x)/a**3 - b*(a**2 - b**2)*log(a*sin(x) + b)/a**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_3_d_cos_pow_n_a_plus_b_csc_pow_m_11():
    f = cos(x)**2/(a + b*csc(x))
    F = -(-a*sin(x) + 2*b)*cos(x)/(2*a**2) + 2*b*sqrt(a**2 - b**2)*atanh((a + b*tan(x/2))/sqrt(a**2 - b**2))/a**3 + x*(a**2 - 2*b**2)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_3_d_cos_pow_n_a_plus_b_csc_pow_m_12():
    f = cos(x)/(a + b*csc(x))
    F = sin(x)/a - b*log(a*sin(x) + b)/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_3_d_cos_pow_n_a_plus_b_csc_pow_m_13():
    f = sec(x)/(a + b*csc(x))
    F = -b*log(a*sin(x) + b)/(a**2 - b**2) - log(1 - sin(x))/(2*a + 2*b) + log(sin(x) + 1)/(2*a - 2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_3_d_cos_pow_n_a_plus_b_csc_pow_m_14():
    f = sec(x)**2/(a + b*csc(x))
    F = 2*a*b*atanh((a + b*tan(x/2))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(3)/2) - (-a*sin(x) + b)*sec(x)/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_3_d_cos_pow_n_a_plus_b_csc_pow_m_15():
    f = sec(x)**3/(a + b*csc(x))
    F = -a**2*b*log(a*sin(x) + b)/(a**2 - b**2)**2 - a*log(1 - sin(x))/(4*(a + b)**2) + a*log(sin(x) + 1)/(4*(a - b)**2) - (-a*sin(x) + b)*sec(x)**2/(2*a**2 - 2*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_3_d_cos_pow_n_a_plus_b_csc_pow_m_16():
    f = sec(x)**4/(a + b*csc(x))
    F = 2*a**3*b*atanh((a + b*tan(x/2))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(5)/2) - (-a*sin(x) + b)*sec(x)**3/(3*a**2 - 3*b**2) - (3*a**2*b - a*(2*a**2 + b**2)*sin(x))*sec(x)/(3*(a**2 - b**2)**2)
    assert integrate(f, x) == F

