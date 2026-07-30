"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.4 Cotangent/4.4.1.3 (d cos)^m (a+b cot)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b = symbols('a b')

def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_3_d_cos_pow_m_a_plus_b_cot_pow_n_1():
    f = cos(x)**4/(cot(x) + I)
    F = -I*x/16 + 3*I/(16*cot(x) + 16*I) + 5/(32*(cot(x) + I)**2) - I/(24*(cot(x) + I)**3) + 1/(32*(-cot(x) + I)**2) + I/(-8*cot(x) + 8*I)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_3_d_cos_pow_m_a_plus_b_cot_pow_n_2():
    f = cos(x)**3/(cot(x) + I)
    F = I*sin(x)**5/5 - I*sin(x)**3/3 - cos(x)**5/5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_3_d_cos_pow_m_a_plus_b_cot_pow_n_3():
    f = cos(x)**2/(cot(x) + I)
    F = -I*x/8 + I/(4*cot(x) + 4*I) + 1/(8*(cot(x) + I)**2) + I/(-8*cot(x) + 8*I)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_3_d_cos_pow_m_a_plus_b_cot_pow_n_4():
    f = cos(x)/(cot(x) + I)
    F = -I*sin(x)**3/3 - cos(x)**3/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_3_d_cos_pow_m_a_plus_b_cot_pow_n_5():
    f = sec(x)/(cot(x) + I)
    F = I*sin(x) - cos(x) - I*atanh(sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_3_d_cos_pow_m_a_plus_b_cot_pow_n_6():
    f = sec(x)**2/(cot(x) + I)
    F = I*x - log(sin(x)) + log(tan(x)) - I*tan(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_3_d_cos_pow_m_a_plus_b_cot_pow_n_7():
    f = sec(x)**3/(cot(x) + I)
    F = -I*tan(x)*sec(x)/2 + I*atanh(sin(x))/2 + sec(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_3_d_cos_pow_m_a_plus_b_cot_pow_n_8():
    f = sec(x)**4/(cot(x) + I)
    F = -I*tan(x)**3/3 + tan(x)**2/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_3_d_cos_pow_m_a_plus_b_cot_pow_n_9():
    f = sec(x)**5/(cot(x) + I)
    F = -I*tan(x)*sec(x)**3/4 + I*tan(x)*sec(x)/8 + I*atanh(sin(x))/8 + sec(x)**3/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_3_d_cos_pow_m_a_plus_b_cot_pow_n_10():
    f = sec(x)**6/(cot(x) + I)
    F = -I*tan(x)**5/5 + tan(x)**4/4 - I*tan(x)**3/3 + tan(x)**2/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_3_d_cos_pow_m_a_plus_b_cot_pow_n_11():
    f = cos(x)**4/(a + b*cot(x))
    F = -a**4*b*log(a*sin(x) + b*cos(x))/(a**2 + b**2)**3 + a*x*(3*a**4 - 6*a**2*b**2 - b**4)/(8*(a**2 + b**2)**3) - (a*cot(x) + b)*sin(x)**4/(4*a**2 + 4*b**2) + (a*(5*a**2 + b**2)*cot(x) + 4*b*(2*a**2 + b**2))*sin(x)**2/(8*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_3_d_cos_pow_m_a_plus_b_cot_pow_n_12():
    f = cos(x)**3/(a + b*cot(x))
    F = a**3*b*atanh((a*cos(x) - b*sin(x))/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(5)/2) - a**2*b*cos(x)/(a**2 + b**2)**2 - a*b**2*sin(x)/(a**2 + b**2)**2 - a*sin(x)**3/(3*a**2 + 3*b**2) + a*sin(x)/(a**2 + b**2) - b*cos(x)**3/(3*a**2 + 3*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_3_d_cos_pow_m_a_plus_b_cot_pow_n_13():
    f = cos(x)**2/(a + b*cot(x))
    F = -a**2*b*log(a*sin(x) + b*cos(x))/(a**2 + b**2)**2 + a*x*(a**2 - b**2)/(2*(a**2 + b**2)**2) + (a*cot(x) + b)*sin(x)**2/(2*a**2 + 2*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_3_d_cos_pow_m_a_plus_b_cot_pow_n_14():
    f = cos(x)/(a + b*cot(x))
    F = a*b*atanh((a*cos(x) - b*sin(x))/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(3)/2) + a*sin(x)/(a**2 + b**2) - b*cos(x)/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_3_d_cos_pow_m_a_plus_b_cot_pow_n_15():
    f = sec(x)/(a + b*cot(x))
    F = b*atanh((a*cos(x) - b*sin(x))/sqrt(a**2 + b**2))/(a*sqrt(a**2 + b**2)) + atanh(sin(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_3_d_cos_pow_m_a_plus_b_cot_pow_n_16():
    f = sec(x)**2/(a + b*cot(x))
    F = tan(x)/a - b*log(a + b*cot(x))/a**2 - b*log(tan(x))/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_3_d_cos_pow_m_a_plus_b_cot_pow_n_17():
    f = sec(x)**3/(a + b*cot(x))
    F = tan(x)*sec(x)/(2*a) + atanh(sin(x))/(2*a) - b*sec(x)/a**2 + b**2*atanh(sin(x))/a**3 + b*sqrt(a**2 + b**2)*atanh((a*cos(x) - b*sin(x))/sqrt(a**2 + b**2))/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_3_d_cos_pow_m_a_plus_b_cot_pow_n_18():
    f = sec(x)**4/(a + b*cot(x))
    F = tan(x)**3/(3*a) - b*tan(x)**2/(2*a**2) + (a**2 + b**2)*tan(x)/a**3 - b*(a**2 + b**2)*log(a + b*cot(x))/a**4 - b*(a**2 + b**2)*log(tan(x))/a**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_3_d_cos_pow_m_a_plus_b_cot_pow_n_19():
    f = sec(x)/(2*cot(x) + 1)
    F = 2*sqrt(5)*atanh(sqrt(5)*(-2*sin(x) + cos(x))/5)/5 + atanh(sin(x))
    assert integrate(f, x) == F

