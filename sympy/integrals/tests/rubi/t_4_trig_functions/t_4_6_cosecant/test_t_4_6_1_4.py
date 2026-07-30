"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.6 Cosecant/4.6.1.4 (d cot)^n (a+b csc)^m.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b = symbols('a b')

def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_4_d_cot_pow_n_a_plus_b_csc_pow_m_1():
    f = tan(x)**4/(a*csc(x) + a)
    F = x/a - (1 - csc(x))*tan(x)**5/(5*a) + (5 - 4*csc(x))*tan(x)**3/(15*a) - (15 - 8*csc(x))*tan(x)/(15*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_4_d_cot_pow_n_a_plus_b_csc_pow_m_2():
    f = tan(x)**3/(a*csc(x) + a)
    F = 5*log(1 - sin(x))/(16*a) + 11*log(sin(x) + 1)/(16*a) + 3/(4*a*(sin(x) + 1)) - 1/(8*a*(sin(x) + 1)**2) + 1/(8*a*(1 - sin(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_4_d_cot_pow_n_a_plus_b_csc_pow_m_3():
    f = tan(x)**2/(a*csc(x) + a)
    F = -x/a - (1 - csc(x))*tan(x)**3/(3*a) + (3 - 2*csc(x))*tan(x)/(3*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_4_d_cot_pow_n_a_plus_b_csc_pow_m_4():
    f = tan(x)/(a*csc(x) + a)
    F = -log(1 - sin(x))/(4*a) - 3*log(sin(x) + 1)/(4*a) - 1/(2*a*(sin(x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_4_d_cot_pow_n_a_plus_b_csc_pow_m_5():
    f = cot(x)/(a*csc(x) + a)
    F = log(sin(x) + 1)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_4_d_cot_pow_n_a_plus_b_csc_pow_m_6():
    f = cot(x)**2/(a*csc(x) + a)
    F = -x/a - atanh(cos(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_4_d_cot_pow_n_a_plus_b_csc_pow_m_7():
    f = cot(x)**3/(a*csc(x) + a)
    F = -log(sin(x))/a - csc(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_4_d_cot_pow_n_a_plus_b_csc_pow_m_8():
    f = cot(x)**4/(a*csc(x) + a)
    F = x/a + (2 - csc(x))*cot(x)/(2*a) + atanh(cos(x))/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_4_d_cot_pow_n_a_plus_b_csc_pow_m_9():
    f = cot(x)**5/(a*csc(x) + a)
    F = log(sin(x))/a - csc(x)**3/(3*a) + csc(x)**2/(2*a) + csc(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_4_d_cot_pow_n_a_plus_b_csc_pow_m_10():
    f = cot(x)**6/(a*csc(x) + a)
    F = -x/a + (4 - 3*csc(x))*cot(x)**3/(12*a) - (8 - 3*csc(x))*cot(x)/(8*a) - 3*atanh(cos(x))/(8*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_4_d_cot_pow_n_a_plus_b_csc_pow_m_11():
    f = cot(x)**7/(a*csc(x) + a)
    F = -log(sin(x))/a - csc(x)**5/(5*a) + csc(x)**4/(4*a) + 2*csc(x)**3/(3*a) - csc(x)**2/a - csc(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_4_d_cot_pow_n_a_plus_b_csc_pow_m_12():
    f = tan(x)**5/(a + b*csc(x))
    F = 1/((16*a - 16*b)*(csc(x) + 1)**2) - (8*a**2 + 21*a*b + 15*b**2)*log(1 - csc(x))/(16*(a + b)**3) + (5*a - 7*b)/(16*(a - b)**2*(csc(x) + 1)) - (8*a**2 - 21*a*b + 15*b**2)*log(csc(x) + 1)/(16*(a - b)**3) + (5*a + 7*b)/(16*(1 - csc(x))*(a + b)**2) + 1/((1 - csc(x))**2*(16*a + 16*b)) + b**6*log(a + b*csc(x))/(a*(a**2 - b**2)**3) - log(sin(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_4_d_cot_pow_n_a_plus_b_csc_pow_m_13():
    f = tan(x)**3/(a + b*csc(x))
    F = -1/((4*a - 4*b)*(csc(x) + 1)) + (2*a + 3*b)*log(1 - csc(x))/(4*(a + b)**2) + (2*a - 3*b)*log(csc(x) + 1)/(4*(a - b)**2) - 1/((1 - csc(x))*(4*a + 4*b)) + b**4*log(a + b*csc(x))/(a*(a**2 - b**2)**2) + log(sin(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_4_d_cot_pow_n_a_plus_b_csc_pow_m_14():
    f = tan(x)/(a + b*csc(x))
    F = -log(1 - csc(x))/(2*a + 2*b) - log(csc(x) + 1)/(2*a - 2*b) + b**2*log(a + b*csc(x))/(a*(a**2 - b**2)) - log(sin(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_4_d_cot_pow_n_a_plus_b_csc_pow_m_15():
    f = cot(x)/(a + b*csc(x))
    F = log(a + b*csc(x))/a + log(sin(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_4_d_cot_pow_n_a_plus_b_csc_pow_m_16():
    f = cot(x)**3/(a + b*csc(x))
    F = -csc(x)/b - (-a**2/b**2 + 1)*log(a + b*csc(x))/a - log(sin(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_4_d_cot_pow_n_a_plus_b_csc_pow_m_17():
    f = cot(x)**5/(a + b*csc(x))
    F = a*csc(x)**2/(2*b**2) - csc(x)**3/(3*b) - (a**2 - 2*b**2)*csc(x)/b**3 + log(sin(x))/a + (a**2 - b**2)**2*log(a + b*csc(x))/(a*b**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_4_d_cot_pow_n_a_plus_b_csc_pow_m_18():
    f = cot(x)**7/(a + b*csc(x))
    F = a*csc(x)**4/(4*b**2) + a*(a**2 - 3*b**2)*csc(x)**2/(2*b**4) - csc(x)**5/(5*b) - (a**2 - 3*b**2)*csc(x)**3/(3*b**3) - (a**4 - 3*a**2*b**2 + 3*b**4)*csc(x)/b**5 - log(sin(x))/a + (a**2 - b**2)**3*log(a + b*csc(x))/(a*b**6)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_4_d_cot_pow_n_a_plus_b_csc_pow_m_19():
    f = cot(x)**2/(a + b*csc(x))
    F = -atanh(cos(x))/b - x/a + 2*sqrt(a**2 - b**2)*atanh((a + b*tan(x/2))/sqrt(a**2 - b**2))/(a*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_4_d_cot_pow_n_a_plus_b_csc_pow_m_20():
    f = cot(x)**4/(a + b*csc(x))
    F = a*cot(x)/b**2 - cot(x)*csc(x)/(2*b) - (2*a**2 - 3*b**2)*atanh(cos(x))/(2*b**3) + x/a + 2*(a**2 - b**2)**(sympy.S(3)/2)*atanh((a + b*tan(x/2))/sqrt(a**2 - b**2))/(a*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_4_d_cot_pow_n_a_plus_b_csc_pow_m_21():
    f = cot(x)**6/(a + b*csc(x))
    F = a*cot(x)**3/(3*b**2) + a*cot(x)/b**2 + a*(a**2 - 3*b**2)*cot(x)/b**4 - cot(x)*csc(x)**3/(4*b) - 3*cot(x)*csc(x)/(8*b) - 3*atanh(cos(x))/(8*b) - (a**2 - 3*b**2)*cot(x)*csc(x)/(2*b**3) - (a**2 - 3*b**2)*atanh(cos(x))/(2*b**3) - (a**4 - 3*a**2*b**2 + 3*b**4)*atanh(cos(x))/b**5 - x/a + 2*(a**2 - b**2)**(sympy.S(5)/2)*atanh((a + b*tan(x/2))/sqrt(a**2 - b**2))/(a*b**5)
    assert integrate(f, x) == F

