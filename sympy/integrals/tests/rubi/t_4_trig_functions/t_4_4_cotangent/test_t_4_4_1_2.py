"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.4 Cotangent/4.4.1.2 (d csc)^m (a+b cot)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, n = symbols('a b n')

def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_1():
    f = sin(x)**4/(cot(x) + I)
    F = -5*I*x/16 + 3*I/(16*cot(x) + 16*I) - 3/(32*(cot(x) + I)**2) - I/(24*(cot(x) + I)**3) + 1/(32*(-cot(x) + I)**2) - I/(-8*cot(x) + 8*I)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_2():
    f = sin(x)**3/(cot(x) + I)
    F = -4*I*cos(x)**3/15 + 4*I*cos(x)/5 + I*sin(x)**3/(5*cot(x) + 5*I)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_3():
    f = sin(x)**2/(cot(x) + I)
    F = -3*I*x/8 + I/(4*cot(x) + 4*I) - 1/(8*(cot(x) + I)**2) - I/(-8*cot(x) + 8*I)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_4():
    f = sin(x)/(cot(x) + I)
    F = 2*I*cos(x)/3 + I*sin(x)/(3*cot(x) + 3*I)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_5():
    f = csc(x)/(cot(x) + I)
    F = I*csc(x)/(cot(x) + I)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_6():
    f = csc(x)**2/(cot(x) + I)
    F = -I*x + log(sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_7():
    f = csc(x)**3/(cot(x) + I)
    F = I*atanh(cos(x)) - csc(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_8():
    f = csc(x)**4/(cot(x) + I)
    F = -cot(x)**2/2 + I*cot(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_9():
    f = csc(x)**5/(cot(x) + I)
    F = I*cot(x)*csc(x)/2 + I*atanh(cos(x))/2 - csc(x)**3/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_10():
    f = csc(x)**6/(cot(x) + I)
    F = -cot(x)**4/4 + I*cot(x)**3/3 - cot(x)**2/2 + I*cot(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_11():
    f = csc(x)**7/(cot(x) + I)
    F = I*cot(x)*csc(x)**3/4 + 3*I*cot(x)*csc(x)/8 + 3*I*atanh(cos(x))/8 - csc(x)**5/5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_12():
    f = csc(x)**6/(a + b*cot(x))
    F = a*cot(x)**3/(3*b**2) + a*(a**2 + 2*b**2)*cot(x)/b**4 - cot(x)**4/(4*b) - (a**2 + 2*b**2)*cot(x)**2/(2*b**3) - (a**2 + b**2)**2*log(a + b*cot(x))/b**5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_13():
    f = csc(x)**4/(a + b*cot(x))
    F = a*cot(x)/b**2 - cot(x)**2/(2*b) - (a**2 + b**2)*log(a + b*cot(x))/b**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_14():
    f = csc(x)**2/(a + b*cot(x))
    F = -log(a + b*cot(x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_15():
    f = sin(x)**2/(a + b*cot(x))
    F = a*x*(a**2 + 3*b**2)/(2*(a**2 + b**2)**2) - b**3*log(a*sin(x) + b*cos(x))/(a**2 + b**2)**2 - (a*cot(x) + b)*sin(x)**2/(2*a**2 + 2*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_16():
    f = sin(x)**4/(a + b*cot(x))
    F = a*x*(3*a**4 + 10*a**2*b**2 + 15*b**4)/(8*(a**2 + b**2)**3) - b**5*log(a*sin(x) + b*cos(x))/(a**2 + b**2)**3 - (a*cot(x) + b)*sin(x)**4/(4*a**2 + 4*b**2) - (a*(3*a**2 + 7*b**2)*cot(x) + 4*b**3)*sin(x)**2/(8*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_17():
    f = csc(x)**5/(a + b*cot(x))
    F = a*cot(x)*csc(x)/(2*b**2) + a*atanh(cos(x))/(2*b**2) + a*(a**2 + b**2)*atanh(cos(x))/b**4 - csc(x)**3/(3*b) - (a**2 + b**2)*csc(x)/b**3 + (a**2 + b**2)**(sympy.S(3)/2)*atanh((-a*cot(x) + b)*sin(x)/sqrt(a**2 + b**2))/b**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_18():
    f = csc(x)**3/(a + b*cot(x))
    F = a*atanh(cos(x))/b**2 - csc(x)/b + sqrt(a**2 + b**2)*atanh((-a*cot(x) + b)*sin(x)/sqrt(a**2 + b**2))/b**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_19():
    f = csc(x)/(a + b*cot(x))
    F = -atanh((a*cot(x) - b)*sin(x)/sqrt(a**2 + b**2))/sqrt(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_20():
    f = sin(x)/(a + b*cot(x))
    F = -a*cos(x)/(a**2 + b**2) + b**2*atanh((-a*cot(x) + b)*sin(x)/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(3)/2) - b*sin(x)/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_21():
    f = sin(x)**3/(a + b*cot(x))
    F = -a*b**2*cos(x)/(a**2 + b**2)**2 + a*cos(x)**3/(3*a**2 + 3*b**2) - a*cos(x)/(a**2 + b**2) + b**4*atanh((-a*cot(x) + b)*sin(x)/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(5)/2) - b**3*sin(x)/(a**2 + b**2)**2 - b*sin(x)**3/(3*a**2 + 3*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_22():
    f = csc(x)**2/(a + b*cot(x))**2
    F = 1/(b*(a + b*cot(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_1_2_d_csc_pow_m_a_plus_b_cot_pow_n_23():
    f = (a + b*cot(x))**n*csc(x)**2
    F = -(a + b*cot(x))**(n + 1)/(b*(n + 1))
    assert integrate(f, x) == F

