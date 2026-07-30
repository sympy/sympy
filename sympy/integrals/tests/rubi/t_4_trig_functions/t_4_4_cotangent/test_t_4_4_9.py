"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.4 Cotangent/4.4.9 trig^m (a+b cot^n+c cot^(2 n))^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e = symbols('a b c d e')

def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_1():
    f = cot(d + e*x)**5/sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)
    F = -b*atanh((b + 2*c*cot(d + e*x))/(2*sqrt(c)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(2*c**(sympy.S(3)/2)*e) + b*(-12*a*c + 5*b**2)*atanh((b + 2*c*cot(d + e*x))/(2*sqrt(c)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(16*c**(sympy.S(7)/2)*e) - sqrt(2)*sqrt(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(a + b*cot(d + e*x) - c - sqrt(a**2 - 2*a*c + b**2 + c**2))/(2*sqrt(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(2*e*sqrt(a**2 - 2*a*c + b**2 + c**2)) + sqrt(2)*sqrt(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(a + b*cot(d + e*x) - c + sqrt(a**2 - 2*a*c + b**2 + c**2))/(2*sqrt(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(2*e*sqrt(a**2 - 2*a*c + b**2 + c**2)) - sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*cot(d + e*x)**2/(3*c*e) + sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)/(c*e) - sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*(-16*a*c + 15*b**2 - 10*b*c*cot(d + e*x))/(24*c**3*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_2():
    f = cot(d + e*x)**3/sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)
    F = b*atanh((b + 2*c*cot(d + e*x))/(2*sqrt(c)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(2*c**(sympy.S(3)/2)*e) + sqrt(2)*sqrt(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(a + b*cot(d + e*x) - c - sqrt(a**2 - 2*a*c + b**2 + c**2))/(2*sqrt(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(2*e*sqrt(a**2 - 2*a*c + b**2 + c**2)) - sqrt(2)*sqrt(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(a + b*cot(d + e*x) - c + sqrt(a**2 - 2*a*c + b**2 + c**2))/(2*sqrt(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(2*e*sqrt(a**2 - 2*a*c + b**2 + c**2)) - sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)/(c*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_3():
    f = cot(d + e*x)/sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)
    F = -sqrt(2)*sqrt(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(a + b*cot(d + e*x) - c - sqrt(a**2 - 2*a*c + b**2 + c**2))/(2*sqrt(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(2*e*sqrt(a**2 - 2*a*c + b**2 + c**2)) + sqrt(2)*sqrt(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(a + b*cot(d + e*x) - c + sqrt(a**2 - 2*a*c + b**2 + c**2))/(2*sqrt(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(2*e*sqrt(a**2 - 2*a*c + b**2 + c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_4():
    f = tan(d + e*x)/sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)
    F = sqrt(2)*sqrt(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(a + b*cot(d + e*x) - c - sqrt(a**2 - 2*a*c + b**2 + c**2))/(2*sqrt(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(2*e*sqrt(a**2 - 2*a*c + b**2 + c**2)) - sqrt(2)*sqrt(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(a + b*cot(d + e*x) - c + sqrt(a**2 - 2*a*c + b**2 + c**2))/(2*sqrt(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(2*e*sqrt(a**2 - 2*a*c + b**2 + c**2)) + atanh((2*a + b*cot(d + e*x))/(2*sqrt(a)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(sqrt(a)*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_5():
    f = tan(d + e*x)**3/sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)
    F = -sqrt(2)*sqrt(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(a + b*cot(d + e*x) - c - sqrt(a**2 - 2*a*c + b**2 + c**2))/(2*sqrt(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(2*e*sqrt(a**2 - 2*a*c + b**2 + c**2)) + sqrt(2)*sqrt(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(a + b*cot(d + e*x) - c + sqrt(a**2 - 2*a*c + b**2 + c**2))/(2*sqrt(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(2*e*sqrt(a**2 - 2*a*c + b**2 + c**2)) + sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*tan(d + e*x)**2/(2*a*e) - 3*b*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*tan(d + e*x)/(4*a**2*e) - atanh((2*a + b*cot(d + e*x))/(2*sqrt(a)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(sqrt(a)*e) + (-4*a*c + 3*b**2)*atanh((2*a + b*cot(d + e*x))/(2*sqrt(a)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(8*a**(sympy.S(5)/2)*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_6():
    f = sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*cot(d + e*x)**5
    F = -b*(b + 2*c*cot(d + e*x))*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)/(8*c**2*e) + b*(b + 2*c*cot(d + e*x))*(-12*a*c + 7*b**2)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)/(128*c**4*e) - b*atanh((b + 2*c*cot(d + e*x))/(2*sqrt(c)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(2*sqrt(c)*e) + b*(-4*a*c + b**2)*atanh((b + 2*c*cot(d + e*x))/(2*sqrt(c)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(16*c**(sympy.S(5)/2)*e) - b*(-12*a*c + 7*b**2)*(-4*a*c + b**2)*atanh((b + 2*c*cot(d + e*x))/(2*sqrt(c)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(256*c**(sympy.S(9)/2)*e) - sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)/e + sqrt(2)*sqrt(a**2 - a*(2*c - sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c - sqrt(a**2 - 2*a*c + b**2 + c**2)))*atanh(sqrt(2)*(b**2 + b*sqrt(a**2 - 2*a*c + b**2 + c**2)*cot(d + e*x) + (a - c)*(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4)*sqrt(a**2 - a*(2*c - sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c - sqrt(a**2 - 2*a*c + b**2 + c**2)))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4)) - sqrt(2)*sqrt(a**2 - a*(2*c + sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c + sqrt(a**2 - 2*a*c + b**2 + c**2)))*atan(sqrt(2)*(b**2 - b*sqrt(a**2 - 2*a*c + b**2 + c**2)*cot(d + e*x) + (a - c)*(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4)*sqrt(a**2 - a*(2*c + sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c + sqrt(a**2 - 2*a*c + b**2 + c**2)))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4)) - (a + b*cot(d + e*x) + c*cot(d + e*x)**2)**(sympy.S(3)/2)*cot(d + e*x)**2/(5*c*e) + (a + b*cot(d + e*x) + c*cot(d + e*x)**2)**(sympy.S(3)/2)/(3*c*e) - (a + b*cot(d + e*x) + c*cot(d + e*x)**2)**(sympy.S(3)/2)*(-32*a*c + 35*b**2 - 42*b*c*cot(d + e*x))/(240*c**3*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_7():
    f = sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*cot(d + e*x)**3
    F = b*(b + 2*c*cot(d + e*x))*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)/(8*c**2*e) + b*atanh((b + 2*c*cot(d + e*x))/(2*sqrt(c)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(2*sqrt(c)*e) - b*(-4*a*c + b**2)*atanh((b + 2*c*cot(d + e*x))/(2*sqrt(c)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(16*c**(sympy.S(5)/2)*e) + sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)/e - sqrt(2)*sqrt(a**2 - a*(2*c - sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c - sqrt(a**2 - 2*a*c + b**2 + c**2)))*atanh(sqrt(2)*(b**2 + b*sqrt(a**2 - 2*a*c + b**2 + c**2)*cot(d + e*x) + (a - c)*(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4)*sqrt(a**2 - a*(2*c - sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c - sqrt(a**2 - 2*a*c + b**2 + c**2)))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4)) + sqrt(2)*sqrt(a**2 - a*(2*c + sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c + sqrt(a**2 - 2*a*c + b**2 + c**2)))*atan(sqrt(2)*(b**2 - b*sqrt(a**2 - 2*a*c + b**2 + c**2)*cot(d + e*x) + (a - c)*(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4)*sqrt(a**2 - a*(2*c + sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c + sqrt(a**2 - 2*a*c + b**2 + c**2)))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4)) - (a + b*cot(d + e*x) + c*cot(d + e*x)**2)**(sympy.S(3)/2)/(3*c*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_8():
    f = sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*cot(d + e*x)
    F = -b*atanh((b + 2*c*cot(d + e*x))/(2*sqrt(c)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(2*sqrt(c)*e) - sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)/e + sqrt(2)*sqrt(a**2 - a*(2*c - sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c - sqrt(a**2 - 2*a*c + b**2 + c**2)))*atanh(sqrt(2)*(b**2 + b*sqrt(a**2 - 2*a*c + b**2 + c**2)*cot(d + e*x) + (a - c)*(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4)*sqrt(a**2 - a*(2*c - sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c - sqrt(a**2 - 2*a*c + b**2 + c**2)))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4)) - sqrt(2)*sqrt(a**2 - a*(2*c + sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c + sqrt(a**2 - 2*a*c + b**2 + c**2)))*atan(sqrt(2)*(b**2 - b*sqrt(a**2 - 2*a*c + b**2 + c**2)*cot(d + e*x) + (a - c)*(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4)*sqrt(a**2 - a*(2*c + sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c + sqrt(a**2 - 2*a*c + b**2 + c**2)))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_9():
    f = sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*tan(d + e*x)
    F = sqrt(a)*atanh((2*a + b*cot(d + e*x))/(2*sqrt(a)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/e - sqrt(2)*sqrt(a**2 - a*(2*c - sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c - sqrt(a**2 - 2*a*c + b**2 + c**2)))*atanh(sqrt(2)*(b**2 + b*sqrt(a**2 - 2*a*c + b**2 + c**2)*cot(d + e*x) + (a - c)*(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4)*sqrt(a**2 - a*(2*c - sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c - sqrt(a**2 - 2*a*c + b**2 + c**2)))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4)) + sqrt(2)*sqrt(a**2 - a*(2*c + sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c + sqrt(a**2 - 2*a*c + b**2 + c**2)))*atan(sqrt(2)*(b**2 - b*sqrt(a**2 - 2*a*c + b**2 + c**2)*cot(d + e*x) + (a - c)*(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4)*sqrt(a**2 - a*(2*c + sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c + sqrt(a**2 - 2*a*c + b**2 + c**2)))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_10():
    f = sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*tan(d + e*x)**3
    F = -sqrt(a)*atanh((2*a + b*cot(d + e*x))/(2*sqrt(a)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/e + sqrt(2)*sqrt(a**2 - a*(2*c - sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c - sqrt(a**2 - 2*a*c + b**2 + c**2)))*atanh(sqrt(2)*(b**2 + b*sqrt(a**2 - 2*a*c + b**2 + c**2)*cot(d + e*x) + (a - c)*(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4)*sqrt(a**2 - a*(2*c - sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c - sqrt(a**2 - 2*a*c + b**2 + c**2)))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4)) - sqrt(2)*sqrt(a**2 - a*(2*c + sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c + sqrt(a**2 - 2*a*c + b**2 + c**2)))*atan(sqrt(2)*(b**2 - b*sqrt(a**2 - 2*a*c + b**2 + c**2)*cot(d + e*x) + (a - c)*(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4)*sqrt(a**2 - a*(2*c + sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c + sqrt(a**2 - 2*a*c + b**2 + c**2)))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4)) + (2*a + b*cot(d + e*x))*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*tan(d + e*x)**2/(4*a*e) - (-4*a*c + b**2)*atanh((2*a + b*cot(d + e*x))/(2*sqrt(a)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(8*a**(sympy.S(3)/2)*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_11():
    f = cot(d + e*x)**7/(a + b*cot(d + e*x) + c*cot(d + e*x)**2)**(sympy.S(3)/2)
    F = 2*b*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*cot(d + e*x)**3/(c*e*(-4*a*c + b**2)) - 3*b*atanh((b + 2*c*cot(d + e*x))/(2*sqrt(c)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(2*c**(sympy.S(5)/2)*e) + 5*b*(-12*a*c + 7*b**2)*atanh((b + 2*c*cot(d + e*x))/(2*sqrt(c)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(16*c**(sympy.S(9)/2)*e) - 2*(2*a + b*cot(d + e*x))*cot(d + e*x)**4/(e*(-4*a*c + b**2)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)) + 2*(2*a + b*cot(d + e*x))*cot(d + e*x)**2/(e*(-4*a*c + b**2)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)) - (4*a + 2*b*cot(d + e*x))/(e*(-4*a*c + b**2)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)) + sqrt(2)*sqrt(2*a - 2*c - sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 + (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(b**2 - b*(2*a - 2*c - sqrt(a**2 - 2*a*c + b**2 + c**2))*cot(d + e*x) - (a - c)*(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*sqrt(2*a - 2*c - sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 + (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(3)/2)) - sqrt(2)*sqrt(2*a - 2*c + sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 - (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(b**2 - b*(2*a - 2*c + sqrt(a**2 - 2*a*c + b**2 + c**2))*cot(d + e*x) - (a - c)*(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*sqrt(2*a - 2*c + sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 - (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(3)/2)) + (2*a*(b**2 - c*(2*a - 2*c)) + 2*b*c*(a + c)*cot(d + e*x))/(e*(b**2 + (a - c)**2)*(-4*a*c + b**2)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)) - (-16*a*c + 7*b**2)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*cot(d + e*x)**2/(3*c**2*e*(-4*a*c + b**2)) + sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*(-8*a*c + 3*b**2 - 2*b*c*cot(d + e*x))/(c**2*e*(-4*a*c + b**2)) - sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*(256*a**2*c**2 - 460*a*b**2*c + 105*b**4 - 2*b*c*(-116*a*c + 35*b**2)*cot(d + e*x))/(24*c**4*e*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_12():
    f = cot(d + e*x)**5/(a + b*cot(d + e*x) + c*cot(d + e*x)**2)**(sympy.S(3)/2)
    F = 3*b*atanh((b + 2*c*cot(d + e*x))/(2*sqrt(c)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(2*c**(sympy.S(5)/2)*e) - 2*(2*a + b*cot(d + e*x))*cot(d + e*x)**2/(e*(-4*a*c + b**2)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)) + (4*a + 2*b*cot(d + e*x))/(e*(-4*a*c + b**2)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)) - sqrt(2)*sqrt(2*a - 2*c - sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 + (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(b**2 - b*(2*a - 2*c - sqrt(a**2 - 2*a*c + b**2 + c**2))*cot(d + e*x) - (a - c)*(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*sqrt(2*a - 2*c - sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 + (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(3)/2)) + sqrt(2)*sqrt(2*a - 2*c + sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 - (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(b**2 - b*(2*a - 2*c + sqrt(a**2 - 2*a*c + b**2 + c**2))*cot(d + e*x) - (a - c)*(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*sqrt(2*a - 2*c + sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 - (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(3)/2)) - (2*a*(b**2 - c*(2*a - 2*c)) + 2*b*c*(a + c)*cot(d + e*x))/(e*(b**2 + (a - c)**2)*(-4*a*c + b**2)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)) - sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*(-8*a*c + 3*b**2 - 2*b*c*cot(d + e*x))/(c**2*e*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_13():
    f = cot(d + e*x)**3/(a + b*cot(d + e*x) + c*cot(d + e*x)**2)**(sympy.S(3)/2)
    F = -(4*a + 2*b*cot(d + e*x))/(e*(-4*a*c + b**2)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)) + sqrt(2)*sqrt(2*a - 2*c - sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 + (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(b**2 - b*(2*a - 2*c - sqrt(a**2 - 2*a*c + b**2 + c**2))*cot(d + e*x) - (a - c)*(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*sqrt(2*a - 2*c - sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 + (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(3)/2)) - sqrt(2)*sqrt(2*a - 2*c + sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 - (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(b**2 - b*(2*a - 2*c + sqrt(a**2 - 2*a*c + b**2 + c**2))*cot(d + e*x) - (a - c)*(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*sqrt(2*a - 2*c + sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 - (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(3)/2)) + (2*a*(b**2 - c*(2*a - 2*c)) + 2*b*c*(a + c)*cot(d + e*x))/(e*(b**2 + (a - c)**2)*(-4*a*c + b**2)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_14():
    f = cot(d + e*x)/(a + b*cot(d + e*x) + c*cot(d + e*x)**2)**(sympy.S(3)/2)
    F = -sqrt(2)*sqrt(2*a - 2*c - sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 + (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(b**2 - b*(2*a - 2*c - sqrt(a**2 - 2*a*c + b**2 + c**2))*cot(d + e*x) - (a - c)*(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*sqrt(2*a - 2*c - sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 + (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(3)/2)) + sqrt(2)*sqrt(2*a - 2*c + sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 - (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(b**2 - b*(2*a - 2*c + sqrt(a**2 - 2*a*c + b**2 + c**2))*cot(d + e*x) - (a - c)*(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*sqrt(2*a - 2*c + sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 - (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(3)/2)) - (2*a*(b**2 - c*(2*a - 2*c)) + 2*b*c*(a + c)*cot(d + e*x))/(e*(b**2 + (a - c)**2)*(-4*a*c + b**2)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_15():
    f = tan(d + e*x)/(a + b*cot(d + e*x) + c*cot(d + e*x)**2)**(sympy.S(3)/2)
    F = sqrt(2)*sqrt(2*a - 2*c - sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 + (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(b**2 - b*(2*a - 2*c - sqrt(a**2 - 2*a*c + b**2 + c**2))*cot(d + e*x) - (a - c)*(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*sqrt(2*a - 2*c - sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 + (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(3)/2)) - sqrt(2)*sqrt(2*a - 2*c + sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 - (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(b**2 - b*(2*a - 2*c + sqrt(a**2 - 2*a*c + b**2 + c**2))*cot(d + e*x) - (a - c)*(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*sqrt(2*a - 2*c + sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 - (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(3)/2)) + (2*a*(b**2 - c*(2*a - 2*c)) + 2*b*c*(a + c)*cot(d + e*x))/(e*(b**2 + (a - c)**2)*(-4*a*c + b**2)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)) - (-4*a*c + 2*b**2 + 2*b*c*cot(d + e*x))/(a*e*(-4*a*c + b**2)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)) + atanh((2*a + b*cot(d + e*x))/(2*sqrt(a)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(a**(sympy.S(3)/2)*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_16():
    f = tan(d + e*x)**3/(a + b*cot(d + e*x) + c*cot(d + e*x)**2)**(sympy.S(3)/2)
    F = -sqrt(2)*sqrt(2*a - 2*c - sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 + (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(b**2 - b*(2*a - 2*c - sqrt(a**2 - 2*a*c + b**2 + c**2))*cot(d + e*x) - (a - c)*(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*sqrt(2*a - 2*c - sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 + (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(3)/2)) + sqrt(2)*sqrt(2*a - 2*c + sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 - (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))*atanh(sqrt(2)*(b**2 - b*(2*a - 2*c + sqrt(a**2 - 2*a*c + b**2 + c**2))*cot(d + e*x) - (a - c)*(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2)))/(2*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*sqrt(2*a - 2*c + sqrt(a**2 - 2*a*c + b**2 + c**2))*sqrt(a**2 - 2*a*c - b**2 + c**2 - (a - c)*sqrt(a**2 - 2*a*c + b**2 + c**2))))/(2*e*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(3)/2)) - (2*a*(b**2 - c*(2*a - 2*c)) + 2*b*c*(a + c)*cot(d + e*x))/(e*(b**2 + (a - c)**2)*(-4*a*c + b**2)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)) - (-4*a*c + 2*b**2 + 2*b*c*cot(d + e*x))*tan(d + e*x)**2/(a*e*(-4*a*c + b**2)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)) + (-4*a*c + 2*b**2 + 2*b*c*cot(d + e*x))/(a*e*(-4*a*c + b**2)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)) + (-12*a*c + 5*b**2)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*tan(d + e*x)**2/(2*a**2*e*(-4*a*c + b**2)) - b*(-52*a*c + 15*b**2)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)*tan(d + e*x)/(4*a**3*e*(-4*a*c + b**2)) - atanh((2*a + b*cot(d + e*x))/(2*sqrt(a)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(a**(sympy.S(3)/2)*e) + (-12*a*c + 15*b**2)*atanh((2*a + b*cot(d + e*x))/(2*sqrt(a)*sqrt(a + b*cot(d + e*x) + c*cot(d + e*x)**2)))/(8*a**(sympy.S(7)/2)*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_17():
    f = cot(d + e*x)**5/sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)
    F = atanh((2*a - b + (b - 2*c)*cot(d + e*x)**2)/(2*sqrt(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*e*sqrt(a - b + c)) - sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)/(2*c*e) + (b + 2*c)*atanh((b + 2*c*cot(d + e*x)**2)/(2*sqrt(c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(4*c**(sympy.S(3)/2)*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_18():
    f = cot(d + e*x)**3/sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)
    F = -atanh((2*a - b + (b - 2*c)*cot(d + e*x)**2)/(2*sqrt(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*e*sqrt(a - b + c)) - atanh((b + 2*c*cot(d + e*x)**2)/(2*sqrt(c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*sqrt(c)*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_19():
    f = cot(d + e*x)/sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)
    F = atanh((2*a - b + (b - 2*c)*cot(d + e*x)**2)/(2*sqrt(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*e*sqrt(a - b + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_20():
    f = tan(d + e*x)/sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)
    F = -atanh((2*a - b + (b - 2*c)*cot(d + e*x)**2)/(2*sqrt(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*e*sqrt(a - b + c)) + atanh((2*a + b*cot(d + e*x)**2)/(2*sqrt(a)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*sqrt(a)*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_21():
    f = tan(d + e*x)**3/sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)
    F = atanh((2*a - b + (b - 2*c)*cot(d + e*x)**2)/(2*sqrt(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*e*sqrt(a - b + c)) + sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)*tan(d + e*x)**2/(2*a*e) - atanh((2*a + b*cot(d + e*x)**2)/(2*sqrt(a)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*sqrt(a)*e) - b*atanh((2*a + b*cot(d + e*x)**2)/(2*sqrt(a)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(4*a**(sympy.S(3)/2)*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_22():
    f = sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)*cot(d + e*x)**5
    F = sqrt(a - b + c)*atanh((2*a - b + (b - 2*c)*cot(d + e*x)**2)/(2*sqrt(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*e) - (a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)**(sympy.S(3)/2)/(6*c*e) + (2*c*(b + 2*c)*cot(d + e*x)**2 + (b - 2*c)*(b + 4*c))*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)/(16*c**2*e) - (b**3 + 2*b**2*c - 4*b*c*(a - 2*c) - 8*c**2*(a + 2*c))*atanh((b + 2*c*cot(d + e*x)**2)/(2*sqrt(c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(32*c**(sympy.S(5)/2)*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_23():
    f = sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)*cot(d + e*x)**3
    F = -sqrt(a - b + c)*atanh((2*a - b + (b - 2*c)*cot(d + e*x)**2)/(2*sqrt(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*e) - sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)*(b + 2*c*cot(d + e*x)**2 - 4*c)/(8*c*e) + (b**2 + 4*b*c - 4*c*(a + 2*c))*atanh((b + 2*c*cot(d + e*x)**2)/(2*sqrt(c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(16*c**(sympy.S(3)/2)*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_24():
    f = sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)*cot(d + e*x)
    F = sqrt(a - b + c)*atanh((2*a - b + (b - 2*c)*cot(d + e*x)**2)/(2*sqrt(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*e) - sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)/(2*e) - (b - 2*c)*atanh((b + 2*c*cot(d + e*x)**2)/(2*sqrt(c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(4*sqrt(c)*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_25():
    f = sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)*tan(d + e*x)
    F = sqrt(a)*atanh((2*a + b*cot(d + e*x)**2)/(2*sqrt(a)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*e) - sqrt(c)*atanh((b + 2*c*cot(d + e*x)**2)/(2*sqrt(c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*e) - sqrt(a - b + c)*atanh((2*a - b + (b - 2*c)*cot(d + e*x)**2)/(2*sqrt(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_26():
    f = sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)*tan(d + e*x)**3
    F = -sqrt(a)*atanh((2*a + b*cot(d + e*x)**2)/(2*sqrt(a)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*e) + b*atanh((b + 2*c*cot(d + e*x)**2)/(2*sqrt(c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(4*sqrt(c)*e) - sqrt(c)*atanh((b + 2*c*cot(d + e*x)**2)/(2*sqrt(c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*e) + sqrt(a - b + c)*atanh((2*a - b + (b - 2*c)*cot(d + e*x)**2)/(2*sqrt(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*e) + sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)*tan(d + e*x)**2/(2*e) - (b - 2*c)*atanh((b + 2*c*cot(d + e*x)**2)/(2*sqrt(c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(4*sqrt(c)*e) + b*atanh((2*a + b*cot(d + e*x)**2)/(2*sqrt(a)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(4*sqrt(a)*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_27():
    f = cot(d + e*x)**7/(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)**(sympy.S(3)/2)
    F = -atanh((2*a - b + (b - 2*c)*cot(d + e*x)**2)/(2*sqrt(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*e*(a - b + c)**(sympy.S(3)/2)) - (a*(-a*(b + 2*c) + b**2) + (2*a**2*c - a*b*(b + 3*c) + b**3)*cot(d + e*x)**2)/(c*e*(-4*a*c + b**2)*(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)) - atanh((b + 2*c*cot(d + e*x)**2)/(2*sqrt(c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*c**(sympy.S(3)/2)*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_28():
    f = cot(d + e*x)**5/(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)**(sympy.S(3)/2)
    F = atanh((2*a - b + (b - 2*c)*cot(d + e*x)**2)/(2*sqrt(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*e*(a - b + c)**(sympy.S(3)/2)) - (a*(2*a - b) + (2*a*c + b*(a - b))*cot(d + e*x)**2)/(e*(-4*a*c + b**2)*(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_29():
    f = cot(d + e*x)**3/(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)**(sympy.S(3)/2)
    F = -atanh((2*a - b + (b - 2*c)*cot(d + e*x)**2)/(2*sqrt(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*e*(a - b + c)**(sympy.S(3)/2)) + (a*(b - 2*c) + c*(2*a - b)*cot(d + e*x)**2)/(e*(-4*a*c + b**2)*(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_30():
    f = cot(d + e*x)/(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)**(sympy.S(3)/2)
    F = atanh((2*a - b + (b - 2*c)*cot(d + e*x)**2)/(2*sqrt(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*e*(a - b + c)**(sympy.S(3)/2)) - (-2*a*c + b**2 - b*c + c*(b - 2*c)*cot(d + e*x)**2)/(e*(-4*a*c + b**2)*(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_31():
    f = tan(d + e*x)/(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)**(sympy.S(3)/2)
    F = -atanh((2*a - b + (b - 2*c)*cot(d + e*x)**2)/(2*sqrt(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*e*(a - b + c)**(sympy.S(3)/2)) + (-2*a*c + b**2 - b*c + c*(b - 2*c)*cot(d + e*x)**2)/(e*(-4*a*c + b**2)*(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)) - (-2*a*c + b**2 + b*c*cot(d + e*x)**2)/(a*e*(-4*a*c + b**2)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)) + atanh((2*a + b*cot(d + e*x)**2)/(2*sqrt(a)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*a**(sympy.S(3)/2)*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_9_trig_pow_m_a_plus_b_cot_pow_n_plus_c_cot_pow_2_n_pow_p_32():
    f = tan(d + e*x)**3/(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)**(sympy.S(3)/2)
    F = atanh((2*a - b + (b - 2*c)*cot(d + e*x)**2)/(2*sqrt(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*e*(a - b + c)**(sympy.S(3)/2)) - (-2*a*c + b**2 - b*c + c*(b - 2*c)*cot(d + e*x)**2)/(e*(-4*a*c + b**2)*(a - b + c)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)) - (-2*a*c + b**2 + b*c*cot(d + e*x)**2)*tan(d + e*x)**2/(a*e*(-4*a*c + b**2)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)) + (-2*a*c + b**2 + b*c*cot(d + e*x)**2)/(a*e*(-4*a*c + b**2)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)) + (-8*a*c + 3*b**2)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)*tan(d + e*x)**2/(2*a**2*e*(-4*a*c + b**2)) - atanh((2*a + b*cot(d + e*x)**2)/(2*sqrt(a)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(2*a**(sympy.S(3)/2)*e) - 3*b*atanh((2*a + b*cot(d + e*x)**2)/(2*sqrt(a)*sqrt(a + b*cot(d + e*x)**2 + c*cot(d + e*x)**4)))/(4*a**(sympy.S(5)/2)*e)
    assert integrate(f, x) == F

