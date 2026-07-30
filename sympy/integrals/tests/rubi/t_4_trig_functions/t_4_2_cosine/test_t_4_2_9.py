"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.2 Cosine/4.2.9 trig^m (a+b cos^n+c cos^(2 n))^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c = symbols('a b c')

def test_integrate_4_Trig_functions_4_2_Cosine_4_2_9_trig_pow_m_a_plus_b_cos_pow_n_plus_c_cos_pow_2_n_pow_p_1():
    f = sin(x)**5/(a + b*cos(x) + c*cos(x)**2)
    F = b*cos(x)**2/(2*c**2) + b*(b**2 - 2*c*(a + c))*log(a + b*cos(x) + c*cos(x)**2)/(2*c**4) - cos(x)**3/(3*c) - (b**2 - c*(a + 2*c))*cos(x)/c**3 + (b**4 - 2*b**2*c*(2*a + c) + 2*c**2*(a + c)**2)*atanh((b + 2*c*cos(x))/sqrt(-4*a*c + b**2))/(c**4*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_9_trig_pow_m_a_plus_b_cos_pow_n_plus_c_cos_pow_2_n_pow_p_2():
    f = sin(x)**3/(a + b*cos(x) + c*cos(x)**2)
    F = -b*log(a + b*cos(x) + c*cos(x)**2)/(2*c**2) + cos(x)/c - (b**2 - 2*c*(a + c))*atanh((b + 2*c*cos(x))/sqrt(-4*a*c + b**2))/(c**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_9_trig_pow_m_a_plus_b_cos_pow_n_plus_c_cos_pow_2_n_pow_p_3():
    f = sin(x)/(a + b*cos(x) + c*cos(x)**2)
    F = 2*atanh((b + 2*c*cos(x))/sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_9_trig_pow_m_a_plus_b_cos_pow_n_plus_c_cos_pow_2_n_pow_p_4():
    f = csc(x)/(a + b*cos(x) + c*cos(x)**2)
    F = b*log(a + b*cos(x) + c*cos(x)**2)/((a + b + c)*(2*a - 2*b + 2*c)) + log(1 - cos(x))/(2*a + 2*b + 2*c) - log(cos(x) + 1)/(2*a - 2*b + 2*c) - (-2*a*c + b**2 - 2*c**2)*atanh((b + 2*c*cos(x))/sqrt(-4*a*c + b**2))/(sqrt(-4*a*c + b**2)*(a - b + c)*(a + b + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_9_trig_pow_m_a_plus_b_cos_pow_n_plus_c_cos_pow_2_n_pow_p_5():
    f = csc(x)**3/(a + b*cos(x) + c*cos(x)**2)
    F = -b*(b**2 - 2*c*(a + c))*log(a + b*cos(x) + c*cos(x)**2)/(2*(a**2 + 2*a*c - b**2 + c**2)**2) + (b - (a + c)*cos(x))*csc(x)**2/((a + b + c)*(2*a - 2*b + 2*c)) - (a - 2*b + 3*c)*log(cos(x) + 1)/(4*(a - b + c)**2) + (a + 2*b + 3*c)*log(1 - cos(x))/(4*(a + b + c)**2) + (b**4 - 2*b**2*c*(2*a + c) + 2*c**2*(a + c)**2)*atanh((b + 2*c*cos(x))/sqrt(-4*a*c + b**2))/(sqrt(-4*a*c + b**2)*(a**2 + 2*a*c - b**2 + c**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_9_trig_pow_m_a_plus_b_cos_pow_n_plus_c_cos_pow_2_n_pow_p_6():
    f = sin(x)**4/(a + b*cos(x) + c*cos(x)**2)
    F = -b*sin(x)/c**2 + x/(2*c) + sin(x)*cos(x)/(2*c) + x*(b**2 - c*(a + 2*c))/c**3 - (2*b*(b**2 - 2*c*(a + c)) - 2*(b**4 - 2*b**2*c*(2*a + c) + 2*c**2*(a + c)**2)/sqrt(-4*a*c + b**2))*atan(sqrt(b - 2*c - sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/(c**3*sqrt(b - 2*c - sqrt(-4*a*c + b**2))*sqrt(b + 2*c - sqrt(-4*a*c + b**2))) - (2*b**4 + 2*b**3*sqrt(-4*a*c + b**2) - 4*b**2*c*(2*a + c) - 4*b*c*(a + c)*sqrt(-4*a*c + b**2) + 4*c**2*(a + c)**2)*atan(sqrt(b - 2*c + sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/(c**3*sqrt(-4*a*c + b**2)*sqrt(b - 2*c + sqrt(-4*a*c + b**2))*sqrt(b + 2*c + sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_9_trig_pow_m_a_plus_b_cos_pow_n_plus_c_cos_pow_2_n_pow_p_7():
    f = sin(x)**2/(a + b*cos(x) + c*cos(x)**2)
    F = -x/c + (2*b - 2*(b**2 - 2*c*(a + c))/sqrt(-4*a*c + b**2))*atan(sqrt(b - 2*c - sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/(c*sqrt(b - 2*c - sqrt(-4*a*c + b**2))*sqrt(b + 2*c - sqrt(-4*a*c + b**2))) + (2*b + 2*(b**2 - 2*c*(a + c))/sqrt(-4*a*c + b**2))*atan(sqrt(b - 2*c + sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/(c*sqrt(b - 2*c + sqrt(-4*a*c + b**2))*sqrt(b + 2*c + sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_9_trig_pow_m_a_plus_b_cos_pow_n_plus_c_cos_pow_2_n_pow_p_8():
    f = csc(x)**2/(a + b*cos(x) + c*cos(x)**2)
    F = -2*b*c*(1 - (b**2 - 2*c*(a + c))/(b*sqrt(-4*a*c + b**2)))*atan(sqrt(b - 2*c + sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/((a - b + c)*(a + b + c)*sqrt(b - 2*c + sqrt(-4*a*c + b**2))*sqrt(b + 2*c + sqrt(-4*a*c + b**2))) - 2*b*c*(1 + (b**2 - 2*c*(a + c))/(b*sqrt(-4*a*c + b**2)))*atan(sqrt(b - 2*c - sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/((a - b + c)*(a + b + c)*sqrt(b - 2*c - sqrt(-4*a*c + b**2))*sqrt(b + 2*c - sqrt(-4*a*c + b**2))) + sin(x)/((cos(x) + 1)*(2*a - 2*b + 2*c)) - sin(x)/((1 - cos(x))*(2*a + 2*b + 2*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_9_trig_pow_m_a_plus_b_cos_pow_n_plus_c_cos_pow_2_n_pow_p_9():
    f = sin(x)/(cos(x)**2 + cos(x) - 2)
    F = -log(1 - cos(x))/3 + log(cos(x) + 2)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_9_trig_pow_m_a_plus_b_cos_pow_n_plus_c_cos_pow_2_n_pow_p_10():
    f = sin(x)/(cos(x)**2 - 5*cos(x) + 4)
    F = log(1 - cos(x))/3 - log(4 - cos(x))/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_9_trig_pow_m_a_plus_b_cos_pow_n_plus_c_cos_pow_2_n_pow_p_11():
    f = sin(x)/(cos(x)**2 - 2*cos(x) + 3)
    F = sqrt(2)*atan(sqrt(2)*(1 - cos(x))/2)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_9_trig_pow_m_a_plus_b_cos_pow_n_plus_c_cos_pow_2_n_pow_p_12():
    f = sin(x)/(cos(x)**2 - 4*cos(x) + 13)**2
    F = (2 - cos(x))/(18*cos(x)**2 - 72*cos(x) + 234) - atan(cos(x)/3 + sympy.S(-2)/3)/54
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_9_trig_pow_m_a_plus_b_cos_pow_n_plus_c_cos_pow_2_n_pow_p_13():
    f = cos(x)**4/(a + b*cos(x) + c*cos(x)**2)
    F = -b*sin(x)/c**2 + x/(2*c) + sin(x)*cos(x)/(2*c) + x*(-a*c + b**2)/c**3 - (-4*a*b*c + 2*b**3 + 2*(2*a**2*c**2 - 4*a*b**2*c + b**4)/sqrt(-4*a*c + b**2))*atan(sqrt(b - 2*c + sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/(c**3*sqrt(b - 2*c + sqrt(-4*a*c + b**2))*sqrt(b + 2*c + sqrt(-4*a*c + b**2))) - (-4*a*b*c + 2*b**3 - 2*(2*a**2*c**2 - 4*a*b**2*c + b**4)/sqrt(-4*a*c + b**2))*atan(sqrt(b - 2*c - sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/(c**3*sqrt(b - 2*c - sqrt(-4*a*c + b**2))*sqrt(b + 2*c - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_9_trig_pow_m_a_plus_b_cos_pow_n_plus_c_cos_pow_2_n_pow_p_14():
    f = cos(x)**3/(a + b*cos(x) + c*cos(x)**2)
    F = -b*x/c**2 + sin(x)/c + (-6*a*b*c/sqrt(-4*a*c + b**2) - 2*a*c + 2*b**3/sqrt(-4*a*c + b**2) + 2*b**2)*atan(sqrt(b - 2*c + sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/(c**2*sqrt(b - 2*c + sqrt(-4*a*c + b**2))*sqrt(b + 2*c + sqrt(-4*a*c + b**2))) + (6*a*b*c/sqrt(-4*a*c + b**2) - 2*a*c - 2*b**3/sqrt(-4*a*c + b**2) + 2*b**2)*atan(sqrt(b - 2*c - sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/(c**2*sqrt(b - 2*c - sqrt(-4*a*c + b**2))*sqrt(b + 2*c - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_9_trig_pow_m_a_plus_b_cos_pow_n_plus_c_cos_pow_2_n_pow_p_15():
    f = cos(x)**2/(a + b*cos(x) + c*cos(x)**2)
    F = x/c - (2*b - 2*(-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(b - 2*c - sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/(c*sqrt(b - 2*c - sqrt(-4*a*c + b**2))*sqrt(b + 2*c - sqrt(-4*a*c + b**2))) - (2*b + 2*(-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(b - 2*c + sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/(c*sqrt(b - 2*c + sqrt(-4*a*c + b**2))*sqrt(b + 2*c + sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_9_trig_pow_m_a_plus_b_cos_pow_n_plus_c_cos_pow_2_n_pow_p_16():
    f = cos(x)/(a + b*cos(x) + c*cos(x)**2)
    F = (-2*b/sqrt(-4*a*c + b**2) + 2)*atan(sqrt(b - 2*c - sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/(sqrt(b - 2*c - sqrt(-4*a*c + b**2))*sqrt(b + 2*c - sqrt(-4*a*c + b**2))) + (2*b/sqrt(-4*a*c + b**2) + 2)*atan(sqrt(b - 2*c + sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/(sqrt(b - 2*c + sqrt(-4*a*c + b**2))*sqrt(b + 2*c + sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_9_trig_pow_m_a_plus_b_cos_pow_n_plus_c_cos_pow_2_n_pow_p_17():
    f = 1/(a + b*cos(x) + c*cos(x)**2)
    F = -4*c*atan(sqrt(b - 2*c + sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/(sqrt(-4*a*c + b**2)*sqrt(b - 2*c + sqrt(-4*a*c + b**2))*sqrt(b + 2*c + sqrt(-4*a*c + b**2))) + 4*c*atan(sqrt(b - 2*c - sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/(sqrt(-4*a*c + b**2)*sqrt(b - 2*c - sqrt(-4*a*c + b**2))*sqrt(b + 2*c - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_9_trig_pow_m_a_plus_b_cos_pow_n_plus_c_cos_pow_2_n_pow_p_18():
    f = sec(x)/(a + b*cos(x) + c*cos(x)**2)
    F = -2*c*(-b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(b - 2*c + sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/(a*sqrt(b - 2*c + sqrt(-4*a*c + b**2))*sqrt(b + 2*c + sqrt(-4*a*c + b**2))) - 2*c*(b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(b - 2*c - sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/(a*sqrt(b - 2*c - sqrt(-4*a*c + b**2))*sqrt(b + 2*c - sqrt(-4*a*c + b**2))) + atanh(sin(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_9_trig_pow_m_a_plus_b_cos_pow_n_plus_c_cos_pow_2_n_pow_p_19():
    f = sec(x)**2/(a + b*cos(x) + c*cos(x)**2)
    F = tan(x)/a + 2*b*c*(1 - (-2*a*c + b**2)/(b*sqrt(-4*a*c + b**2)))*atan(sqrt(b - 2*c + sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/(a**2*sqrt(b - 2*c + sqrt(-4*a*c + b**2))*sqrt(b + 2*c + sqrt(-4*a*c + b**2))) + 2*b*c*(1 + (-2*a*c + b**2)/(b*sqrt(-4*a*c + b**2)))*atan(sqrt(b - 2*c - sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/(a**2*sqrt(b - 2*c - sqrt(-4*a*c + b**2))*sqrt(b + 2*c - sqrt(-4*a*c + b**2))) - b*atanh(sin(x))/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_9_trig_pow_m_a_plus_b_cos_pow_n_plus_c_cos_pow_2_n_pow_p_20():
    f = sec(x)**3/(a + b*cos(x) + c*cos(x)**2)
    F = tan(x)*sec(x)/(2*a) + atanh(sin(x))/(2*a) - b*tan(x)/a**2 + 2*c*(-3*a*b*c + b**3 - sqrt(-4*a*c + b**2)*(-a*c + b**2))*atan(sqrt(b - 2*c + sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/(a**3*sqrt(-4*a*c + b**2)*sqrt(b - 2*c + sqrt(-4*a*c + b**2))*sqrt(b + 2*c + sqrt(-4*a*c + b**2))) - 2*c*(-3*a*b*c + b**3 + sqrt(-4*a*c + b**2)*(-a*c + b**2))*atan(sqrt(b - 2*c - sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/(a**3*sqrt(-4*a*c + b**2)*sqrt(b - 2*c - sqrt(-4*a*c + b**2))*sqrt(b + 2*c - sqrt(-4*a*c + b**2))) + (-a*c + b**2)*atanh(sin(x))/a**3
    assert integrate(f, x) == F

