"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.1 Sine/4.1.9 trig^m (a+b sin^n+c sin^(2 n))^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c = symbols('a b c')

def test_integrate_4_Trig_functions_4_1_Sine_4_1_9_trig_pow_m_a_plus_b_sin_pow_n_plus_c_sin_pow_2_n_pow_p_1():
    f = sin(x)**4/(a + b*sin(x) + c*sin(x)**2)
    F = b*cos(x)/c**2 + x/(2*c) - sin(x)*cos(x)/(2*c) + x*(-a*c + b**2)/c**3 - sqrt(2)*(-2*a*b*c + b**3 + (2*a**2*c**2 - 4*a*b**2*c + b**4)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*(2*c + (b + sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/(c**3*sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c))) - sqrt(2)*(-2*a*b*c + b**3 - (2*a**2*c**2 - 4*a*b**2*c + b**4)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*(2*c + (b - sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/(c**3*sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_9_trig_pow_m_a_plus_b_sin_pow_n_plus_c_sin_pow_2_n_pow_p_2():
    f = sin(x)**3/(a + b*sin(x) + c*sin(x)**2)
    F = -b*x/c**2 + sqrt(2)*b*(-3*a*c/sqrt(-4*a*c + b**2) - a*c/b + b**2/sqrt(-4*a*c + b**2) + b)*atan(sqrt(2)*(2*c + (b + sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/(c**2*sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c))) + sqrt(2)*b*(3*a*c/sqrt(-4*a*c + b**2) - a*c/b - b**2/sqrt(-4*a*c + b**2) + b)*atan(sqrt(2)*(2*c + (b - sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/(c**2*sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c))) - cos(x)/c
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_9_trig_pow_m_a_plus_b_sin_pow_n_plus_c_sin_pow_2_n_pow_p_3():
    f = sin(x)**2/(a + b*sin(x) + c*sin(x)**2)
    F = x/c - sqrt(2)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*(2*c + (b - sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/(c*sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c))) - sqrt(2)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*(2*c + (b + sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/(c*sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_9_trig_pow_m_a_plus_b_sin_pow_n_plus_c_sin_pow_2_n_pow_p_4():
    f = sin(x)/(a + b*sin(x) + c*sin(x)**2)
    F = sqrt(2)*(-b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(2)*(2*c + (b - sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c)) + sqrt(2)*(b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(2)*(2*c + (b + sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_9_trig_pow_m_a_plus_b_sin_pow_n_plus_c_sin_pow_2_n_pow_p_5():
    f = 1/(a + b*sin(x) + c*sin(x)**2)
    F = -2*sqrt(2)*c*atan(sqrt(2)*(2*c + (b + sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/(sqrt(-4*a*c + b**2)*sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c))) + 2*sqrt(2)*c*atan(sqrt(2)*(2*c + (b - sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/(sqrt(-4*a*c + b**2)*sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_9_trig_pow_m_a_plus_b_sin_pow_n_plus_c_sin_pow_2_n_pow_p_6():
    f = csc(x)/(a + b*sin(x) + c*sin(x)**2)
    F = -sqrt(2)*c*(-b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(2)*(2*c + (b + sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/(a*sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c))) - sqrt(2)*c*(b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(2)*(2*c + (b - sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/(a*sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c))) - atanh(cos(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_9_trig_pow_m_a_plus_b_sin_pow_n_plus_c_sin_pow_2_n_pow_p_7():
    f = csc(x)**2/(a + b*sin(x) + c*sin(x)**2)
    F = -cot(x)/a + sqrt(2)*b*c*(1 - (-2*a*c + b**2)/(b*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*(2*c + (b + sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/(a**2*sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c))) + sqrt(2)*b*c*(1 + (-2*a*c + b**2)/(b*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*(2*c + (b - sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/(a**2*sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c))) + b*atanh(cos(x))/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_9_trig_pow_m_a_plus_b_sin_pow_n_plus_c_sin_pow_2_n_pow_p_8():
    f = csc(x)**3/(a + b*sin(x) + c*sin(x)**2)
    F = -cot(x)*csc(x)/(2*a) - atanh(cos(x))/(2*a) + b*cot(x)/a**2 + sqrt(2)*c*(-3*a*b*c + b**3 - sqrt(-4*a*c + b**2)*(-a*c + b**2))*atan(sqrt(2)*(2*c + (b + sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/(a**3*sqrt(-4*a*c + b**2)*sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c))) - sqrt(2)*c*(-3*a*b*c + b**3 + sqrt(-4*a*c + b**2)*(-a*c + b**2))*atan(sqrt(2)*(2*c + (b - sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/(a**3*sqrt(-4*a*c + b**2)*sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c))) - (-a*c + b**2)*atanh(cos(x))/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_9_trig_pow_m_a_plus_b_sin_pow_n_plus_c_sin_pow_2_n_pow_p_9():
    f = cos(x)**3/(a + b*sin(x) + c*sin(x)**2)
    F = b*log(a + b*sin(x) + c*sin(x)**2)/(2*c**2) - sin(x)/c + (b**2 - 2*c*(a + c))*atanh((b + 2*c*sin(x))/sqrt(-4*a*c + b**2))/(c**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_9_trig_pow_m_a_plus_b_sin_pow_n_plus_c_sin_pow_2_n_pow_p_10():
    f = cos(x)**2/(a + b*sin(x) + c*sin(x)**2)
    F = -x/c - sqrt(2)*sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c))*atan(sqrt(2)*(2*c + (b - sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/(c*sqrt(-4*a*c + b**2)) + sqrt(2)*sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c))*atan(sqrt(2)*(2*c + (b + sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/(c*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_9_trig_pow_m_a_plus_b_sin_pow_n_plus_c_sin_pow_2_n_pow_p_11():
    f = cos(x)/(a + b*sin(x) + c*sin(x)**2)
    F = -2*atanh((b + 2*c*sin(x))/sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_9_trig_pow_m_a_plus_b_sin_pow_n_plus_c_sin_pow_2_n_pow_p_12():
    f = sec(x)/(a + b*sin(x) + c*sin(x)**2)
    F = -b*log(a + b*sin(x) + c*sin(x)**2)/((a + b + c)*(2*a - 2*b + 2*c)) - log(1 - sin(x))/(2*a + 2*b + 2*c) + log(sin(x) + 1)/(2*a - 2*b + 2*c) + (-2*a*c + b**2 - 2*c**2)*atanh((b + 2*c*sin(x))/sqrt(-4*a*c + b**2))/(sqrt(-4*a*c + b**2)*(a - b + c)*(a + b + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_9_trig_pow_m_a_plus_b_sin_pow_n_plus_c_sin_pow_2_n_pow_p_13():
    f = sec(x)**2/(a + b*sin(x) + c*sin(x)**2)
    F = -sqrt(2)*b*c*(1 - (b**2 - 2*c*(a + c))/(b*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*(2*c + (b + sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/((a - b + c)*(a + b + c)*sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c))) - sqrt(2)*b*c*(1 + (b**2 - 2*c*(a + c))/(b*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*(2*c + (b - sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/((a - b + c)*(a + b + c)*sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c))) - cos(x)/((sin(x) + 1)*(2*a - 2*b + 2*c)) + cos(x)/((1 - sin(x))*(2*a + 2*b + 2*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_9_trig_pow_m_a_plus_b_sin_pow_n_plus_c_sin_pow_2_n_pow_p_14():
    f = sec(x)**3/(a + b*sin(x) + c*sin(x)**2)
    F = b*(b**2 - 2*c*(a + c))*log(a + b*sin(x) + c*sin(x)**2)/(2*(a**2 + 2*a*c - b**2 + c**2)**2) - (b - (a + c)*sin(x))*sec(x)**2/((a + b + c)*(2*a - 2*b + 2*c)) + (a - 2*b + 3*c)*log(sin(x) + 1)/(4*(a - b + c)**2) - (a + 2*b + 3*c)*log(1 - sin(x))/(4*(a + b + c)**2) - (b**4 - 2*b**2*c*(2*a + c) + 2*c**2*(a + c)**2)*atanh((b + 2*c*sin(x))/sqrt(-4*a*c + b**2))/(sqrt(-4*a*c + b**2)*(a**2 + 2*a*c - b**2 + c**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_9_trig_pow_m_a_plus_b_sin_pow_n_plus_c_sin_pow_2_n_pow_p_15():
    f = cos(x)/(sin(x)**2 + sin(x) - 6)
    F = log(2 - sin(x))/5 - log(sin(x) + 3)/5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_9_trig_pow_m_a_plus_b_sin_pow_n_plus_c_sin_pow_2_n_pow_p_16():
    f = cos(x)/(sin(x)**2 - 3*sin(x) + 2)
    F = -log(1 - sin(x)) + log(2 - sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_9_trig_pow_m_a_plus_b_sin_pow_n_plus_c_sin_pow_2_n_pow_p_17():
    f = cos(x)/(sin(x)**2 + 4*sin(x) - 5)
    F = log(1 - sin(x))/6 - log(sin(x) + 5)/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_9_trig_pow_m_a_plus_b_sin_pow_n_plus_c_sin_pow_2_n_pow_p_18():
    f = cos(x)/(sin(x)**2 - 6*sin(x) + 10)
    F = atan(sin(x) - 3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_9_trig_pow_m_a_plus_b_sin_pow_n_plus_c_sin_pow_2_n_pow_p_19():
    f = cos(x)/(sin(x)**2 + 2*sin(x) + 2)
    F = atan(sin(x) + 1)
    assert integrate(f, x) == F

