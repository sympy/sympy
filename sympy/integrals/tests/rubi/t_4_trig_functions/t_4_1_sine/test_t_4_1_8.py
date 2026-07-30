"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.1 Sine/4.1.8 (a+b sin)^m (c+d trig)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, a, b, c = symbols('A B a b c')

def test_integrate_4_Trig_functions_4_1_Sine_4_1_8_a_plus_b_sin_pow_m_c_plus_d_trig_pow_n_1():
    f = (A + B*cos(x))/(a + b*sin(x))
    F = 2*A*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/sqrt(a**2 - b**2) + B*log(a + b*sin(x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_8_a_plus_b_sin_pow_m_c_plus_d_trig_pow_n_2():
    f = (A + B*cos(x))/(sin(x) + 1)
    F = -A*cos(x)/(sin(x) + 1) + B*log(sin(x) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_8_a_plus_b_sin_pow_m_c_plus_d_trig_pow_n_3():
    f = (A + B*cos(x))/(1 - sin(x))
    F = A*cos(x)/(1 - sin(x)) - B*log(1 - sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_8_a_plus_b_sin_pow_m_c_plus_d_trig_pow_n_4():
    f = (b + c + cos(x))/(a + b*sin(x))
    F = (2*b + 2*c)*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/sqrt(a**2 - b**2) + log(a + b*sin(x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_8_a_plus_b_sin_pow_m_c_plus_d_trig_pow_n_5():
    f = (b + c + cos(x))/(a - b*sin(x))
    F = -(2*b + 2*c)*atan((-a*tan(x/2) + b)/sqrt(a**2 - b**2))/sqrt(a**2 - b**2) - log(a - b*sin(x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_8_a_plus_b_sin_pow_m_c_plus_d_trig_pow_n_6():
    f = (A + B*tan(x))/(a + b*sin(x))
    F = 2*A*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/sqrt(a**2 - b**2) + B*a*log(a + b*sin(x))/(a**2 - b**2) - B*log(1 - sin(x))/(2*a + 2*b) - B*log(sin(x) + 1)/(2*a - 2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_8_a_plus_b_sin_pow_m_c_plus_d_trig_pow_n_7():
    f = (A + B*cot(x))/(a + b*sin(x))
    F = 2*A*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/sqrt(a**2 - b**2) - B*log(a + b*sin(x))/a + B*log(sin(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_8_a_plus_b_sin_pow_m_c_plus_d_trig_pow_n_8():
    f = (A + B*sec(x))/(a + b*sin(x))
    F = 2*A*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/sqrt(a**2 - b**2) - B*b*log(a + b*sin(x))/(a**2 - b**2) - B*log(1 - sin(x))/(2*a + 2*b) + B*log(sin(x) + 1)/(2*a - 2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_8_a_plus_b_sin_pow_m_c_plus_d_trig_pow_n_9():
    f = (A + B*csc(x))/(a + b*sin(x))
    F = -B*atanh(cos(x))/a + (2*A*a - 2*B*b)*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(a*sqrt(a**2 - b**2))
    assert integrate(f, x) == F

