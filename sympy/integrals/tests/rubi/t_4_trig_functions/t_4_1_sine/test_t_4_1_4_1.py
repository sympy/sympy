"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.1 Sine/4.1.4.1 (a+b sin)^m (A+B sin+C sin^2).m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, C, a, b, e, f, m = symbols('A B C a b e f m')

def test_integrate_4_Trig_functions_4_1_Sine_4_1_4_1_a_plus_b_sin_pow_m_A_plus_B_sin_plus_C_sin_pow_2_1():
    f = (m - (m + 2)*sin(e + f*x)**2 + 1)*sin(e + f*x)**m
    F = sin(e + f*x)**(m + 1)*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_4_1_a_plus_b_sin_pow_m_A_plus_B_sin_plus_C_sin_pow_2_2():
    f = (6 - 7*sin(e + f*x)**2)*sin(e + f*x)**5
    F = sin(e + f*x)**6*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_4_1_a_plus_b_sin_pow_m_A_plus_B_sin_plus_C_sin_pow_2_3():
    f = (5 - 6*sin(e + f*x)**2)*sin(e + f*x)**4
    F = sin(e + f*x)**5*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_4_1_a_plus_b_sin_pow_m_A_plus_B_sin_plus_C_sin_pow_2_4():
    f = (4 - 5*sin(e + f*x)**2)*sin(e + f*x)**3
    F = sin(e + f*x)**4*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_4_1_a_plus_b_sin_pow_m_A_plus_B_sin_plus_C_sin_pow_2_5():
    f = (3 - 4*sin(e + f*x)**2)*sin(e + f*x)**2
    F = sin(e + f*x)**3*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_4_1_a_plus_b_sin_pow_m_A_plus_B_sin_plus_C_sin_pow_2_6():
    f = (2 - 3*sin(e + f*x)**2)*sin(e + f*x)
    F = sin(e + f*x)**2*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_4_1_a_plus_b_sin_pow_m_A_plus_B_sin_plus_C_sin_pow_2_7():
    f = 1 - 2*sin(e + f*x)**2
    F = sin(e + f*x)*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_4_1_a_plus_b_sin_pow_m_A_plus_B_sin_plus_C_sin_pow_2_8():
    f = -sin(e + f*x)**2*csc(e + f*x)
    F = cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_4_1_a_plus_b_sin_pow_m_A_plus_B_sin_plus_C_sin_pow_2_9():
    f = -csc(e + f*x)**2
    F = cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_4_1_a_plus_b_sin_pow_m_A_plus_B_sin_plus_C_sin_pow_2_10():
    f = (sin(e + f*x)**2 - 2)*csc(e + f*x)**3
    F = cot(e + f*x)*csc(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_4_1_a_plus_b_sin_pow_m_A_plus_B_sin_plus_C_sin_pow_2_11():
    f = (2*sin(e + f*x)**2 - 3)*csc(e + f*x)**4
    F = cot(e + f*x)*csc(e + f*x)**2/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_4_1_a_plus_b_sin_pow_m_A_plus_B_sin_plus_C_sin_pow_2_12():
    f = (3*sin(e + f*x)**2 - 4)*csc(e + f*x)**5
    F = cot(e + f*x)*csc(e + f*x)**3/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_4_1_a_plus_b_sin_pow_m_A_plus_B_sin_plus_C_sin_pow_2_13():
    f = (A + C*sin(e + f*x)**2)*(a*sin(e + f*x) + a)**m
    F = -2**(m + sympy.S.Half)*(A*(m**2 + 3*m + 2) + C*(m**2 + m + 1))*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(-m + sympy.S(-1)/2)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - sin(e + f*x)/2)/(f*(m + 1)*(m + 2)) + C*(a*sin(e + f*x) + a)**m*cos(e + f*x)/(f*(m**2 + 3*m + 2)) - C*(a*sin(e + f*x) + a)**(m + 1)*cos(e + f*x)/(a*f*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_4_1_a_plus_b_sin_pow_m_A_plus_B_sin_plus_C_sin_pow_2_14():
    f = (a + b*sin(e + f*x))**m*(-A*sin(e + f*x)**2 + A)
    F = 4*sqrt(2)*A*(a + b*sin(e + f*x))**m*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S(-3)/2, -m, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, b*(1 - sin(e + f*x))/(a + b))/(f*((a + b*sin(e + f*x))/(a + b))**m*sqrt(sin(e + f*x) + 1)) - 4*sqrt(2)*A*(a + b*sin(e + f*x))**m*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S(-1)/2, -m, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, b*(1 - sin(e + f*x))/(a + b))/(f*((a + b*sin(e + f*x))/(a + b))**m*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_4_1_a_plus_b_sin_pow_m_A_plus_B_sin_plus_C_sin_pow_2_15():
    f = (A + C*sin(e + f*x)**2)*(a + b*sin(e + f*x))**m
    F = sqrt(2)*C*a*(a + b)*(a + b*sin(e + f*x))**m*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m - 1, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, b*(1 - sin(e + f*x))/(a + b))/(b**2*f*((a + b*sin(e + f*x))/(a + b))**m*(m + 2)*sqrt(sin(e + f*x) + 1)) - C*(a + b*sin(e + f*x))**(m + 1)*cos(e + f*x)/(b*f*(m + 2)) - sqrt(2)*(a + b*sin(e + f*x))**m*(C*a**2 + b**2*(A*(m + 2) + C*(m + 1)))*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, b*(1 - sin(e + f*x))/(a + b))/(b**2*f*((a + b*sin(e + f*x))/(a + b))**m*(m + 2)*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_4_1_a_plus_b_sin_pow_m_A_plus_B_sin_plus_C_sin_pow_2_16():
    f = (A + C*sin(e + f*x)**2)*sin(e + f*x)**5
    F = C*cos(e + f*x)**7/(7*f) - (A + C)*cos(e + f*x)/f - (A + 3*C)*cos(e + f*x)**5/(5*f) + (2*A + 3*C)*cos(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_4_1_a_plus_b_sin_pow_m_A_plus_B_sin_plus_C_sin_pow_2_17():
    f = (a*sin(e + f*x) + a)**m*(A + B*sin(e + f*x) + C*sin(e + f*x)**2)
    F = -2**(m + sympy.S.Half)*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(-m + sympy.S(-1)/2)*(A*(m**2 + 3*m + 2) + B*m*(m + 2) + C*(m**2 + m + 1))*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - sin(e + f*x)/2)/(f*(m + 1)*(m + 2)) - C*(a*sin(e + f*x) + a)**(m + 1)*cos(e + f*x)/(a*f*(m + 2)) + (-B*(m + 2) + C)*(a*sin(e + f*x) + a)**m*cos(e + f*x)/(f*(m + 1)*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_4_1_a_plus_b_sin_pow_m_A_plus_B_sin_plus_C_sin_pow_2_18():
    f = (a + b*sin(e + f*x))**m*(A + C*sin(e + f*x)**2 + (A + C)*sin(e + f*x))
    F = -4*sqrt(2)*C*(a + b*sin(e + f*x))**m*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S(-3)/2, -m, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, b*(1 - sin(e + f*x))/(a + b))/(f*((a + b*sin(e + f*x))/(a + b))**m*sqrt(sin(e + f*x) + 1)) - 2*sqrt(2)*(A - C)*(a + b*sin(e + f*x))**m*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S(-1)/2, -m, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, b*(1 - sin(e + f*x))/(a + b))/(f*((a + b*sin(e + f*x))/(a + b))**m*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_4_1_a_plus_b_sin_pow_m_A_plus_B_sin_plus_C_sin_pow_2_19():
    f = (a + b*sin(e + f*x))**m*(A + B*sin(e + f*x) + C*sin(e + f*x)**2)
    F = -C*(a + b*sin(e + f*x))**(m + 1)*cos(e + f*x)/(b*f*(m + 2)) + sqrt(2)*(a + b)*(a + b*sin(e + f*x))**m*(-B*b*(m + 2) + C*a)*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m - 1, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, b*(1 - sin(e + f*x))/(a + b))/(b**2*f*((a + b*sin(e + f*x))/(a + b))**m*(m + 2)*sqrt(sin(e + f*x) + 1)) - sqrt(2)*(a + b*sin(e + f*x))**m*(A*b**2*(m + 2) - B*a*b*(m + 2) + C*a**2 + C*b**2*(m + 1))*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, b*(1 - sin(e + f*x))/(a + b))/(b**2*f*((a + b*sin(e + f*x))/(a + b))**m*(m + 2)*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F

