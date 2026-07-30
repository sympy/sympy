"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.1 Sine/4.1.2.3 (g sin)^p (a+b sin)^m (c+d sin)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, a, b, c, d, e, f, g, m, n, p = symbols('A B a b c d e f g m n p')

def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_1():
    f = (a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)*sin(e + f*x)**3
    F = a**2*c*x/16 + a**2*c*sin(e + f*x)**5*cos(e + f*x)/(6*f) - a**2*c*sin(e + f*x)**3*cos(e + f*x)/(24*f) - a**2*c*sin(e + f*x)*cos(e + f*x)/(16*f) + a**2*c*cos(e + f*x)**5/(5*f) - a**2*c*cos(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_2():
    f = (a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)*sin(e + f*x)**2
    F = a**2*c*x/8 + a**2*c*sin(e + f*x)**3*cos(e + f*x)/(4*f) - a**2*c*sin(e + f*x)*cos(e + f*x)/(8*f) + a**2*c*cos(e + f*x)**5/(5*f) - a**2*c*cos(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_3():
    f = (a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)*sin(e + f*x)
    F = a**2*c*x/8 + a**2*c*sin(e + f*x)**3*cos(e + f*x)/(4*f) - a**2*c*sin(e + f*x)*cos(e + f*x)/(8*f) - a**2*c*cos(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_4():
    f = (a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)
    F = a**2*c*x/2 + a**2*c*sin(e + f*x)*cos(e + f*x)/(2*f) - a**2*c*cos(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_5():
    f = (a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)*csc(e + f*x)
    F = a**2*c*x/2 + a**2*c*sin(e + f*x)*cos(e + f*x)/(2*f) + a**2*c*cos(e + f*x)/f - a**2*c*atanh(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_6():
    f = (a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)*csc(e + f*x)**2
    F = -a**2*c*x + a**2*c*cos(e + f*x)/f - a**2*c*cot(e + f*x)/f - a**2*c*atanh(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_7():
    f = (a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)*csc(e + f*x)**3
    F = -a**2*c*x - a**2*c*cot(e + f*x)*csc(e + f*x)/(2*f) - a**2*c*cot(e + f*x)/f + a**2*c*atanh(cos(e + f*x))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_8():
    f = (a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)*csc(e + f*x)**4
    F = -a**2*c*cot(e + f*x)**3/(3*f) - a**2*c*cot(e + f*x)*csc(e + f*x)/(2*f) + a**2*c*atanh(cos(e + f*x))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_9():
    f = (a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)*csc(e + f*x)**5
    F = -a**2*c*cot(e + f*x)**3/(3*f) - a**2*c*cot(e + f*x)*csc(e + f*x)**3/(4*f) + a**2*c*cot(e + f*x)*csc(e + f*x)/(8*f) + a**2*c*atanh(cos(e + f*x))/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_10():
    f = (a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)*csc(e + f*x)**6
    F = -a**2*c*cot(e + f*x)**5/(5*f) - a**2*c*cot(e + f*x)**3/(3*f) - a**2*c*cot(e + f*x)*csc(e + f*x)**3/(4*f) + a**2*c*cot(e + f*x)*csc(e + f*x)/(8*f) + a**2*c*atanh(cos(e + f*x))/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_11():
    f = (a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)*csc(e + f*x)**7
    F = -a**2*c*cot(e + f*x)**5/(5*f) - a**2*c*cot(e + f*x)**3/(3*f) - a**2*c*cot(e + f*x)*csc(e + f*x)**5/(6*f) + a**2*c*cot(e + f*x)*csc(e + f*x)**3/(24*f) + a**2*c*cot(e + f*x)*csc(e + f*x)/(16*f) + a**2*c*atanh(cos(e + f*x))/(16*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_12():
    f = sqrt(a*sin(e + f*x) + a)/((-c*sin(e + f*x) + c)*sin(e + f*x))
    F = -2*sqrt(a)*atanh(sqrt(a)*cos(e + f*x)/sqrt(a*sin(e + f*x) + a))/(c*f) + 2*sqrt(a*sin(e + f*x) + a)*sec(e + f*x)/(c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_13():
    f = 1/(sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)*sin(e + f*x))
    F = sqrt(a*sin(e + f*x) + a)*sec(e + f*x)/(a*c*f) - 2*atanh(sqrt(a)*cos(e + f*x)/sqrt(a*sin(e + f*x) + a))/(sqrt(a)*c*f) + sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(2*sqrt(a)*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_14():
    f = sqrt(g*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)
    F = 2*sqrt(a)*sqrt(g)*atan(sqrt(a)*sqrt(g)*cos(e + f*x)/(sqrt(g*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(c*f) + 2*sqrt(g*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)*sec(e + f*x)/(c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_15():
    f = sqrt(a*sin(e + f*x) + a)/(sqrt(g*sin(e + f*x))*(-c*sin(e + f*x) + c))
    F = 2*sqrt(g*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)*sec(e + f*x)/(c*f*g)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_16():
    f = sqrt(g*sin(e + f*x))/(sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c))
    F = sqrt(g*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)*sec(e + f*x)/(a*c*f) + sqrt(2)*sqrt(g)*atan(sqrt(2)*sqrt(a)*sqrt(g)*cos(e + f*x)/(2*sqrt(g*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(2*sqrt(a)*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_17():
    f = 1/(sqrt(g*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c))
    F = sqrt(g*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)*sec(e + f*x)/(a*c*f*g) - sqrt(2)*atan(sqrt(2)*sqrt(a)*sqrt(g)*cos(e + f*x)/(2*sqrt(g*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(2*sqrt(a)*c*f*sqrt(g))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_18():
    f = sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)/sin(e + f*x)
    F = sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)*log(sin(e + f*x))*sec(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_19():
    f = sqrt(a*sin(e + f*x) + a)/(sqrt(-c*sin(e + f*x) + c)*sin(e + f*x))
    F = -a*log(1 - sin(e + f*x))*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) + sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)*log(sin(e + f*x))*sec(e + f*x)/(c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_20():
    f = sqrt(-c*sin(e + f*x) + c)/(sqrt(a*sin(e + f*x) + a)*sin(e + f*x))
    F = -c*log(sin(e + f*x) + 1)*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) + sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)*log(sin(e + f*x))*sec(e + f*x)/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_21():
    f = 1/(sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)*sin(e + f*x))
    F = log(tan(e + f*x))*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_22():
    f = sqrt(a*sin(e + f*x) + a)/((c + d*sin(e + f*x))*sin(e + f*x))
    F = 2*sqrt(a)*sqrt(d)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(c*f*sqrt(c + d)) - 2*sqrt(a)*atanh(sqrt(a)*cos(e + f*x)/sqrt(a*sin(e + f*x) + a))/(c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_23():
    f = 1/((c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)*sin(e + f*x))
    F = sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f*(c - d)) - 2*d**(sympy.S(3)/2)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*c*f*(c - d)*sqrt(c + d)) - 2*atanh(sqrt(a)*cos(e + f*x)/sqrt(a*sin(e + f*x) + a))/(sqrt(a)*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_24():
    f = sqrt(g*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)/(c + d*sin(e + f*x))
    F = 2*sqrt(a)*sqrt(c)*sqrt(g)*atan(sqrt(a)*sqrt(c)*sqrt(g)*cos(e + f*x)/(sqrt(g*sin(e + f*x))*sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(d*f*sqrt(c + d)) - 2*sqrt(a)*sqrt(g)*atan(sqrt(a)*sqrt(g)*cos(e + f*x)/(sqrt(g*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_25():
    f = sqrt(a*sin(e + f*x) + a)/(sqrt(g*sin(e + f*x))*(c + d*sin(e + f*x)))
    F = -2*sqrt(a)*atan(sqrt(a)*sqrt(c)*sqrt(g)*cos(e + f*x)/(sqrt(g*sin(e + f*x))*sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(sqrt(c)*f*sqrt(g)*sqrt(c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_26():
    f = sqrt(g*sin(e + f*x))/((c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a))
    F = -2*sqrt(c)*sqrt(g)*atan(sqrt(a)*sqrt(c)*sqrt(g)*cos(e + f*x)/(sqrt(g*sin(e + f*x))*sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f*(c - d)*sqrt(c + d)) + sqrt(2)*sqrt(g)*atan(sqrt(2)*sqrt(a)*sqrt(g)*cos(e + f*x)/(2*sqrt(g*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f*(c - d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_27():
    f = 1/(sqrt(g*sin(e + f*x))*(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a))
    F = -sqrt(2)*atan(sqrt(2)*sqrt(a)*sqrt(g)*cos(e + f*x)/(2*sqrt(g*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f*sqrt(g)*(c - d)) + 2*d*atan(sqrt(a)*sqrt(c)*sqrt(g)*cos(e + f*x)/(sqrt(g*sin(e + f*x))*sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*sqrt(c)*f*sqrt(g)*(c - d)*sqrt(c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_28():
    f = sqrt(a + b*sin(e + f*x))/((c*sin(e + f*x) + c)*sin(e + f*x))
    F = ((sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('c') * Symbol('f') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('c') * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('c') * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('f') * (Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_29():
    f = 1/(sqrt(a + b*sin(e + f*x))*(c*sin(e + f*x) + c)*sin(e + f*x))
    F = ((sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('c') * Symbol('f') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('c') * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((Integer(2) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('c') * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('f') * (Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_30():
    f = sqrt(g*sin(e + f*x))*sqrt(a + b*sin(e + f*x))/(c*sin(e + f*x) + c)
    F = ((Integer(2) * sympy.sqrt(Symbol('g')) * sympy.Function('EllipticPi')((Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('g') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((sympy.sqrt(Symbol('g')) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('c') * Symbol('f')))**(Integer(-1))) + ((Symbol('g') * sympy.elliptic_e(sympy.asin((sympy.cos((Symbol('e') + (Symbol('f') * x))) * ((Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x)))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * sympy.sqrt((sympy.sin((Symbol('e') + (Symbol('f') * x))) * ((Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('c') * Symbol('f') * sympy.sqrt((Symbol('g') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + Symbol('b')) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_31():
    f = sqrt(a + b*sin(e + f*x))/(sqrt(g*sin(e + f*x))*(c*sin(e + f*x) + c))
    F = -sqrt(sin(e + f*x)/(sin(e + f*x) + 1))*sqrt(a + b*sin(e + f*x))*elliptic_e(asin(cos(e + f*x)/(sin(e + f*x) + 1)), -(a - b)/(a + b))/(c*f*sqrt(g*sin(e + f*x))*sqrt((a + b*sin(e + f*x))/((a + b)*(sin(e + f*x) + 1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_32():
    f = sqrt(g*sin(e + f*x))/(sqrt(a + b*sin(e + f*x))*(c*sin(e + f*x) + c))
    F = -2*sqrt(g)*sqrt(a*(1 - csc(e + f*x))/(a + b))*sqrt(a*(csc(e + f*x) + 1)/(a - b))*sqrt(a + b)*tan(e + f*x)*elliptic_f(asin(sqrt(g)*sqrt(a + b*sin(e + f*x))/(sqrt(g*sin(e + f*x))*sqrt(a + b))), -(a + b)/(a - b))/(c*f*(a - b)) + g*sqrt(sin(e + f*x)/(sin(e + f*x) + 1))*sqrt(a + b*sin(e + f*x))*elliptic_e(asin(cos(e + f*x)/(sin(e + f*x) + 1)), -(a - b)/(a + b))/(c*f*sqrt(g*sin(e + f*x))*sqrt((a + b*sin(e + f*x))/((a + b)*(sin(e + f*x) + 1)))*(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_33():
    f = 1/(sqrt(g*sin(e + f*x))*sqrt(a + b*sin(e + f*x))*(c*sin(e + f*x) + c))
    F = -sqrt(sin(e + f*x)/(sin(e + f*x) + 1))*sqrt(a + b*sin(e + f*x))*elliptic_e(asin(cos(e + f*x)/(sin(e + f*x) + 1)), -(a - b)/(a + b))/(c*f*sqrt(g*sin(e + f*x))*sqrt((a + b*sin(e + f*x))/((a + b)*(sin(e + f*x) + 1)))*(a - b)) + 2*b*sqrt(a*(1 - csc(e + f*x))/(a + b))*sqrt(a*(csc(e + f*x) + 1)/(a - b))*sqrt(a + b)*tan(e + f*x)*elliptic_f(asin(sqrt(g)*sqrt(a + b*sin(e + f*x))/(sqrt(g*sin(e + f*x))*sqrt(a + b))), -(a + b)/(a - b))/(a*c*f*sqrt(g)*(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_34():
    f = sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)/sin(e + f*x)
    F = -2*sqrt(a)*sqrt(c)*atanh(sqrt(a)*sqrt(c)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/f - 2*sqrt(a)*sqrt(d)*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_35():
    f = sqrt(a*sin(e + f*x) + a)/(sqrt(c + d*sin(e + f*x))*sin(e + f*x))
    F = -2*sqrt(a)*atanh(sqrt(a)*sqrt(c)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_36():
    f = sqrt(c + d*sin(e + f*x))/(sqrt(a*sin(e + f*x) + a)*sin(e + f*x))
    F = -2*sqrt(c)*atanh(sqrt(a)*sqrt(c)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f) + sqrt(2)*sqrt(c - d)*atanh(sqrt(2)*sqrt(a)*sqrt(c - d)*cos(e + f*x)/(2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_37():
    f = 1/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)*sin(e + f*x))
    F = sqrt(2)*atanh(sqrt(2)*sqrt(a)*sqrt(c - d)*cos(e + f*x)/(2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f*sqrt(c - d)) - 2*atanh(sqrt(a)*sqrt(c)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_38():
    f = sin(e + f*x)**2/((a + b*sin(e + f*x))**2*(c + d*sin(e + f*x)))
    F = a**2*cos(e + f*x)/(f*(a + b*sin(e + f*x))*(a**2 - b**2)*(-a*d + b*c)) - 2*a*(a**2*c + a*b*d - 2*b**2*c)*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(f*(a**2 - b**2)**(sympy.S(3)/2)*(-a*d + b*c)**2) + 2*c**2*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*sqrt(c**2 - d**2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_39():
    f = sqrt(c + d*sin(e + f*x))/((a + b*sin(e + f*x))*sin(e + f*x))
    F = ((Integer(2) * Symbol('c') * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Symbol('a') * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Symbol('a') * (Symbol('a') + Symbol('b')) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_40():
    f = 1/((a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))*sin(e + f*x))
    F = ((Integer(2) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Symbol('a') * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Symbol('a') * (Symbol('a') + Symbol('b')) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_41():
    f = sqrt(g*sin(e + f*x))*sqrt(a + b*sin(e + f*x))/(c + d*sin(e + f*x))
    F = ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(Symbol('g')) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.csc((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.csc((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin(((sympy.sqrt(Symbol('g')) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('g') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('d') * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(-1) * (sympy.cot((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.csc((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('c')) * ((Symbol('c') + Symbol('d')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.csc((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('d') * (Symbol('c') + Symbol('d')) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_42():
    f = sqrt(a + b*sin(e + f*x))/(sqrt(g*sin(e + f*x))*(c + d*sin(e + f*x)))
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.csc((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.csc((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt(Symbol('g')) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('g') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('c') * Symbol('f') * sympy.sqrt(Symbol('g'))))**(Integer(-1)))) + ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(-1) * (sympy.cot((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.csc((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('c')) * ((Symbol('c') + Symbol('d')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.csc((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('c') * (Symbol('c') + Symbol('d')) * Symbol('f') * Symbol('g') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_43():
    f = sqrt(g*sin(e + f*x))/(sqrt(a + b*sin(e + f*x))*(c + d*sin(e + f*x)))
    F = (Integer(2) * sympy.sqrt((Integer(-1) * (sympy.cot((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.csc((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('c')) * ((Symbol('c') + Symbol('d')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.csc((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * (((Symbol('c') + Symbol('d')) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_44():
    f = 1/(sqrt(g*sin(e + f*x))*sqrt(a + b*sin(e + f*x))*(c + d*sin(e + f*x)))
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.csc((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.csc((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt(Symbol('g')) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('g') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') * Symbol('c') * Symbol('f') * sympy.sqrt(Symbol('g'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * sympy.sqrt((Integer(-1) * (sympy.cot((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.csc((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('c')) * ((Symbol('c') + Symbol('d')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.csc((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('c') * (Symbol('c') + Symbol('d')) * Symbol('f') * Symbol('g') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_45():
    f = sqrt(g*sin(e + f*x))*sqrt(c + d*sin(e + f*x))/(a + b*sin(e + f*x))
    F = ((Integer(2) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt(Symbol('g')) * sympy.sqrt(((Symbol('c') * (Integer(1) + (Integer(-1) * sympy.csc((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') * (Integer(1) + sympy.csc((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Symbol('c') + Symbol('d')) * (Symbol('d'))**(Integer(-1))), sympy.asin(((sympy.sqrt(Symbol('g')) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('g') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('c') + Symbol('d')) * ((Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('b') * Symbol('f')))**(Integer(-1))) + ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(-1) * (sympy.cot((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))) * sympy.sqrt(((Symbol('d') + (Symbol('c') * sympy.csc((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.csc((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('c')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('b') * (Symbol('a') + Symbol('b')) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_46():
    f = sqrt(g*sin(e + f*x))/((a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x)))
    F = (Integer(2) * sympy.sqrt((Integer(-1) * (sympy.cot((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))) * sympy.sqrt(((Symbol('d') + (Symbol('c') * sympy.csc((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.csc((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('c')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * (((Symbol('a') + Symbol('b')) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_47():
    f = sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))/sin(e + f*x)
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.Function('EllipticPi')(((Symbol('a') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('c')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('f')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_48():
    f = sqrt(a + b*sin(e + f*x))/(sqrt(c + d*sin(e + f*x))*sin(e + f*x))
    F = Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.Function('EllipticPi')(((Symbol('a') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('c')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('c') * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_49():
    f = 1/(sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))*sin(e + f*x))
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.Function('EllipticPi')(((Symbol('a') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('c')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('c') * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') * sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_3_g_sin_pow_p_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_50():
    f = (A + B*sin(e + f*x))**p*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**n
    F = 2**(n + sympy.S.Half)*(1 - sin(e + f*x))**(sympy.S.Half - n)*(A + B*sin(e + f*x))**p*(a*sin(e + f*x) + a)**(m + 1)*(-c*sin(e + f*x) + c)**n*appellf1(m + sympy.S.Half, -p, sympy.S.Half - n, m + sympy.S(3)/2, -B*(sin(e + f*x) + 1)/(A - B), sin(e + f*x)/2 + sympy.S.Half)*sec(e + f*x)/(a*f*((A + B*sin(e + f*x))/(A - B))**p*(2*m + 1))
    assert integrate(f, x) == F

