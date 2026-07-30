"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.2 Cosine/4.2.1.3 (g tan)^p (a+b cos)^m.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, f, g, m, p = symbols('a b c d e f g m p')

def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_1():
    f = tan(x)**4/(a*cos(x) + a)
    F = tan(x)**3/(3*a) - tan(x)*sec(x)/(2*a) + atanh(sin(x))/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_2():
    f = tan(x)**3/(a*cos(x) + a)
    F = sec(x)**2/(2*a) - sec(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_3():
    f = tan(x)**2/(a*cos(x) + a)
    F = tan(x)/a - atanh(sin(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_4():
    f = tan(x)/(a*cos(x) + a)
    F = log(cos(x) + 1)/a - log(cos(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_5():
    f = cot(x)/(a*cos(x) + a)
    F = cot(x)*csc(x)/(2*a) - atanh(cos(x))/(2*a) - csc(x)**2/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_6():
    f = cot(x)**2/(a*cos(x) + a)
    F = -cot(x)**3/(3*a) + csc(x)**3/(3*a) - csc(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_7():
    f = cot(x)**3/(a*cos(x) + a)
    F = -cot(x)**4/(4*a) + cot(x)**3*csc(x)/(4*a) - 3*cot(x)*csc(x)/(8*a) + 3*atanh(cos(x))/(8*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_8():
    f = cot(x)**4/(a*cos(x) + a)
    F = -cot(x)**5/(5*a) + csc(x)**5/(5*a) - 2*csc(x)**3/(3*a) + csc(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_9():
    f = tan(3*x)/(cos(3*x) + 1)**2
    F = log(cos(3*x) + 1)/3 - log(cos(3*x))/3 - 1/(3*cos(3*x) + 3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_10():
    f = tan(x)**4/(a + b*cos(x))
    F = tan(x)*sec(x)**2/(3*a) - b*tan(x)*sec(x)/(2*a**2) - (4*a**2 - 3*b**2)*tan(x)/(3*a**3) + b*(3*a**2 - 2*b**2)*atanh(sin(x))/(2*a**4) + 2*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)*atan(sqrt(a - b)*tan(x/2)/sqrt(a + b))/a**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_11():
    f = tan(x)**3/(a + b*cos(x))
    F = sec(x)**2/(2*a) - b*sec(x)/a**2 - (a**2 - b**2)*log(a + b*cos(x))/a**3 + (a**2 - b**2)*log(cos(x))/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_12():
    f = tan(x)**2/(a + b*cos(x))
    F = tan(x)/a - b*atanh(sin(x))/a**2 - 2*sqrt(a - b)*sqrt(a + b)*atan(sqrt(a - b)*tan(x/2)/sqrt(a + b))/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_13():
    f = tan(x)/(a + b*cos(x))
    F = log(a + b*cos(x))/a - log(cos(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_14():
    f = cot(x)/(a + b*cos(x))
    F = -a*log(a + b*cos(x))/(a**2 - b**2) + log(1 - cos(x))/(2*a + 2*b) + log(cos(x) + 1)/(2*a - 2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_15():
    f = cot(x)**2/(a + b*cos(x))
    F = -2*a**2*atan(sqrt(a - b)*tan(x/2)/sqrt(a + b))/((a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) - a*cot(x)/(a**2 - b**2) + b*csc(x)/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_16():
    f = cot(x)**3/(a + b*cos(x))
    F = a**3*log(a + b*cos(x))/(a**2 - b**2)**2 - (a - b*cos(x))*csc(x)**2/(2*a**2 - 2*b**2) - (2*a + b)*log(1 - cos(x))/(4*(a + b)**2) - (2*a - b)*log(cos(x) + 1)/(4*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_17():
    f = cot(x)**4/(a + b*cos(x))
    F = 2*a**4*atan(sqrt(a - b)*tan(x/2)/sqrt(a + b))/((a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + a**3*cot(x)/(a**2 - b**2)**2 - a**2*b*csc(x)/(a**2 - b**2)**2 - a*cot(x)**3/(3*a**2 - 3*b**2) + b*csc(x)**3/(3*a**2 - 3*b**2) - b*csc(x)/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_18():
    f = cot(x)/sqrt(3 - cos(x))
    F = -sqrt(2)*atanh(sqrt(2)*sqrt(3 - cos(x))/2)/2 - atanh(sqrt(3 - cos(x))/2)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_19():
    f = sqrt(a + b*cos(x))*tan(x)
    F = 2*sqrt(a)*atanh(sqrt(a + b*cos(x))/sqrt(a)) - 2*sqrt(a + b*cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_20():
    f = tan(x)/sqrt(a + b*cos(x))
    F = 2*atanh(sqrt(a + b*cos(x))/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_21():
    f = sqrt(e*tan(c + d*x))/(a + b*cos(c + d*x))
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')((Integer(-1) * (sympy.sqrt(((Integer(-1) * Symbol('a')) + Symbol('b'))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), sympy.asin((sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt((Integer(1) + sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), Integer(-1)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt(((Integer(-1) * Symbol('a')) + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')((sympy.sqrt(((Integer(-1) * Symbol('a')) + Symbol('b'))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1))), sympy.asin((sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt((Integer(1) + sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), Integer(-1)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt(((Integer(-1) * Symbol('a')) + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_1_3_g_tan_pow_p_a_plus_b_cos_pow_m_22():
    f = (g*tan(e + f*x))**p*(a + b*cos(e + f*x))**m
    F = ((Symbol('g') * sympy.cot((Symbol('e') + (Symbol('f') * x)))))**(Symbol('p')) * ((Symbol('g') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))**(Symbol('p')) * sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m')) * (((Symbol('g') * sympy.cot((Symbol('e') + (Symbol('f') * x)))))**(Symbol('p')))**(Integer(-1))), x)
    assert integrate(f, x) == F

