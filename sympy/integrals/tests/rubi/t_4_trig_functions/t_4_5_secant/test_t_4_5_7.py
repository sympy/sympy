"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.5 Secant/4.5.7 (d trig)^m (a+b (c sec)^n)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n, p = symbols('a b c d e f m n p')

def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_1():
    f = (a + b*sec(e + f*x)**2)*sin(e + f*x)**7
    F = a*cos(e + f*x)**7/(7*f) + b*sec(e + f*x)/f - (a - 3*b)*cos(e + f*x)/f + (a - b)*cos(e + f*x)**3/f - (3*a - b)*cos(e + f*x)**5/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_2():
    f = (a + b*sec(e + f*x)**2)*sin(e + f*x)**5
    F = -a*cos(e + f*x)**5/(5*f) + b*sec(e + f*x)/f - (a - 2*b)*cos(e + f*x)/f + (2*a - b)*cos(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_3():
    f = (a + b*sec(e + f*x)**2)*sin(e + f*x)**3
    F = a*cos(e + f*x)**3/(3*f) + b*sec(e + f*x)/f - (a - b)*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_4():
    f = (a + b*sec(e + f*x)**2)*sin(e + f*x)
    F = -a*cos(e + f*x)/f + b*sec(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_5():
    f = (a + b*sec(e + f*x)**2)*csc(e + f*x)
    F = b*sec(e + f*x)/f - (a + b)*atanh(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_6():
    f = (a + b*sec(e + f*x)**2)*csc(e + f*x)**3
    F = b*sec(e + f*x)/f - (a + b)*cot(e + f*x)*csc(e + f*x)/(2*f) - (a + 3*b)*atanh(cos(e + f*x))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_7():
    f = (a + b*sec(e + f*x)**2)*csc(e + f*x)**5
    F = b*sec(e + f*x)/f + (-3*a - 15*b)*atanh(cos(e + f*x))/(8*f) - (a + b)*cot(e + f*x)*csc(e + f*x)**3/(4*f) - (3*a + 7*b)*cot(e + f*x)*csc(e + f*x)/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_8():
    f = (a + b*sec(e + f*x)**2)*sin(e + f*x)**6
    F = -a*sin(e + f*x)*cos(e + f*x)**5/(6*f) + b*tan(e + f*x)/f + x*(5*a - 30*b)/16 - (11*a - 18*b)*sin(e + f*x)*cos(e + f*x)/(16*f) + (13*a - 6*b)*sin(e + f*x)*cos(e + f*x)**3/(24*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_9():
    f = (a + b*sec(e + f*x)**2)*sin(e + f*x)**4
    F = a*sin(e + f*x)*cos(e + f*x)**3/(4*f) + b*tan(e + f*x)/f + x*(3*a - 12*b)/8 - (5*a - 4*b)*sin(e + f*x)*cos(e + f*x)/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_10():
    f = (a + b*sec(e + f*x)**2)*sin(e + f*x)**2
    F = -a*sin(e + f*x)*cos(e + f*x)/(2*f) + b*tan(e + f*x)/f + x*(a - 2*b)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_11():
    f = a + b*sec(e + f*x)**2
    F = a*x + b*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_12():
    f = (a + b*sec(e + f*x)**2)*csc(e + f*x)**2
    F = b*tan(e + f*x)/f - (a + b)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_13():
    f = (a + b*sec(e + f*x)**2)*csc(e + f*x)**4
    F = b*tan(e + f*x)/f - (a + b)*cot(e + f*x)**3/(3*f) - (a + 2*b)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_14():
    f = (a + b*sec(e + f*x)**2)*csc(e + f*x)**6
    F = b*tan(e + f*x)/f - (a + b)*cot(e + f*x)**5/(5*f) - (a + 3*b)*cot(e + f*x)/f - (2*a + 3*b)*cot(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_15():
    f = (a + b*sec(e + f*x)**2)**2*sin(e + f*x)**5
    F = -a**2*cos(e + f*x)**5/(5*f) + 2*a*(a - b)*cos(e + f*x)**3/(3*f) + b**2*sec(e + f*x)**3/(3*f) + b*(2*a - 2*b)*sec(e + f*x)/f - (a**2 - 4*a*b + b**2)*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_16():
    f = (a + b*sec(e + f*x)**2)**2*sin(e + f*x)**3
    F = a**2*cos(e + f*x)**3/(3*f) - a*(a - 2*b)*cos(e + f*x)/f + b**2*sec(e + f*x)**3/(3*f) + b*(2*a - b)*sec(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_17():
    f = (a + b*sec(e + f*x)**2)**2*sin(e + f*x)
    F = -a**2*cos(e + f*x)/f + 2*a*b*sec(e + f*x)/f + b**2*sec(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_18():
    f = (a + b*sec(e + f*x)**2)**2*csc(e + f*x)
    F = b**2*sec(e + f*x)**3/(3*f) + b*(2*a + b)*sec(e + f*x)/f - (a + b)**2*atanh(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_19():
    f = (a + b*sec(e + f*x)**2)**2*csc(e + f*x)**3
    F = b**2*csc(e + f*x)**2*sec(e + f*x)**3/(3*f) + b*(6*a + 5*b)*sec(e + f*x)/(3*f) - (a + b)*(a + 5*b)*atanh(cos(e + f*x))/(2*f) - (3*a**2 + 6*a*b + 5*b**2)*cot(e + f*x)*csc(e + f*x)/(6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_20():
    f = (a + b*sec(e + f*x)**2)**2*csc(e + f*x)**5
    F = b**2*csc(e + f*x)**4*sec(e + f*x)**3/(3*f) + b*(6*a + 7*b)*sec(e + f*x)/(3*f) - (3*a + 7*b)**2*cot(e + f*x)*csc(e + f*x)/(24*f) - (3*a**2 + 6*a*b + 7*b**2)*cot(e + f*x)*csc(e + f*x)**3/(12*f) - (3*a**2 + 30*a*b + 35*b**2)*atanh(cos(e + f*x))/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_21():
    f = (a + b*sec(e + f*x)**2)**2*sin(e + f*x)**6
    F = a**2*sin(e + f*x)**6*tan(e + f*x)/(6*f) + a*(a - 12*b)*sin(e + f*x)*cos(e + f*x)**3/(24*f) + b**2*tan(e + f*x)**3/(3*f) + x*(5*a**2 - 60*a*b + 40*b**2)/16 - (a**2 - 12*a*b + 12*b**2)*tan(e + f*x)/(6*f) - (3*a**2 - 36*a*b + 8*b**2)*sin(e + f*x)*cos(e + f*x)/(16*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_22():
    f = (a + b*sec(e + f*x)**2)**2*sin(e + f*x)**4
    F = a**2*sin(e + f*x)**4*tan(e + f*x)/(4*f) - a*(a - 8*b)*sin(e + f*x)*cos(e + f*x)/(8*f) + b**2*tan(e + f*x)**3/(3*f) + x*(3*a**2 - 24*a*b + 8*b**2)/8 - (a**2 - 8*a*b + 4*b**2)*tan(e + f*x)/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_23():
    f = (a + b*sec(e + f*x)**2)**2*sin(e + f*x)**2
    F = a**2*sin(e + f*x)**2*tan(e + f*x)/(2*f) + a*x*(a - 4*b)/2 - a*(a - 4*b)*tan(e + f*x)/(2*f) + b**2*tan(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_24():
    f = (a + b*sec(e + f*x)**2)**2
    F = a**2*x + b**2*tan(e + f*x)**3/(3*f) + b*(2*a + b)*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_25():
    f = (a + b*sec(e + f*x)**2)**2*csc(e + f*x)**2
    F = b**2*tan(e + f*x)**3/(3*f) + 2*b*(a + b)*tan(e + f*x)/f - (a + b)**2*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_26():
    f = (a + b*sec(e + f*x)**2)**2*csc(e + f*x)**4
    F = b**2*tan(e + f*x)**3/(3*f) + b*(2*a + 3*b)*tan(e + f*x)/f - (a + b)**2*cot(e + f*x)**3/(3*f) - (a + b)*(a + 3*b)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_27():
    f = (a + b*sec(e + f*x)**2)**2*csc(e + f*x)**6
    F = b**2*tan(e + f*x)**3/(3*f) + 2*b*(a + 2*b)*tan(e + f*x)/f - (a + b)**2*cot(e + f*x)**5/(5*f) - (a + 2*b)*(2*a + 2*b)*cot(e + f*x)**3/(3*f) - (a**2 + 6*a*b + 6*b**2)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_28():
    f = sin(e + f*x)**5/(a + b*sec(e + f*x)**2)
    F = -cos(e + f*x)**5/(5*a*f) + (2*a + b)*cos(e + f*x)**3/(3*a**2*f) - (a + b)**2*cos(e + f*x)/(a**3*f) + sqrt(b)*(a + b)**2*atan(sqrt(a)*cos(e + f*x)/sqrt(b))/(a**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_29():
    f = sin(e + f*x)**3/(a + b*sec(e + f*x)**2)
    F = cos(e + f*x)**3/(3*a*f) - (a + b)*cos(e + f*x)/(a**2*f) + sqrt(b)*(a + b)*atan(sqrt(a)*cos(e + f*x)/sqrt(b))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_30():
    f = sin(e + f*x)/(a + b*sec(e + f*x)**2)
    F = -cos(e + f*x)/(a*f) + sqrt(b)*atan(sqrt(a)*cos(e + f*x)/sqrt(b))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_31():
    f = csc(e + f*x)/(a + b*sec(e + f*x)**2)
    F = -atanh(cos(e + f*x))/(f*(a + b)) + sqrt(b)*atan(sqrt(a)*cos(e + f*x)/sqrt(b))/(sqrt(a)*f*(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_32():
    f = csc(e + f*x)**3/(a + b*sec(e + f*x)**2)
    F = sqrt(a)*sqrt(b)*atan(sqrt(a)*cos(e + f*x)/sqrt(b))/(f*(a + b)**2) - (a - b)*atanh(cos(e + f*x))/(2*f*(a + b)**2) - cot(e + f*x)*csc(e + f*x)/(f*(2*a + 2*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_33():
    f = csc(e + f*x)**5/(a + b*sec(e + f*x)**2)
    F = a**(sympy.S(3)/2)*sqrt(b)*atan(sqrt(a)*cos(e + f*x)/sqrt(b))/(f*(a + b)**3) - cot(e + f*x)*csc(e + f*x)**3/(f*(4*a + 4*b)) - (3*a - b)*cot(e + f*x)*csc(e + f*x)/(8*f*(a + b)**2) - (3*a**2 - 6*a*b - b**2)*atanh(cos(e + f*x))/(8*f*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_34():
    f = sin(e + f*x)**6/(a + b*sec(e + f*x)**2)
    F = sin(e + f*x)**3*cos(e + f*x)**3/(6*a*f) + (3*a + 2*b)*sin(e + f*x)*cos(e + f*x)**3/(8*a**2*f) - (11*a**2 + 18*a*b + 8*b**2)*sin(e + f*x)*cos(e + f*x)/(16*a**3*f) - sqrt(b)*(a + b)**(sympy.S(5)/2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(a**4*f) + x*(5*a**3 + 30*a**2*b + 40*a*b**2 + 16*b**3)/(16*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_35():
    f = sin(e + f*x)**4/(a + b*sec(e + f*x)**2)
    F = sin(e + f*x)*cos(e + f*x)**3/(4*a*f) - (5*a + 4*b)*sin(e + f*x)*cos(e + f*x)/(8*a**2*f) - sqrt(b)*(a + b)**(sympy.S(3)/2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(a**3*f) + x*(3*a**2 + 12*a*b + 8*b**2)/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_36():
    f = sin(e + f*x)**2/(a + b*sec(e + f*x)**2)
    F = -sin(e + f*x)*cos(e + f*x)/(2*a*f) - sqrt(b)*sqrt(a + b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(a**2*f) + x*(a + 2*b)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_37():
    f = 1/(a + b*sec(e + f*x)**2)
    F = sqrt(b)*atan(sqrt(a + b)*cot(e + f*x)/sqrt(b))/(a*f*sqrt(a + b)) + x/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_38():
    f = csc(e + f*x)**2/(a + b*sec(e + f*x)**2)
    F = -sqrt(b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(f*(a + b)**(sympy.S(3)/2)) - cot(e + f*x)/(f*(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_39():
    f = csc(e + f*x)**4/(a + b*sec(e + f*x)**2)
    F = -a*sqrt(b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(f*(a + b)**(sympy.S(5)/2)) - a*cot(e + f*x)/(f*(a + b)**2) - cot(e + f*x)**3/(f*(3*a + 3*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_40():
    f = csc(e + f*x)**6/(a + b*sec(e + f*x)**2)
    F = -a**2*sqrt(b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(f*(a + b)**(sympy.S(7)/2)) - a**2*cot(e + f*x)/(f*(a + b)**3) - cot(e + f*x)**5/(f*(5*a + 5*b)) - (2*a + b)*cot(e + f*x)**3/(3*f*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_41():
    f = sin(e + f*x)**5/(a + b*sec(e + f*x)**2)**2
    F = -cos(e + f*x)**5/(5*a**2*f) - (a + b)**2*cos(e + f*x)**5/(2*a**2*b*f*(a*cos(e + f*x)**2 + b)) + (a + b)*(3*a + 7*b)*cos(e + f*x)**3/(6*a**3*b*f) - (a + b)*(3*a + 7*b)*cos(e + f*x)/(2*a**4*f) + sqrt(b)*(a + b)*(3*a + 7*b)*atan(sqrt(a)*cos(e + f*x)/sqrt(b))/(2*a**(sympy.S(9)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_42():
    f = sin(e + f*x)**3/(a + b*sec(e + f*x)**2)**2
    F = cos(e + f*x)**3/(3*a**2*f) - b*(a + b)*cos(e + f*x)/(2*a**3*f*(a*cos(e + f*x)**2 + b)) - (a + 2*b)*cos(e + f*x)/(a**3*f) + sqrt(b)*(3*a + 5*b)*atan(sqrt(a)*cos(e + f*x)/sqrt(b))/(2*a**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_43():
    f = sin(e + f*x)/(a + b*sec(e + f*x)**2)**2
    F = cos(e + f*x)**3/(2*a*f*(a*cos(e + f*x)**2 + b)) - 3*cos(e + f*x)/(2*a**2*f) + 3*sqrt(b)*atan(sqrt(a)*cos(e + f*x)/sqrt(b))/(2*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_44():
    f = csc(e + f*x)/(a + b*sec(e + f*x)**2)**2
    F = -atanh(cos(e + f*x))/(f*(a + b)**2) - b*cos(e + f*x)/(2*a*f*(a + b)*(a*cos(e + f*x)**2 + b)) + sqrt(b)*(3*a + b)*atan(sqrt(a)*cos(e + f*x)/sqrt(b))/(2*a**(sympy.S(3)/2)*f*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_45():
    f = csc(e + f*x)**3/(a + b*sec(e + f*x)**2)**2
    F = -(a - 3*b)*atanh(cos(e + f*x))/(2*f*(a + b)**3) + (a - b)*cos(e + f*x)/(2*f*(a + b)**2*(a*cos(e + f*x)**2 + b)) - cot(e + f*x)*csc(e + f*x)/(f*(2*a + 2*b)*(a*cos(e + f*x)**2 + b)) + sqrt(b)*(3*a - b)*atan(sqrt(a)*cos(e + f*x)/sqrt(b))/(2*sqrt(a)*f*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_46():
    f = csc(e + f*x)**5/(a + b*sec(e + f*x)**2)**2
    F = 3*sqrt(a)*sqrt(b)*(a - b)*atan(sqrt(a)*cos(e + f*x)/sqrt(b))/(2*f*(a + b)**4) + 3*a*(a - 3*b)*cos(e + f*x)/(8*f*(a + b)**3*(a*cos(e + f*x)**2 + b)) - (a - 5*b)*cot(e + f*x)*csc(e + f*x)/(8*f*(a + b)**2*(a*cos(e + f*x)**2 + b)) - cot(e + f*x)*csc(e + f*x)**3/(f*(4*a + 4*b)*(a*cos(e + f*x)**2 + b)) - (3*a**2 - 18*a*b + 3*b**2)*atanh(cos(e + f*x))/(8*f*(a + b)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_47():
    f = sin(e + f*x)**6/(a + b*sec(e + f*x)**2)**2
    F = sin(e + f*x)**3*cos(e + f*x)**3/(6*a*f*(a + b*tan(e + f*x)**2 + b)) + (9*a + 8*b)*sin(e + f*x)*cos(e + f*x)**3/(24*a**2*f*(a + b*tan(e + f*x)**2 + b)) - (33*a**2 + 82*a*b + 48*b**2)*sin(e + f*x)*cos(e + f*x)/(48*a**3*f*(a + b*tan(e + f*x)**2 + b)) - b*(19*a**2 + 52*a*b + 32*b**2)*tan(e + f*x)/(16*a**4*f*(a + b*tan(e + f*x)**2 + b)) - sqrt(b)*(a + b)**(sympy.S(3)/2)*(3*a + 8*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(2*a**5*f) + x*(5*a**3 + 60*a**2*b + 120*a*b**2 + 64*b**3)/(16*a**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_48():
    f = sin(e + f*x)**4/(a + b*sec(e + f*x)**2)**2
    F = sin(e + f*x)*cos(e + f*x)**3/(4*a*f*(a + b*tan(e + f*x)**2 + b)) - (5*a + 6*b)*sin(e + f*x)*cos(e + f*x)/(8*a**2*f*(a + b*tan(e + f*x)**2 + b)) - 3*b*(3*a + 4*b)*tan(e + f*x)/(8*a**3*f*(a + b*tan(e + f*x)**2 + b)) - 3*sqrt(b)*sqrt(a + b)*(a + 2*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(2*a**4*f) + x*(3*a**2 + 24*a*b + 24*b**2)/(8*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_49():
    f = sin(e + f*x)**2/(a + b*sec(e + f*x)**2)**2
    F = -sin(e + f*x)*cos(e + f*x)/(2*a*f*(a + b*tan(e + f*x)**2 + b)) - b*tan(e + f*x)/(a**2*f*(a + b*tan(e + f*x)**2 + b)) - sqrt(b)*(3*a + 4*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(2*a**3*f*sqrt(a + b)) + x*(a + 4*b)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_50():
    f = (a + b*sec(e + f*x)**2)**(-2)
    F = -b*tan(e + f*x)/(2*a*f*(a + b)*(a + b*tan(e + f*x)**2 + b)) - sqrt(b)*(3*a + 2*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(2*a**2*f*(a + b)**(sympy.S(3)/2)) + x/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_51():
    f = csc(e + f*x)**2/(a + b*sec(e + f*x)**2)**2
    F = -3*sqrt(b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(2*f*(a + b)**(sympy.S(5)/2)) + cot(e + f*x)/(f*(2*a + 2*b)*(a + b*tan(e + f*x)**2 + b)) - 3*cot(e + f*x)/(2*f*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_52():
    f = csc(e + f*x)**4/(a + b*sec(e + f*x)**2)**2
    F = -a*b*tan(e + f*x)/(2*f*(a + b)**3*(a + b*tan(e + f*x)**2 + b)) - sqrt(b)*(3*a - 2*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(2*f*(a + b)**(sympy.S(7)/2)) - (a - b)*cot(e + f*x)/(f*(a + b)**3) - cot(e + f*x)**3/(3*f*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_53():
    f = csc(e + f*x)**6/(a + b*sec(e + f*x)**2)**2
    F = -a*sqrt(b)*(3*a - 4*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(2*f*(a + b)**(sympy.S(9)/2)) - b*(5*a**2 + 2*b**2)*tan(e + f*x)/(10*f*(a + b)**4*(a + b*tan(e + f*x)**2 + b)) - cot(e + f*x)**5/(f*(5*a + 5*b)*(a + b*tan(e + f*x)**2 + b)) - (10*a + 3*b)*cot(e + f*x)**3/(15*f*(a + b)**3) - (5*a**2 - 10*a*b - b**2)*cot(e + f*x)/(5*f*(a + b)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_54():
    f = sin(e + f*x)**5/(a + b*sec(e + f*x)**2)**3
    F = -(a + b)**2*cos(e + f*x)**7/(4*a**2*b*f*(a*cos(e + f*x)**2 + b)**2) - cos(e + f*x)**5/(5*a**3*f) + (a + 3*b)*(3*a + 5*b)*cos(e + f*x)**3/(12*a**4*b*f) - b*(a + b)*(3*a + 11*b)*cos(e + f*x)/(8*a**5*f*(a*cos(e + f*x)**2 + b)) - (3*a**2 + 14*a*b + 13*b**2)*cos(e + f*x)/(2*a**5*f) + sqrt(b)*(15*a**2 + 70*a*b + 63*b**2)*atan(sqrt(a)*cos(e + f*x)/sqrt(b))/(8*a**(sympy.S(11)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_55():
    f = sin(e + f*x)**3/(a + b*sec(e + f*x)**2)**3
    F = cos(e + f*x)**3/(3*a**3*f) + b**2*(a + b)*cos(e + f*x)/(4*a**4*f*(a*cos(e + f*x)**2 + b)**2) - b*(9*a + 13*b)*cos(e + f*x)/(8*a**4*f*(a*cos(e + f*x)**2 + b)) - (a + 3*b)*cos(e + f*x)/(a**4*f) + 5*sqrt(b)*(3*a + 7*b)*atan(sqrt(a)*cos(e + f*x)/sqrt(b))/(8*a**(sympy.S(9)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_56():
    f = sin(e + f*x)/(a + b*sec(e + f*x)**2)**3
    F = cos(e + f*x)**5/(4*a*f*(a*cos(e + f*x)**2 + b)**2) + 5*cos(e + f*x)**3/(8*a**2*f*(a*cos(e + f*x)**2 + b)) - 15*cos(e + f*x)/(8*a**3*f) + 15*sqrt(b)*atan(sqrt(a)*cos(e + f*x)/sqrt(b))/(8*a**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_57():
    f = csc(e + f*x)/(a + b*sec(e + f*x)**2)**3
    F = -atanh(cos(e + f*x))/(f*(a + b)**3) - b*cos(e + f*x)**3/(4*a*f*(a + b)*(a*cos(e + f*x)**2 + b)**2) - b*(7*a + 3*b)*cos(e + f*x)/(8*a**2*f*(a + b)**2*(a*cos(e + f*x)**2 + b)) + sqrt(b)*(15*a**2 + 10*a*b + 3*b**2)*atan(sqrt(a)*cos(e + f*x)/sqrt(b))/(8*a**(sympy.S(5)/2)*f*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_58():
    f = csc(e + f*x)**3/(a + b*sec(e + f*x)**2)**3
    F = -(a - 5*b)*atanh(cos(e + f*x))/(2*f*(a + b)**4) - cos(e + f*x)*cot(e + f*x)**2/(f*(2*a + 2*b)*(a*cos(e + f*x)**2 + b)**2) - b*(2*a - b)*cos(e + f*x)/(4*a*f*(a + b)**2*(a*cos(e + f*x)**2 + b)**2) + (4*a**2 - 9*a*b - b**2)*cos(e + f*x)/(8*a*f*(a + b)**3*(a*cos(e + f*x)**2 + b)) + sqrt(b)*(15*a**2 - 10*a*b - b**2)*atan(sqrt(a)*cos(e + f*x)/sqrt(b))/(8*a**(sympy.S(3)/2)*f*(a + b)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_59():
    f = csc(e + f*x)**5/(a + b*sec(e + f*x)**2)**3
    F = -(a - 7*b)*cot(e + f*x)*csc(e + f*x)/(8*f*(a + b)**2*(a*cos(e + f*x)**2 + b)**2) - cot(e + f*x)**3*csc(e + f*x)/(f*(4*a + 4*b)*(a*cos(e + f*x)**2 + b)**2) + (a**2 - 9*a*b + 2*b**2)*cos(e + f*x)/(8*f*(a + b)**3*(a*cos(e + f*x)**2 + b)**2) + (3*a**2 - 18*a*b + 3*b**2)*cos(e + f*x)/(8*f*(a + b)**4*(a*cos(e + f*x)**2 + b)) - (3*a**2 - 30*a*b + 15*b**2)*atanh(cos(e + f*x))/(8*f*(a + b)**5) + 3*sqrt(b)*(5*a**2 - 10*a*b + b**2)*atan(sqrt(a)*cos(e + f*x)/sqrt(b))/(8*sqrt(a)*f*(a + b)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_60():
    f = sin(e + f*x)**6/(a + b*sec(e + f*x)**2)**3
    F = sin(e + f*x)**3*cos(e + f*x)**3/(6*a*f*(a + b*tan(e + f*x)**2 + b)**2) + (9*a + 10*b)*sin(e + f*x)*cos(e + f*x)**3/(24*a**2*f*(a + b*tan(e + f*x)**2 + b)**2) - (33*a**2 + 110*a*b + 80*b**2)*sin(e + f*x)*cos(e + f*x)/(48*a**3*f*(a + b*tan(e + f*x)**2 + b)**2) - 5*b*(9*a**2 + 32*a*b + 24*b**2)*tan(e + f*x)/(48*a**4*f*(a + b*tan(e + f*x)**2 + b)**2) - 5*b*(5*a**2 + 20*a*b + 16*b**2)*tan(e + f*x)/(16*a**5*f*(a + b*tan(e + f*x)**2 + b)) - 5*sqrt(b)*sqrt(a + b)*(a + 4*b)*(3*a + 4*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(8*a**6*f) + x*(5*a + 10*b)*(a**2 + 16*a*b + 16*b**2)/(16*a**6)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_61():
    f = sin(e + f*x)**4/(a + b*sec(e + f*x)**2)**3
    F = sin(e + f*x)*cos(e + f*x)**3/(4*a*f*(a + b*tan(e + f*x)**2 + b)**2) - (5*a + 8*b)*sin(e + f*x)*cos(e + f*x)/(8*a**2*f*(a + b*tan(e + f*x)**2 + b)**2) - b*(7*a + 12*b)*tan(e + f*x)/(8*a**3*f*(a + b*tan(e + f*x)**2 + b)**2) - 3*b*(a + 2*b)*tan(e + f*x)/(2*a**4*f*(a + b*tan(e + f*x)**2 + b)) - 3*sqrt(b)*(5*a**2 + 20*a*b + 16*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(8*a**5*f*sqrt(a + b)) + x*(3*a**2 + 36*a*b + 48*b**2)/(8*a**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_62():
    f = sin(e + f*x)**2/(a + b*sec(e + f*x)**2)**3
    F = -sin(e + f*x)*cos(e + f*x)/(2*a*f*(a + b*tan(e + f*x)**2 + b)**2) - 3*b*tan(e + f*x)/(4*a**2*f*(a + b*tan(e + f*x)**2 + b)**2) - b*(11*a + 12*b)*tan(e + f*x)/(8*a**3*f*(a + b)*(a + b*tan(e + f*x)**2 + b)) - sqrt(b)*(15*a**2 + 40*a*b + 24*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(8*a**4*f*(a + b)**(sympy.S(3)/2)) + x*(a + 6*b)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_63():
    f = (a + b*sec(e + f*x)**2)**(-3)
    F = -b*tan(e + f*x)/(4*a*f*(a + b)*(a + b*tan(e + f*x)**2 + b)**2) - b*(7*a + 4*b)*tan(e + f*x)/(8*a**2*f*(a + b)**2*(a + b*tan(e + f*x)**2 + b)) - sqrt(b)*(15*a**2 + 20*a*b + 8*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(8*a**3*f*(a + b)**(sympy.S(5)/2)) + x/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_64():
    f = csc(e + f*x)**2/(a + b*sec(e + f*x)**2)**3
    F = -15*sqrt(b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(8*f*(a + b)**(sympy.S(7)/2)) + cot(e + f*x)/(f*(4*a + 4*b)*(a + b*tan(e + f*x)**2 + b)**2) + 5*cot(e + f*x)/(8*f*(a + b)**2*(a + b*tan(e + f*x)**2 + b)) - 15*cot(e + f*x)/(8*f*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_65():
    f = csc(e + f*x)**4/(a + b*sec(e + f*x)**2)**3
    F = -a*b*tan(e + f*x)/(4*f*(a + b)**3*(a + b*tan(e + f*x)**2 + b)**2) + sqrt(b)*(-15*a + 20*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(8*f*(a + b)**(sympy.S(9)/2)) - b*(7*a - 4*b)*tan(e + f*x)/(8*f*(a + b)**4*(a + b*tan(e + f*x)**2 + b)) - (a - 2*b)*cot(e + f*x)/(f*(a + b)**4) - cot(e + f*x)**3/(3*f*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_66():
    f = csc(e + f*x)**6/(a + b*sec(e + f*x)**2)**3
    F = -sqrt(b)*(15*a**2 - 40*a*b + 8*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(8*f*(a + b)**(sympy.S(11)/2)) - b*(5*a**2 + 4*b**2)*tan(e + f*x)/(20*f*(a + b)**4*(a + b*tan(e + f*x)**2 + b)**2) - b*(35*a**2 - 40*a*b + 24*b**2)*tan(e + f*x)/(40*f*(a + b)**5*(a + b*tan(e + f*x)**2 + b)) - cot(e + f*x)**5/(f*(5*a + 5*b)*(a + b*tan(e + f*x)**2 + b)**2) - (10*a + b)*cot(e + f*x)**3/(15*f*(a + b)**4) - (5*a**2 - 20*a*b + 2*b**2)*cot(e + f*x)/(5*f*(a + b)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_67():
    f = sqrt(a + b*sec(e + f*x)**2)*sin(e + f*x)**5
    F = sqrt(b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/f - sqrt(a + b*sec(e + f*x)**2)*cos(e + f*x)/f - (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*cos(e + f*x)**5/(5*a*f) + (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*(10*a + 2*b)*cos(e + f*x)**3/(15*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_68():
    f = sqrt(a + b*sec(e + f*x)**2)*sin(e + f*x)**3
    F = sqrt(b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/f - sqrt(a + b*sec(e + f*x)**2)*cos(e + f*x)/f + (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*cos(e + f*x)**3/(3*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_69():
    f = sqrt(a + b*sec(e + f*x)**2)*sin(e + f*x)
    F = sqrt(b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/f - sqrt(a + b*sec(e + f*x)**2)*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_70():
    f = sqrt(a + b*sec(e + f*x)**2)*csc(e + f*x)
    F = sqrt(b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/f - sqrt(a + b)*atanh(sqrt(a + b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_71():
    f = sqrt(a + b*sec(e + f*x)**2)*csc(e + f*x)**3
    F = sqrt(b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/f - sqrt(a + b*sec(e + f*x)**2)*cot(e + f*x)*csc(e + f*x)/(2*f) - (a + 2*b)*atanh(sqrt(a + b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/(2*f*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_72():
    f = sqrt(a + b*sec(e + f*x)**2)*csc(e + f*x)**5
    F = sqrt(b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/f - sqrt(a + b*sec(e + f*x)**2)*(3*a + 4*b)*cot(e + f*x)*csc(e + f*x)/(f*(8*a + 8*b)) - sqrt(a + b*sec(e + f*x)**2)*cot(e + f*x)*csc(e + f*x)**3/(4*f) - (3*a**2 + 12*a*b + 8*b**2)*atanh(sqrt(a + b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/(8*f*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_73():
    f = sqrt(a + b*sec(e + f*x)**2)*sin(e + f*x)**6
    F = sqrt(b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f - sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)**5*cos(e + f*x)/(6*f) - (5*a - b)*sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)**3*cos(e + f*x)/(24*a*f) - (a - b)*(5*a + b)*sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)/(16*a**2*f) + (5*a**3 - 15*a**2*b - 5*a*b**2 - b**3)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(16*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_74():
    f = sqrt(a + b*sec(e + f*x)**2)*sin(e + f*x)**4
    F = sqrt(b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f - sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)**3*cos(e + f*x)/(4*f) - (3*a - b)*sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)/(8*a*f) + (3*a**2 - 6*a*b - b**2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(8*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_75():
    f = sqrt(a + b*sec(e + f*x)**2)*sin(e + f*x)**2
    F = sqrt(b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f - sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)/(2*f) + (a - b)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_76():
    f = sqrt(a + b*sec(e + f*x)**2)
    F = sqrt(a)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f + sqrt(b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_77():
    f = sqrt(a + b*sec(e + f*x)**2)*csc(e + f*x)**2
    F = sqrt(b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f - sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_78():
    f = sqrt(a + b*sec(e + f*x)**2)*csc(e + f*x)**4
    F = sqrt(b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f - sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)/f - (a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)*cot(e + f*x)**3/(f*(3*a + 3*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_79():
    f = sqrt(a + b*sec(e + f*x)**2)*csc(e + f*x)**6
    F = sqrt(b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f - sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)/f - (a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)*cot(e + f*x)**5/(f*(5*a + 5*b)) - (10*a + 8*b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)*cot(e + f*x)**3/(15*f*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_80():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*sin(e + f*x)**5
    F = sqrt(b)*(3*a - 4*b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/(2*f) + b*sqrt(a + b*sec(e + f*x)**2)*(3*a - 4*b)*sec(e + f*x)/(2*a*f) - (a + b*sec(e + f*x)**2)**(sympy.S(5)/2)*cos(e + f*x)**5/(5*a*f) + 2*(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)*cos(e + f*x)**3/(3*a*f) - (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*(3*a - 4*b)*cos(e + f*x)/(3*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_81():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*sin(e + f*x)**3
    F = sqrt(b)*(3*a - 2*b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/(2*f) + b*sqrt(a + b*sec(e + f*x)**2)*(3*a - 2*b)*sec(e + f*x)/(2*a*f) + (a + b*sec(e + f*x)**2)**(sympy.S(5)/2)*cos(e + f*x)**3/(3*a*f) - (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*(3*a - 2*b)*cos(e + f*x)/(3*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_82():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*sin(e + f*x)
    F = 3*a*sqrt(b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/(2*f) + 3*b*sqrt(a + b*sec(e + f*x)**2)*sec(e + f*x)/(2*f) - (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_83():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*csc(e + f*x)
    F = sqrt(b)*(3*a + 2*b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/(2*f) + b*sqrt(a + b*sec(e + f*x)**2)*sec(e + f*x)/(2*f) - (a + b)**(sympy.S(3)/2)*atanh(sqrt(a + b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_84():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*csc(e + f*x)**3
    F = sqrt(b)*(3*a + 4*b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/(2*f) + b*sqrt(a + b*sec(e + f*x)**2)*sec(e + f*x)/f - sqrt(a + b)*(a + 4*b)*atanh(sqrt(a + b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/(2*f) - (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)*csc(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_85():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*csc(e + f*x)**5
    F = 3*sqrt(b)*(a + 2*b)*atanh(sqrt(b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/(2*f) - (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)*csc(e + f*x)**3/(4*f) - sqrt(a + b*sec(e + f*x)**2)*(3*a + 6*b)*csc(e + f*x)**2*sec(e + f*x)/(8*f) + sqrt(a + b*sec(e + f*x)**2)*(3*a + 12*b)*sec(e + f*x)/(8*f) - (3*a**2 + 24*a*b + 24*b**2)*atanh(sqrt(a + b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/(8*f*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_86():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*sin(e + f*x)**6
    F = sqrt(b)*(3*a - 5*b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*f) + (5*a - 3*b)*sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)**4*tan(e + f*x)/(24*f) - (a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)*sin(e + f*x)**5*cos(e + f*x)/(6*f) + sqrt(a + b*tan(e + f*x)**2 + b)*(5*a**2 - 40*a*b + 3*b**2)*sin(e + f*x)**2*tan(e + f*x)/(48*a*f) - sqrt(a + b*tan(e + f*x)**2 + b)*(5*a**2 - 26*a*b + b**2)*tan(e + f*x)/(16*a*f) + (5*a**3 - 45*a**2*b + 15*a*b**2 + b**3)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(16*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_87():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*sin(e + f*x)**4
    F = sqrt(b)*(3*a - 3*b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*f) - (3*a - 9*b)*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/(8*f) + (3*a - 3*b)*sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)**2*tan(e + f*x)/(8*f) - (a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)*sin(e + f*x)**3*cos(e + f*x)/(4*f) + (3*a**2 - 18*a*b + 3*b**2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(8*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_88():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*sin(e + f*x)**2
    F = sqrt(a)*(a - 3*b)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*f) + sqrt(b)*(3*a - b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*f) + b*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/f - (a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)*sin(e + f*x)*cos(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_89():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = a**(sympy.S(3)/2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f + sqrt(b)*(3*a + b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*f) + b*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_90():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*csc(e + f*x)**2
    F = 3*sqrt(b)*(a + b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*f) + 3*b*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/(2*f) - (a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_91():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*csc(e + f*x)**4
    F = sqrt(b)*(3*a + 5*b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*f) + b*(3*a + 5*b)*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/(f*(2*a + 2*b)) - (3*a + 5*b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)*cot(e + f*x)/(f*(3*a + 3*b)) - (a + b*tan(e + f*x)**2 + b)**(sympy.S(5)/2)*cot(e + f*x)**3/(f*(3*a + 3*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_92():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*csc(e + f*x)**6
    F = sqrt(b)*(3*a + 7*b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*f) + b*(3*a + 7*b)*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/(f*(2*a + 2*b)) - (a + b*tan(e + f*x)**2 + b)**(sympy.S(5)/2)*cot(e + f*x)**5/(f*(5*a + 5*b)) - (3*a + 7*b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)*cot(e + f*x)/(f*(3*a + 3*b)) - 2*(a + b*tan(e + f*x)**2 + b)**(sympy.S(5)/2)*cot(e + f*x)**3/(f*(3*a + 3*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_93():
    f = sin(e + f*x)**5/sqrt(a + b*sec(e + f*x)**2)
    F = -sqrt(a + b*sec(e + f*x)**2)*cos(e + f*x)**5/(5*a*f) + sqrt(a + b*sec(e + f*x)**2)*(10*a + 4*b)*cos(e + f*x)**3/(15*a**2*f) - sqrt(a + b*sec(e + f*x)**2)*(15*a**2 + 20*a*b + 8*b**2)*cos(e + f*x)/(15*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_94():
    f = sin(e + f*x)**3/sqrt(a + b*sec(e + f*x)**2)
    F = sqrt(a + b*sec(e + f*x)**2)*cos(e + f*x)**3/(3*a*f) - sqrt(a + b*sec(e + f*x)**2)*(3*a + 2*b)*cos(e + f*x)/(3*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_95():
    f = sin(e + f*x)/sqrt(a + b*sec(e + f*x)**2)
    F = -sqrt(a + b*sec(e + f*x)**2)*cos(e + f*x)/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_96():
    f = csc(e + f*x)/sqrt(a + b*sec(e + f*x)**2)
    F = -atanh(sqrt(a + b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/(f*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_97():
    f = csc(e + f*x)**3/sqrt(a + b*sec(e + f*x)**2)
    F = -a*atanh(sqrt(a + b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/(2*f*(a + b)**(sympy.S(3)/2)) - sqrt(a + b*sec(e + f*x)**2)*cot(e + f*x)*csc(e + f*x)/(f*(2*a + 2*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_98():
    f = csc(e + f*x)**5/sqrt(a + b*sec(e + f*x)**2)
    F = -3*a**2*atanh(sqrt(a + b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/(8*f*(a + b)**(sympy.S(5)/2)) - sqrt(a + b*sec(e + f*x)**2)*cot(e + f*x)**3*csc(e + f*x)/(f*(4*a + 4*b)) - sqrt(a + b*sec(e + f*x)**2)*(5*a + 2*b)*cot(e + f*x)*csc(e + f*x)/(8*f*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_99():
    f = sin(e + f*x)**6/sqrt(a + b*sec(e + f*x)**2)
    F = sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)**3*cos(e + f*x)**3/(6*a*f) + (9*a + 5*b)*sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)**3/(24*a**2*f) - sqrt(a + b*tan(e + f*x)**2 + b)*(33*a**2 + 40*a*b + 15*b**2)*sin(e + f*x)*cos(e + f*x)/(48*a**3*f) + 5*(a + b)**3*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(16*a**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_100():
    f = sin(e + f*x)**4/sqrt(a + b*sec(e + f*x)**2)
    F = sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)**3/(4*a*f) - (5*a + 3*b)*sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)/(8*a**2*f) + 3*(a + b)**2*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(8*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_101():
    f = sin(e + f*x)**2/sqrt(a + b*sec(e + f*x)**2)
    F = -sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)/(2*a*f) + (a + b)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_102():
    f = 1/sqrt(a + b*sec(e + f*x)**2)
    F = atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_103():
    f = csc(e + f*x)**2/sqrt(a + b*sec(e + f*x)**2)
    F = -sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)/(f*(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_104():
    f = csc(e + f*x)**4/sqrt(a + b*sec(e + f*x)**2)
    F = -sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)**3/(f*(3*a + 3*b)) - (3*a + b)*sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)/(3*f*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_105():
    f = csc(e + f*x)**6/sqrt(a + b*sec(e + f*x)**2)
    F = -sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)**5/(f*(5*a + 5*b)) - (10*a + 6*b)*sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)**3/(15*f*(a + b)**2) - sqrt(a + b*tan(e + f*x)**2 + b)*(15*a**2 + 10*a*b + 3*b**2)*cot(e + f*x)/(15*f*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_106():
    f = sin(e + f*x)**5/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -cos(e + f*x)**5/(5*a*f*sqrt(a + b*sec(e + f*x)**2)) + (10*a + 6*b)*cos(e + f*x)**3/(15*a**2*f*sqrt(a + b*sec(e + f*x)**2)) - (15*a**2 + 40*a*b + 24*b**2)*cos(e + f*x)/(15*a**3*f*sqrt(a + b*sec(e + f*x)**2)) - 2*b*(15*a**2 + 40*a*b + 24*b**2)*sec(e + f*x)/(15*a**4*f*sqrt(a + b*sec(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_107():
    f = sin(e + f*x)**3/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = cos(e + f*x)**3/(3*a*f*sqrt(a + b*sec(e + f*x)**2)) - (3*a + 4*b)*cos(e + f*x)/(3*a**2*f*sqrt(a + b*sec(e + f*x)**2)) - 2*b*(3*a + 4*b)*sec(e + f*x)/(3*a**3*f*sqrt(a + b*sec(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_108():
    f = sin(e + f*x)/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -cos(e + f*x)/(a*f*sqrt(a + b*sec(e + f*x)**2)) - 2*b*sec(e + f*x)/(a**2*f*sqrt(a + b*sec(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_109():
    f = csc(e + f*x)/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -atanh(sqrt(a + b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/(f*(a + b)**(sympy.S(3)/2)) - b*sec(e + f*x)/(a*f*(a + b)*sqrt(a + b*sec(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_110():
    f = csc(e + f*x)**3/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -3*b*sec(e + f*x)/(2*f*(a + b)**2*sqrt(a + b*sec(e + f*x)**2)) - (a - 2*b)*atanh(sqrt(a + b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/(2*f*(a + b)**(sympy.S(5)/2)) - cot(e + f*x)*csc(e + f*x)/(f*sqrt(a + b*sec(e + f*x)**2)*(2*a + 2*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_111():
    f = csc(e + f*x)**5/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -3*a*(a - 4*b)*atanh(sqrt(a + b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/(8*f*(a + b)**(sympy.S(7)/2)) - 5*a*cot(e + f*x)*csc(e + f*x)/(8*f*(a + b)**2*sqrt(a + b*sec(e + f*x)**2)) - b*(13*a - 2*b)*sec(e + f*x)/(8*f*(a + b)**3*sqrt(a + b*sec(e + f*x)**2)) - cot(e + f*x)**3*csc(e + f*x)/(f*sqrt(a + b*sec(e + f*x)**2)*(4*a + 4*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_112():
    f = sin(e + f*x)**6/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = sin(e + f*x)**3*cos(e + f*x)**3/(6*a*f*sqrt(a + b*tan(e + f*x)**2 + b)) + (9*a + 7*b)*sin(e + f*x)*cos(e + f*x)**3/(24*a**2*f*sqrt(a + b*tan(e + f*x)**2 + b)) - (a + b)*(33*a + 35*b)*sin(e + f*x)*cos(e + f*x)/(48*a**3*f*sqrt(a + b*tan(e + f*x)**2 + b)) - b*(81*a**2 + 190*a*b + 105*b**2)*tan(e + f*x)/(48*a**4*f*sqrt(a + b*tan(e + f*x)**2 + b)) + 5*(a + b)**2*(a + 7*b)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(16*a**(sympy.S(9)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_113():
    f = sin(e + f*x)**4/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = sin(e + f*x)*cos(e + f*x)**3/(4*a*f*sqrt(a + b*tan(e + f*x)**2 + b)) - (5*a + 5*b)*sin(e + f*x)*cos(e + f*x)/(8*a**2*f*sqrt(a + b*tan(e + f*x)**2 + b)) - b*(13*a + 15*b)*tan(e + f*x)/(8*a**3*f*sqrt(a + b*tan(e + f*x)**2 + b)) + (a + 5*b)*(3*a + 3*b)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(8*a**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_114():
    f = sin(e + f*x)**2/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -sin(e + f*x)*cos(e + f*x)/(2*a*f*sqrt(a + b*tan(e + f*x)**2 + b)) - 3*b*tan(e + f*x)/(2*a**2*f*sqrt(a + b*tan(e + f*x)**2 + b)) + (a + 3*b)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_115():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(-3)/2)
    F = -b*tan(e + f*x)/(a*f*(a + b)*sqrt(a + b*tan(e + f*x)**2 + b)) + atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_116():
    f = csc(e + f*x)**2/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -2*b*tan(e + f*x)/(f*(a + b)**2*sqrt(a + b*tan(e + f*x)**2 + b)) - cot(e + f*x)/(f*(a + b)*sqrt(a + b*tan(e + f*x)**2 + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_117():
    f = csc(e + f*x)**4/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -b*(6*a - 2*b)*tan(e + f*x)/(3*f*(a + b)**3*sqrt(a + b*tan(e + f*x)**2 + b)) - cot(e + f*x)**3/(f*(3*a + 3*b)*sqrt(a + b*tan(e + f*x)**2 + b)) - (3*a - b)*cot(e + f*x)/(3*f*(a + b)**2*sqrt(a + b*tan(e + f*x)**2 + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_118():
    f = csc(e + f*x)**6/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -2*b*(15*a**2 - 10*a*b - b**2)*tan(e + f*x)/(15*f*(a + b)**4*sqrt(a + b*tan(e + f*x)**2 + b)) - cot(e + f*x)**5/(f*(5*a + 5*b)*sqrt(a + b*tan(e + f*x)**2 + b)) - (10*a + 4*b)*cot(e + f*x)**3/(15*f*(a + b)**2*sqrt(a + b*tan(e + f*x)**2 + b)) - (15*a**2 - 10*a*b - b**2)*cot(e + f*x)/(15*f*(a + b)**3*sqrt(a + b*tan(e + f*x)**2 + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_119():
    f = sin(e + f*x)**5/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -cos(e + f*x)**5/(5*a*f*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)) + (10*a + 8*b)*cos(e + f*x)**3/(15*a**2*f*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)) - (5*a**2 + 20*a*b + 16*b**2)*cos(e + f*x)/(5*a**3*f*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)) - 4*b*(5*a**2 + 20*a*b + 16*b**2)*sec(e + f*x)/(15*a**4*f*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)) - 8*b*(5*a**2 + 20*a*b + 16*b**2)*sec(e + f*x)/(15*a**5*f*sqrt(a + b*sec(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_120():
    f = sin(e + f*x)**3/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = cos(e + f*x)**3/(3*a*f*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)) - (a + 2*b)*cos(e + f*x)/(a**2*f*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)) - 4*b*(a + 2*b)*sec(e + f*x)/(3*a**3*f*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)) - 8*b*(a + 2*b)*sec(e + f*x)/(3*a**4*f*sqrt(a + b*sec(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_121():
    f = sin(e + f*x)/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -cos(e + f*x)/(a*f*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)) - 4*b*sec(e + f*x)/(3*a**2*f*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)) - 8*b*sec(e + f*x)/(3*a**3*f*sqrt(a + b*sec(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_122():
    f = csc(e + f*x)/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -atanh(sqrt(a + b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/(f*(a + b)**(sympy.S(5)/2)) - b*sec(e + f*x)/(3*a*f*(a + b)*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)) - b*(5*a + 2*b)*sec(e + f*x)/(3*a**2*f*(a + b)**2*sqrt(a + b*sec(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_123():
    f = csc(e + f*x)**3/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -5*b*sec(e + f*x)/(6*f*(a + b)**2*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)) - (a - 4*b)*atanh(sqrt(a + b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/(2*f*(a + b)**(sympy.S(7)/2)) - cot(e + f*x)*csc(e + f*x)/(f*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*(2*a + 2*b)) - b*(13*a - 2*b)*sec(e + f*x)/(6*a*f*(a + b)**3*sqrt(a + b*sec(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_124():
    f = csc(e + f*x)**5/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -b*(23*a - 12*b)*sec(e + f*x)/(24*f*(a + b)**3*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)) - b*(55*a - 50*b)*sec(e + f*x)/(24*f*(a + b)**4*sqrt(a + b*sec(e + f*x)**2)) - cot(e + f*x)**3*csc(e + f*x)/(f*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*(4*a + 4*b)) - (5*a - 2*b)*cot(e + f*x)*csc(e + f*x)/(8*f*(a + b)**2*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)) - (3*a**2 - 24*a*b + 8*b**2)*atanh(sqrt(a + b)*sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2))/(8*f*(a + b)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_125():
    f = sin(e + f*x)**6/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = sin(e + f*x)**3*cos(e + f*x)**3/(6*a*f*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) + (3*a + 3*b)*sin(e + f*x)*cos(e + f*x)**3/(8*a**2*f*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - (a + b)*(11*a + 21*b)*sin(e + f*x)*cos(e + f*x)/(16*a**3*f*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - 7*b*(a + b)*(7*a + 15*b)*tan(e + f*x)/(48*a**4*f*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - b*(113*a**2 + 420*a*b + 315*b**2)*tan(e + f*x)/(48*a**5*f*sqrt(a + b*tan(e + f*x)**2 + b)) + (5*a + 5*b)*(a**2 + 14*a*b + 21*b**2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(16*a**(sympy.S(11)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_126():
    f = sin(e + f*x)**4/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = sin(e + f*x)*cos(e + f*x)**3/(4*a*f*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - (5*a + 7*b)*sin(e + f*x)*cos(e + f*x)/(8*a**2*f*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - b*(23*a + 35*b)*tan(e + f*x)/(24*a**3*f*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - 5*b*(11*a + 21*b)*tan(e + f*x)/(24*a**4*f*sqrt(a + b*tan(e + f*x)**2 + b)) + (3*a**2 + 30*a*b + 35*b**2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(8*a**(sympy.S(9)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_127():
    f = sin(e + f*x)**2/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -sin(e + f*x)*cos(e + f*x)/(2*a*f*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - 5*b*tan(e + f*x)/(6*a**2*f*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - b*(13*a + 15*b)*tan(e + f*x)/(6*a**3*f*(a + b)*sqrt(a + b*tan(e + f*x)**2 + b)) + (a + 5*b)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*a**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_128():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(-5)/2)
    F = -b*tan(e + f*x)/(3*a*f*(a + b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - b*(5*a + 3*b)*tan(e + f*x)/(3*a**2*f*(a + b)**2*sqrt(a + b*tan(e + f*x)**2 + b)) + atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_129():
    f = csc(e + f*x)**2/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -4*b*tan(e + f*x)/(3*f*(a + b)**2*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - 8*b*tan(e + f*x)/(3*f*(a + b)**3*sqrt(a + b*tan(e + f*x)**2 + b)) - cot(e + f*x)/(f*(a + b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_130():
    f = csc(e + f*x)**4/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -b*(4*a - 4*b)*tan(e + f*x)/(3*f*(a + b)**3*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - b*(8*a - 8*b)*tan(e + f*x)/(3*f*(a + b)**4*sqrt(a + b*tan(e + f*x)**2 + b)) - (a - b)*cot(e + f*x)/(f*(a + b)**2*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - cot(e + f*x)**3/(f*(3*a + 3*b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_131():
    f = csc(e + f*x)**6/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -4*b*(5*a**2 - 10*a*b + b**2)*tan(e + f*x)/(15*f*(a + b)**4*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - 8*b*(5*a**2 - 10*a*b + b**2)*tan(e + f*x)/(15*f*(a + b)**5*sqrt(a + b*tan(e + f*x)**2 + b)) - cot(e + f*x)**5/(f*(5*a + 5*b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - (10*a + 2*b)*cot(e + f*x)**3/(15*f*(a + b)**2*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - (5*a**2 - 10*a*b + b**2)*cot(e + f*x)/(5*f*(a + b)**3*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_132():
    f = (d*sin(e + f*x))**m*(a + b*sec(e + f*x)**2)**p
    F = (d*sin(e + f*x))**m*(a + b*sec(e + f*x)**2)**p*(cos(e + f*x)**2)**(p + sympy.S.Half)*tan(e + f*x)*appellf1(m/2 + sympy.S.Half, -p, p + sympy.S.Half, m/2 + sympy.S(3)/2, a*sin(e + f*x)**2/(a + b), sin(e + f*x)**2)/(f*((-a*sin(e + f*x)**2 + a + b)/(a + b))**p*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_133():
    f = (a + b*sec(e + f*x)**2)**p*sin(e + f*x)**5
    F = -(a + b*sec(e + f*x)**2)**(p + 1)*cos(e + f*x)**5/(5*a*f) + (a + b*sec(e + f*x)**2)**(p + 1)*(10*a + b*(3 - 2*p))*cos(e + f*x)**3/(15*a**2*f) - (a + b*sec(e + f*x)**2)**p*(15*a**2 + 10*a*b*(1 - 2*p) + b**2*(4*p**2 - 8*p + 3))*cos(e + f*x)*hyper((sympy.S(-1)/2, -p), (sympy.S.Half,), -b*sec(e + f*x)**2/a)/(15*a**2*f*(1 + b*sec(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_134():
    f = (a + b*sec(e + f*x)**2)**p*sin(e + f*x)**3
    F = (a + b*sec(e + f*x)**2)**(p + 1)*cos(e + f*x)**3/(3*a*f) - (a + b*sec(e + f*x)**2)**p*(3*a - 2*b*p + b)*cos(e + f*x)*hyper((sympy.S(-1)/2, -p), (sympy.S.Half,), -b*sec(e + f*x)**2/a)/(3*a*f*(1 + b*sec(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_135():
    f = (a + b*sec(e + f*x)**2)**p*sin(e + f*x)
    F = -(a + b*sec(e + f*x)**2)**p*cos(e + f*x)*hyper((sympy.S(-1)/2, -p), (sympy.S.Half,), -b*sec(e + f*x)**2/a)/(f*(1 + b*sec(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_136():
    f = (a + b*sec(e + f*x)**2)**p*csc(e + f*x)
    F = -(a + b*sec(e + f*x)**2)**p*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, sec(e + f*x)**2, -b*sec(e + f*x)**2/a)*sec(e + f*x)/(f*(1 + b*sec(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_137():
    f = (a + b*sec(e + f*x)**2)**p*csc(e + f*x)**3
    F = (a + b*sec(e + f*x)**2)**p*appellf1(sympy.S(3)/2, 2, -p, sympy.S(5)/2, sec(e + f*x)**2, -b*sec(e + f*x)**2/a)*sec(e + f*x)**3/(3*f*(1 + b*sec(e + f*x)**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_138():
    f = (a + b*sec(e + f*x)**2)**p*sin(e + f*x)**4
    F = (a + b*tan(e + f*x)**2 + b)**p*tan(e + f*x)**5*appellf1(sympy.S(5)/2, 3, -p, sympy.S(7)/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/(a + b))/(5*f*(b*tan(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_139():
    f = (a + b*sec(e + f*x)**2)**p*sin(e + f*x)**2
    F = (a + b*tan(e + f*x)**2 + b)**p*tan(e + f*x)**3*appellf1(sympy.S(3)/2, 2, -p, sympy.S(5)/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/(a + b))/(3*f*(b*tan(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_140():
    f = (a + b*sec(e + f*x)**2)**p
    F = (a + b*tan(e + f*x)**2 + b)**p*tan(e + f*x)*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/(a + b))/(f*(b*tan(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_141():
    f = (a + b*sec(e + f*x)**2)**p*csc(e + f*x)**2
    F = -(a + b*tan(e + f*x)**2 + b)**p*cot(e + f*x)*hyper((sympy.S(-1)/2, -p), (sympy.S.Half,), -b*tan(e + f*x)**2/(a + b))/(f*(b*tan(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_142():
    f = (a + b*sec(e + f*x)**2)**p*csc(e + f*x)**4
    F = -(3*a + 2*b*(p + 1))*(a + b*tan(e + f*x)**2 + b)**p*cot(e + f*x)*hyper((sympy.S(-1)/2, -p), (sympy.S.Half,), -b*tan(e + f*x)**2/(a + b))/(f*(3*a + 3*b)*(b*tan(e + f*x)**2/(a + b) + 1)**p) - (a + b*tan(e + f*x)**2 + b)**(p + 1)*cot(e + f*x)**3/(f*(3*a + 3*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_143():
    f = (a + b*sec(e + f*x)**2)**p*csc(e + f*x)**6
    F = -(a + b*tan(e + f*x)**2 + b)**(p + 1)*cot(e + f*x)**5/(f*(5*a + 5*b)) - (10*a + b*(2*p + 7))*(a + b*tan(e + f*x)**2 + b)**(p + 1)*cot(e + f*x)**3/(15*f*(a + b)**2) - (a + b*tan(e + f*x)**2 + b)**p*(15*a**2 + 20*a*b*(p + 1) + 4*b**2*(p**2 + 3*p + 2))*cot(e + f*x)*hyper((sympy.S(-1)/2, -p), (sympy.S.Half,), -b*tan(e + f*x)**2/(a + b))/(15*f*(a + b)**2*(b*tan(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_144():
    f = (-a*sec(c + d*x)**2 + a)**4
    F = a**4*x + a**4*tan(c + d*x)**7/(7*d) - a**4*tan(c + d*x)**5/(5*d) + a**4*tan(c + d*x)**3/(3*d) - a**4*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_145():
    f = (-a*sec(c + d*x)**2 + a)**3
    F = a**3*x - a**3*tan(c + d*x)**5/(5*d) + a**3*tan(c + d*x)**3/(3*d) - a**3*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_146():
    f = (-a*sec(c + d*x)**2 + a)**2
    F = a**2*x + a**2*tan(c + d*x)**3/(3*d) - a**2*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_147():
    f = -a*sec(c + d*x)**2 + a
    F = a*x - a*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_148():
    f = 1/(-a*sec(c + d*x)**2 + a)
    F = x/a + cot(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_149():
    f = (-a*sec(c + d*x)**2 + a)**(-2)
    F = x/a**2 - cot(c + d*x)**3/(3*a**2*d) + cot(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_150():
    f = (-a*sec(c + d*x)**2 + a)**(-3)
    F = x/a**3 + cot(c + d*x)**5/(5*a**3*d) - cot(c + d*x)**3/(3*a**3*d) + cot(c + d*x)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_151():
    f = (-a*sec(c + d*x)**2 + a)**(-4)
    F = x/a**4 - cot(c + d*x)**7/(7*a**4*d) + cot(c + d*x)**5/(5*a**4*d) - cot(c + d*x)**3/(3*a**4*d) + cot(c + d*x)/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_152():
    f = (a + b*sec(e + f*x)**2)*sec(e + f*x)**5
    F = b*tan(e + f*x)*sec(e + f*x)**5/(6*f) + (6*a + 5*b)*tan(e + f*x)*sec(e + f*x)**3/(24*f) + (6*a + 5*b)*tan(e + f*x)*sec(e + f*x)/(16*f) + (6*a + 5*b)*atanh(sin(e + f*x))/(16*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_153():
    f = (a + b*sec(e + f*x)**2)*sec(e + f*x)**3
    F = b*tan(e + f*x)*sec(e + f*x)**3/(4*f) + (4*a + 3*b)*tan(e + f*x)*sec(e + f*x)/(8*f) + (4*a + 3*b)*atanh(sin(e + f*x))/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_154():
    f = (a + b*sec(e + f*x)**2)*sec(e + f*x)
    F = b*tan(e + f*x)*sec(e + f*x)/(2*f) + (2*a + b)*atanh(sin(e + f*x))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_155():
    f = (a + b*sec(e + f*x)**2)*cos(e + f*x)
    F = a*sin(e + f*x)/f + b*atanh(sin(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_156():
    f = (a + b*sec(e + f*x)**2)*cos(e + f*x)**3
    F = -a*sin(e + f*x)**3/(3*f) + (a + b)*sin(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_157():
    f = (a + b*sec(e + f*x)**2)*cos(e + f*x)**5
    F = a*sin(e + f*x)**5/(5*f) + (a + b)*sin(e + f*x)/f - (2*a + b)*sin(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_158():
    f = (a + b*sec(e + f*x)**2)*sec(e + f*x)**6
    F = b*tan(e + f*x)*sec(e + f*x)**6/(7*f) + (7*a + 6*b)*tan(e + f*x)**5/(35*f) + (7*a + 6*b)*tan(e + f*x)/(7*f) + (14*a + 12*b)*tan(e + f*x)**3/(21*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_159():
    f = (a + b*sec(e + f*x)**2)*sec(e + f*x)**4
    F = b*tan(e + f*x)*sec(e + f*x)**4/(5*f) + (5*a + 4*b)*tan(e + f*x)**3/(15*f) + (5*a + 4*b)*tan(e + f*x)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_160():
    f = (a + b*sec(e + f*x)**2)*sec(e + f*x)**2
    F = b*tan(e + f*x)*sec(e + f*x)**2/(3*f) + (3*a + 2*b)*tan(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_161():
    f = a + b*sec(e + f*x)**2
    F = a*x + b*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_162():
    f = (a + b*sec(e + f*x)**2)*cos(e + f*x)**2
    F = a*sin(e + f*x)*cos(e + f*x)/(2*f) + x*(a + 2*b)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_163():
    f = (a + b*sec(e + f*x)**2)*cos(e + f*x)**4
    F = a*sin(e + f*x)*cos(e + f*x)**3/(4*f) + x*(3*a + 4*b)/8 + (3*a + 4*b)*sin(e + f*x)*cos(e + f*x)/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_164():
    f = (a + b*sec(e + f*x)**2)*cos(e + f*x)**6
    F = a*sin(e + f*x)*cos(e + f*x)**5/(6*f) + x*(5*a + 6*b)/16 + (5*a + 6*b)*sin(e + f*x)*cos(e + f*x)**3/(24*f) + (5*a + 6*b)*sin(e + f*x)*cos(e + f*x)/(16*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_165():
    f = (a + b*sec(e + f*x)**2)**2*sec(e + f*x)**5
    F = b*(10*a + 7*b)*tan(e + f*x)*sec(e + f*x)**5/(48*f) + b*(-a*sin(e + f*x)**2 + a + b)*tan(e + f*x)*sec(e + f*x)**7/(8*f) + (48*a**2 + 80*a*b + 35*b**2)*tan(e + f*x)*sec(e + f*x)**3/(192*f) + (48*a**2 + 80*a*b + 35*b**2)*tan(e + f*x)*sec(e + f*x)/(128*f) + (48*a**2 + 80*a*b + 35*b**2)*atanh(sin(e + f*x))/(128*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_166():
    f = (a + b*sec(e + f*x)**2)**2*sec(e + f*x)**3
    F = b*(8*a + 5*b)*tan(e + f*x)*sec(e + f*x)**3/(24*f) + b*(-a*sin(e + f*x)**2 + a + b)*tan(e + f*x)*sec(e + f*x)**5/(6*f) + (8*a**2 + 12*a*b + 5*b**2)*tan(e + f*x)*sec(e + f*x)/(16*f) + (8*a**2 + 12*a*b + 5*b**2)*atanh(sin(e + f*x))/(16*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_167():
    f = (a + b*sec(e + f*x)**2)**2*sec(e + f*x)
    F = 3*b*(2*a + b)*tan(e + f*x)*sec(e + f*x)/(8*f) + b*(-a*sin(e + f*x)**2 + a + b)*tan(e + f*x)*sec(e + f*x)**3/(4*f) + (8*a**2 + 8*a*b + 3*b**2)*atanh(sin(e + f*x))/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_168():
    f = (a + b*sec(e + f*x)**2)**2*cos(e + f*x)
    F = a**2*sin(e + f*x)/f + b**2*tan(e + f*x)*sec(e + f*x)/(2*f) + b*(4*a + b)*atanh(sin(e + f*x))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_169():
    f = (a + b*sec(e + f*x)**2)**2*cos(e + f*x)**3
    F = -a**2*sin(e + f*x)**3/(3*f) + a*(a + 2*b)*sin(e + f*x)/f + b**2*atanh(sin(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_170():
    f = (a + b*sec(e + f*x)**2)**2*cos(e + f*x)**5
    F = a**2*sin(e + f*x)**5/(5*f) - 2*a*(a + b)*sin(e + f*x)**3/(3*f) + (a + b)**2*sin(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_171():
    f = (a + b*sec(e + f*x)**2)**2*sec(e + f*x)**6
    F = b**2*tan(e + f*x)**9/(9*f) + 2*b*(a + 2*b)*tan(e + f*x)**7/(7*f) + (a + b)**2*tan(e + f*x)/f + (a + 2*b)*(2*a + 2*b)*tan(e + f*x)**3/(3*f) + (a**2 + 6*a*b + 6*b**2)*tan(e + f*x)**5/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_172():
    f = (a + b*sec(e + f*x)**2)**2*sec(e + f*x)**4
    F = b**2*tan(e + f*x)**7/(7*f) + b*(2*a + 3*b)*tan(e + f*x)**5/(5*f) + (a + b)**2*tan(e + f*x)/f + (a + b)*(a + 3*b)*tan(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_173():
    f = (a + b*sec(e + f*x)**2)**2*sec(e + f*x)**2
    F = b**2*tan(e + f*x)**5/(5*f) + 2*b*(a + b)*tan(e + f*x)**3/(3*f) + (a + b)**2*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_174():
    f = (a + b*sec(e + f*x)**2)**2
    F = a**2*x + b**2*tan(e + f*x)**3/(3*f) + b*(2*a + b)*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_175():
    f = (a + b*sec(e + f*x)**2)**2*cos(e + f*x)**2
    F = a**2*sin(e + f*x)*cos(e + f*x)/(2*f) + a*x*(a + 4*b)/2 + b**2*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_176():
    f = (a + b*sec(e + f*x)**2)**2*cos(e + f*x)**4
    F = 3*a*(a + 2*b)*sin(e + f*x)*cos(e + f*x)/(8*f) + a*(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)**3/(4*f) + x*(3*a**2 + 8*a*b + 8*b**2)/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_177():
    f = (a + b*sec(e + f*x)**2)**2*cos(e + f*x)**6
    F = a*(5*a + 8*b)*sin(e + f*x)*cos(e + f*x)**3/(24*f) + a*(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)**5/(6*f) + x*(5*a**2 + 12*a*b + 8*b**2)/16 + (5*a**2 + 12*a*b + 8*b**2)*sin(e + f*x)*cos(e + f*x)/(16*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_178():
    f = (a + b*sec(c + d*x)**2)**3
    F = a**3*x + b**3*tan(c + d*x)**5/(5*d) + b**2*(3*a + 2*b)*tan(c + d*x)**3/(3*d) + b*(3*a**2 + 3*a*b + b**2)*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_179():
    f = (a + b*sec(c + d*x)**2)**4
    F = a**4*x + b**4*tan(c + d*x)**7/(7*d) + b**3*(4*a + 3*b)*tan(c + d*x)**5/(5*d) + b**2*(6*a**2 + 8*a*b + 3*b**2)*tan(c + d*x)**3/(3*d) + b*(2*a + b)*(2*a**2 + 2*a*b + b**2)*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_180():
    f = sec(e + f*x)**5/(a + b*sec(e + f*x)**2)
    F = a**(sympy.S(3)/2)*atanh(sqrt(a)*sin(e + f*x)/sqrt(a + b))/(b**2*f*sqrt(a + b)) + tan(e + f*x)*sec(e + f*x)/(2*b*f) - (2*a - b)*atanh(sin(e + f*x))/(2*b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_181():
    f = sec(e + f*x)**3/(a + b*sec(e + f*x)**2)
    F = -sqrt(a)*atanh(sqrt(a)*sin(e + f*x)/sqrt(a + b))/(b*f*sqrt(a + b)) + atanh(sin(e + f*x))/(b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_182():
    f = sec(e + f*x)/(a + b*sec(e + f*x)**2)
    F = atanh(sqrt(a)*sin(e + f*x)/sqrt(a + b))/(sqrt(a)*f*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_183():
    f = cos(e + f*x)/(a + b*sec(e + f*x)**2)
    F = sin(e + f*x)/(a*f) - b*atanh(sqrt(a)*sin(e + f*x)/sqrt(a + b))/(a**(sympy.S(3)/2)*f*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_184():
    f = cos(e + f*x)**3/(a + b*sec(e + f*x)**2)
    F = -sin(e + f*x)**3/(3*a*f) + (a - b)*sin(e + f*x)/(a**2*f) + b**2*atanh(sqrt(a)*sin(e + f*x)/sqrt(a + b))/(a**(sympy.S(5)/2)*f*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_185():
    f = cos(e + f*x)**5/(a + b*sec(e + f*x)**2)
    F = sin(e + f*x)**5/(5*a*f) - (2*a - b)*sin(e + f*x)**3/(3*a**2*f) + (a**2 - a*b + b**2)*sin(e + f*x)/(a**3*f) - b**3*atanh(sqrt(a)*sin(e + f*x)/sqrt(a + b))/(a**(sympy.S(7)/2)*f*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_186():
    f = sec(e + f*x)**6/(a + b*sec(e + f*x)**2)
    F = a**2*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(b**(sympy.S(5)/2)*f*sqrt(a + b)) + tan(e + f*x)**3/(3*b*f) - (a - b)*tan(e + f*x)/(b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_187():
    f = sec(e + f*x)**4/(a + b*sec(e + f*x)**2)
    F = -a*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(b**(sympy.S(3)/2)*f*sqrt(a + b)) + tan(e + f*x)/(b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_188():
    f = sec(e + f*x)**2/(a + b*sec(e + f*x)**2)
    F = atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(sqrt(b)*f*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_189():
    f = 1/(a + b*sec(e + f*x)**2)
    F = sqrt(b)*atan(sqrt(a + b)*cot(e + f*x)/sqrt(b))/(a*f*sqrt(a + b)) + x/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_190():
    f = cos(e + f*x)**2/(a + b*sec(e + f*x)**2)
    F = sin(e + f*x)*cos(e + f*x)/(2*a*f) + b**(sympy.S(3)/2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(a**2*f*sqrt(a + b)) + x*(a - 2*b)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_191():
    f = cos(e + f*x)**4/(a + b*sec(e + f*x)**2)
    F = sin(e + f*x)*cos(e + f*x)**3/(4*a*f) + (3*a - 4*b)*sin(e + f*x)*cos(e + f*x)/(8*a**2*f) - b**(sympy.S(5)/2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(a**3*f*sqrt(a + b)) + x*(3*a**2 - 4*a*b + 8*b**2)/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_192():
    f = cos(e + f*x)**6/(a + b*sec(e + f*x)**2)
    F = sin(e + f*x)*cos(e + f*x)**5/(6*a*f) + (5*a - 6*b)*sin(e + f*x)*cos(e + f*x)**3/(24*a**2*f) + (5*a**2 - 6*a*b + 8*b**2)*sin(e + f*x)*cos(e + f*x)/(16*a**3*f) + b**(sympy.S(7)/2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(a**4*f*sqrt(a + b)) + x*(5*a**3 - 6*a**2*b + 8*a*b**2 - 16*b**3)/(16*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_193():
    f = sec(e + f*x)**5/(a + b*sec(e + f*x)**2)**2
    F = -sqrt(a)*(2*a + 3*b)*atanh(sqrt(a)*sin(e + f*x)/sqrt(a + b))/(2*b**2*f*(a + b)**(sympy.S(3)/2)) - a*sin(e + f*x)/(2*b*f*(a + b)*(-a*sin(e + f*x)**2 + a + b)) + atanh(sin(e + f*x))/(b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_194():
    f = sec(e + f*x)**3/(a + b*sec(e + f*x)**2)**2
    F = sin(e + f*x)/(f*(2*a + 2*b)*(-a*sin(e + f*x)**2 + a + b)) + atanh(sqrt(a)*sin(e + f*x)/sqrt(a + b))/(2*sqrt(a)*f*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_195():
    f = sec(e + f*x)/(a + b*sec(e + f*x)**2)**2
    F = -b*sin(e + f*x)/(2*a*f*(a + b)*(-a*sin(e + f*x)**2 + a + b)) + (2*a + b)*atanh(sqrt(a)*sin(e + f*x)/sqrt(a + b))/(2*a**(sympy.S(3)/2)*f*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_196():
    f = cos(e + f*x)/(a + b*sec(e + f*x)**2)**2
    F = b**2*sin(e + f*x)/(2*a**2*f*(a + b)*(-a*sin(e + f*x)**2 + a + b)) + sin(e + f*x)/(a**2*f) - b*(4*a + 3*b)*atanh(sqrt(a)*sin(e + f*x)/sqrt(a + b))/(2*a**(sympy.S(5)/2)*f*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_197():
    f = cos(e + f*x)**3/(a + b*sec(e + f*x)**2)**2
    F = -sin(e + f*x)**3/(3*a**2*f) - b**3*sin(e + f*x)/(2*a**3*f*(a + b)*(-a*sin(e + f*x)**2 + a + b)) + (a - 2*b)*sin(e + f*x)/(a**3*f) + b**2*(6*a + 5*b)*atanh(sqrt(a)*sin(e + f*x)/sqrt(a + b))/(2*a**(sympy.S(7)/2)*f*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_198():
    f = cos(e + f*x)**5/(a + b*sec(e + f*x)**2)**2
    F = sin(e + f*x)**5/(5*a**2*f) - (2*a - 2*b)*sin(e + f*x)**3/(3*a**3*f) + b**4*sin(e + f*x)/(2*a**4*f*(a + b)*(-a*sin(e + f*x)**2 + a + b)) + (a**2 - 2*a*b + 3*b**2)*sin(e + f*x)/(a**4*f) - b**3*(8*a + 7*b)*atanh(sqrt(a)*sin(e + f*x)/sqrt(a + b))/(2*a**(sympy.S(9)/2)*f*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_199():
    f = sec(e + f*x)**6/(a + b*sec(e + f*x)**2)**2
    F = a**2*tan(e + f*x)/(2*b**2*f*(a + b)*(a + b*tan(e + f*x)**2 + b)) - a*(3*a + 4*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(2*b**(sympy.S(5)/2)*f*(a + b)**(sympy.S(3)/2)) + tan(e + f*x)/(b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_200():
    f = sec(e + f*x)**4/(a + b*sec(e + f*x)**2)**2
    F = -a*tan(e + f*x)/(2*b*f*(a + b)*(a + b*tan(e + f*x)**2 + b)) + (a + 2*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(2*b**(sympy.S(3)/2)*f*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_201():
    f = sec(e + f*x)**2/(a + b*sec(e + f*x)**2)**2
    F = tan(e + f*x)/(f*(2*a + 2*b)*(a + b*tan(e + f*x)**2 + b)) + atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(2*sqrt(b)*f*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_202():
    f = (a + b*sec(e + f*x)**2)**(-2)
    F = -b*tan(e + f*x)/(2*a*f*(a + b)*(a + b*tan(e + f*x)**2 + b)) - sqrt(b)*(3*a + 2*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(2*a**2*f*(a + b)**(sympy.S(3)/2)) + x/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_203():
    f = cos(e + f*x)**2/(a + b*sec(e + f*x)**2)**2
    F = sin(e + f*x)*cos(e + f*x)/(2*a*f*(a + b*tan(e + f*x)**2 + b)) + b*(a + 2*b)*tan(e + f*x)/(2*a**2*f*(a + b)*(a + b*tan(e + f*x)**2 + b)) + b**(sympy.S(3)/2)*(5*a + 4*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(2*a**3*f*(a + b)**(sympy.S(3)/2)) + x*(a - 4*b)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_204():
    f = cos(e + f*x)**4/(a + b*sec(e + f*x)**2)**2
    F = sin(e + f*x)*cos(e + f*x)**3/(4*a*f*(a + b*tan(e + f*x)**2 + b)) + (3*a - 6*b)*sin(e + f*x)*cos(e + f*x)/(8*a**2*f*(a + b*tan(e + f*x)**2 + b)) + b*(a - 3*b)*(3*a + 4*b)*tan(e + f*x)/(8*a**3*f*(a + b)*(a + b*tan(e + f*x)**2 + b)) - b**(sympy.S(5)/2)*(7*a + 6*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(2*a**4*f*(a + b)**(sympy.S(3)/2)) + x*(3*a**2 - 8*a*b + 24*b**2)/(8*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_205():
    f = cos(e + f*x)**6/(a + b*sec(e + f*x)**2)**2
    F = sin(e + f*x)*cos(e + f*x)**5/(6*a*f*(a + b*tan(e + f*x)**2 + b)) + (5*a - 8*b)*sin(e + f*x)*cos(e + f*x)**3/(24*a**2*f*(a + b*tan(e + f*x)**2 + b)) + (15*a**2 - 26*a*b + 48*b**2)*sin(e + f*x)*cos(e + f*x)/(48*a**3*f*(a + b*tan(e + f*x)**2 + b)) + b*(5*a**3 - 7*a**2*b + 12*a*b**2 + 32*b**3)*tan(e + f*x)/(16*a**4*f*(a + b)*(a + b*tan(e + f*x)**2 + b)) + b**(sympy.S(7)/2)*(9*a + 8*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(2*a**5*f*(a + b)**(sympy.S(3)/2)) + x*(5*a**3 - 12*a**2*b + 24*a*b**2 - 64*b**3)/(16*a**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_206():
    f = sec(e + f*x)**5/(a + b*sec(e + f*x)**2)**3
    F = sin(e + f*x)/(f*(4*a + 4*b)*(-a*sin(e + f*x)**2 + a + b)**2) + 3*sin(e + f*x)/(8*f*(a + b)**2*(-a*sin(e + f*x)**2 + a + b)) + 3*atanh(sqrt(a)*sin(e + f*x)/sqrt(a + b))/(8*sqrt(a)*f*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_207():
    f = sec(e + f*x)**3/(a + b*sec(e + f*x)**2)**3
    F = -b*sin(e + f*x)/(4*a*f*(a + b)*(-a*sin(e + f*x)**2 + a + b)**2) + (4*a + b)*sin(e + f*x)/(8*a*f*(a + b)**2*(-a*sin(e + f*x)**2 + a + b)) + (4*a + b)*atanh(sqrt(a)*sin(e + f*x)/sqrt(a + b))/(8*a**(sympy.S(3)/2)*f*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_208():
    f = sec(e + f*x)/(a + b*sec(e + f*x)**2)**3
    F = -b*sin(e + f*x)*cos(e + f*x)**2/(4*a*f*(a + b)*(-a*sin(e + f*x)**2 + a + b)**2) - 3*b*(2*a + b)*sin(e + f*x)/(8*a**2*f*(a + b)**2*(-a*sin(e + f*x)**2 + a + b)) + (8*a**2 + 8*a*b + 3*b**2)*atanh(sqrt(a)*sin(e + f*x)/sqrt(a + b))/(8*a**(sympy.S(5)/2)*f*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_209():
    f = cos(e + f*x)/(a + b*sec(e + f*x)**2)**3
    F = -b**3*sin(e + f*x)/(4*a**3*f*(a + b)*(-a*sin(e + f*x)**2 + a + b)**2) + 3*b**2*(4*a + 3*b)*sin(e + f*x)/(8*a**3*f*(a + b)**2*(-a*sin(e + f*x)**2 + a + b)) + sin(e + f*x)/(a**3*f) - 3*b*(4*(a + b)**2 + (2*a + b)**2)*atanh(sqrt(a)*sin(e + f*x)/sqrt(a + b))/(8*a**(sympy.S(7)/2)*f*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_210():
    f = cos(e + f*x)**3/(a + b*sec(e + f*x)**2)**3
    F = -sin(e + f*x)**3/(3*a**3*f) + b**4*sin(e + f*x)/(4*a**4*f*(a + b)*(-a*sin(e + f*x)**2 + a + b)**2) - b**3*(16*a + 13*b)*sin(e + f*x)/(8*a**4*f*(a + b)**2*(-a*sin(e + f*x)**2 + a + b)) + (a - 3*b)*sin(e + f*x)/(a**4*f) + b**2*(48*a**2 + 80*a*b + 35*b**2)*atanh(sqrt(a)*sin(e + f*x)/sqrt(a + b))/(8*a**(sympy.S(9)/2)*f*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_211():
    f = cos(e + f*x)**5/(a + b*sec(e + f*x)**2)**3
    F = sin(e + f*x)**5/(5*a**3*f) - (2*a - 3*b)*sin(e + f*x)**3/(3*a**4*f) - b**5*sin(e + f*x)/(4*a**5*f*(a + b)*(-a*sin(e + f*x)**2 + a + b)**2) + b**4*(20*a + 17*b)*sin(e + f*x)/(8*a**5*f*(a + b)**2*(-a*sin(e + f*x)**2 + a + b)) + (a**2 - 3*a*b + 6*b**2)*sin(e + f*x)/(a**5*f) - b**3*(80*a**2 + 140*a*b + 63*b**2)*atanh(sqrt(a)*sin(e + f*x)/sqrt(a + b))/(8*a**(sympy.S(11)/2)*f*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_212():
    f = sec(e + f*x)**6/(a + b*sec(e + f*x)**2)**3
    F = -a*tan(e + f*x)*sec(e + f*x)**2/(4*b*f*(a + b)*(a + b*tan(e + f*x)**2 + b)**2) - 3*a*(a + 2*b)*tan(e + f*x)/(8*b**2*f*(a + b)**2*(a + b*tan(e + f*x)**2 + b)) + (3*a**2 + 8*a*b + 8*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(8*b**(sympy.S(5)/2)*f*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_213():
    f = sec(e + f*x)**4/(a + b*sec(e + f*x)**2)**3
    F = -a*tan(e + f*x)/(4*b*f*(a + b)*(a + b*tan(e + f*x)**2 + b)**2) + (a + 4*b)*tan(e + f*x)/(8*b*f*(a + b)**2*(a + b*tan(e + f*x)**2 + b)) + (a + 4*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(8*b**(sympy.S(3)/2)*f*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_214():
    f = sec(e + f*x)**2/(a + b*sec(e + f*x)**2)**3
    F = tan(e + f*x)/(f*(4*a + 4*b)*(a + b*tan(e + f*x)**2 + b)**2) + 3*tan(e + f*x)/(8*f*(a + b)**2*(a + b*tan(e + f*x)**2 + b)) + 3*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(8*sqrt(b)*f*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_215():
    f = (a + b*sec(e + f*x)**2)**(-3)
    F = -b*tan(e + f*x)/(4*a*f*(a + b)*(a + b*tan(e + f*x)**2 + b)**2) - b*(7*a + 4*b)*tan(e + f*x)/(8*a**2*f*(a + b)**2*(a + b*tan(e + f*x)**2 + b)) - sqrt(b)*(15*a**2 + 20*a*b + 8*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(8*a**3*f*(a + b)**(sympy.S(5)/2)) + x/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_216():
    f = cos(e + f*x)**2/(a + b*sec(e + f*x)**2)**3
    F = sin(e + f*x)*cos(e + f*x)/(2*a*f*(a + b*tan(e + f*x)**2 + b)**2) + b*(2*a + 3*b)*tan(e + f*x)/(4*a**2*f*(a + b)*(a + b*tan(e + f*x)**2 + b)**2) + b*(a + 4*b)*(4*a + 3*b)*tan(e + f*x)/(8*a**3*f*(a + b)**2*(a + b*tan(e + f*x)**2 + b)) + b**(sympy.S(3)/2)*(35*a**2 + 56*a*b + 24*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(8*a**4*f*(a + b)**(sympy.S(5)/2)) + x*(a - 6*b)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_217():
    f = cos(e + f*x)**4/(a + b*sec(e + f*x)**2)**3
    F = sin(e + f*x)*cos(e + f*x)**3/(4*a*f*(a + b*tan(e + f*x)**2 + b)**2) + (3*a - 8*b)*sin(e + f*x)*cos(e + f*x)/(8*a**2*f*(a + b*tan(e + f*x)**2 + b)**2) + b*(3*a**2 - 7*a*b - 12*b**2)*tan(e + f*x)/(8*a**3*f*(a + b)*(a + b*tan(e + f*x)**2 + b)**2) + 3*b*(a + 2*b)*(a**2 - 4*a*b - 4*b**2)*tan(e + f*x)/(8*a**4*f*(a + b)**2*(a + b*tan(e + f*x)**2 + b)) - 3*b**(sympy.S(5)/2)*(21*a**2 + 36*a*b + 16*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(8*a**5*f*(a + b)**(sympy.S(5)/2)) + x*(3*a**2 - 12*a*b + 48*b**2)/(8*a**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_218():
    f = cos(e + f*x)**6/(a + b*sec(e + f*x)**2)**3
    F = sin(e + f*x)*cos(e + f*x)**5/(6*a*f*(a + b*tan(e + f*x)**2 + b)**2) + (5*a - 10*b)*sin(e + f*x)*cos(e + f*x)**3/(24*a**2*f*(a + b*tan(e + f*x)**2 + b)**2) + (15*a**2 - 34*a*b + 80*b**2)*sin(e + f*x)*cos(e + f*x)/(48*a**3*f*(a + b*tan(e + f*x)**2 + b)**2) + b*(15*a**3 - 29*a**2*b + 64*a*b**2 + 120*b**3)*tan(e + f*x)/(48*a**4*f*(a + b)*(a + b*tan(e + f*x)**2 + b)**2) + b*(5*a**4 - 8*a**3*b + 17*a**2*b**2 + 116*a*b**3 + 80*b**4)*tan(e + f*x)/(16*a**5*f*(a + b)**2*(a + b*tan(e + f*x)**2 + b)) + b**(sympy.S(7)/2)*(99*a**2 + 176*a*b + 80*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(8*a**6*f*(a + b)**(sympy.S(5)/2)) + x*(5*a**3 - 18*a**2*b + 48*a*b**2 - 160*b**3)/(16*a**6)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_219():
    f = (a + b*sec(c + d*x)**2)**(-4)
    F = -b*tan(c + d*x)/(6*a*d*(a + b)*(a + b*tan(c + d*x)**2 + b)**3) - b*(11*a + 6*b)*tan(c + d*x)/(24*a**2*d*(a + b)**2*(a + b*tan(c + d*x)**2 + b)**2) - b*(19*a**2 + 22*a*b + 8*b**2)*tan(c + d*x)/(16*a**3*d*(a + b)**3*(a + b*tan(c + d*x)**2 + b)) - sqrt(b)*(35*a**3 + 70*a**2*b + 56*a*b**2 + 16*b**3)*atan(sqrt(b)*tan(c + d*x)/sqrt(a + b))/(16*a**4*d*(a + b)**(sympy.S(7)/2)) + x/a**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_220():
    f = (-a*sec(c + d*x)**2 + a)**(sympy.S(7)/2)
    F = -a**3*sqrt(-a*tan(c + d*x)**2)*log(cos(c + d*x))*cot(c + d*x)/d - a**3*sqrt(-a*tan(c + d*x)**2)*tan(c + d*x)**5/(6*d) + a**3*sqrt(-a*tan(c + d*x)**2)*tan(c + d*x)**3/(4*d) - a**3*sqrt(-a*tan(c + d*x)**2)*tan(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_221():
    f = (-a*sec(c + d*x)**2 + a)**(sympy.S(5)/2)
    F = -a**2*sqrt(-a*tan(c + d*x)**2)*log(cos(c + d*x))*cot(c + d*x)/d + a**2*sqrt(-a*tan(c + d*x)**2)*tan(c + d*x)**3/(4*d) - a**2*sqrt(-a*tan(c + d*x)**2)*tan(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_222():
    f = (-a*sec(c + d*x)**2 + a)**(sympy.S(3)/2)
    F = -a*sqrt(-a*tan(c + d*x)**2)*log(cos(c + d*x))*cot(c + d*x)/d - a*sqrt(-a*tan(c + d*x)**2)*tan(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_223():
    f = sqrt(-a*sec(c + d*x)**2 + a)
    F = -sqrt(-a*tan(c + d*x)**2)*log(cos(c + d*x))*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_224():
    f = 1/sqrt(-a*sec(c + d*x)**2 + a)
    F = log(sin(c + d*x))*tan(c + d*x)/(d*sqrt(-a*tan(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_225():
    f = (-a*sec(c + d*x)**2 + a)**(sympy.S(-3)/2)
    F = log(sin(c + d*x))*tan(c + d*x)/(a*d*sqrt(-a*tan(c + d*x)**2)) + cot(c + d*x)/(2*a*d*sqrt(-a*tan(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_226():
    f = (-a*sec(c + d*x)**2 + a)**(sympy.S(-5)/2)
    F = log(sin(c + d*x))*tan(c + d*x)/(a**2*d*sqrt(-a*tan(c + d*x)**2)) - cot(c + d*x)**3/(4*a**2*d*sqrt(-a*tan(c + d*x)**2)) + cot(c + d*x)/(2*a**2*d*sqrt(-a*tan(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_227():
    f = (-a*sec(c + d*x)**2 + a)**(sympy.S(-7)/2)
    F = log(sin(c + d*x))*tan(c + d*x)/(a**3*d*sqrt(-a*tan(c + d*x)**2)) + cot(c + d*x)**5/(6*a**3*d*sqrt(-a*tan(c + d*x)**2)) - cot(c + d*x)**3/(4*a**3*d*sqrt(-a*tan(c + d*x)**2)) + cot(c + d*x)/(2*a**3*d*sqrt(-a*tan(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_228():
    f = sqrt(a + b*sec(e + f*x)**2)*sec(e + f*x)**5
    F = sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*tan(e + f*x)*sec(e + f*x)**3/(5*f*sqrt(a*cos(e + f*x)**2 + b)) - (a - 8*b)*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(15*b*f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)) + (a + 4*b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*tan(e + f*x)*sec(e + f*x)/(15*b*f*sqrt(a*cos(e + f*x)**2 + b)) - sqrt(a + b*sec(e + f*x)**2)*(2*a**2 - 3*a*b - 8*b**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sin(e + f*x)/(15*b**2*f*sqrt(a*cos(e + f*x)**2 + b)) + sqrt(a + b*sec(e + f*x)**2)*(2*a**2 - 3*a*b - 8*b**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(15*b**2*f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_229():
    f = sqrt(a + b*sec(e + f*x)**2)*sec(e + f*x)**3
    F = sqrt(a + b*sec(e + f*x)**2)*(2*a + 2*b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(3*f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)) + sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*tan(e + f*x)*sec(e + f*x)/(3*f*sqrt(a*cos(e + f*x)**2 + b)) + (a + 2*b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sin(e + f*x)/(3*b*f*sqrt(a*cos(e + f*x)**2 + b)) - (a + 2*b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(3*b*f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_230():
    f = sqrt(a + b*sec(e + f*x)**2)*sec(e + f*x)
    F = (a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)) + sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sin(e + f*x)/(f*sqrt(a*cos(e + f*x)**2 + b)) - sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_231():
    f = sqrt(a + b*sec(e + f*x)**2)*cos(e + f*x)
    F = sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_232():
    f = sqrt(a + b*sec(e + f*x)**2)*cos(e + f*x)**3
    F = sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sin(e + f*x)*cos(e + f*x)**2/(3*f*sqrt(a*cos(e + f*x)**2 + b)) - b*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(3*a*f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)) + sqrt(a + b*sec(e + f*x)**2)*(2*a + b)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(3*a*f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_233():
    f = sqrt(a + b*sec(e + f*x)**2)*cos(e + f*x)**5
    F = sqrt(a + b*sec(e + f*x)**2)*(4*a - 2*b)*sqrt(-a*sin(e + f*x)**2 + a + b)*sin(e + f*x)*cos(e + f*x)**2/(15*a*f*sqrt(a*cos(e + f*x)**2 + b)) + sqrt(a + b*sec(e + f*x)**2)*(-a*sin(e + f*x)**2 + a + b)**(sympy.S(3)/2)*sin(e + f*x)*cos(e + f*x)**2/(5*a*f*sqrt(a*cos(e + f*x)**2 + b)) - b*(a + b)*sqrt(a + b*sec(e + f*x)**2)*(4*a - 2*b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(15*a**2*f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)) + sqrt(a + b*sec(e + f*x)**2)*(8*a**2 + 3*a*b - 2*b**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(15*a**2*f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_234():
    f = sqrt(a + b*sec(e + f*x)**2)*sec(e + f*x)**6
    F = (a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)*tan(e + f*x)*sec(e + f*x)**2/(6*b*f) - (3*a - 5*b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)*tan(e + f*x)/(24*b**2*f) + sqrt(a + b*tan(e + f*x)**2 + b)*(a**2 - 2*a*b + 5*b**2)*tan(e + f*x)/(16*b**2*f) + (a + b)*(a**2 - 2*a*b + 5*b**2)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(16*b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_235():
    f = sqrt(a + b*sec(e + f*x)**2)*sec(e + f*x)**4
    F = -(a - 3*b)*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/(8*b*f) + (a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)*tan(e + f*x)/(4*b*f) - (a - 3*b)*(a + b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(8*b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_236():
    f = sqrt(a + b*sec(e + f*x)**2)*sec(e + f*x)**2
    F = sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/(2*f) + (a + b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_237():
    f = sqrt(a + b*sec(e + f*x)**2)
    F = sqrt(a)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f + sqrt(b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_238():
    f = sqrt(a + b*sec(e + f*x)**2)*cos(e + f*x)**2
    F = sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)/(2*f) + (a + b)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_239():
    f = sqrt(a + b*sec(e + f*x)**2)*cos(e + f*x)**4
    F = (3*a - b)*sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)/(8*a*f) + (a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)*sin(e + f*x)*cos(e + f*x)**3/(4*a*f) + (a + b)*(3*a - b)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(8*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_240():
    f = sqrt(a + b*sec(e + f*x)**2)*cos(e + f*x)**6
    F = sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)**5/(6*f) + (5*a + b)*sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)**3/(24*a*f) + (3*a - b)*(5*a + 3*b)*sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)/(48*a**2*f) + (a + b)*(5*a**2 - 2*a*b + b**2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(16*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_241():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*sec(e + f*x)**5
    F = b*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*tan(e + f*x)*sec(e + f*x)**5/(7*f*sqrt(a*cos(e + f*x)**2 + b)) + sqrt(a + b*sec(e + f*x)**2)*(8*a + 6*b)*sqrt(-a*sin(e + f*x)**2 + a + b)*tan(e + f*x)*sec(e + f*x)**3/(35*f*sqrt(a*cos(e + f*x)**2 + b)) - (a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*(a**2 - 16*a*b - 16*b**2)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(35*b*f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)) + sqrt(a + b*sec(e + f*x)**2)*(a**2 + 11*a*b + 8*b**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*tan(e + f*x)*sec(e + f*x)/(35*b*f*sqrt(a*cos(e + f*x)**2 + b)) + (-2*a - 4*b)*sqrt(a + b*sec(e + f*x)**2)*(a**2 - 4*a*b - 4*b**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sin(e + f*x)/(35*b**2*f*sqrt(a*cos(e + f*x)**2 + b)) + sqrt(a + b*sec(e + f*x)**2)*(2*a + 4*b)*(a**2 - 4*a*b - 4*b**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(35*b**2*f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_242():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*sec(e + f*x)**3
    F = b*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*tan(e + f*x)*sec(e + f*x)**3/(5*f*sqrt(a*cos(e + f*x)**2 + b)) + (a + b)*sqrt(a + b*sec(e + f*x)**2)*(9*a + 8*b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(15*f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)) + sqrt(a + b*sec(e + f*x)**2)*(6*a + 4*b)*sqrt(-a*sin(e + f*x)**2 + a + b)*tan(e + f*x)*sec(e + f*x)/(15*f*sqrt(a*cos(e + f*x)**2 + b)) + sqrt(a + b*sec(e + f*x)**2)*(3*a**2 + 13*a*b + 8*b**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sin(e + f*x)/(15*b*f*sqrt(a*cos(e + f*x)**2 + b)) - sqrt(a + b*sec(e + f*x)**2)*(3*a**2 + 13*a*b + 8*b**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(15*b*f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_243():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*sec(e + f*x)
    F = b*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*tan(e + f*x)*sec(e + f*x)/(3*f*sqrt(a*cos(e + f*x)**2 + b)) + (a + b)*sqrt(a + b*sec(e + f*x)**2)*(3*a + 2*b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(3*f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)) + sqrt(a + b*sec(e + f*x)**2)*(4*a + 2*b)*sqrt(-a*sin(e + f*x)**2 + a + b)*sin(e + f*x)/(3*f*sqrt(a*cos(e + f*x)**2 + b)) - sqrt(a + b*sec(e + f*x)**2)*(4*a + 2*b)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(3*f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_244():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*cos(e + f*x)
    F = b*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)) + b*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sin(e + f*x)/(f*sqrt(a*cos(e + f*x)**2 + b)) + (a - b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_245():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*cos(e + f*x)**3
    F = a*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sin(e + f*x)*cos(e + f*x)**2/(3*f*sqrt(a*cos(e + f*x)**2 + b)) - b*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(3*f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)) + sqrt(a + b*sec(e + f*x)**2)*(2*a + 4*b)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(3*f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_246():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*cos(e + f*x)**5
    F = a*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sin(e + f*x)*cos(e + f*x)**4/(5*f*sqrt(a*cos(e + f*x)**2 + b)) + sqrt(a + b*sec(e + f*x)**2)*(4*a + 6*b)*sqrt(-a*sin(e + f*x)**2 + a + b)*sin(e + f*x)*cos(e + f*x)**2/(15*f*sqrt(a*cos(e + f*x)**2 + b)) - b*(a + b)*sqrt(a + b*sec(e + f*x)**2)*(4*a + 3*b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(15*a*f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)) + sqrt(a + b*sec(e + f*x)**2)*(8*a**2 + 13*a*b + 3*b**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(15*a*f*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_247():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*sec(e + f*x)**6
    F = (a + b*tan(e + f*x)**2 + b)**(sympy.S(5)/2)*tan(e + f*x)*sec(e + f*x)**2/(8*b*f) + (a + b)*sqrt(a + b*tan(e + f*x)**2 + b)*(3*a**2 - 10*a*b + 35*b**2)*tan(e + f*x)/(128*b**2*f) - (3*a - 7*b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(5)/2)*tan(e + f*x)/(48*b**2*f) + (a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)*(3*a**2 - 10*a*b + 35*b**2)*tan(e + f*x)/(192*b**2*f) + (a + b)**2*(3*a**2 - 10*a*b + 35*b**2)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(128*b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_248():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*sec(e + f*x)**4
    F = -(a - 5*b)*(a + b)*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/(16*b*f) - (a - 5*b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)*tan(e + f*x)/(24*b*f) + (a + b*tan(e + f*x)**2 + b)**(sympy.S(5)/2)*tan(e + f*x)/(6*b*f) - (a - 5*b)*(a + b)**2*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(16*b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_249():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*sec(e + f*x)**2
    F = (3*a + 3*b)*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/(8*f) + (a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)*tan(e + f*x)/(4*f) + 3*(a + b)**2*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(8*sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_250():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = a**(sympy.S(3)/2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f + sqrt(b)*(3*a + b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*f) + b*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_251():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*cos(e + f*x)**2
    F = sqrt(a)*(a + 3*b)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*f) + a*sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)/(2*f) + b**(sympy.S(3)/2)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_252():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*cos(e + f*x)**4
    F = (3*a + 3*b)*sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)/(8*f) + (a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)*sin(e + f*x)*cos(e + f*x)**3/(4*f) + 3*(a + b)**2*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(8*sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_253():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*cos(e + f*x)**6
    F = (a + b)*(5*a - b)*sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)/(16*a*f) + (5*a - b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)*sin(e + f*x)*cos(e + f*x)**3/(24*a*f) + (a + b*tan(e + f*x)**2 + b)**(sympy.S(5)/2)*sin(e + f*x)*cos(e + f*x)**5/(6*a*f) + (a + b)**2*(5*a - b)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(16*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_254():
    f = (a + b*sec(c + d*x)**2)**(sympy.S(5)/2)
    F = a**(sympy.S(5)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a + b*tan(c + d*x)**2 + b))/d + sqrt(b)*(15*a**2 + 10*a*b + 3*b**2)*atanh(sqrt(b)*tan(c + d*x)/sqrt(a + b*tan(c + d*x)**2 + b))/(8*d) + b*(7*a + 3*b)*sqrt(a + b*tan(c + d*x)**2 + b)*tan(c + d*x)/(8*d) + b*(a + b*tan(c + d*x)**2 + b)**(sympy.S(3)/2)*tan(c + d*x)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_255():
    f = (sec(x)**2 + 1)**(sympy.S(3)/2)
    F = sqrt(tan(x)**2 + 2)*tan(x)/2 + 2*asinh(sqrt(2)*tan(x)/2) + atan(tan(x)/sqrt(tan(x)**2 + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_256():
    f = sqrt(sec(x)**2 + 1)
    F = asinh(sqrt(2)*tan(x)/2) + atan(tan(x)/sqrt(tan(x)**2 + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_257():
    f = sec(e + f*x)**5/sqrt(a + b*sec(e + f*x)**2)
    F = -(a - 2*b)*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(3*b*f*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)) + sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)*tan(e + f*x)*sec(e + f*x)**3/(3*b*f*sqrt(a + b*sec(e + f*x)**2)) - (2*a - 2*b)*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)*tan(e + f*x)*sec(e + f*x)/(3*b**2*f*sqrt(a + b*sec(e + f*x)**2)) + (2*a - 2*b)*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(3*b**2*f*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_258():
    f = sec(e + f*x)**3/sqrt(a + b*sec(e + f*x)**2)
    F = -sqrt(a)*sqrt(a + b)*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*elliptic_e(asin(sqrt(a)*sin(e + f*x)/sqrt(a + b)), (a + b)/a)/(b*f*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)) + sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)*tan(e + f*x)*sec(e + f*x)/(b*f*sqrt(a + b*sec(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_259():
    f = sec(e + f*x)/sqrt(a + b*sec(e + f*x)**2)
    F = sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(f*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_260():
    f = cos(e + f*x)/sqrt(a + b*sec(e + f*x)**2)
    F = sqrt(a + b)*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*elliptic_e(asin(sqrt(a)*sin(e + f*x)/sqrt(a + b)), (a + b)/a)/(sqrt(a)*f*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_261():
    f = cos(e + f*x)**3/sqrt(a + b*sec(e + f*x)**2)
    F = sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)*sin(e + f*x)/(3*a*f*sqrt(a + b*sec(e + f*x)**2)) - b*(a - 2*b)*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(3*a**2*f*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)) + (2*a - 2*b)*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(3*a**2*f*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_262():
    f = cos(e + f*x)**5/sqrt(a + b*sec(e + f*x)**2)
    F = sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)*sin(e + f*x)*cos(e + f*x)**2/(5*a*f*sqrt(a + b*sec(e + f*x)**2)) + (4*a - 4*b)*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)*sin(e + f*x)/(15*a**2*f*sqrt(a + b*sec(e + f*x)**2)) - b*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*(4*a**2 - 3*a*b + 8*b**2)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(15*a**3*f*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)) + sqrt(a*cos(e + f*x)**2 + b)*(8*a**2 - 7*a*b + 8*b**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(15*a**3*f*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_263():
    f = sec(e + f*x)**6/sqrt(a + b*sec(e + f*x)**2)
    F = sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)*sec(e + f*x)**2/(4*b*f) - (3*a - 3*b)*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/(8*b**2*f) + (3*a**2 - 2*a*b + 3*b**2)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(8*b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_264():
    f = sec(e + f*x)**4/sqrt(a + b*sec(e + f*x)**2)
    F = sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/(2*b*f) - (a - b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_265():
    f = sec(e + f*x)**2/sqrt(a + b*sec(e + f*x)**2)
    F = atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_266():
    f = 1/sqrt(a + b*sec(e + f*x)**2)
    F = atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_267():
    f = cos(e + f*x)**2/sqrt(a + b*sec(e + f*x)**2)
    F = sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)/(2*a*f) + (a - b)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_268():
    f = cos(e + f*x)**4/sqrt(a + b*sec(e + f*x)**2)
    F = sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)**3/(4*a*f) + (3*a - 3*b)*sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)/(8*a**2*f) + (3*a**2 - 2*a*b + 3*b**2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(8*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_269():
    f = cos(e + f*x)**6/sqrt(a + b*sec(e + f*x)**2)
    F = sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)**5/(6*a*f) + (5*a - 5*b)*sqrt(a + b*tan(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)**3/(24*a**2*f) + sqrt(a + b*tan(e + f*x)**2 + b)*(15*a**2 - 14*a*b + 15*b**2)*sin(e + f*x)*cos(e + f*x)/(48*a**3*f) + (a - b)*(5*a**2 + 2*a*b + 5*b**2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(16*a**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_270():
    f = sec(e + f*x)**5/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = a*(2*a + b)*sqrt(a*cos(e + f*x)**2 + b)*sin(e + f*x)/(b**2*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)) + sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(b*f*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)) + sqrt(a*cos(e + f*x)**2 + b)*tan(e + f*x)*sec(e + f*x)/(b*f*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)) - (2*a + b)*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(b**2*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_271():
    f = sec(e + f*x)**3/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -a*sqrt(a*cos(e + f*x)**2 + b)*sin(e + f*x)/(b*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)) + sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(b*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_272():
    f = sec(e + f*x)/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = sqrt(a*cos(e + f*x)**2 + b)*sin(e + f*x)/(f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)) + sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(a*f*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)) - sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(a*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_273():
    f = cos(e + f*x)/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -b*sqrt(a*cos(e + f*x)**2 + b)*sin(e + f*x)/(a*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)) - 2*b*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(a**2*f*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)) + (a + 2*b)*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(a**2*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_274():
    f = cos(e + f*x)**3/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -b*sqrt(a*cos(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)**2/(a*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)) + (a + 4*b)*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)*sin(e + f*x)/(3*a**2*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)) - b*(a - 8*b)*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(3*a**3*f*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)) + sqrt(a*cos(e + f*x)**2 + b)*(2*a**2 - 3*a*b - 8*b**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(3*a**3*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_275():
    f = cos(e + f*x)**5/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -b*sqrt(a*cos(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)**4/(a*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)) + (a + 6*b)*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)*sin(e + f*x)*cos(e + f*x)**2/(5*a**2*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)) + sqrt(a*cos(e + f*x)**2 + b)*(4*a**2 - 5*a*b - 24*b**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sin(e + f*x)/(15*a**3*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)) - 4*b*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*(a**2 - 2*a*b + 12*b**2)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(15*a**4*f*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)) + sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)*(8*a**3 - 9*a**2*b + 16*a*b**2 + 48*b**3)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(15*a**4*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_276():
    f = sec(e + f*x)**6/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -a*tan(e + f*x)*sec(e + f*x)**2/(b*f*(a + b)*sqrt(a + b*tan(e + f*x)**2 + b)) + (3*a + b)*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/(2*b**2*f*(a + b)) - (3*a - b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_277():
    f = sec(e + f*x)**4/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -a*tan(e + f*x)/(b*f*(a + b)*sqrt(a + b*tan(e + f*x)**2 + b)) + atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_278():
    f = sec(e + f*x)**2/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = tan(e + f*x)/(f*(a + b)*sqrt(a + b*tan(e + f*x)**2 + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_279():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(-3)/2)
    F = -b*tan(e + f*x)/(a*f*(a + b)*sqrt(a + b*tan(e + f*x)**2 + b)) + atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_280():
    f = cos(e + f*x)**2/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = sin(e + f*x)*cos(e + f*x)/(2*a*f*sqrt(a + b*tan(e + f*x)**2 + b)) + b*(a + 3*b)*tan(e + f*x)/(2*a**2*f*(a + b)*sqrt(a + b*tan(e + f*x)**2 + b)) + (a - 3*b)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_281():
    f = cos(e + f*x)**4/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = sin(e + f*x)*cos(e + f*x)**3/(4*a*f*sqrt(a + b*tan(e + f*x)**2 + b)) + (3*a - 5*b)*sin(e + f*x)*cos(e + f*x)/(8*a**2*f*sqrt(a + b*tan(e + f*x)**2 + b)) + b*(a - 3*b)*(3*a + 5*b)*tan(e + f*x)/(8*a**3*f*(a + b)*sqrt(a + b*tan(e + f*x)**2 + b)) + (3*a**2 - 6*a*b + 15*b**2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(8*a**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_282():
    f = cos(e + f*x)**6/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = sin(e + f*x)*cos(e + f*x)**5/(6*a*f*sqrt(a + b*tan(e + f*x)**2 + b)) + (5*a - 7*b)*sin(e + f*x)*cos(e + f*x)**3/(24*a**2*f*sqrt(a + b*tan(e + f*x)**2 + b)) + (15*a**2 - 22*a*b + 35*b**2)*sin(e + f*x)*cos(e + f*x)/(48*a**3*f*sqrt(a + b*tan(e + f*x)**2 + b)) + b*(15*a**3 - 17*a**2*b + 25*a*b**2 + 105*b**3)*tan(e + f*x)/(48*a**4*f*(a + b)*sqrt(a + b*tan(e + f*x)**2 + b)) + (5*a**3 - 9*a**2*b + 15*a*b**2 - 35*b**3)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(16*a**(sympy.S(9)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_283():
    f = sec(e + f*x)**5/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -a*sqrt(a*cos(e + f*x)**2 + b)*sin(e + f*x)/(3*b*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*(-a*sin(e + f*x)**2 + a + b)**(sympy.S(3)/2)) - 2*a*(a + 2*b)*sqrt(a*cos(e + f*x)**2 + b)*sin(e + f*x)/(3*b**2*f*(a + b)**2*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)) - sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(3*b*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)) + (2*a + 4*b)*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(3*b**2*f*(a + b)**2*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_284():
    f = sec(e + f*x)**3/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = sqrt(a*cos(e + f*x)**2 + b)*sin(e + f*x)/(f*sqrt(a + b*sec(e + f*x)**2)*(3*a + 3*b)*(-a*sin(e + f*x)**2 + a + b)**(sympy.S(3)/2)) - (a - b)*sqrt(a*cos(e + f*x)**2 + b)*sin(e + f*x)/(3*b*f*(a + b)**2*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)) + sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(3*a*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)) + (a - b)*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(3*a*b*f*(a + b)**2*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_285():
    f = sec(e + f*x)/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -b*sqrt(a*cos(e + f*x)**2 + b)*sin(e + f*x)/(3*a*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*(-a*sin(e + f*x)**2 + a + b)**(sympy.S(3)/2)) + (4*a + 2*b)*sqrt(a*cos(e + f*x)**2 + b)*sin(e + f*x)/(3*a*f*(a + b)**2*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)) + (3*a + 2*b)*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(3*a**2*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)) - (4*a + 2*b)*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(3*a**2*f*(a + b)**2*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_286():
    f = cos(e + f*x)/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -b*sqrt(a*cos(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)**2/(3*a*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*(-a*sin(e + f*x)**2 + a + b)**(sympy.S(3)/2)) - 2*b*(3*a + 2*b)*sqrt(a*cos(e + f*x)**2 + b)*sin(e + f*x)/(3*a**2*f*(a + b)**2*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)) - b*(9*a + 8*b)*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(3*a**3*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)) + sqrt(a*cos(e + f*x)**2 + b)*(3*a**2 + 13*a*b + 8*b**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(3*a**3*f*(a + b)**2*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_287():
    f = cos(e + f*x)**3/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -b*sqrt(a*cos(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)**4/(3*a*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*(-a*sin(e + f*x)**2 + a + b)**(sympy.S(3)/2)) - 2*b*(4*a + 3*b)*sqrt(a*cos(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)**2/(3*a**2*f*(a + b)**2*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)) + sqrt(a*cos(e + f*x)**2 + b)*(a**2 + 11*a*b + 8*b**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sin(e + f*x)/(3*a**3*f*(a + b)**2*sqrt(a + b*sec(e + f*x)**2)) - b*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*(a**2 - 16*a*b - 16*b**2)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(3*a**4*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)) + (2*a + 4*b)*sqrt(a*cos(e + f*x)**2 + b)*(a**2 - 4*a*b - 4*b**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(3*a**4*f*(a + b)**2*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_288():
    f = cos(e + f*x)**5/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -b*sqrt(a*cos(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)**6/(3*a*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*(-a*sin(e + f*x)**2 + a + b)**(sympy.S(3)/2)) - 2*b*(5*a + 4*b)*sqrt(a*cos(e + f*x)**2 + b)*sin(e + f*x)*cos(e + f*x)**4/(3*a**2*f*(a + b)**2*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)) + sqrt(a*cos(e + f*x)**2 + b)*(3*a**2 + 61*a*b + 48*b**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sin(e + f*x)*cos(e + f*x)**2/(15*a**3*f*(a + b)**2*sqrt(a + b*sec(e + f*x)**2)) + sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)*(4*a**3 - 6*a**2*b - 84*a*b**2 - 64*b**3)*sin(e + f*x)/(15*a**4*f*(a + b)**2*sqrt(a + b*sec(e + f*x)**2)) - b*sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*(4*a**3 - 9*a**2*b + 120*a*b**2 + 128*b**3)*elliptic_f(asin(sin(e + f*x)), a/(a + b))/(15*a**5*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2 + a + b)*sqrt(cos(e + f*x)**2)) + sqrt(a*cos(e + f*x)**2 + b)*sqrt(-a*sin(e + f*x)**2 + a + b)*(8*a**4 - 11*a**3*b + 27*a**2*b**2 + 184*a*b**3 + 128*b**4)*elliptic_e(asin(sin(e + f*x)), a/(a + b))/(15*a**5*f*(a + b)**2*sqrt(a + b*sec(e + f*x)**2)*sqrt(-a*sin(e + f*x)**2/(a + b) + 1)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_289():
    f = sec(e + f*x)**6/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -a*tan(e + f*x)*sec(e + f*x)**2/(3*b*f*(a + b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - a*(3*a + 5*b)*tan(e + f*x)/(3*b**2*f*(a + b)**2*sqrt(a + b*tan(e + f*x)**2 + b)) + atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_290():
    f = sec(e + f*x)**4/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = tan(e + f*x)*sec(e + f*x)**2/(f*(3*a + 3*b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) + 2*tan(e + f*x)/(3*f*(a + b)**2*sqrt(a + b*tan(e + f*x)**2 + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_291():
    f = sec(e + f*x)**2/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = tan(e + f*x)/(f*(3*a + 3*b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) + 2*tan(e + f*x)/(3*f*(a + b)**2*sqrt(a + b*tan(e + f*x)**2 + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_292():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(-5)/2)
    F = -b*tan(e + f*x)/(3*a*f*(a + b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - b*(5*a + 3*b)*tan(e + f*x)/(3*a**2*f*(a + b)**2*sqrt(a + b*tan(e + f*x)**2 + b)) + atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_293():
    f = cos(e + f*x)**2/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = sin(e + f*x)*cos(e + f*x)/(2*a*f*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) + b*(3*a + 5*b)*tan(e + f*x)/(6*a**2*f*(a + b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) + b*(3*a**2 + 22*a*b + 15*b**2)*tan(e + f*x)/(6*a**3*f*(a + b)**2*sqrt(a + b*tan(e + f*x)**2 + b)) + (a - 5*b)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*a**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_294():
    f = cos(e + f*x)**4/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = sin(e + f*x)*cos(e + f*x)**3/(4*a*f*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) + (3*a - 7*b)*sin(e + f*x)*cos(e + f*x)/(8*a**2*f*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) + b*(9*a**2 - 18*a*b - 35*b**2)*tan(e + f*x)/(24*a**3*f*(a + b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) + b*(9*a**3 - 15*a**2*b - 145*a*b**2 - 105*b**3)*tan(e + f*x)/(24*a**4*f*(a + b)**2*sqrt(a + b*tan(e + f*x)**2 + b)) + (3*a**2 - 10*a*b + 35*b**2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(8*a**(sympy.S(9)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_295():
    f = cos(e + f*x)**6/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = sin(e + f*x)*cos(e + f*x)**5/(6*a*f*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) + (5*a - 9*b)*sin(e + f*x)*cos(e + f*x)**3/(24*a**2*f*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) + (5*a**2 - 10*a*b + 21*b**2)*sin(e + f*x)*cos(e + f*x)/(16*a**3*f*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) + b*(15*a**3 - 25*a**2*b + 49*a*b**2 + 105*b**3)*tan(e + f*x)/(48*a**4*f*(a + b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) + b*(15*a**4 - 20*a**3*b + 38*a**2*b**2 + 420*a*b**3 + 315*b**4)*tan(e + f*x)/(48*a**5*f*(a + b)**2*sqrt(a + b*tan(e + f*x)**2 + b)) + (5*a - 15*b)*(a**2 + 7*b**2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(16*a**(sympy.S(11)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_296():
    f = (a + b*sec(c + d*x)**2)**(sympy.S(-7)/2)
    F = -b*tan(c + d*x)/(5*a*d*(a + b)*(a + b*tan(c + d*x)**2 + b)**(sympy.S(5)/2)) - b*(9*a + 5*b)*tan(c + d*x)/(15*a**2*d*(a + b)**2*(a + b*tan(c + d*x)**2 + b)**(sympy.S(3)/2)) - b*(33*a**2 + 40*a*b + 15*b**2)*tan(c + d*x)/(15*a**3*d*(a + b)**3*sqrt(a + b*tan(c + d*x)**2 + b)) + atan(sqrt(a)*tan(c + d*x)/sqrt(a + b*tan(c + d*x)**2 + b))/(a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_297():
    f = 1/sqrt(sec(x)**2 + 1)
    F = atan(tan(x)/sqrt(tan(x)**2 + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_298():
    f = (d*sec(e + f*x))**m*(a + b*sec(e + f*x)**2)**p
    F = sqrt(-tan(e + f*x)**2)*(d*sec(e + f*x))**m*(a + b*sec(e + f*x)**2)**p*cos(e + f*x)*appellf1(m/2, sympy.S.Half, -p, m/2 + 1, sec(e + f*x)**2, -b*sec(e + f*x)**2/a)/(f*m*(1 + b*sec(e + f*x)**2/a)**p*sin(e + f*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_299():
    f = (a + b*sec(e + f*x)**2)**p*sec(e + f*x)**5
    F = (a + b*sec(e + f*x)**2)**p*(-a*sin(e + f*x)**2 + a + b)**p*(cos(e + f*x)**2)**p*sin(e + f*x)*appellf1(sympy.S.Half, -p, p + 3, sympy.S(3)/2, a*sin(e + f*x)**2/(a + b), sin(e + f*x)**2)/(f*(a*cos(e + f*x)**2 + b)**p*(-a*sin(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_300():
    f = (a + b*sec(e + f*x)**2)**p*sec(e + f*x)**3
    F = (a + b*sec(e + f*x)**2)**p*(-a*sin(e + f*x)**2 + a + b)**p*(cos(e + f*x)**2)**p*sin(e + f*x)*appellf1(sympy.S.Half, -p, p + 2, sympy.S(3)/2, a*sin(e + f*x)**2/(a + b), sin(e + f*x)**2)/(f*(a*cos(e + f*x)**2 + b)**p*(-a*sin(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_301():
    f = (a + b*sec(e + f*x)**2)**p*sec(e + f*x)
    F = (a + b*sec(e + f*x)**2)**p*(-a*sin(e + f*x)**2 + a + b)**p*(cos(e + f*x)**2)**p*sin(e + f*x)*appellf1(sympy.S.Half, -p, p + 1, sympy.S(3)/2, a*sin(e + f*x)**2/(a + b), sin(e + f*x)**2)/(f*(a*cos(e + f*x)**2 + b)**p*(-a*sin(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_302():
    f = (a + b*sec(e + f*x)**2)**p*cos(e + f*x)
    F = (a + b*sec(e + f*x)**2)**p*(-a*sin(e + f*x)**2 + a + b)**p*(cos(e + f*x)**2)**p*sin(e + f*x)*appellf1(sympy.S.Half, -p, p, sympy.S(3)/2, a*sin(e + f*x)**2/(a + b), sin(e + f*x)**2)/(f*(a*cos(e + f*x)**2 + b)**p*(-a*sin(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_303():
    f = (a + b*sec(e + f*x)**2)**p*cos(e + f*x)**3
    F = (a + b*sec(e + f*x)**2)**p*(-a*sin(e + f*x)**2 + a + b)**p*(cos(e + f*x)**2)**p*sin(e + f*x)*appellf1(sympy.S.Half, -p, p - 1, sympy.S(3)/2, a*sin(e + f*x)**2/(a + b), sin(e + f*x)**2)/(f*(a*cos(e + f*x)**2 + b)**p*(-a*sin(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_304():
    f = (a + b*sec(e + f*x)**2)**p*cos(e + f*x)**5
    F = (a + b*sec(e + f*x)**2)**p*(-a*sin(e + f*x)**2 + a + b)**p*(cos(e + f*x)**2)**p*sin(e + f*x)*appellf1(sympy.S.Half, -p, p - 2, sympy.S(3)/2, a*sin(e + f*x)**2/(a + b), sin(e + f*x)**2)/(f*(a*cos(e + f*x)**2 + b)**p*(-a*sin(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_305():
    f = (a + b*sec(e + f*x)**2)**p*sec(e + f*x)**6
    F = (a + b*tan(e + f*x)**2 + b)**(p + 1)*tan(e + f*x)*sec(e + f*x)**2/(b*f*(2*p + 5)) - (3*a - 2*b*(p + 2))*(a + b*tan(e + f*x)**2 + b)**(p + 1)*tan(e + f*x)/(b**2*f*(2*p + 3)*(2*p + 5)) + (a + b*tan(e + f*x)**2 + b)**p*(3*a**2 - 4*a*b*(p + 1) + 4*b**2*(p**2 + 3*p + 2))*tan(e + f*x)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*tan(e + f*x)**2/(a + b))/(b**2*f*(2*p + 3)*(2*p + 5)*(b*tan(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_306():
    f = (a + b*sec(e + f*x)**2)**p*sec(e + f*x)**4
    F = -(a - 2*b*(p + 1))*(a + b*tan(e + f*x)**2 + b)**p*tan(e + f*x)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*tan(e + f*x)**2/(a + b))/(b*f*(2*p + 3)*(b*tan(e + f*x)**2/(a + b) + 1)**p) + (a + b*tan(e + f*x)**2 + b)**(p + 1)*tan(e + f*x)/(b*f*(2*p + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_307():
    f = (a + b*sec(e + f*x)**2)**p*sec(e + f*x)**2
    F = (a + b*tan(e + f*x)**2 + b)**p*tan(e + f*x)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*tan(e + f*x)**2/(a + b))/(f*(b*tan(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_308():
    f = (a + b*sec(e + f*x)**2)**p
    F = (a + b*tan(e + f*x)**2 + b)**p*tan(e + f*x)*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/(a + b))/(f*(b*tan(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_309():
    f = (a + b*sec(e + f*x)**2)**p*cos(e + f*x)**2
    F = (a + b*tan(e + f*x)**2 + b)**p*tan(e + f*x)*appellf1(sympy.S.Half, 2, -p, sympy.S(3)/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/(a + b))/(f*(b*tan(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_310():
    f = (a + b*sec(e + f*x)**2)**p*cos(e + f*x)**4
    F = (a + b*tan(e + f*x)**2 + b)**p*tan(e + f*x)*appellf1(sympy.S.Half, 3, -p, sympy.S(3)/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/(a + b))/(f*(b*tan(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_311():
    f = (a + b*sec(e + f*x)**2)**p*cos(e + f*x)**6
    F = (a + b*tan(e + f*x)**2 + b)**p*tan(e + f*x)*appellf1(sympy.S.Half, 4, -p, sympy.S(3)/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/(a + b))/(f*(b*tan(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_312():
    f = (a + b*sec(e + f*x)**2)*tan(e + f*x)**5
    F = -a*log(cos(e + f*x))/f + b*sec(e + f*x)**6/(6*f) + (a - 2*b)*sec(e + f*x)**4/(4*f) - (2*a - b)*sec(e + f*x)**2/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_313():
    f = (a + b*sec(e + f*x)**2)*tan(e + f*x)**3
    F = a*log(cos(e + f*x))/f + b*sec(e + f*x)**4/(4*f) + (a - b)*sec(e + f*x)**2/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_314():
    f = (a + b*sec(e + f*x)**2)*tan(e + f*x)
    F = -a*log(cos(e + f*x))/f + b*sec(e + f*x)**2/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_315():
    f = (a + b*sec(e + f*x)**2)*cot(e + f*x)
    F = -b*log(cos(e + f*x))/f + (a + b)*log(sin(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_316():
    f = (a + b*sec(e + f*x)**2)*cot(e + f*x)**3
    F = -a*log(sin(e + f*x))/f - (a + b)*csc(e + f*x)**2/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_317():
    f = (a + b*sec(e + f*x)**2)*cot(e + f*x)**5
    F = a*log(sin(e + f*x))/f - (a + b)*csc(e + f*x)**4/(4*f) + (2*a + b)*csc(e + f*x)**2/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_318():
    f = (a + b*sec(e + f*x)**2)*tan(e + f*x)**6
    F = -a*x + a*tan(e + f*x)**5/(5*f) - a*tan(e + f*x)**3/(3*f) + a*tan(e + f*x)/f + b*tan(e + f*x)**7/(7*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_319():
    f = (a + b*sec(e + f*x)**2)*tan(e + f*x)**4
    F = a*x + a*tan(e + f*x)**3/(3*f) - a*tan(e + f*x)/f + b*tan(e + f*x)**5/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_320():
    f = (a + b*sec(e + f*x)**2)*tan(e + f*x)**2
    F = -a*x + a*tan(e + f*x)/f + b*tan(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_321():
    f = a + b*sec(e + f*x)**2
    F = a*x + b*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_322():
    f = (a + b*sec(e + f*x)**2)*cot(e + f*x)**2
    F = -a*x - (a + b)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_323():
    f = (a + b*sec(e + f*x)**2)*cot(e + f*x)**4
    F = a*x + a*cot(e + f*x)/f - (a + b)*cot(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_324():
    f = (a + b*sec(e + f*x)**2)*cot(e + f*x)**6
    F = -a*x + a*cot(e + f*x)**3/(3*f) - a*cot(e + f*x)/f - (a + b)*cot(e + f*x)**5/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_325():
    f = (a + b*sec(e + f*x)**2)**2*tan(e + f*x)**5
    F = -a**2*log(cos(e + f*x))/f - a*(a - b)*sec(e + f*x)**2/f + b**2*sec(e + f*x)**8/(8*f) + b*(a - b)*sec(e + f*x)**6/(3*f) + (a**2 - 4*a*b + b**2)*sec(e + f*x)**4/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_326():
    f = (a + b*sec(e + f*x)**2)**2*tan(e + f*x)**3
    F = a**2*log(cos(e + f*x))/f + a*(a - 2*b)*sec(e + f*x)**2/(2*f) + b**2*sec(e + f*x)**6/(6*f) + b*(2*a - b)*sec(e + f*x)**4/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_327():
    f = (a + b*sec(e + f*x)**2)**2*tan(e + f*x)
    F = -a**2*log(cos(e + f*x))/f + a*b*sec(e + f*x)**2/f + b**2*sec(e + f*x)**4/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_328():
    f = (a + b*sec(e + f*x)**2)**2*cot(e + f*x)
    F = b**2*sec(e + f*x)**2/(2*f) - b*(2*a + b)*log(cos(e + f*x))/f + (a + b)**2*log(sin(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_329():
    f = (a + b*sec(e + f*x)**2)**2*cot(e + f*x)**3
    F = -b**2*log(cos(e + f*x))/f - (a + b)**2*csc(e + f*x)**2/(2*f) - (a**2 - b**2)*log(sin(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_330():
    f = (a + b*sec(e + f*x)**2)**2*cot(e + f*x)**5
    F = a**2*log(sin(e + f*x))/f + a*(a + b)*csc(e + f*x)**2/f - (a + b)**2*csc(e + f*x)**4/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_331():
    f = (a + b*sec(e + f*x)**2)**2*tan(e + f*x)**6
    F = -a**2*x + a**2*tan(e + f*x)**5/(5*f) - a**2*tan(e + f*x)**3/(3*f) + a**2*tan(e + f*x)/f + b**2*tan(e + f*x)**9/(9*f) + b*(2*a + b)*tan(e + f*x)**7/(7*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_332():
    f = (a + b*sec(e + f*x)**2)**2*tan(e + f*x)**4
    F = a**2*x + a**2*tan(e + f*x)**3/(3*f) - a**2*tan(e + f*x)/f + b**2*tan(e + f*x)**7/(7*f) + b*(2*a + b)*tan(e + f*x)**5/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_333():
    f = (a + b*sec(e + f*x)**2)**2*tan(e + f*x)**2
    F = -a**2*x + a**2*tan(e + f*x)/f + b**2*tan(e + f*x)**5/(5*f) + b*(2*a + b)*tan(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_334():
    f = (a + b*sec(e + f*x)**2)**2
    F = a**2*x + b**2*tan(e + f*x)**3/(3*f) + b*(2*a + b)*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_335():
    f = (a + b*sec(e + f*x)**2)**2*cot(e + f*x)**2
    F = -a**2*x + b**2*tan(e + f*x)/f - (a + b)**2*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_336():
    f = (a + b*sec(e + f*x)**2)**2*cot(e + f*x)**4
    F = a**2*x - (a + b)**2*cot(e + f*x)**3/(3*f) + (a**2 - b**2)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_337():
    f = (a + b*sec(e + f*x)**2)**2*cot(e + f*x)**6
    F = -a**2*x - a**2*cot(e + f*x)/f - (a + b)**2*cot(e + f*x)**5/(5*f) + (a**2 - b**2)*cot(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_338():
    f = tan(e + f*x)**5/(a + b*sec(e + f*x)**2)
    F = sec(e + f*x)**2/(2*b*f) + (a + 2*b)*log(cos(e + f*x))/(b**2*f) - (a + b)**2*log(a*cos(e + f*x)**2 + b)/(2*a*b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_339():
    f = tan(e + f*x)**3/(a + b*sec(e + f*x)**2)
    F = -log(cos(e + f*x))/(b*f) + (a + b)*log(a*cos(e + f*x)**2 + b)/(2*a*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_340():
    f = tan(e + f*x)/(a + b*sec(e + f*x)**2)
    F = -log(a*cos(e + f*x)**2 + b)/(2*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_341():
    f = cot(e + f*x)/(a + b*sec(e + f*x)**2)
    F = log(sin(e + f*x))/(f*(a + b)) + b*log(a*cos(e + f*x)**2 + b)/(2*a*f*(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_342():
    f = cot(e + f*x)**3/(a + b*sec(e + f*x)**2)
    F = -csc(e + f*x)**2/(f*(2*a + 2*b)) - (a + 2*b)*log(sin(e + f*x))/(f*(a + b)**2) - b**2*log(a*cos(e + f*x)**2 + b)/(2*a*f*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_343():
    f = cot(e + f*x)**5/(a + b*sec(e + f*x)**2)
    F = -csc(e + f*x)**4/(f*(4*a + 4*b)) + (2*a + 3*b)*csc(e + f*x)**2/(2*f*(a + b)**2) + (a**2 + 3*a*b + 3*b**2)*log(sin(e + f*x))/(f*(a + b)**3) + b**3*log(a*cos(e + f*x)**2 + b)/(2*a*f*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_344():
    f = tan(e + f*x)**6/(a + b*sec(e + f*x)**2)
    F = tan(e + f*x)**3/(3*b*f) - (a + 2*b)*tan(e + f*x)/(b**2*f) - x/a + (a + b)**(sympy.S(5)/2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(a*b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_345():
    f = tan(e + f*x)**4/(a + b*sec(e + f*x)**2)
    F = tan(e + f*x)/(b*f) + x/a - (a + b)**(sympy.S(3)/2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(a*b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_346():
    f = tan(e + f*x)**2/(a + b*sec(e + f*x)**2)
    F = -x/a + sqrt(a + b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(a*sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_347():
    f = 1/(a + b*sec(e + f*x)**2)
    F = sqrt(b)*atan(sqrt(a + b)*cot(e + f*x)/sqrt(b))/(a*f*sqrt(a + b)) + x/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_348():
    f = cot(e + f*x)**2/(a + b*sec(e + f*x)**2)
    F = -cot(e + f*x)/(f*(a + b)) + b**(sympy.S(3)/2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(a*f*(a + b)**(sympy.S(3)/2)) - x/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_349():
    f = cot(e + f*x)**4/(a + b*sec(e + f*x)**2)
    F = -cot(e + f*x)**3/(f*(3*a + 3*b)) + (a + 2*b)*cot(e + f*x)/(f*(a + b)**2) - b**(sympy.S(5)/2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(a*f*(a + b)**(sympy.S(5)/2)) + x/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_350():
    f = cot(e + f*x)**6/(a + b*sec(e + f*x)**2)
    F = -cot(e + f*x)**5/(f*(5*a + 5*b)) + (a + 2*b)*cot(e + f*x)**3/(3*f*(a + b)**2) - (a**2 + 3*a*b + 3*b**2)*cot(e + f*x)/(f*(a + b)**3) + b**(sympy.S(7)/2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(a*f*(a + b)**(sympy.S(7)/2)) - x/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_351():
    f = tan(e + f*x)**5/(a + b*sec(e + f*x)**2)**2
    F = -(-1/b**2 + a**(-2))*log(a*cos(e + f*x)**2 + b)/(2*f) - log(cos(e + f*x))/(b**2*f) - (a + b)**2/(2*a**2*b*f*(a*cos(e + f*x)**2 + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_352():
    f = tan(e + f*x)**3/(a + b*sec(e + f*x)**2)**2
    F = (a + b)/(2*a**2*f*(a*cos(e + f*x)**2 + b)) + log(a*cos(e + f*x)**2 + b)/(2*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_353():
    f = tan(e + f*x)/(a + b*sec(e + f*x)**2)**2
    F = -b/(2*a**2*f*(a*cos(e + f*x)**2 + b)) - log(a*cos(e + f*x)**2 + b)/(2*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_354():
    f = cot(e + f*x)/(a + b*sec(e + f*x)**2)**2
    F = log(sin(e + f*x))/(f*(a + b)**2) + b**2/(2*a**2*f*(a + b)*(a*cos(e + f*x)**2 + b)) + b*(2*a + b)*log(a*cos(e + f*x)**2 + b)/(2*a**2*f*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_355():
    f = cot(e + f*x)**3/(a + b*sec(e + f*x)**2)**2
    F = -csc(e + f*x)**2/(2*f*(a + b)**2) - (a + 3*b)*log(sin(e + f*x))/(f*(a + b)**3) - b**3/(2*a**2*f*(a + b)**2*(a*cos(e + f*x)**2 + b)) - b**2*(3*a + b)*log(a*cos(e + f*x)**2 + b)/(2*a**2*f*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_356():
    f = cot(e + f*x)**5/(a + b*sec(e + f*x)**2)**2
    F = -csc(e + f*x)**4/(4*f*(a + b)**2) + (a + 2*b)*csc(e + f*x)**2/(f*(a + b)**3) + (a**2 + 4*a*b + 6*b**2)*log(sin(e + f*x))/(f*(a + b)**4) + b**4/(2*a**2*f*(a + b)**3*(a*cos(e + f*x)**2 + b)) + b**3*(4*a + b)*log(a*cos(e + f*x)**2 + b)/(2*a**2*f*(a + b)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_357():
    f = tan(e + f*x)**6/(a + b*sec(e + f*x)**2)**2
    F = -(a + b)*tan(e + f*x)**3/(2*a*b*f*(a + b*tan(e + f*x)**2 + b)) + (3*a + b)*tan(e + f*x)/(2*a*b**2*f) - x/a**2 - (a + b)**(sympy.S(3)/2)*(3*a - 2*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(2*a**2*b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_358():
    f = tan(e + f*x)**4/(a + b*sec(e + f*x)**2)**2
    F = -(a + b)*tan(e + f*x)/(2*a*b*f*(a + b*tan(e + f*x)**2 + b)) + x/a**2 + (a - 2*b)*sqrt(a + b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(2*a**2*b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_359():
    f = tan(e + f*x)**2/(a + b*sec(e + f*x)**2)**2
    F = tan(e + f*x)/(2*a*f*(a + b*tan(e + f*x)**2 + b)) - x/a**2 + (a + 2*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(2*a**2*sqrt(b)*f*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_360():
    f = (a + b*sec(e + f*x)**2)**(-2)
    F = -b*tan(e + f*x)/(2*a*f*(a + b)*(a + b*tan(e + f*x)**2 + b)) - sqrt(b)*(3*a + 2*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(2*a**2*f*(a + b)**(sympy.S(3)/2)) + x/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_361():
    f = cot(e + f*x)**2/(a + b*sec(e + f*x)**2)**2
    F = -b*cot(e + f*x)/(2*a*f*(a + b)*(a + b*tan(e + f*x)**2 + b)) - (2*a - b)*cot(e + f*x)/(2*a*f*(a + b)**2) + b**(sympy.S(3)/2)*(5*a + 2*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(2*a**2*f*(a + b)**(sympy.S(5)/2)) - x/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_362():
    f = cot(e + f*x)**4/(a + b*sec(e + f*x)**2)**2
    F = -b*cot(e + f*x)**3/(2*a*f*(a + b)*(a + b*tan(e + f*x)**2 + b)) - (2*a - 3*b)*cot(e + f*x)**3/(6*a*f*(a + b)**2) + (2*a**2 + 6*a*b - b**2)*cot(e + f*x)/(2*a*f*(a + b)**3) - b**(sympy.S(5)/2)*(7*a + 2*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(2*a**2*f*(a + b)**(sympy.S(7)/2)) + x/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_363():
    f = cot(e + f*x)**6/(a + b*sec(e + f*x)**2)**2
    F = -b*cot(e + f*x)**5/(2*a*f*(a + b)*(a + b*tan(e + f*x)**2 + b)) - (2*a - 5*b)*cot(e + f*x)**5/(10*a*f*(a + b)**2) + (2*a**2 + 6*a*b - 3*b**2)*cot(e + f*x)**3/(6*a*f*(a + b)**3) - (2*a**3 + 8*a**2*b + 12*a*b**2 - b**3)*cot(e + f*x)/(2*a*f*(a + b)**4) + b**(sympy.S(7)/2)*(9*a + 2*b)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(2*a**2*f*(a + b)**(sympy.S(9)/2)) - x/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_364():
    f = tan(e + f*x)**5/(a + b*sec(e + f*x)**2)**3
    F = (a + b)**2/(4*a**3*f*(a*cos(e + f*x)**2 + b)**2) - (a + b)/(a**3*f*(a*cos(e + f*x)**2 + b)) - log(a*cos(e + f*x)**2 + b)/(2*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_365():
    f = tan(e + f*x)**3/(a + b*sec(e + f*x)**2)**3
    F = -b*(a + b)/(4*a**3*f*(a*cos(e + f*x)**2 + b)**2) + (a + 2*b)/(2*a**3*f*(a*cos(e + f*x)**2 + b)) + log(a*cos(e + f*x)**2 + b)/(2*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_366():
    f = tan(e + f*x)/(a + b*sec(e + f*x)**2)**3
    F = b**2/(4*a**3*f*(a*cos(e + f*x)**2 + b)**2) - b/(a**3*f*(a*cos(e + f*x)**2 + b)) - log(a*cos(e + f*x)**2 + b)/(2*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_367():
    f = cot(e + f*x)/(a + b*sec(e + f*x)**2)**3
    F = log(sin(e + f*x))/(f*(a + b)**3) - b**3/(4*a**3*f*(a + b)*(a*cos(e + f*x)**2 + b)**2) + b**2*(3*a + 2*b)/(2*a**3*f*(a + b)**2*(a*cos(e + f*x)**2 + b)) + b*(3*a**2 + 3*a*b + b**2)*log(a*cos(e + f*x)**2 + b)/(2*a**3*f*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_368():
    f = cot(e + f*x)**3/(a + b*sec(e + f*x)**2)**3
    F = -csc(e + f*x)**2/(2*f*(a + b)**3) - (a + 4*b)*log(sin(e + f*x))/(f*(a + b)**4) + b**4/(4*a**3*f*(a + b)**2*(a*cos(e + f*x)**2 + b)**2) - b**3*(2*a + b)/(a**3*f*(a + b)**3*(a*cos(e + f*x)**2 + b)) - b**2*(6*a**2 + 4*a*b + b**2)*log(a*cos(e + f*x)**2 + b)/(2*a**3*f*(a + b)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_369():
    f = cot(e + f*x)**5/(a + b*sec(e + f*x)**2)**3
    F = -csc(e + f*x)**4/(4*f*(a + b)**3) + (2*a + 5*b)*csc(e + f*x)**2/(2*f*(a + b)**4) + (a**2 + 5*a*b + 10*b**2)*log(sin(e + f*x))/(f*(a + b)**5) - b**5/(4*a**3*f*(a + b)**3*(a*cos(e + f*x)**2 + b)**2) + b**4*(5*a + 2*b)/(2*a**3*f*(a + b)**4*(a*cos(e + f*x)**2 + b)) + b**3*(10*a**2 + 5*a*b + b**2)*log(a*cos(e + f*x)**2 + b)/(2*a**3*f*(a + b)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_370():
    f = tan(e + f*x)**6/(a + b*sec(e + f*x)**2)**3
    F = -(a + b)*tan(e + f*x)**3/(4*a*b*f*(a + b*tan(e + f*x)**2 + b)**2) - (a + b)*(3*a - 4*b)*tan(e + f*x)/(8*a**2*b**2*f*(a + b*tan(e + f*x)**2 + b)) - x/a**3 + sqrt(a + b)*(3*a**2 - 4*a*b + 8*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(8*a**3*b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_371():
    f = tan(e + f*x)**4/(a + b*sec(e + f*x)**2)**3
    F = -(a + b)*tan(e + f*x)/(4*a*b*f*(a + b*tan(e + f*x)**2 + b)**2) + (a - 4*b)*tan(e + f*x)/(8*a**2*b*f*(a + b*tan(e + f*x)**2 + b)) + x/a**3 + (a**2 - 4*a*b - 8*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(8*a**3*b**(sympy.S(3)/2)*f*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_372():
    f = tan(e + f*x)**2/(a + b*sec(e + f*x)**2)**3
    F = tan(e + f*x)/(4*a*f*(a + b*tan(e + f*x)**2 + b)**2) + (3*a + 4*b)*tan(e + f*x)/(8*a**2*f*(a + b)*(a + b*tan(e + f*x)**2 + b)) - x/a**3 + (3*a**2 + 12*a*b + 8*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(8*a**3*sqrt(b)*f*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_373():
    f = (a + b*sec(e + f*x)**2)**(-3)
    F = -b*tan(e + f*x)/(4*a*f*(a + b)*(a + b*tan(e + f*x)**2 + b)**2) - b*(7*a + 4*b)*tan(e + f*x)/(8*a**2*f*(a + b)**2*(a + b*tan(e + f*x)**2 + b)) - sqrt(b)*(15*a**2 + 20*a*b + 8*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(8*a**3*f*(a + b)**(sympy.S(5)/2)) + x/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_374():
    f = cot(e + f*x)**2/(a + b*sec(e + f*x)**2)**3
    F = -b*cot(e + f*x)/(4*a*f*(a + b)*(a + b*tan(e + f*x)**2 + b)**2) - b*(9*a + 4*b)*cot(e + f*x)/(8*a**2*f*(a + b)**2*(a + b*tan(e + f*x)**2 + b)) - (8*a**2 - 11*a*b - 4*b**2)*cot(e + f*x)/(8*a**2*f*(a + b)**3) + b**(sympy.S(3)/2)*(35*a**2 + 28*a*b + 8*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(8*a**3*f*(a + b)**(sympy.S(7)/2)) - x/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_375():
    f = cot(e + f*x)**4/(a + b*sec(e + f*x)**2)**3
    F = -b*cot(e + f*x)**3/(4*a*f*(a + b)*(a + b*tan(e + f*x)**2 + b)**2) - b*(11*a + 4*b)*cot(e + f*x)**3/(8*a**2*f*(a + b)**2*(a + b*tan(e + f*x)**2 + b)) - (8*a**2 - 39*a*b - 12*b**2)*cot(e + f*x)**3/(24*a**2*f*(a + b)**3) + (8*a**3 + 32*a**2*b - 15*a*b**2 - 4*b**3)*cot(e + f*x)/(8*a**2*f*(a + b)**4) - b**(sympy.S(5)/2)*(63*a**2 + 36*a*b + 8*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(8*a**3*f*(a + b)**(sympy.S(9)/2)) + x/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_376():
    f = cot(e + f*x)**6/(a + b*sec(e + f*x)**2)**3
    F = -b*cot(e + f*x)**5/(4*a*f*(a + b)*(a + b*tan(e + f*x)**2 + b)**2) - b*(13*a + 4*b)*cot(e + f*x)**5/(8*a**2*f*(a + b)**2*(a + b*tan(e + f*x)**2 + b)) - (8*a**2 - 75*a*b - 20*b**2)*cot(e + f*x)**5/(40*a**2*f*(a + b)**3) + (8*a**3 + 32*a**2*b - 51*a*b**2 - 12*b**3)*cot(e + f*x)**3/(24*a**2*f*(a + b)**4) - (8*a**4 + 40*a**3*b + 80*a**2*b**2 - 19*a*b**3 - 4*b**4)*cot(e + f*x)/(8*a**2*f*(a + b)**5) + b**(sympy.S(7)/2)*(99*a**2 + 44*a*b + 8*b**2)*atan(sqrt(b)*tan(e + f*x)/sqrt(a + b))/(8*a**3*f*(a + b)**(sympy.S(11)/2)) - x/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_377():
    f = sqrt(a + b*sec(e + f*x)**2)*tan(e + f*x)**5
    F = -sqrt(a)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/f + sqrt(a + b*sec(e + f*x)**2)/f - (a + 2*b)*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)/(3*b**2*f) + (a + b*sec(e + f*x)**2)**(sympy.S(5)/2)/(5*b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_378():
    f = sqrt(a + b*sec(e + f*x)**2)*tan(e + f*x)**3
    F = sqrt(a)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/f - sqrt(a + b*sec(e + f*x)**2)/f + (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)/(3*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_379():
    f = sqrt(a + b*sec(e + f*x)**2)*tan(e + f*x)
    F = -sqrt(a)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/f + sqrt(a + b*sec(e + f*x)**2)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_380():
    f = sqrt(a + b*sec(e + f*x)**2)*cot(e + f*x)
    F = sqrt(a)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/f - sqrt(a + b)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a + b))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_381():
    f = sqrt(a + b*sec(e + f*x)**2)*cot(e + f*x)**3
    F = -sqrt(a)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/f - sqrt(a + b*sec(e + f*x)**2)*cot(e + f*x)**2/(2*f) + (2*a + b)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a + b))/(2*f*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_382():
    f = sqrt(a + b*sec(e + f*x)**2)*cot(e + f*x)**5
    F = sqrt(a)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/f + sqrt(a + b*sec(e + f*x)**2)*(4*a + 3*b)*cot(e + f*x)**2/(f*(8*a + 8*b)) - sqrt(a + b*sec(e + f*x)**2)*cot(e + f*x)**4/(4*f) - (8*a**2 + 12*a*b + 3*b**2)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a + b))/(8*f*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_383():
    f = sqrt(a + b*sec(e + f*x)**2)*tan(e + f*x)**6
    F = -sqrt(a)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f + sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)**5/(6*f) + (a - 5*b)*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)**3/(24*b*f) - (a - b)*(a + 5*b)*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/(16*b**2*f) + (a**3 + 5*a**2*b + 15*a*b**2 - 5*b**3)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(16*b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_384():
    f = sqrt(a + b*sec(e + f*x)**2)*tan(e + f*x)**4
    F = sqrt(a)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f + sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)**3/(4*f) + (a - 3*b)*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/(8*b*f) - (a**2 + 6*a*b - 3*b**2)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(8*b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_385():
    f = sqrt(a + b*sec(e + f*x)**2)*tan(e + f*x)**2
    F = -sqrt(a)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f + sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/(2*f) + (a - b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_386():
    f = sqrt(a + b*sec(e + f*x)**2)
    F = sqrt(a)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f + sqrt(b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_387():
    f = sqrt(a + b*sec(e + f*x)**2)*cot(e + f*x)**2
    F = -sqrt(a)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f - sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_388():
    f = sqrt(a + b*sec(e + f*x)**2)*cot(e + f*x)**4
    F = sqrt(a)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f + (3*a + 2*b)*sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)/(f*(3*a + 3*b)) - sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_389():
    f = sqrt(a + b*sec(e + f*x)**2)*cot(e + f*x)**6
    F = -sqrt(a)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f - (-5*a - 4*b)*sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)**3/(f*(15*a + 15*b)) - sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)**5/(5*f) - sqrt(a + b*tan(e + f*x)**2 + b)*(15*a**2 + 25*a*b + 8*b**2)*cot(e + f*x)/(15*f*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_390():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*tan(e + f*x)**5
    F = -a**(sympy.S(3)/2)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/f + a*sqrt(a + b*sec(e + f*x)**2)/f + (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)/(3*f) - (a + 2*b)*(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)/(5*b**2*f) + (a + b*sec(e + f*x)**2)**(sympy.S(7)/2)/(7*b**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_391():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*tan(e + f*x)**3
    F = a**(sympy.S(3)/2)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/f - a*sqrt(a + b*sec(e + f*x)**2)/f - (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)/(3*f) + (a + b*sec(e + f*x)**2)**(sympy.S(5)/2)/(5*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_392():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*tan(e + f*x)
    F = -a**(sympy.S(3)/2)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/f + a*sqrt(a + b*sec(e + f*x)**2)/f + (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_393():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)
    F = a**(sympy.S(3)/2)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/f + b*sqrt(a + b*sec(e + f*x)**2)/f - (a + b)**(sympy.S(3)/2)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a + b))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_394():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)**3
    F = -a**(sympy.S(3)/2)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/f + sqrt(a + b)*(2*a - b)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a + b))/(2*f) - (a + b)*sqrt(a + b*sec(e + f*x)**2)*cot(e + f*x)**2/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_395():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)**5
    F = a**(sympy.S(3)/2)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/f - (a + b)*sqrt(a + b*sec(e + f*x)**2)*cot(e + f*x)**4/(4*f) + sqrt(a + b*sec(e + f*x)**2)*(4*a - b)*cot(e + f*x)**2/(8*f) - (8*a**2 + 4*a*b - b**2)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a + b))/(8*f*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_396():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*tan(e + f*x)**6
    F = -a**(sympy.S(3)/2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f + b*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)**7/(8*f) + (9*a + b)*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)**5/(48*f) + sqrt(a + b*tan(e + f*x)**2 + b)*(3*a**2 - 50*a*b - 5*b**2)*tan(e + f*x)**3/(192*b*f) - sqrt(a + b*tan(e + f*x)**2 + b)*(3*a**3 + 17*a**2*b - 55*a*b**2 - 5*b**3)*tan(e + f*x)/(128*b**2*f) + (3*a**4 + 20*a**3*b + 90*a**2*b**2 - 60*a*b**3 - 5*b**4)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(128*b**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_397():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*tan(e + f*x)**4
    F = a**(sympy.S(3)/2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f + b*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)**5/(6*f) + (7*a + b)*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)**3/(24*f) + sqrt(a + b*tan(e + f*x)**2 + b)*(a**2 - 8*a*b - b**2)*tan(e + f*x)/(16*b*f) - (a - b)*(a**2 + 10*a*b + b**2)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(16*b**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_398():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*tan(e + f*x)**2
    F = -a**(sympy.S(3)/2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f + b*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)**3/(4*f) + (5*a + b)*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/(8*f) + (3*a**2 - 6*a*b - b**2)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(8*sqrt(b)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_399():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = a**(sympy.S(3)/2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f + sqrt(b)*(3*a + b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*f) + b*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_400():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)**2
    F = -a**(sympy.S(3)/2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f + b**(sympy.S(3)/2)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f - (a + b)*sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_401():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)**4
    F = a**(sympy.S(3)/2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f - (a + b)*sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)**3/(3*f) + (3*a - b)*sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_402():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*cot(e + f*x)**6
    F = -a**(sympy.S(3)/2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/f - (a + b)*sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)**5/(5*f) + (5*a - b)*sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)**3/(15*f) - sqrt(a + b*tan(e + f*x)**2 + b)*(15*a**2 + 10*a*b - 2*b**2)*cot(e + f*x)/(f*(15*a + 15*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_403():
    f = tan(e + f*x)**5/sqrt(a + b*sec(e + f*x)**2)
    F = -(a + 2*b)*sqrt(a + b*sec(e + f*x)**2)/(b**2*f) + (a + b*sec(e + f*x)**2)**(sympy.S(3)/2)/(3*b**2*f) - atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_404():
    f = tan(e + f*x)**3/sqrt(a + b*sec(e + f*x)**2)
    F = sqrt(a + b*sec(e + f*x)**2)/(b*f) + atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_405():
    f = tan(e + f*x)/sqrt(a + b*sec(e + f*x)**2)
    F = -atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_406():
    f = cot(e + f*x)/sqrt(a + b*sec(e + f*x)**2)
    F = -atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a + b))/(f*sqrt(a + b)) + atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_407():
    f = cot(e + f*x)**3/sqrt(a + b*sec(e + f*x)**2)
    F = -sqrt(a + b*sec(e + f*x)**2)*cot(e + f*x)**2/(f*(2*a + 2*b)) + (2*a + 3*b)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a + b))/(2*f*(a + b)**(sympy.S(3)/2)) - atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_408():
    f = cot(e + f*x)**5/sqrt(a + b*sec(e + f*x)**2)
    F = -sqrt(a + b*sec(e + f*x)**2)*cot(e + f*x)**4/(f*(4*a + 4*b)) + sqrt(a + b*sec(e + f*x)**2)*(4*a + 7*b)*cot(e + f*x)**2/(8*f*(a + b)**2) - (8*a**2 + 20*a*b + 15*b**2)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a + b))/(8*f*(a + b)**(sympy.S(5)/2)) + atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_409():
    f = tan(e + f*x)**6/sqrt(a + b*sec(e + f*x)**2)
    F = sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)**3/(4*b*f) - (3*a + 7*b)*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/(8*b**2*f) + (3*a**2 + 10*a*b + 15*b**2)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(8*b**(sympy.S(5)/2)*f) - atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_410():
    f = tan(e + f*x)**4/sqrt(a + b*sec(e + f*x)**2)
    F = sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/(2*b*f) - (a + 3*b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*b**(sympy.S(3)/2)*f) + atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_411():
    f = tan(e + f*x)**2/sqrt(a + b*sec(e + f*x)**2)
    F = atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(sqrt(b)*f) - atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_412():
    f = 1/sqrt(a + b*sec(e + f*x)**2)
    F = atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_413():
    f = cot(e + f*x)**2/sqrt(a + b*sec(e + f*x)**2)
    F = -sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)/(f*(a + b)) - atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_414():
    f = cot(e + f*x)**4/sqrt(a + b*sec(e + f*x)**2)
    F = -sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)**3/(f*(3*a + 3*b)) + (3*a + 5*b)*sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)/(3*f*(a + b)**2) + atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_415():
    f = cot(e + f*x)**6/sqrt(a + b*sec(e + f*x)**2)
    F = -sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)**5/(f*(5*a + 5*b)) + (5*a + 9*b)*sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)**3/(15*f*(a + b)**2) - sqrt(a + b*tan(e + f*x)**2 + b)*(15*a**2 + 40*a*b + 33*b**2)*cot(e + f*x)/(15*f*(a + b)**3) - atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_416():
    f = tan(e + f*x)**5/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = sqrt(a + b*sec(e + f*x)**2)/(b**2*f) + (a + b)**2/(a*b**2*f*sqrt(a + b*sec(e + f*x)**2)) - atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_417():
    f = tan(e + f*x)**3/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -(a + b)/(a*b*f*sqrt(a + b*sec(e + f*x)**2)) + atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_418():
    f = tan(e + f*x)/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = 1/(a*f*sqrt(a + b*sec(e + f*x)**2)) - atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_419():
    f = cot(e + f*x)/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a + b))/(f*(a + b)**(sympy.S(3)/2)) - b/(a*f*(a + b)*sqrt(a + b*sec(e + f*x)**2)) + atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_420():
    f = cot(e + f*x)**3/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -cot(e + f*x)**2/(f*sqrt(a + b*sec(e + f*x)**2)*(2*a + 2*b)) + (2*a + 5*b)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a + b))/(2*f*(a + b)**(sympy.S(5)/2)) - b*(a - 2*b)/(2*a*f*(a + b)**2*sqrt(a + b*sec(e + f*x)**2)) - atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_421():
    f = cot(e + f*x)**5/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -cot(e + f*x)**4/(f*sqrt(a + b*sec(e + f*x)**2)*(4*a + 4*b)) + (4*a + 9*b)*cot(e + f*x)**2/(8*f*(a + b)**2*sqrt(a + b*sec(e + f*x)**2)) - (8*a**2 + 28*a*b + 35*b**2)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a + b))/(8*f*(a + b)**(sympy.S(7)/2)) + b*(4*a**2 + 11*a*b - 8*b**2)/(8*a*f*(a + b)**3*sqrt(a + b*sec(e + f*x)**2)) + atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_422():
    f = tan(e + f*x)**6/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -(3*a + 5*b)*atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(2*b**(sympy.S(5)/2)*f) - (a + b)*tan(e + f*x)**3/(a*b*f*sqrt(a + b*tan(e + f*x)**2 + b)) + (3*a + 2*b)*sqrt(a + b*tan(e + f*x)**2 + b)*tan(e + f*x)/(2*a*b**2*f) - atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_423():
    f = tan(e + f*x)**4/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(b**(sympy.S(3)/2)*f) - (a + b)*tan(e + f*x)/(a*b*f*sqrt(a + b*tan(e + f*x)**2 + b)) + atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_424():
    f = tan(e + f*x)**2/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = tan(e + f*x)/(a*f*sqrt(a + b*tan(e + f*x)**2 + b)) - atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_425():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(-3)/2)
    F = -b*tan(e + f*x)/(a*f*(a + b)*sqrt(a + b*tan(e + f*x)**2 + b)) + atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_426():
    f = cot(e + f*x)**2/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -b*cot(e + f*x)/(a*f*(a + b)*sqrt(a + b*tan(e + f*x)**2 + b)) - (a - b)*sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)/(a*f*(a + b)**2) - atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_427():
    f = cot(e + f*x)**4/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -b*cot(e + f*x)**3/(a*f*(a + b)*sqrt(a + b*tan(e + f*x)**2 + b)) - (a - 3*b)*sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)**3/(3*a*f*(a + b)**2) + (a + 3*b)*(3*a - b)*sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)/(3*a*f*(a + b)**3) + atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_428():
    f = cot(e + f*x)**6/(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)
    F = -b*cot(e + f*x)**5/(a*f*(a + b)*sqrt(a + b*tan(e + f*x)**2 + b)) - (a - 5*b)*sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)**5/(5*a*f*(a + b)**2) + sqrt(a + b*tan(e + f*x)**2 + b)*(5*a**2 + 14*a*b - 15*b**2)*cot(e + f*x)**3/(15*a*f*(a + b)**3) - sqrt(a + b*tan(e + f*x)**2 + b)*(15*a**3 + 55*a**2*b + 73*a*b**2 - 15*b**3)*cot(e + f*x)/(15*a*f*(a + b)**4) - atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_429():
    f = tan(e + f*x)**5/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = (-1/b**2 + a**(-2))/(f*sqrt(a + b*sec(e + f*x)**2)) + (a + b)**2/(3*a*b**2*f*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)) - atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_430():
    f = tan(e + f*x)**3/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -(a + b)/(3*a*b*f*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)) - 1/(a**2*f*sqrt(a + b*sec(e + f*x)**2)) + atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_431():
    f = tan(e + f*x)/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = 1/(3*a*f*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)) + 1/(a**2*f*sqrt(a + b*sec(e + f*x)**2)) - atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_432():
    f = cot(e + f*x)/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a + b))/(f*(a + b)**(sympy.S(5)/2)) - b/(3*a*f*(a + b)*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)) - b*(2*a + b)/(a**2*f*(a + b)**2*sqrt(a + b*sec(e + f*x)**2)) + atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_433():
    f = cot(e + f*x)**3/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -cot(e + f*x)**2/(f*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*(2*a + 2*b)) + (2*a + 7*b)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a + b))/(2*f*(a + b)**(sympy.S(7)/2)) - b*(3*a - 2*b)/(6*a*f*(a + b)**2*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)) - b*(a**2 - 6*a*b - 2*b**2)/(2*a**2*f*(a + b)**3*sqrt(a + b*sec(e + f*x)**2)) - atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_434():
    f = cot(e + f*x)**5/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -cot(e + f*x)**4/(f*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)*(4*a + 4*b)) + (4*a + 11*b)*cot(e + f*x)**2/(8*f*(a + b)**2*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)) - (8*a**2 + 36*a*b + 63*b**2)*atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a + b))/(8*f*(a + b)**(sympy.S(9)/2)) + b*(12*a**2 + 39*a*b - 8*b**2)/(24*a*f*(a + b)**3*(a + b*sec(e + f*x)**2)**(sympy.S(3)/2)) + b*(4*a**3 + 15*a**2*b - 32*a*b**2 - 8*b**3)/(8*a**2*f*(a + b)**4*sqrt(a + b*sec(e + f*x)**2)) + atanh(sqrt(a + b*sec(e + f*x)**2)/sqrt(a))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_435():
    f = tan(e + f*x)**6/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = (-1/b**2 + a**(-2))*tan(e + f*x)/(f*sqrt(a + b*tan(e + f*x)**2 + b)) + atanh(sqrt(b)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(b**(sympy.S(5)/2)*f) - (a + b)*tan(e + f*x)**3/(3*a*b*f*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_436():
    f = tan(e + f*x)**4/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -(a + b)*tan(e + f*x)/(3*a*b*f*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) + (a - 3*b)*tan(e + f*x)/(3*a**2*b*f*sqrt(a + b*tan(e + f*x)**2 + b)) + atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_437():
    f = tan(e + f*x)**2/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = tan(e + f*x)/(3*a*f*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) + (2*a + 3*b)*tan(e + f*x)/(3*a**2*f*(a + b)*sqrt(a + b*tan(e + f*x)**2 + b)) - atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_438():
    f = (a + b*sec(e + f*x)**2)**(sympy.S(-5)/2)
    F = -b*tan(e + f*x)/(3*a*f*(a + b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - b*(5*a + 3*b)*tan(e + f*x)/(3*a**2*f*(a + b)**2*sqrt(a + b*tan(e + f*x)**2 + b)) + atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_439():
    f = cot(e + f*x)**2/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -b*cot(e + f*x)/(3*a*f*(a + b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - b*(7*a + 3*b)*cot(e + f*x)/(3*a**2*f*(a + b)**2*sqrt(a + b*tan(e + f*x)**2 + b)) - (a - 3*b)*(3*a + b)*sqrt(a + b*tan(e + f*x)**2 + b)*cot(e + f*x)/(3*a**2*f*(a + b)**3) - atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_440():
    f = cot(e + f*x)**4/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -b*cot(e + f*x)**3/(3*a*f*(a + b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - b*(3*a + b)*cot(e + f*x)**3/(a**2*f*(a + b)**2*sqrt(a + b*tan(e + f*x)**2 + b)) + (a - b)*sqrt(a + b*tan(e + f*x)**2 + b)*(3*a**2 + 14*a*b + 3*b**2)*cot(e + f*x)/(3*a**2*f*(a + b)**4) - sqrt(a + b*tan(e + f*x)**2 + b)*(a**2 - 10*a*b - 3*b**2)*cot(e + f*x)**3/(3*a**2*f*(a + b)**3) + atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_441():
    f = cot(e + f*x)**6/(a + b*sec(e + f*x)**2)**(sympy.S(5)/2)
    F = -b*cot(e + f*x)**5/(3*a*f*(a + b)*(a + b*tan(e + f*x)**2 + b)**(sympy.S(3)/2)) - b*(11*a + 3*b)*cot(e + f*x)**5/(3*a**2*f*(a + b)**2*sqrt(a + b*tan(e + f*x)**2 + b)) - sqrt(a + b*tan(e + f*x)**2 + b)*(a**2 - 20*a*b - 5*b**2)*cot(e + f*x)**5/(5*a**2*f*(a + b)**3) + sqrt(a + b*tan(e + f*x)**2 + b)*(5*a**3 + 19*a**2*b - 65*a*b**2 - 15*b**3)*cot(e + f*x)**3/(15*a**2*f*(a + b)**4) - sqrt(a + b*tan(e + f*x)**2 + b)*(15*a**4 + 70*a**3*b + 128*a**2*b**2 - 70*a*b**3 - 15*b**4)*cot(e + f*x)/(15*a**2*f*(a + b)**5) - atan(sqrt(a)*tan(e + f*x)/sqrt(a + b*tan(e + f*x)**2 + b))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_442():
    f = (d*tan(e + f*x))**m*(a + b*sec(e + f*x)**2)**p
    F = (d*tan(e + f*x))**(m + 1)*(a + b*tan(e + f*x)**2 + b)**p*appellf1(m/2 + sympy.S.Half, 1, -p, m/2 + sympy.S(3)/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/(a + b))/(d*f*(m + 1)*(b*tan(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_443():
    f = (a + b*sec(e + f*x)**2)**p*tan(e + f*x)**5
    F = -(a + 2*b)*(a + b*sec(e + f*x)**2)**(p + 1)/(2*b**2*f*(p + 1)) + (a + b*sec(e + f*x)**2)**(p + 2)/(2*b**2*f*(p + 2)) - (a + b*sec(e + f*x)**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*sec(e + f*x)**2/a)/(2*a*f*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_444():
    f = (a + b*sec(e + f*x)**2)**p*tan(e + f*x)**3
    F = (a + b*sec(e + f*x)**2)**(p + 1)/(2*b*f*(p + 1)) + (a + b*sec(e + f*x)**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*sec(e + f*x)**2/a)/(2*a*f*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_445():
    f = (a + b*sec(e + f*x)**2)**p*tan(e + f*x)
    F = -(a + b*sec(e + f*x)**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*sec(e + f*x)**2/a)/(2*a*f*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_446():
    f = (a + b*sec(e + f*x)**2)**p*cot(e + f*x)
    F = -(a + b*sec(e + f*x)**2)**(p + 1)*hyper((1, p + 1), (p + 2,), (a + b*sec(e + f*x)**2)/(a + b))/(f*(2*a + 2*b)*(p + 1)) + (a + b*sec(e + f*x)**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*sec(e + f*x)**2/a)/(2*a*f*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_447():
    f = (a + b*sec(e + f*x)**2)**p*cot(e + f*x)**3
    F = -(a + b*sec(e + f*x)**2)**(p + 1)*cot(e + f*x)**2/(f*(2*a + 2*b)) + (a + b*sec(e + f*x)**2)**(p + 1)*(a - b*p + b)*hyper((1, p + 1), (p + 2,), (a + b*sec(e + f*x)**2)/(a + b))/(2*f*(a + b)**2*(p + 1)) - (a + b*sec(e + f*x)**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*sec(e + f*x)**2/a)/(2*a*f*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_448():
    f = (a + b*sec(e + f*x)**2)**p*tan(e + f*x)**4
    F = (a + b*tan(e + f*x)**2 + b)**p*tan(e + f*x)**5*appellf1(sympy.S(5)/2, 1, -p, sympy.S(7)/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/(a + b))/(5*f*(b*tan(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_449():
    f = (a + b*sec(e + f*x)**2)**p*tan(e + f*x)**2
    F = (a + b*tan(e + f*x)**2 + b)**p*tan(e + f*x)**3*appellf1(sympy.S(3)/2, 1, -p, sympy.S(5)/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/(a + b))/(3*f*(b*tan(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_450():
    f = (a + b*sec(e + f*x)**2)**p
    F = (a + b*tan(e + f*x)**2 + b)**p*tan(e + f*x)*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/(a + b))/(f*(b*tan(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_451():
    f = (a + b*sec(e + f*x)**2)**p*cot(e + f*x)**2
    F = -(a + b*tan(e + f*x)**2 + b)**p*cot(e + f*x)*appellf1(sympy.S(-1)/2, 1, -p, sympy.S.Half, -tan(e + f*x)**2, -b*tan(e + f*x)**2/(a + b))/(f*(b*tan(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_452():
    f = (a + b*sec(e + f*x)**2)**p*cot(e + f*x)**4
    F = -(a + b*tan(e + f*x)**2 + b)**p*cot(e + f*x)**3*appellf1(sympy.S(-3)/2, 1, -p, sympy.S(-1)/2, -tan(e + f*x)**2, -b*tan(e + f*x)**2/(a + b))/(3*f*(b*tan(e + f*x)**2/(a + b) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_453():
    f = (a + b*sec(e + f*x)**3)*tan(e + f*x)**5
    F = -a*log(cos(e + f*x))/f + a*sec(e + f*x)**4/(4*f) - a*sec(e + f*x)**2/f + b*sec(e + f*x)**7/(7*f) - 2*b*sec(e + f*x)**5/(5*f) + b*sec(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_454():
    f = (a + b*sec(e + f*x)**3)*tan(e + f*x)**3
    F = a*log(cos(e + f*x))/f + a*sec(e + f*x)**2/(2*f) + b*sec(e + f*x)**5/(5*f) - b*sec(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_455():
    f = (a + b*sec(e + f*x)**3)*tan(e + f*x)
    F = -a*log(cos(e + f*x))/f + b*sec(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_456():
    f = (a + b*sec(e + f*x)**3)*cot(e + f*x)
    F = b*sec(e + f*x)/f + (a - b)*log(cos(e + f*x) + 1)/(2*f) + (a + b)*log(1 - cos(e + f*x))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_457():
    f = (a + b*sec(e + f*x)**3)*cot(e + f*x)**3
    F = -(a + b*cos(e + f*x))*csc(e + f*x)**2/(2*f) - (2*a - b)*log(1 - cos(e + f*x))/(4*f) - (2*a + b)*log(cos(e + f*x) + 1)/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_458():
    f = tan(e + f*x)**5/(a + b*sec(e + f*x)**3)
    F = sec(e + f*x)/(b*f) - log(a*cos(e + f*x)**3 + b)/(3*a*f) - (a**(sympy.S(2)/3) - 2*b**(sympy.S(2)/3))*log(a**(sympy.S(1)/3)*cos(e + f*x) + b**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)*b**(sympy.S(4)/3)*f) + (a**(sympy.S(2)/3) - 2*b**(sympy.S(2)/3))*log(a**(sympy.S(2)/3)*cos(e + f*x)**2 - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*cos(e + f*x) + b**(sympy.S(2)/3))/(6*a**(sympy.S(1)/3)*b**(sympy.S(4)/3)*f) - sqrt(3)*(a**(sympy.S(2)/3) + 2*b**(sympy.S(2)/3))*atan(sqrt(3)*(-2*a**(sympy.S(1)/3)*cos(e + f*x) + b**(sympy.S(1)/3))/(3*b**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3)*b**(sympy.S(4)/3)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_459():
    f = tan(e + f*x)**3/(a + b*sec(e + f*x)**3)
    F = log(a*cos(e + f*x)**3 + b)/(3*a*f) - log(a**(sympy.S(1)/3)*cos(e + f*x) + b**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*f) + log(a**(sympy.S(2)/3)*cos(e + f*x)**2 - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*cos(e + f*x) + b**(sympy.S(2)/3))/(6*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*f) + sqrt(3)*atan(sqrt(3)*(-2*a**(sympy.S(1)/3)*cos(e + f*x) + b**(sympy.S(1)/3))/(3*b**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_460():
    f = tan(e + f*x)/(a + b*sec(e + f*x)**3)
    F = -log(a*cos(e + f*x)**3 + b)/(3*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_461():
    f = cot(e + f*x)/(a + b*sec(e + f*x)**3)
    F = log(1 - cos(e + f*x))/(f*(2*a + 2*b)) + log(cos(e + f*x) + 1)/(f*(2*a - 2*b)) - b**2*log(a*cos(e + f*x)**3 + b)/(3*a*f*(a**2 - b**2)) - b**(sympy.S(2)/3)*(a**(sympy.S(2)/3) + b**(sympy.S(2)/3))*log(a**(sympy.S(1)/3)*cos(e + f*x) + b**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)*f*(a**2 - b**2)) + b**(sympy.S(2)/3)*(a**(sympy.S(2)/3) + b**(sympy.S(2)/3))*log(a**(sympy.S(2)/3)*cos(e + f*x)**2 - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*cos(e + f*x) + b**(sympy.S(2)/3))/(6*a**(sympy.S(1)/3)*f*(a**2 - b**2)) - sqrt(3)*b**(sympy.S(2)/3)*atan(sqrt(3)*(-2*a**(sympy.S(1)/3)*cos(e + f*x) + b**(sympy.S(1)/3))/(3*b**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3)*f*(a**(sympy.S(4)/3) + a**(sympy.S(2)/3)*b**(sympy.S(2)/3) + b**(sympy.S(4)/3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_462():
    f = cot(e + f*x)**3/(a + b*sec(e + f*x)**3)
    F = -1/(f*(4*a - 4*b)*(cos(e + f*x) + 1)) - (2*a + 5*b)*log(1 - cos(e + f*x))/(4*f*(a + b)**2) - (2*a - 5*b)*log(cos(e + f*x) + 1)/(4*f*(a - b)**2) - 1/(f*(1 - cos(e + f*x))*(4*a + 4*b)) - b**2*(2*a**2 + b**2)*log(a*cos(e + f*x)**3 + b)/(3*a*f*(a**2 - b**2)**2) + sqrt(3)*b**(sympy.S(4)/3)*(-3*a**(sympy.S(2)/3)*b**(sympy.S(4)/3) + a**2 + 2*b**2)*atan(sqrt(3)*(-2*a**(sympy.S(1)/3)*cos(e + f*x) + b**(sympy.S(1)/3))/(3*b**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3)*f*(a**2 - b**2)**2) - b**(sympy.S(4)/3)*(3*a**(sympy.S(2)/3)*b**(sympy.S(4)/3) + a**2 + 2*b**2)*log(a**(sympy.S(1)/3)*cos(e + f*x) + b**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)*f*(a**2 - b**2)**2) + b**(sympy.S(4)/3)*(3*a**(sympy.S(2)/3)*b**(sympy.S(4)/3) + a**2 + 2*b**2)*log(a**(sympy.S(2)/3)*cos(e + f*x)**2 - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*cos(e + f*x) + b**(sympy.S(2)/3))/(6*a**(sympy.S(1)/3)*f*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_463():
    f = (d*tan(e + f*x))**m*(a + b*(c*sec(e + f*x))**n)**p
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * ((Symbol('c') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))**(Symbol('p')) * ((Symbol('d') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_464():
    f = (a + b*(c*sec(e + f*x))**n)**p*tan(e + f*x)**5
    F = -(a + b*(c*sec(e + f*x))**n)**p*hyper((-p, 2/n), ((n + 2)/n,), -b*(c*sec(e + f*x))**n/a)*sec(e + f*x)**2/(f*(1 + b*(c*sec(e + f*x))**n/a)**p) + (a + b*(c*sec(e + f*x))**n)**p*hyper((-p, 4/n), ((n + 4)/n,), -b*(c*sec(e + f*x))**n/a)*sec(e + f*x)**4/(4*f*(1 + b*(c*sec(e + f*x))**n/a)**p) - (a + b*(c*sec(e + f*x))**n)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*(c*sec(e + f*x))**n/a)/(a*f*n*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_465():
    f = (a + b*(c*sec(e + f*x))**n)**p*tan(e + f*x)**3
    F = (a + b*(c*sec(e + f*x))**n)**p*hyper((-p, 2/n), ((n + 2)/n,), -b*(c*sec(e + f*x))**n/a)*sec(e + f*x)**2/(2*f*(1 + b*(c*sec(e + f*x))**n/a)**p) + (a + b*(c*sec(e + f*x))**n)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*(c*sec(e + f*x))**n/a)/(a*f*n*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_466():
    f = (a + b*(c*sec(e + f*x))**n)**p*tan(e + f*x)
    F = -(a + b*(c*sec(e + f*x))**n)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*(c*sec(e + f*x))**n/a)/(a*f*n*(p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_467():
    f = (a + b*(c*sec(e + f*x))**n)**p*cot(e + f*x)
    F = sympy.Function('Unintegrable')((sympy.cot((Symbol('e') + (Symbol('f') * x))) * ((Symbol('a') + (Symbol('b') * ((Symbol('c') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_468():
    f = (a + b*(c*sec(e + f*x))**n)**p*cot(e + f*x)**3
    F = sympy.Function('Unintegrable')(((sympy.cot((Symbol('e') + (Symbol('f') * x))))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * ((Symbol('c') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_469():
    f = (a + b*(c*sec(e + f*x))**n)**p*tan(e + f*x)**2
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * ((Symbol('c') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))**(Symbol('p')) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_470():
    f = (a + b*(c*sec(e + f*x))**n)**p
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * ((Symbol('c') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))**(Symbol('p')), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_7_d_trig_pow_m_a_plus_b_c_sec_pow_n_pow_p_471():
    f = (a + b*(c*sec(e + f*x))**n)**p*cot(e + f*x)**2
    F = sympy.Function('Unintegrable')(((sympy.cot((Symbol('e') + (Symbol('f') * x))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * ((Symbol('c') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')))))**(Symbol('p'))), x)
    assert integrate(f, x) == F

