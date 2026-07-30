"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.2 Cosine/4.2.7 (d trig)^m (a+b (c cos)^n)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, n = symbols('a b n')

def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_1():
    f = sin(x)**6/(-a*cos(x)**2 + a)
    F = 3*x/(8*a) - sin(x)**3*cos(x)/(4*a) - 3*sin(x)*cos(x)/(8*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_2():
    f = sin(x)**5/(-a*cos(x)**2 + a)
    F = cos(x)**3/(3*a) - cos(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_3():
    f = sin(x)**4/(-a*cos(x)**2 + a)
    F = x/(2*a) - sin(x)*cos(x)/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_4():
    f = sin(x)**3/(-a*cos(x)**2 + a)
    F = -cos(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_5():
    f = sin(x)**2/(-a*cos(x)**2 + a)
    F = x/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_6():
    f = sin(x)/(-a*cos(x)**2 + a)
    F = -atanh(cos(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_7():
    f = csc(x)/(-a*cos(x)**2 + a)
    F = -cot(x)*csc(x)/(2*a) - atanh(cos(x))/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_8():
    f = csc(x)**2/(-a*cos(x)**2 + a)
    F = -cot(x)**3/(3*a) - cot(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_9():
    f = csc(x)**3/(-a*cos(x)**2 + a)
    F = -cot(x)*csc(x)**3/(4*a) - 3*cot(x)*csc(x)/(8*a) - 3*atanh(cos(x))/(8*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_10():
    f = sin(x)**7/(a + b*cos(x)**2)
    F = cos(x)**5/(5*b) - (a + 3*b)*cos(x)**3/(3*b**2) + (a**2 + 3*a*b + 3*b**2)*cos(x)/b**3 - (a + b)**3*atan(sqrt(b)*cos(x)/sqrt(a))/(sqrt(a)*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_11():
    f = sin(x)**5/(a + b*cos(x)**2)
    F = -cos(x)**3/(3*b) + (a + 2*b)*cos(x)/b**2 - (a + b)**2*atan(sqrt(b)*cos(x)/sqrt(a))/(sqrt(a)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_12():
    f = sin(x)**3/(a + b*cos(x)**2)
    F = cos(x)/b - (a + b)*atan(sqrt(b)*cos(x)/sqrt(a))/(sqrt(a)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_13():
    f = sin(x)/(a + b*cos(x)**2)
    F = -atan(sqrt(b)*cos(x)/sqrt(a))/(sqrt(a)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_14():
    f = csc(x)/(a + b*cos(x)**2)
    F = -atanh(cos(x))/(a + b) - sqrt(b)*atan(sqrt(b)*cos(x)/sqrt(a))/(sqrt(a)*(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_15():
    f = csc(x)**3/(a + b*cos(x)**2)
    F = -cot(x)*csc(x)/(2*a + 2*b) - (a + 3*b)*atanh(cos(x))/(2*(a + b)**2) - b**(sympy.S(3)/2)*atan(sqrt(b)*cos(x)/sqrt(a))/(sqrt(a)*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_16():
    f = csc(x)**5/(a + b*cos(x)**2)
    F = -cot(x)*csc(x)**3/(4*a + 4*b) - (3*a + 7*b)*cot(x)*csc(x)/(8*(a + b)**2) - (3*a**2 + 10*a*b + 15*b**2)*atanh(cos(x))/(8*(a + b)**3) - b**(sympy.S(5)/2)*atan(sqrt(b)*cos(x)/sqrt(a))/(sqrt(a)*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_17():
    f = sin(x)**6/(a + b*cos(x)**2)
    F = sin(x)**3*cos(x)/(4*b) + (4*a + 7*b)*sin(x)*cos(x)/(8*b**2) - x*(8*a**2 + 20*a*b + 15*b**2)/(8*b**3) - (a + b)**(sympy.S(5)/2)*atan(sqrt(a + b)*cot(x)/sqrt(a))/(sqrt(a)*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_18():
    f = sin(x)**4/(a + b*cos(x)**2)
    F = sin(x)*cos(x)/(2*b) - x*(2*a + 3*b)/(2*b**2) - (a + b)**(sympy.S(3)/2)*atan(sqrt(a + b)*cot(x)/sqrt(a))/(sqrt(a)*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_19():
    f = sin(x)**2/(a + b*cos(x)**2)
    F = -x/b - sqrt(a + b)*atan(sqrt(a + b)*cot(x)/sqrt(a))/(sqrt(a)*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_20():
    f = 1/(a + b*cos(x)**2)
    F = -atan(sqrt(a + b)*cot(x)/sqrt(a))/(sqrt(a)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_21():
    f = csc(x)**2/(a + b*cos(x)**2)
    F = -cot(x)/(a + b) - b*atan(sqrt(a + b)*cot(x)/sqrt(a))/(sqrt(a)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_22():
    f = csc(x)**4/(a + b*cos(x)**2)
    F = -cot(x)**3/(3*a + 3*b) - (a + 2*b)*cot(x)/(a + b)**2 - b**2*atan(sqrt(a + b)*cot(x)/sqrt(a))/(sqrt(a)*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_23():
    f = csc(x)**6/(a + b*cos(x)**2)
    F = -cot(x)**5/(5*a + 5*b) - (2*a + 3*b)*cot(x)**3/(3*(a + b)**2) - (a**2 + 3*a*b + 3*b**2)*cot(x)/(a + b)**3 - b**3*atan(sqrt(a + b)*cot(x)/sqrt(a))/(sqrt(a)*(a + b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_24():
    f = sin(x)/(4 - 3*cos(x)**3)
    F = 6**(sympy.S(2)/3)*log(-3**(sympy.S(1)/3)*cos(x) + 2**(sympy.S(2)/3))/36 - 6**(sympy.S(2)/3)*log(3**(sympy.S(2)/3)*cos(x)**2 + 2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*cos(x) + 2*2**(sympy.S(1)/3))/72 - 2**(sympy.S(2)/3)*3**(sympy.S(1)/6)*atan(sqrt(3)*(6**(sympy.S(1)/3)*cos(x) + 1)/3)/12
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_25():
    f = 1/(1 - cos(x)**2)
    F = -cot(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_26():
    f = (1 - cos(x)**2)**(-2)
    F = -cot(x)**3/3 - cot(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_27():
    f = (1 - cos(x)**2)**(-3)
    F = -cot(x)**5/5 - 2*cot(x)**3/3 - cot(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_28():
    f = cos(x)**7/(a + b*cos(x)**2)
    F = -a**3*atanh(sqrt(b)*sin(x)/sqrt(a + b))/(b**(sympy.S(7)/2)*sqrt(a + b)) + sin(x)**5/(5*b) + (a - 2*b)*sin(x)**3/(3*b**2) + (a**2 - a*b + b**2)*sin(x)/b**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_29():
    f = cos(x)**5/(a + b*cos(x)**2)
    F = a**2*atanh(sqrt(b)*sin(x)/sqrt(a + b))/(b**(sympy.S(5)/2)*sqrt(a + b)) - sin(x)**3/(3*b) - (a - b)*sin(x)/b**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_30():
    f = cos(x)**3/(a + b*cos(x)**2)
    F = -a*atanh(sqrt(b)*sin(x)/sqrt(a + b))/(b**(sympy.S(3)/2)*sqrt(a + b)) + sin(x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_31():
    f = cos(x)/(a + b*cos(x)**2)
    F = atanh(sqrt(b)*sin(x)/sqrt(a + b))/(sqrt(b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_32():
    f = sec(x)/(a + b*cos(x)**2)
    F = -sqrt(b)*atanh(sqrt(b)*sin(x)/sqrt(a + b))/(a*sqrt(a + b)) + atanh(sin(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_33():
    f = sec(x)**3/(a + b*cos(x)**2)
    F = tan(x)*sec(x)/(2*a) + b**(sympy.S(3)/2)*atanh(sqrt(b)*sin(x)/sqrt(a + b))/(a**2*sqrt(a + b)) + (a - 2*b)*atanh(sin(x))/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_34():
    f = sec(x)**5/(a + b*cos(x)**2)
    F = tan(x)*sec(x)**3/(4*a) + (3*a - 4*b)*tan(x)*sec(x)/(8*a**2) - b**(sympy.S(5)/2)*atanh(sqrt(b)*sin(x)/sqrt(a + b))/(a**3*sqrt(a + b)) + (3*a**2 - 4*a*b + 8*b**2)*atanh(sin(x))/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_35():
    f = cos(x)**6/(a + b*cos(x)**2)
    F = a**(sympy.S(5)/2)*atan(sqrt(a + b)*cot(x)/sqrt(a))/(b**3*sqrt(a + b)) + sin(x)*cos(x)**3/(4*b) - (4*a - 3*b)*sin(x)*cos(x)/(8*b**2) + x*(8*a**2 - 4*a*b + 3*b**2)/(8*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_36():
    f = cos(x)**4/(a + b*cos(x)**2)
    F = -a**(sympy.S(3)/2)*atan(sqrt(a + b)*cot(x)/sqrt(a))/(b**2*sqrt(a + b)) + sin(x)*cos(x)/(2*b) - x*(2*a - b)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_37():
    f = cos(x)**2/(a + b*cos(x)**2)
    F = sqrt(a)*atan(sqrt(a + b)*cot(x)/sqrt(a))/(b*sqrt(a + b)) + x/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_38():
    f = 1/(a + b*cos(x)**2)
    F = -atan(sqrt(a + b)*cot(x)/sqrt(a))/(sqrt(a)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_39():
    f = sec(x)**2/(a + b*cos(x)**2)
    F = tan(x)/a + b*atan(sqrt(a + b)*cot(x)/sqrt(a))/(a**(sympy.S(3)/2)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_40():
    f = sec(x)**4/(a + b*cos(x)**2)
    F = tan(x)**3/(3*a) + (a - b)*tan(x)/a**2 - b**2*atan(sqrt(a + b)*cot(x)/sqrt(a))/(a**(sympy.S(5)/2)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_41():
    f = sec(x)**6/(a + b*cos(x)**2)
    F = tan(x)**5/(5*a) + (2*a - b)*tan(x)**3/(3*a**2) + (a**2 - a*b + b**2)*tan(x)/a**3 + b**3*atan(sqrt(a + b)*cot(x)/sqrt(a))/(a**(sympy.S(7)/2)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_42():
    f = (a + b*cos(x)**2)**(-2)
    F = -b*sin(x)*cos(x)/(2*a*(a + b)*(a + b*cos(x)**2)) - (2*a + b)*atan(sqrt(a + b)*cot(x)/sqrt(a))/(2*a**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_43():
    f = (a + b*cos(x)**2)**(-3)
    F = -b*sin(x)*cos(x)/(4*a*(a + b)*(a + b*cos(x)**2)**2) - 3*b*(2*a + b)*sin(x)*cos(x)/(8*a**2*(a + b)**2*(a + b*cos(x)**2)) - (8*a**2 + 8*a*b + 3*b**2)*atan(sqrt(a + b)*cot(x)/sqrt(a))/(8*a**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_44():
    f = 1/(cos(x)**2 + 1)
    F = sqrt(2)*x/2 - sqrt(2)*atan(sin(x)*cos(x)/(cos(x)**2 + 1 + sqrt(2)))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_45():
    f = (cos(x)**2 + 1)**(-2)
    F = 3*sqrt(2)*x/8 - 3*sqrt(2)*atan(sin(x)*cos(x)/(cos(x)**2 + 1 + sqrt(2)))/8 - sin(x)*cos(x)/(4*cos(x)**2 + 4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_46():
    f = (cos(x)**2 + 1)**(-3)
    F = 19*sqrt(2)*x/64 - 19*sqrt(2)*atan(sin(x)*cos(x)/(cos(x)**2 + 1 + sqrt(2)))/64 - 9*sin(x)*cos(x)/(32*cos(x)**2 + 32) - sin(x)*cos(x)/(8*(cos(x)**2 + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_47():
    f = sqrt(1 - cos(x)**2)
    F = -sqrt(sin(x)**2)*cot(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_48():
    f = sqrt(cos(x)**2 - 1)
    F = -sqrt(-sin(x)**2)*cot(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_49():
    f = (1 - cos(x)**2)**(sympy.S(3)/2)
    F = -(sin(x)**2)**(sympy.S(3)/2)*cot(x)/3 - 2*sqrt(sin(x)**2)*cot(x)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_50():
    f = (cos(x)**2 - 1)**(sympy.S(3)/2)
    F = -(-sin(x)**2)**(sympy.S(3)/2)*cot(x)/3 + 2*sqrt(-sin(x)**2)*cot(x)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_51():
    f = 1/sqrt(1 - cos(x)**2)
    F = -sin(x)*atanh(cos(x))/sqrt(sin(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_52():
    f = 1/sqrt(cos(x)**2 - 1)
    F = -sin(x)*atanh(cos(x))/sqrt(-sin(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_53():
    f = (1 - cos(x)**2)**(sympy.S(-3)/2)
    F = -sin(x)*atanh(cos(x))/(2*sqrt(sin(x)**2)) - cot(x)/(2*sqrt(sin(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_54():
    f = (cos(x)**2 - 1)**(sympy.S(-3)/2)
    F = sin(x)*atanh(cos(x))/(2*sqrt(-sin(x)**2)) + cot(x)/(2*sqrt(-sin(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_55():
    f = sqrt(cos(x)**2 + 1)
    F = elliptic_e(x + pi/2, -1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_56():
    f = sqrt(-cos(x)**2 - 1)
    F = sqrt(-cos(x)**2 - 1)*elliptic_e(x + pi/2, -1)/sqrt(cos(x)**2 + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_57():
    f = sqrt(a + b*cos(x)**2)
    F = sqrt(a + b*cos(x)**2)*elliptic_e(x + pi/2, -b/a)/sqrt(1 + b*cos(x)**2/a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_58():
    f = (cos(x)**2 + 1)**(sympy.S(3)/2)
    F = sqrt(cos(x)**2 + 1)*sin(x)*cos(x)/3 + 2*elliptic_e(x + pi/2, -1) - 2*elliptic_f(x + pi/2, -1)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_59():
    f = (-cos(x)**2 - 1)**(sympy.S(3)/2)
    F = -sqrt(-cos(x)**2 - 1)*sin(x)*cos(x)/3 - 2*sqrt(-cos(x)**2 - 1)*elliptic_e(x + pi/2, -1)/sqrt(cos(x)**2 + 1) - 2*sqrt(cos(x)**2 + 1)*elliptic_f(x + pi/2, -1)/(3*sqrt(-cos(x)**2 - 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_60():
    f = (a + b*cos(x)**2)**(sympy.S(3)/2)
    F = -a*sqrt(1 + b*cos(x)**2/a)*(a + b)*elliptic_f(x + pi/2, -b/a)/(3*sqrt(a + b*cos(x)**2)) + b*sqrt(a + b*cos(x)**2)*sin(x)*cos(x)/3 + sqrt(a + b*cos(x)**2)*(4*a + 2*b)*elliptic_e(x + pi/2, -b/a)/(3*sqrt(1 + b*cos(x)**2/a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_61():
    f = 1/sqrt(cos(x)**2 + 1)
    F = elliptic_f(x + pi/2, -1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_62():
    f = 1/sqrt(-cos(x)**2 - 1)
    F = sqrt(cos(x)**2 + 1)*elliptic_f(x + pi/2, -1)/sqrt(-cos(x)**2 - 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_63():
    f = 1/sqrt(a + b*cos(x)**2)
    F = sqrt(1 + b*cos(x)**2/a)*elliptic_f(x + pi/2, -b/a)/sqrt(a + b*cos(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_64():
    f = (cos(x)**2 + 1)**(sympy.S(-3)/2)
    F = elliptic_e(x + pi/2, -1)/2 - sin(x)*cos(x)/(2*sqrt(cos(x)**2 + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_65():
    f = (-cos(x)**2 - 1)**(sympy.S(-3)/2)
    F = sqrt(-cos(x)**2 - 1)*elliptic_e(x + pi/2, -1)/(2*sqrt(cos(x)**2 + 1)) + sin(x)*cos(x)/(2*sqrt(-cos(x)**2 - 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_66():
    f = (a + b*cos(x)**2)**(sympy.S(-3)/2)
    F = -b*sin(x)*cos(x)/(a*(a + b)*sqrt(a + b*cos(x)**2)) + sqrt(a + b*cos(x)**2)*elliptic_e(x + pi/2, -b/a)/(a*sqrt(1 + b*cos(x)**2/a)*(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_67():
    f = cos(x)/sqrt(cos(x)**2 + 1)
    F = asin(sqrt(2)*sin(x)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_68():
    f = cos(3*x + 5)/sqrt(cos(3*x + 5)**2 + 3)
    F = asin(sin(3*x + 5)/2)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_69():
    f = cos(x)/sqrt(4 - cos(x)**2)
    F = asinh(sqrt(3)*sin(x)/3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_70():
    f = 1/(a + b*cos(x)**4)
    F = -sqrt(2)*(sqrt(a) - sqrt(a + b))*log(-sqrt(2)*a**(sympy.S(1)/4)*sqrt(-sqrt(a)*sqrt(a + b) + a + b)*cot(x) + sqrt(a)*(a + b)**(sympy.S(1)/4) + (a + b)**(sympy.S(3)/4)*cot(x)**2)/(8*a**(sympy.S(3)/4)*(a + b)**(sympy.S(1)/4)*sqrt(-sqrt(a)*sqrt(a + b) + a + b)) + sqrt(2)*(sqrt(a) - sqrt(a + b))*log(sqrt(2)*a**(sympy.S(1)/4)*sqrt(-sqrt(a)*sqrt(a + b) + a + b)*cot(x) + sqrt(a)*(a + b)**(sympy.S(1)/4) + (a + b)**(sympy.S(3)/4)*cot(x)**2)/(8*a**(sympy.S(3)/4)*(a + b)**(sympy.S(1)/4)*sqrt(-sqrt(a)*sqrt(a + b) + a + b)) + sqrt(2)*(sqrt(a) + sqrt(a + b))*atan((a**(sympy.S(1)/4)*sqrt(-sqrt(a)*sqrt(a + b) + a + b) - sqrt(2)*(a + b)**(sympy.S(3)/4)*cot(x))/(a**(sympy.S(1)/4)*sqrt(sqrt(a)*sqrt(a + b) + a + b)))/(4*a**(sympy.S(3)/4)*(a + b)**(sympy.S(1)/4)*sqrt(sqrt(a)*sqrt(a + b) + a + b)) - sqrt(2)*(sqrt(a) + sqrt(a + b))*atan((a**(sympy.S(1)/4)*sqrt(-sqrt(a)*sqrt(a + b) + a + b) + sqrt(2)*(a + b)**(sympy.S(3)/4)*cot(x))/(a**(sympy.S(1)/4)*sqrt(sqrt(a)*sqrt(a + b) + a + b)))/(4*a**(sympy.S(3)/4)*(a + b)**(sympy.S(1)/4)*sqrt(sqrt(a)*sqrt(a + b) + a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_71():
    f = 1/(a - b*cos(x)**4)
    F = -atan(sqrt(sqrt(a) + sqrt(b))*cot(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*sqrt(sqrt(a) + sqrt(b))) - atan(sqrt(sqrt(a) - sqrt(b))*cot(x)/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*sqrt(sqrt(a) - sqrt(b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_72():
    f = 1/(cos(x)**4 + 1)
    F = x/(2*sqrt(-1 + sqrt(2))) - sqrt(-1 + sqrt(2))*log(sqrt(2)*cot(x)**2 + sqrt(-2 + 2*sqrt(2))*cot(x) + 1)/8 + sqrt(-1 + sqrt(2))*log(2*cot(x)**2 - 2*sqrt(-1 + sqrt(2))*cot(x) + sqrt(2))/8 + atan((sqrt(-1 + sqrt(2))*(1 - 2*sin(x)**2) + (-2 + sqrt(2))*sin(x)*cos(x))/((-2 + sqrt(2))*sin(x)**2 + 2*sqrt(-1 + sqrt(2))*sin(x)*cos(x) + sqrt(1 + sqrt(2)) + 2))/(4*sqrt(-1 + sqrt(2))) + atan((sqrt(-1 + sqrt(2))*(2*sin(x)**2 - 1) + (-2 + sqrt(2))*sin(x)*cos(x))/((-2 + sqrt(2))*sin(x)**2 - 2*sqrt(-1 + sqrt(2))*sin(x)*cos(x) + sqrt(1 + sqrt(2)) + 2))/(4*sqrt(-1 + sqrt(2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_73():
    f = 1/(1 - cos(x)**4)
    F = sqrt(2)*x/4 - cot(x)/2 - sqrt(2)*atan(sin(x)*cos(x)/(cos(x)**2 + 1 + sqrt(2)))/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_74():
    f = 1/(a + b*cos(x)**5)
    F = 2*atan(sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(4)/5)*b**(sympy.S(1)/5))*tan(x/2)/sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(4)/5)*b**(sympy.S(1)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(4)/5)*b**(sympy.S(1)/5))*sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(4)/5)*b**(sympy.S(1)/5))) + 2*atan(sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(3)/5)*b**(sympy.S(1)/5))*tan(x/2)/sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(3)/5)*b**(sympy.S(1)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(3)/5)*b**(sympy.S(1)/5))*sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(3)/5)*b**(sympy.S(1)/5))) + 2*atan(sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(2)/5)*b**(sympy.S(1)/5))*tan(x/2)/sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(2)/5)*b**(sympy.S(1)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(2)/5)*b**(sympy.S(1)/5))*sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(2)/5)*b**(sympy.S(1)/5))) + 2*atan(sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(1)/5)*b**(sympy.S(1)/5))*tan(x/2)/sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(1)/5)*b**(sympy.S(1)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(1)/5)*b**(sympy.S(1)/5))*sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(1)/5)*b**(sympy.S(1)/5))) + 2*atan(sqrt(a**(sympy.S(1)/5) - b**(sympy.S(1)/5))*tan(x/2)/sqrt(a**(sympy.S(1)/5) + b**(sympy.S(1)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(1)/5) - b**(sympy.S(1)/5))*sqrt(a**(sympy.S(1)/5) + b**(sympy.S(1)/5)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_75():
    f = 1/(a + b*cos(x)**6)
    F = -atan(sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))*cot(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))) - atan(sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))*cot(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))) - atan(sqrt(a**(sympy.S(1)/3) + b**(sympy.S(1)/3))*cot(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_76():
    f = 1/(a + b*cos(x)**8)
    F = atan(sqrt(I*b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))*cot(x)/(-a)**(sympy.S(1)/8))/(4*(-a)**(sympy.S(7)/8)*sqrt(I*b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))) + atan(sqrt(-I*b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))*cot(x)/(-a)**(sympy.S(1)/8))/(4*(-a)**(sympy.S(7)/8)*sqrt(-I*b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))) + atan(sqrt(b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))*cot(x)/(-a)**(sympy.S(1)/8))/(4*(-a)**(sympy.S(7)/8)*sqrt(b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))) + atan(sqrt(-b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))*cot(x)/(-a)**(sympy.S(1)/8))/(4*(-a)**(sympy.S(7)/8)*sqrt(-b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_77():
    f = 1/(a - b*cos(x)**5)
    F = 2*atan(sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(4)/5)*b**(sympy.S(1)/5))*tan(x/2)/sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(4)/5)*b**(sympy.S(1)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(4)/5)*b**(sympy.S(1)/5))*sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(4)/5)*b**(sympy.S(1)/5))) + 2*atan(sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(3)/5)*b**(sympy.S(1)/5))*tan(x/2)/sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(3)/5)*b**(sympy.S(1)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(3)/5)*b**(sympy.S(1)/5))*sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(3)/5)*b**(sympy.S(1)/5))) + 2*atan(sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(2)/5)*b**(sympy.S(1)/5))*tan(x/2)/sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(2)/5)*b**(sympy.S(1)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(2)/5)*b**(sympy.S(1)/5))*sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(2)/5)*b**(sympy.S(1)/5))) + 2*atan(sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(1)/5)*b**(sympy.S(1)/5))*tan(x/2)/sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(1)/5)*b**(sympy.S(1)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(1)/5) - (-1)**(sympy.S(1)/5)*b**(sympy.S(1)/5))*sqrt(a**(sympy.S(1)/5) + (-1)**(sympy.S(1)/5)*b**(sympy.S(1)/5))) + 2*atan(sqrt(a**(sympy.S(1)/5) + b**(sympy.S(1)/5))*tan(x/2)/sqrt(a**(sympy.S(1)/5) - b**(sympy.S(1)/5)))/(5*a**(sympy.S(4)/5)*sqrt(a**(sympy.S(1)/5) - b**(sympy.S(1)/5))*sqrt(a**(sympy.S(1)/5) + b**(sympy.S(1)/5)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_78():
    f = 1/(a - b*cos(x)**6)
    F = -atan(sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))*cot(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) - (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3))) - atan(sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))*cot(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3))) - atan(sqrt(a**(sympy.S(1)/3) - b**(sympy.S(1)/3))*cot(x)/a**(sympy.S(1)/6))/(3*a**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_79():
    f = 1/(a - b*cos(x)**8)
    F = -atan(sqrt(a**(sympy.S(1)/4) + I*b**(sympy.S(1)/4))*cot(x)/a**(sympy.S(1)/8))/(4*a**(sympy.S(7)/8)*sqrt(a**(sympy.S(1)/4) + I*b**(sympy.S(1)/4))) - atan(sqrt(a**(sympy.S(1)/4) - I*b**(sympy.S(1)/4))*cot(x)/a**(sympy.S(1)/8))/(4*a**(sympy.S(7)/8)*sqrt(a**(sympy.S(1)/4) - I*b**(sympy.S(1)/4))) - atan(sqrt(a**(sympy.S(1)/4) + b**(sympy.S(1)/4))*cot(x)/a**(sympy.S(1)/8))/(4*a**(sympy.S(7)/8)*sqrt(a**(sympy.S(1)/4) + b**(sympy.S(1)/4))) - atan(sqrt(a**(sympy.S(1)/4) - b**(sympy.S(1)/4))*cot(x)/a**(sympy.S(1)/8))/(4*a**(sympy.S(7)/8)*sqrt(a**(sympy.S(1)/4) - b**(sympy.S(1)/4)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_80():
    f = 1/(cos(x)**5 + 1)
    F = 2*atan(sqrt((1 - (-1)**(sympy.S(2)/5))/(1 + (-1)**(sympy.S(2)/5)))*tan(x/2))/(5*sqrt(1 - (-1)**(sympy.S(4)/5))) + 2*atan(sqrt((1 - (-1)**(sympy.S(4)/5))/(1 + (-1)**(sympy.S(4)/5)))*tan(x/2))/(5*sqrt(1 + (-1)**(sympy.S(3)/5))) - 2*atanh(tan(x/2)/sqrt(-(1 - (-1)**(sympy.S(1)/5))/(1 + (-1)**(sympy.S(1)/5))))/(5*sqrt(-1 + (-1)**(sympy.S(2)/5))) - 2*sqrt(-(1 + (-1)**(sympy.S(3)/5))/(1 - (-1)**(sympy.S(3)/5)))*atanh(sqrt(-(1 + (-1)**(sympy.S(3)/5))/(1 - (-1)**(sympy.S(3)/5)))*tan(x/2))/(5 + 5*(-1)**(sympy.S(3)/5)) + sin(x)/(5*cos(x) + 5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_81():
    f = 1/(cos(x)**8 + 1)
    F = -atan(sqrt(1 - (-1)**(sympy.S(1)/4))*cot(x))/(4*sqrt(1 - (-1)**(sympy.S(1)/4))) - atan(sqrt(1 + (-1)**(sympy.S(1)/4))*cot(x))/(4*sqrt(1 + (-1)**(sympy.S(1)/4))) - atan(sqrt(1 - (-1)**(sympy.S(3)/4))*cot(x))/(4*sqrt(1 - (-1)**(sympy.S(3)/4))) - atan(sqrt(1 + (-1)**(sympy.S(3)/4))*cot(x))/(4*sqrt(1 + (-1)**(sympy.S(3)/4)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_82():
    f = 1/(1 - cos(x)**5)
    F = 2*atan(sqrt((1 - (-1)**(sympy.S(1)/5))/(1 + (-1)**(sympy.S(1)/5)))*tan(x/2))/(5*sqrt(1 - (-1)**(sympy.S(2)/5))) + 2*atan(sqrt((1 - (-1)**(sympy.S(3)/5))/(1 + (-1)**(sympy.S(3)/5)))*tan(x/2))/(5*sqrt(1 + (-1)**(sympy.S(1)/5))) - 2*atanh(tan(x/2)/sqrt(-(1 - (-1)**(sympy.S(2)/5))/(1 + (-1)**(sympy.S(2)/5))))/(5*sqrt(-1 + (-1)**(sympy.S(4)/5))) + 2*atanh(sqrt(-(1 + (-1)**(sympy.S(4)/5))/(1 - (-1)**(sympy.S(4)/5)))*tan(x/2))/(5*sqrt(-1 - (-1)**(sympy.S(3)/5))) - sin(x)/(5 - 5*cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_83():
    f = 1/(1 - cos(x)**6)
    F = -cot(x)/3 - atan(sqrt(1 + (-1)**(sympy.S(1)/3))*cot(x))/(3*sqrt(1 + (-1)**(sympy.S(1)/3))) - atan(sqrt(1 - (-1)**(sympy.S(2)/3))*cot(x))/(3*sqrt(1 - (-1)**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_84():
    f = 1/(1 - cos(x)**8)
    F = sqrt(2)*x/8 - cot(x)/4 - atan(sqrt(1 - I)*cot(x))/(4*sqrt(1 - I)) - atan(sqrt(1 + I)*cot(x))/(4*sqrt(1 + I)) - sqrt(2)*atan(sin(x)*cos(x)/(cos(x)**2 + 1 + sqrt(2)))/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_85():
    f = tan(x)/(cos(x)**2 + 1)
    F = log(cos(x)**2 + 1)/2 - log(cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_86():
    f = sqrt(a + b*cos(x)**2)*tan(x)
    F = sqrt(a)*atanh(sqrt(a + b*cos(x)**2)/sqrt(a)) - sqrt(a + b*cos(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_87():
    f = sqrt(1 - cos(x)**2)*tan(x)
    F = -sqrt(sin(x)**2) + atanh(sqrt(sin(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_88():
    f = tan(x)/sqrt(a + b*cos(x)**2)
    F = atanh(sqrt(a + b*cos(x)**2)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_89():
    f = tan(x)/sqrt(cos(x)**2 + 1)
    F = atanh(sqrt(cos(x)**2 + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_90():
    f = tan(x)/sqrt(1 - cos(x)**2)
    F = atanh(sqrt(sin(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_91():
    f = tan(x)**3/(a + b*cos(x)**3)
    F = -log(a + b*cos(x)**3)/(3*a) + log(cos(x))/a + sec(x)**2/(2*a) + b**(sympy.S(2)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*cos(x))/(3*a**(sympy.S(5)/3)) - b**(sympy.S(2)/3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*cos(x) + b**(sympy.S(2)/3)*cos(x)**2)/(6*a**(sympy.S(5)/3)) - sqrt(3)*b**(sympy.S(2)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*cos(x))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_92():
    f = sqrt(a + b*cos(x)**3)*tan(x)
    F = 2*sqrt(a)*atanh(sqrt(a + b*cos(x)**3)/sqrt(a))/3 - 2*sqrt(a + b*cos(x)**3)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_93():
    f = tan(x)/sqrt(a + b*cos(x)**3)
    F = 2*atanh(sqrt(a + b*cos(x)**3)/sqrt(a))/(3*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_94():
    f = sqrt(a + b*cos(x)**4)*tan(x)
    F = sqrt(a)*atanh(sqrt(a + b*cos(x)**4)/sqrt(a))/2 - sqrt(a + b*cos(x)**4)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_95():
    f = tan(x)/sqrt(a + b*cos(x)**4)
    F = atanh(sqrt(a + b*cos(x)**4)/sqrt(a))/(2*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_96():
    f = sqrt(a + b*cos(x)**n)*tan(x)
    F = 2*sqrt(a)*atanh(sqrt(a + b*cos(x)**n)/sqrt(a))/n - 2*sqrt(a + b*cos(x)**n)/n
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_7_d_trig_pow_m_a_plus_b_c_cos_pow_n_pow_p_97():
    f = tan(x)/sqrt(a + b*cos(x)**n)
    F = 2*atanh(sqrt(a + b*cos(x)**n)/sqrt(a))/(sqrt(a)*n)
    assert integrate(f, x) == F

