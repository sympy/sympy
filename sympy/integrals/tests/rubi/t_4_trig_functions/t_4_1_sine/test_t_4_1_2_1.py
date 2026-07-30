"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.1 Sine/4.1.2.1 (a+b sin)^m (c+d sin)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, a, b, c, d, e, f, m, n, p = symbols('A B a b c d e f m n p')

def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_1():
    f = (a*sin(e + f*x) + a)**2*sin(e + f*x)**3
    F = 3*a**2*x/4 - a**2*sin(e + f*x)**3*cos(e + f*x)/(2*f) - 3*a**2*sin(e + f*x)*cos(e + f*x)/(4*f) - a**2*cos(e + f*x)**5/(5*f) + a**2*cos(e + f*x)**3/f - 2*a**2*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_2():
    f = (a*sin(e + f*x) + a)**3*sin(e + f*x)**3
    F = 23*a**3*x/16 - a**3*sin(e + f*x)**5*cos(e + f*x)/(6*f) - 23*a**3*sin(e + f*x)**3*cos(e + f*x)/(24*f) - 23*a**3*sin(e + f*x)*cos(e + f*x)/(16*f) - 3*a**3*cos(e + f*x)**5/(5*f) + 7*a**3*cos(e + f*x)**3/(3*f) - 4*a**3*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_3():
    f = sin(x)**4/(a*sin(x) + a)
    F = sin(x)**3*cos(x)/(a*sin(x) + a) - 3*x/(2*a) + 3*sin(x)*cos(x)/(2*a) + 4*cos(x)**3/(3*a) - 4*cos(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_4():
    f = sin(x)**3/(a*sin(x) + a)
    F = sin(x)**2*cos(x)/(a*sin(x) + a) + 3*x/(2*a) - 3*sin(x)*cos(x)/(2*a) + 2*cos(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_5():
    f = sin(x)**2/(a*sin(x) + a)
    F = -x/a - cos(x)/a - cos(x)/(a*(sin(x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_6():
    f = sin(x)/(a*sin(x) + a)
    F = cos(x)/(a*sin(x) + a) + x/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_7():
    f = 1/(a*sin(x) + a)
    F = -cos(x)/(a*sin(x) + a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_8():
    f = csc(x)/(a*sin(x) + a)
    F = cos(x)/(a*sin(x) + a) - atanh(cos(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_9():
    f = csc(x)**2/(a*sin(x) + a)
    F = cot(x)/(a*sin(x) + a) - 2*cot(x)/a + atanh(cos(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_10():
    f = csc(x)**3/(a*sin(x) + a)
    F = cot(x)*csc(x)/(a*sin(x) + a) - 3*cot(x)*csc(x)/(2*a) + 2*cot(x)/a - 3*atanh(cos(x))/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_11():
    f = csc(x)**4/(a*sin(x) + a)
    F = cot(x)*csc(x)**2/(a*sin(x) + a) - 4*cot(x)**3/(3*a) + 3*cot(x)*csc(x)/(2*a) - 4*cot(x)/a + 3*atanh(cos(x))/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_12():
    f = sin(x)**4/(a*sin(x) + a)**2
    F = sin(x)**3*cos(x)/(3*(a*sin(x) + a)**2) + 7*x/(2*a**2) - 7*sin(x)*cos(x)/(2*a**2) + 16*cos(x)/(3*a**2) + 8*sin(x)**2*cos(x)/(3*a**2*(sin(x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_13():
    f = sin(x)**3/(a*sin(x) + a)**2
    F = sin(x)**2*cos(x)/(3*(a*sin(x) + a)**2) - 2*x/a**2 - 4*cos(x)/(3*a**2) - 2*cos(x)/(a**2*(sin(x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_14():
    f = sin(x)**2/(a*sin(x) + a)**2
    F = -cos(x)/(3*(a*sin(x) + a)**2) + x/a**2 + 5*cos(x)/(3*a**2*(sin(x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_15():
    f = sin(x)/(a*sin(x) + a)**2
    F = -2*cos(x)/(3*a**2*sin(x) + 3*a**2) + cos(x)/(3*(a*sin(x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_16():
    f = (a*sin(x) + a)**(-2)
    F = -cos(x)/(3*a**2*sin(x) + 3*a**2) - cos(x)/(3*(a*sin(x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_17():
    f = csc(x)/(a*sin(x) + a)**2
    F = cos(x)/(3*(a*sin(x) + a)**2) - atanh(cos(x))/a**2 + 4*cos(x)/(3*a**2*(sin(x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_18():
    f = csc(x)**2/(a*sin(x) + a)**2
    F = cot(x)/(3*(a*sin(x) + a)**2) - 10*cot(x)/(3*a**2) + 2*atanh(cos(x))/a**2 + 2*cot(x)/(a**2*(sin(x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_19():
    f = csc(x)**3/(a*sin(x) + a)**2
    F = cot(x)*csc(x)/(3*(a*sin(x) + a)**2) - 7*cot(x)*csc(x)/(2*a**2) + 16*cot(x)/(3*a**2) - 7*atanh(cos(x))/(2*a**2) + 8*cot(x)*csc(x)/(3*a**2*(sin(x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_20():
    f = sin(x)**6/(a*sin(x) + a)**3
    F = 23*sin(x)**3*cos(x)/(3*a**3*sin(x) + 3*a**3) + sin(x)**5*cos(x)/(5*(a*sin(x) + a)**3) + 13*sin(x)**4*cos(x)/(15*a*(a*sin(x) + a)**2) - 23*x/(2*a**3) + 23*sin(x)*cos(x)/(2*a**3) + 136*cos(x)**3/(15*a**3) - 136*cos(x)/(5*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_21():
    f = sin(x)**5/(a*sin(x) + a)**3
    F = 76*sin(x)**2*cos(x)/(15*a**3*sin(x) + 15*a**3) + sin(x)**4*cos(x)/(5*(a*sin(x) + a)**3) + 11*sin(x)**3*cos(x)/(15*a*(a*sin(x) + a)**2) + 13*x/(2*a**3) - 13*sin(x)*cos(x)/(2*a**3) + 152*cos(x)/(15*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_22():
    f = sin(x)**4/(a*sin(x) + a)**3
    F = -3*cos(x)/(a**3*sin(x) + a**3) + sin(x)**3*cos(x)/(5*(a*sin(x) + a)**3) + 3*sin(x)**2*cos(x)/(5*a*(a*sin(x) + a)**2) - 3*x/a**3 - 9*cos(x)/(5*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_23():
    f = sin(x)**3/(a*sin(x) + a)**3
    F = 29*cos(x)/(15*a**3*sin(x) + 15*a**3) + sin(x)**2*cos(x)/(5*(a*sin(x) + a)**3) - 7*cos(x)/(15*a*(a*sin(x) + a)**2) + x/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_24():
    f = sin(x)**2/(a*sin(x) + a)**3
    F = -7*cos(x)/(15*a**3*sin(x) + 15*a**3) - cos(x)/(5*(a*sin(x) + a)**3) + 8*cos(x)/(15*a*(a*sin(x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_25():
    f = sin(x)/(a*sin(x) + a)**3
    F = -cos(x)/(5*a**3*sin(x) + 5*a**3) + cos(x)/(5*(a*sin(x) + a)**3) - cos(x)/(5*a*(a*sin(x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_26():
    f = (a*sin(x) + a)**(-3)
    F = -2*cos(x)/(15*a**3*sin(x) + 15*a**3) - cos(x)/(5*(a*sin(x) + a)**3) - 2*cos(x)/(15*a*(a*sin(x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_27():
    f = csc(x)/(a*sin(x) + a)**3
    F = 22*cos(x)/(15*a**3*sin(x) + 15*a**3) + cos(x)/(5*(a*sin(x) + a)**3) + 7*cos(x)/(15*a*(a*sin(x) + a)**2) - atanh(cos(x))/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_28():
    f = csc(x)**2/(a*sin(x) + a)**3
    F = 3*cot(x)/(a**3*sin(x) + a**3) + cot(x)/(5*(a*sin(x) + a)**3) + 3*cot(x)/(5*a*(a*sin(x) + a)**2) - 24*cot(x)/(5*a**3) + 3*atanh(cos(x))/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_29():
    f = csc(x)**3/(a*sin(x) + a)**3
    F = 76*cot(x)*csc(x)/(15*a**3*sin(x) + 15*a**3) + cot(x)*csc(x)/(5*(a*sin(x) + a)**3) + 11*cot(x)*csc(x)/(15*a*(a*sin(x) + a)**2) - 13*cot(x)*csc(x)/(2*a**3) + 152*cot(x)/(15*a**3) - 13*atanh(cos(x))/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_30():
    f = csc(x)**4/(a*sin(x) + a)**3
    F = 23*cot(x)*csc(x)**2/(3*a**3*sin(x) + 3*a**3) + cot(x)*csc(x)**2/(5*(a*sin(x) + a)**3) + 13*cot(x)*csc(x)**2/(15*a*(a*sin(x) + a)**2) - 136*cot(x)**3/(15*a**3) + 23*cot(x)*csc(x)/(2*a**3) - 136*cot(x)/(5*a**3) + 23*atanh(cos(x))/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_31():
    f = sqrt(a*sin(c + d*x) + a)*sin(c + d*x)**4
    F = -2*a*sin(c + d*x)**4*cos(c + d*x)/(9*d*sqrt(a*sin(c + d*x) + a)) - 16*a*sin(c + d*x)**3*cos(c + d*x)/(63*d*sqrt(a*sin(c + d*x) + a)) - 32*a*cos(c + d*x)/(45*d*sqrt(a*sin(c + d*x) + a)) + 64*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(315*d) - 32*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)/(105*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_32():
    f = sqrt(a*sin(c + d*x) + a)*sin(c + d*x)**3
    F = -2*a*sin(c + d*x)**3*cos(c + d*x)/(7*d*sqrt(a*sin(c + d*x) + a)) - 4*a*cos(c + d*x)/(5*d*sqrt(a*sin(c + d*x) + a)) + 8*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(35*d) - 12*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)/(35*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_33():
    f = sqrt(a*sin(c + d*x) + a)*sin(c + d*x)**2
    F = -14*a*cos(c + d*x)/(15*d*sqrt(a*sin(c + d*x) + a)) + 4*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(15*d) - 2*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_34():
    f = sqrt(a*sin(c + d*x) + a)*sin(c + d*x)
    F = -2*a*cos(c + d*x)/(3*d*sqrt(a*sin(c + d*x) + a)) - 2*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_35():
    f = sqrt(a*sin(c + d*x) + a)
    F = -2*a*cos(c + d*x)/(d*sqrt(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_36():
    f = sqrt(a*sin(c + d*x) + a)*csc(c + d*x)
    F = -2*sqrt(a)*atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_37():
    f = sqrt(a*sin(c + d*x) + a)*csc(c + d*x)**2
    F = -sqrt(a)*atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/d - a*cot(c + d*x)/(d*sqrt(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_38():
    f = sqrt(a*sin(c + d*x) + a)*csc(c + d*x)**3
    F = -3*sqrt(a)*atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/(4*d) - a*cot(c + d*x)*csc(c + d*x)/(2*d*sqrt(a*sin(c + d*x) + a)) - 3*a*cot(c + d*x)/(4*d*sqrt(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_39():
    f = sqrt(a*sin(c + d*x) + a)*csc(c + d*x)**4
    F = -5*sqrt(a)*atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/(8*d) - a*cot(c + d*x)*csc(c + d*x)**2/(3*d*sqrt(a*sin(c + d*x) + a)) - 5*a*cot(c + d*x)*csc(c + d*x)/(12*d*sqrt(a*sin(c + d*x) + a)) - 5*a*cot(c + d*x)/(8*d*sqrt(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_40():
    f = sqrt(-a*sin(c + d*x) + a)*csc(c + d*x)
    F = -2*sqrt(a)*atanh(sqrt(a)*cos(c + d*x)/sqrt(-a*sin(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_41():
    f = sqrt(a*sin(c + d*x) - a)*csc(c + d*x)
    F = 2*sqrt(a)*atan(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) - a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_42():
    f = sqrt(-a*sin(c + d*x) - a)*csc(c + d*x)
    F = 2*sqrt(a)*atan(sqrt(a)*cos(c + d*x)/sqrt(-a*sin(c + d*x) - a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_43():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)**3
    F = -2*a**2*sin(c + d*x)**4*cos(c + d*x)/(9*d*sqrt(a*sin(c + d*x) + a)) - 34*a**2*sin(c + d*x)**3*cos(c + d*x)/(63*d*sqrt(a*sin(c + d*x) + a)) - 68*a**2*cos(c + d*x)/(45*d*sqrt(a*sin(c + d*x) + a)) + 136*a*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(315*d) - 68*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)/(105*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_44():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)**2
    F = -152*a**2*cos(c + d*x)/(105*d*sqrt(a*sin(c + d*x) + a)) - 38*a*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(105*d) + 4*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)/(35*d) - 2*(a*sin(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)/(7*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_45():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)
    F = -8*a**2*cos(c + d*x)/(5*d*sqrt(a*sin(c + d*x) + a)) - 2*a*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(5*d) - 2*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_46():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = -8*a**2*cos(c + d*x)/(3*d*sqrt(a*sin(c + d*x) + a)) - 2*a*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_47():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)*csc(c + d*x)
    F = -2*a**(sympy.S(3)/2)*atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/d - 2*a**2*cos(c + d*x)/(d*sqrt(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_48():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)*csc(c + d*x)**2
    F = -3*a**(sympy.S(3)/2)*atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/d - a**2*cot(c + d*x)/(d*sqrt(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_49():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)*csc(c + d*x)**3
    F = -7*a**(sympy.S(3)/2)*atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/(4*d) - a**2*cot(c + d*x)*csc(c + d*x)/(2*d*sqrt(a*sin(c + d*x) + a)) - 7*a**2*cot(c + d*x)/(4*d*sqrt(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_50():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)*csc(c + d*x)**4
    F = -11*a**(sympy.S(3)/2)*atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/(8*d) - a**2*cot(c + d*x)*csc(c + d*x)**2/(3*d*sqrt(a*sin(c + d*x) + a)) - 11*a**2*cot(c + d*x)*csc(c + d*x)/(12*d*sqrt(a*sin(c + d*x) + a)) - 11*a**2*cot(c + d*x)/(8*d*sqrt(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_51():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)*sin(c + d*x)**3
    F = -46*a**3*sin(c + d*x)**4*cos(c + d*x)/(99*d*sqrt(a*sin(c + d*x) + a)) - 710*a**3*sin(c + d*x)**3*cos(c + d*x)/(693*d*sqrt(a*sin(c + d*x) + a)) - 284*a**3*cos(c + d*x)/(99*d*sqrt(a*sin(c + d*x) + a)) - 2*a**2*sqrt(a*sin(c + d*x) + a)*sin(c + d*x)**4*cos(c + d*x)/(11*d) + 568*a**2*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(693*d) - 284*a*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)/(231*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_52():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)*sin(c + d*x)**2
    F = -832*a**3*cos(c + d*x)/(315*d*sqrt(a*sin(c + d*x) + a)) - 208*a**2*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(315*d) - 26*a*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)/(105*d) + 4*(a*sin(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)/(63*d) - 2*(a*sin(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)/(9*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_53():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)*sin(c + d*x)
    F = -64*a**3*cos(c + d*x)/(21*d*sqrt(a*sin(c + d*x) + a)) - 16*a**2*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(21*d) - 2*a*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)/(7*d) - 2*(a*sin(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_54():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = -64*a**3*cos(c + d*x)/(15*d*sqrt(a*sin(c + d*x) + a)) - 16*a**2*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(15*d) - 2*a*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_55():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)*csc(c + d*x)
    F = -2*a**(sympy.S(5)/2)*atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/d - 14*a**3*cos(c + d*x)/(3*d*sqrt(a*sin(c + d*x) + a)) - 2*a**2*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_56():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)*csc(c + d*x)**2
    F = -5*a**(sympy.S(5)/2)*atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/d - a**3*cos(c + d*x)/(d*sqrt(a*sin(c + d*x) + a)) - a**2*sqrt(a*sin(c + d*x) + a)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_57():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)*csc(c + d*x)**3
    F = -19*a**(sympy.S(5)/2)*atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/(4*d) - 9*a**3*cot(c + d*x)/(4*d*sqrt(a*sin(c + d*x) + a)) - a**2*sqrt(a*sin(c + d*x) + a)*cot(c + d*x)*csc(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_58():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)*csc(c + d*x)**4
    F = -25*a**(sympy.S(5)/2)*atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/(8*d) - 13*a**3*cot(c + d*x)*csc(c + d*x)/(12*d*sqrt(a*sin(c + d*x) + a)) - 25*a**3*cot(c + d*x)/(8*d*sqrt(a*sin(c + d*x) + a)) - a**2*sqrt(a*sin(c + d*x) + a)*cot(c + d*x)*csc(c + d*x)**2/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_59():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)*csc(c + d*x)**5
    F = -163*a**(sympy.S(5)/2)*atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/(64*d) - 17*a**3*cot(c + d*x)*csc(c + d*x)**2/(24*d*sqrt(a*sin(c + d*x) + a)) - 163*a**3*cot(c + d*x)*csc(c + d*x)/(96*d*sqrt(a*sin(c + d*x) + a)) - 163*a**3*cot(c + d*x)/(64*d*sqrt(a*sin(c + d*x) + a)) - a**2*sqrt(a*sin(c + d*x) + a)*cot(c + d*x)*csc(c + d*x)**3/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_60():
    f = sin(c + d*x)**3/sqrt(a*sin(c + d*x) + a)
    F = -2*sin(c + d*x)**2*cos(c + d*x)/(5*d*sqrt(a*sin(c + d*x) + a)) - 28*cos(c + d*x)/(15*d*sqrt(a*sin(c + d*x) + a)) + 2*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(15*a*d) + sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_61():
    f = sin(c + d*x)**2/sqrt(a*sin(c + d*x) + a)
    F = 4*cos(c + d*x)/(3*d*sqrt(a*sin(c + d*x) + a)) - 2*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(3*a*d) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_62():
    f = sin(c + d*x)/sqrt(a*sin(c + d*x) + a)
    F = -2*cos(c + d*x)/(d*sqrt(a*sin(c + d*x) + a)) + sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_63():
    f = 1/sqrt(a*sin(c + d*x) + a)
    F = -sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_64():
    f = csc(c + d*x)/sqrt(a*sin(c + d*x) + a)
    F = -2*atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/(sqrt(a)*d) + sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_65():
    f = csc(c + d*x)**2/sqrt(a*sin(c + d*x) + a)
    F = -cot(c + d*x)/(d*sqrt(a*sin(c + d*x) + a)) + atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/(sqrt(a)*d) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_66():
    f = csc(c + d*x)**3/sqrt(a*sin(c + d*x) + a)
    F = -cot(c + d*x)*csc(c + d*x)/(2*d*sqrt(a*sin(c + d*x) + a)) + cot(c + d*x)/(4*d*sqrt(a*sin(c + d*x) + a)) - 7*atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/(4*sqrt(a)*d) + sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_67():
    f = sin(c + d*x)**4/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = sin(c + d*x)**3*cos(c + d*x)/(2*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 9*sin(c + d*x)**2*cos(c + d*x)/(10*a*d*sqrt(a*sin(c + d*x) + a)) - 31*cos(c + d*x)/(5*a*d*sqrt(a*sin(c + d*x) + a)) + 13*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(10*a**2*d) + 15*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_68():
    f = sin(c + d*x)**3/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = sin(c + d*x)**2*cos(c + d*x)/(2*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + 13*cos(c + d*x)/(3*a*d*sqrt(a*sin(c + d*x) + a)) - 7*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(6*a**2*d) - 11*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_69():
    f = sin(c + d*x)**2/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = -cos(c + d*x)/(2*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 2*cos(c + d*x)/(a*d*sqrt(a*sin(c + d*x) + a)) + 7*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_70():
    f = sin(c + d*x)/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = cos(c + d*x)/(2*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 3*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_71():
    f = (a*sin(c + d*x) + a)**(sympy.S(-3)/2)
    F = -cos(c + d*x)/(2*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_72():
    f = csc(c + d*x)/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = cos(c + d*x)/(2*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 2*atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/(a**(sympy.S(3)/2)*d) + 5*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_73():
    f = csc(c + d*x)**2/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = cot(c + d*x)/(2*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 3*cot(c + d*x)/(2*a*d*sqrt(a*sin(c + d*x) + a)) + 3*atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/(a**(sympy.S(3)/2)*d) - 9*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_74():
    f = csc(c + d*x)**3/(a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = cot(c + d*x)*csc(c + d*x)/(2*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - cot(c + d*x)*csc(c + d*x)/(a*d*sqrt(a*sin(c + d*x) + a)) + 7*cot(c + d*x)/(4*a*d*sqrt(a*sin(c + d*x) + a)) - 19*atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/(4*a**(sympy.S(3)/2)*d) + 13*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_75():
    f = sin(c + d*x)**5/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = sin(c + d*x)**4*cos(c + d*x)/(4*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) + 21*sin(c + d*x)**3*cos(c + d*x)/(16*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 157*sin(c + d*x)**2*cos(c + d*x)/(80*a**2*d*sqrt(a*sin(c + d*x) + a)) - 1729*cos(c + d*x)/(120*a**2*d*sqrt(a*sin(c + d*x) + a)) + 787*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(240*a**3*d) + 283*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_76():
    f = sin(c + d*x)**4/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = sin(c + d*x)**3*cos(c + d*x)/(4*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) + 17*sin(c + d*x)**2*cos(c + d*x)/(16*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) + 197*cos(c + d*x)/(24*a**2*d*sqrt(a*sin(c + d*x) + a)) - 95*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(48*a**3*d) - 163*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_77():
    f = sin(c + d*x)**3/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = sin(c + d*x)**2*cos(c + d*x)/(4*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) - 13*cos(c + d*x)/(16*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 9*cos(c + d*x)/(4*a**2*d*sqrt(a*sin(c + d*x) + a)) + 75*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_78():
    f = sin(c + d*x)**2/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = -cos(c + d*x)/(4*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) + 13*cos(c + d*x)/(16*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 19*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_79():
    f = sin(c + d*x)/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = cos(c + d*x)/(4*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) - 5*cos(c + d*x)/(16*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 5*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_80():
    f = (a*sin(c + d*x) + a)**(sympy.S(-5)/2)
    F = -cos(c + d*x)/(4*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) - 3*cos(c + d*x)/(16*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 3*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_81():
    f = csc(c + d*x)/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = cos(c + d*x)/(4*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) + 11*cos(c + d*x)/(16*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 2*atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/(a**(sympy.S(5)/2)*d) + 43*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_82():
    f = csc(c + d*x)**2/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = cot(c + d*x)/(4*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) + 15*cot(c + d*x)/(16*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 35*cot(c + d*x)/(16*a**2*d*sqrt(a*sin(c + d*x) + a)) + 5*atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/(a**(sympy.S(5)/2)*d) - 115*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_83():
    f = csc(c + d*x)**3/(a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = cot(c + d*x)*csc(c + d*x)/(4*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) + 19*cot(c + d*x)*csc(c + d*x)/(16*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 31*cot(c + d*x)*csc(c + d*x)/(16*a**2*d*sqrt(a*sin(c + d*x) + a)) + 63*cot(c + d*x)/(16*a**2*d*sqrt(a*sin(c + d*x) + a)) - 39*atanh(sqrt(a)*cos(c + d*x)/sqrt(a*sin(c + d*x) + a))/(4*a**(sympy.S(5)/2)*d) + 219*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_84():
    f = sqrt(a*sin(e + f*x) + a)/sqrt(sin(e + f*x))
    F = -2*sqrt(a)*asin(sqrt(a)*cos(e + f*x)/sqrt(a*sin(e + f*x) + a))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_85():
    f = sqrt(-a*sin(e + f*x) + a)/sqrt(-sin(e + f*x))
    F = 2*sqrt(a)*asin(sqrt(a)*cos(e + f*x)/sqrt(-a*sin(e + f*x) + a))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_86():
    f = 1/(sqrt(sin(x) + 1)*sqrt(sin(x)))
    F = -sqrt(2)*asin(cos(x)/(sin(x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_87():
    f = 1/(sqrt(a*sin(x) + a)*sqrt(sin(x)))
    F = -sqrt(2)*atan(sqrt(2)*sqrt(a)*cos(x)/(2*sqrt(a*sin(x) + a)*sqrt(sin(x))))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_88():
    f = 1/(sqrt(1 - sin(x))*sqrt(sin(x)))
    F = sqrt(2)*atanh(sqrt(2)*cos(x)/(2*sqrt(1 - sin(x))*sqrt(sin(x))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_89():
    f = 1/(sqrt(-a*sin(x) + a)*sqrt(sin(x)))
    F = sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(x)/(2*sqrt(-a*sin(x) + a)*sqrt(sin(x))))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_90():
    f = sin(c + d*x)**(sympy.S(1)/3)/(a*sin(c + d*x) + a)**2
    F = -sin(c + d*x)**(sympy.S(1)/3)*cos(c + d*x)/(3*d*(a*sin(c + d*x) + a)**2) - sin(c + d*x)**(sympy.S(4)/3)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), sin(c + d*x)**2)/(36*a**2*d*sqrt(cos(c + d*x)**2)) + 4*sin(c + d*x)**(sympy.S(1)/3)*cos(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), sin(c + d*x)**2)/(9*a**2*d*sqrt(cos(c + d*x)**2)) - sin(c + d*x)**(sympy.S(1)/3)*cos(c + d*x)/(9*a**2*d*(sin(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_91():
    f = (a*sin(c + d*x) + a)**(sympy.S(2)/3)*sin(c + d*x)**3
    F = -3*(a*sin(c + d*x) + a)**(sympy.S(2)/3)*sin(c + d*x)**2*cos(c + d*x)/(11*d) - 63*(a*sin(c + d*x) + a)**(sympy.S(2)/3)*cos(c + d*x)/(220*d) - 67*2**(sympy.S(1)/6)*(a*sin(c + d*x) + a)**(sympy.S(2)/3)*cos(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(110*d*(sin(c + d*x) + 1)**(sympy.S(7)/6)) - 3*(a*sin(c + d*x) + a)**(sympy.S(5)/3)*cos(c + d*x)/(44*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_92():
    f = (a*sin(c + d*x) + a)**(sympy.S(2)/3)*sin(c + d*x)**2
    F = 9*(a*sin(c + d*x) + a)**(sympy.S(2)/3)*cos(c + d*x)/(40*d) - 19*2**(sympy.S(1)/6)*(a*sin(c + d*x) + a)**(sympy.S(2)/3)*cos(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(20*d*(sin(c + d*x) + 1)**(sympy.S(7)/6)) - 3*(a*sin(c + d*x) + a)**(sympy.S(5)/3)*cos(c + d*x)/(8*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_93():
    f = (a*sin(c + d*x) + a)**(sympy.S(2)/3)*sin(c + d*x)
    F = -3*(a*sin(c + d*x) + a)**(sympy.S(2)/3)*cos(c + d*x)/(5*d) - 4*2**(sympy.S(1)/6)*(a*sin(c + d*x) + a)**(sympy.S(2)/3)*cos(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(5*d*(sin(c + d*x) + 1)**(sympy.S(7)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_94():
    f = (a*sin(c + d*x) + a)**(sympy.S(2)/3)
    F = -2*2**(sympy.S(1)/6)*(a*sin(c + d*x) + a)**(sympy.S(2)/3)*cos(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*(sin(c + d*x) + 1)**(sympy.S(7)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_95():
    f = (a*sin(c + d*x) + a)**(sympy.S(2)/3)*csc(c + d*x)
    F = -2*2**(sympy.S(1)/6)*(a*sin(c + d*x) + a)**(sympy.S(2)/3)*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S(-1)/6, 1, sympy.S(3)/2, sympy.S.Half - sin(c + d*x)/2, 1 - sin(c + d*x))/(d*(sin(c + d*x) + 1)**(sympy.S(7)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_96():
    f = (a*sin(c + d*x) + a)**(sympy.S(2)/3)*csc(c + d*x)**2
    F = -2*2**(sympy.S(1)/6)*(a*sin(c + d*x) + a)**(sympy.S(2)/3)*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S(-1)/6, 2, sympy.S(3)/2, sympy.S.Half - sin(c + d*x)/2, 1 - sin(c + d*x))/(d*(sin(c + d*x) + 1)**(sympy.S(7)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_97():
    f = (a*sin(c + d*x) + a)**(sympy.S(4)/3)*sin(c + d*x)**3
    F = -388*2**(sympy.S(5)/6)*a*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*cos(c + d*x)*hyper((sympy.S(-5)/6, sympy.S.Half), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(455*d*(sin(c + d*x) + 1)**(sympy.S(5)/6)) - 3*(a*sin(c + d*x) + a)**(sympy.S(4)/3)*sin(c + d*x)**2*cos(c + d*x)/(13*d) - 72*(a*sin(c + d*x) + a)**(sympy.S(4)/3)*cos(c + d*x)/(455*d) - 6*(a*sin(c + d*x) + a)**(sympy.S(7)/3)*cos(c + d*x)/(65*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_98():
    f = (a*sin(c + d*x) + a)**(sympy.S(4)/3)*sin(c + d*x)**2
    F = -37*2**(sympy.S(5)/6)*a*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*cos(c + d*x)*hyper((sympy.S(-5)/6, sympy.S.Half), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(35*d*(sin(c + d*x) + 1)**(sympy.S(5)/6)) + 9*(a*sin(c + d*x) + a)**(sympy.S(4)/3)*cos(c + d*x)/(70*d) - 3*(a*sin(c + d*x) + a)**(sympy.S(7)/3)*cos(c + d*x)/(10*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_99():
    f = (a*sin(c + d*x) + a)**(sympy.S(4)/3)*sin(c + d*x)
    F = -8*2**(sympy.S(5)/6)*a*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*cos(c + d*x)*hyper((sympy.S(-5)/6, sympy.S.Half), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(7*d*(sin(c + d*x) + 1)**(sympy.S(5)/6)) - 3*(a*sin(c + d*x) + a)**(sympy.S(4)/3)*cos(c + d*x)/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_100():
    f = (a*sin(c + d*x) + a)**(sympy.S(4)/3)
    F = -2*2**(sympy.S(5)/6)*a*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*cos(c + d*x)*hyper((sympy.S(-5)/6, sympy.S.Half), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*(sin(c + d*x) + 1)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_101():
    f = (a*sin(c + d*x) + a)**(sympy.S(4)/3)*csc(c + d*x)
    F = -2*2**(sympy.S(5)/6)*a*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S(-5)/6, 1, sympy.S(3)/2, sympy.S.Half - sin(c + d*x)/2, 1 - sin(c + d*x))/(d*(sin(c + d*x) + 1)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_102():
    f = (a*sin(c + d*x) + a)**(sympy.S(4)/3)*csc(c + d*x)**2
    F = -2*2**(sympy.S(5)/6)*a*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S(-5)/6, 2, sympy.S(3)/2, sympy.S.Half - sin(c + d*x)/2, 1 - sin(c + d*x))/(d*(sin(c + d*x) + 1)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_103():
    f = sin(c + d*x)**3/(a*sin(c + d*x) + a)**(sympy.S(1)/3)
    F = -3*sin(c + d*x)**2*cos(c + d*x)/(8*d*(a*sin(c + d*x) + a)**(sympy.S(1)/3)) - 99*cos(c + d*x)/(80*d*(a*sin(c + d*x) + a)**(sympy.S(1)/3)) + 37*2**(sympy.S(1)/6)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(80*d*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*(sin(c + d*x) + 1)**(sympy.S(1)/6)) + 3*(a*sin(c + d*x) + a)**(sympy.S(2)/3)*cos(c + d*x)/(40*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_104():
    f = sin(c + d*x)**2/(a*sin(c + d*x) + a)**(sympy.S(1)/3)
    F = 9*cos(c + d*x)/(10*d*(a*sin(c + d*x) + a)**(sympy.S(1)/3)) - 7*2**(sympy.S(1)/6)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(10*d*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*(sin(c + d*x) + 1)**(sympy.S(1)/6)) - 3*(a*sin(c + d*x) + a)**(sympy.S(2)/3)*cos(c + d*x)/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_105():
    f = sin(c + d*x)/(a*sin(c + d*x) + a)**(sympy.S(1)/3)
    F = -3*cos(c + d*x)/(2*d*(a*sin(c + d*x) + a)**(sympy.S(1)/3)) + 2**(sympy.S(1)/6)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(2*d*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*(sin(c + d*x) + 1)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_106():
    f = (a*sin(c + d*x) + a)**(sympy.S(-1)/3)
    F = -2**(sympy.S(1)/6)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*(sin(c + d*x) + 1)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_107():
    f = csc(c + d*x)/(a*sin(c + d*x) + a)**(sympy.S(1)/3)
    F = -2**(sympy.S(1)/6)*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S(5)/6, 1, sympy.S(3)/2, sympy.S.Half - sin(c + d*x)/2, 1 - sin(c + d*x))/(d*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*(sin(c + d*x) + 1)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_108():
    f = csc(c + d*x)**2/(a*sin(c + d*x) + a)**(sympy.S(1)/3)
    F = -2**(sympy.S(1)/6)*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S(5)/6, 2, sympy.S(3)/2, sympy.S.Half - sin(c + d*x)/2, 1 - sin(c + d*x))/(d*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*(sin(c + d*x) + 1)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_109():
    f = sin(c + d*x)**3/(a*sin(c + d*x) + a)**(sympy.S(4)/3)
    F = -3*sin(c + d*x)**2*cos(c + d*x)/(5*d*(a*sin(c + d*x) + a)**(sympy.S(4)/3)) + 6*cos(c + d*x)/(5*d*(a*sin(c + d*x) + a)**(sympy.S(4)/3)) + 6*cos(c + d*x)/(5*a*d*(a*sin(c + d*x) + a)**(sympy.S(1)/3)) - 2*2**(sympy.S(1)/6)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(a*d*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*(sin(c + d*x) + 1)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_110():
    f = sin(c + d*x)**2/(a*sin(c + d*x) + a)**(sympy.S(4)/3)
    F = -3*cos(c + d*x)/(5*d*(a*sin(c + d*x) + a)**(sympy.S(4)/3)) - 3*cos(c + d*x)/(2*a*d*(a*sin(c + d*x) + a)**(sympy.S(1)/3)) + 13*2**(sympy.S(1)/6)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(10*a*d*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*(sin(c + d*x) + 1)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_111():
    f = sin(c + d*x)/(a*sin(c + d*x) + a)**(sympy.S(4)/3)
    F = 3*cos(c + d*x)/(5*d*(a*sin(c + d*x) + a)**(sympy.S(4)/3)) - 4*2**(sympy.S(1)/6)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(5*a*d*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*(sin(c + d*x) + 1)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_112():
    f = (a*sin(c + d*x) + a)**(sympy.S(-4)/3)
    F = -2**(sympy.S(1)/6)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S(11)/6), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(2*a*d*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*(sin(c + d*x) + 1)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_113():
    f = csc(c + d*x)/(a*sin(c + d*x) + a)**(sympy.S(4)/3)
    F = -2**(sympy.S(1)/6)*cos(c + d*x)*appellf1(sympy.S.Half, 1, sympy.S(11)/6, sympy.S(3)/2, 1 - sin(c + d*x), sympy.S.Half - sin(c + d*x)/2)/(2*a*d*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*(sin(c + d*x) + 1)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_114():
    f = csc(c + d*x)**2/(a*sin(c + d*x) + a)**(sympy.S(4)/3)
    F = -2**(sympy.S(1)/6)*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S(11)/6, 2, sympy.S(3)/2, sympy.S.Half - sin(c + d*x)/2, 1 - sin(c + d*x))/(2*a*d*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*(sin(c + d*x) + 1)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_115():
    f = (sin(e + f*x) + 1)**(sympy.S(3)/2)*sin(e + f*x)**n
    F = -(8*n + 10)*cos(e + f*x)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), 1 - sin(e + f*x))/(f*(2*n + 3)*sqrt(sin(e + f*x) + 1)) - 2*sin(e + f*x)**(n + 1)*cos(e + f*x)/(f*(2*n + 3)*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_116():
    f = sqrt(sin(e + f*x) + 1)*sin(e + f*x)**n
    F = -2*cos(e + f*x)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), 1 - sin(e + f*x))/(f*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_117():
    f = sin(e + f*x)**n/sqrt(sin(e + f*x) + 1)
    F = -cos(e + f*x)*appellf1(sympy.S.Half, 1, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, 1 - sin(e + f*x))/(f*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_118():
    f = sin(e + f*x)**n/(sin(e + f*x) + 1)**(sympy.S(3)/2)
    F = -cos(e + f*x)*appellf1(sympy.S.Half, 2, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, 1 - sin(e + f*x))/(2*f*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_119():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)*sin(e + f*x)**n
    F = -2*a**2*(4*n + 5)*cos(e + f*x)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), 1 - sin(e + f*x))/(f*(2*n + 3)*sqrt(a*sin(e + f*x) + a)) - 2*a**2*sin(e + f*x)**(n + 1)*cos(e + f*x)/(f*(2*n + 3)*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_120():
    f = sqrt(a*sin(e + f*x) + a)*sin(e + f*x)**n
    F = -2*a*cos(e + f*x)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), 1 - sin(e + f*x))/(f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_121():
    f = sin(e + f*x)**n/sqrt(a*sin(e + f*x) + a)
    F = -cos(e + f*x)*appellf1(sympy.S.Half, 1, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, 1 - sin(e + f*x))/(f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_122():
    f = sin(e + f*x)**n/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -cos(e + f*x)*appellf1(sympy.S.Half, 2, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, 1 - sin(e + f*x))/(2*a*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_123():
    f = (d*sin(e + f*x))**n*(sin(e + f*x) + 1)**(sympy.S(3)/2)
    F = -2*(d*sin(e + f*x))**(n + 1)*cos(e + f*x)/(d*f*(2*n + 3)*sqrt(sin(e + f*x) + 1)) + (d*sin(e + f*x))**(n + 1)*(4*n + 5)*cos(e + f*x)*hyper((sympy.S.Half, n + 1), (n + 2,), sin(e + f*x))/(d*f*sqrt(1 - sin(e + f*x))*(n + 1)*(2*n + 3)*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_124():
    f = (d*sin(e + f*x))**n*sqrt(sin(e + f*x) + 1)
    F = (d*sin(e + f*x))**(n + 1)*cos(e + f*x)*hyper((sympy.S.Half, n + 1), (n + 2,), sin(e + f*x))/(d*f*sqrt(1 - sin(e + f*x))*(n + 1)*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_125():
    f = (d*sin(e + f*x))**n/sqrt(sin(e + f*x) + 1)
    F = -(d*sin(e + f*x))**n*cos(e + f*x)*appellf1(sympy.S.Half, 1, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, 1 - sin(e + f*x))/(f*sqrt(sin(e + f*x) + 1)*sin(e + f*x)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_126():
    f = (d*sin(e + f*x))**n/(sin(e + f*x) + 1)**(sympy.S(3)/2)
    F = -(d*sin(e + f*x))**n*cos(e + f*x)*appellf1(sympy.S.Half, 2, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, 1 - sin(e + f*x))/(2*f*sqrt(sin(e + f*x) + 1)*sin(e + f*x)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_127():
    f = (d*sin(e + f*x))**n*(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -2*a**2*(d*sin(e + f*x))**n*(4*n + 5)*cos(e + f*x)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), 1 - sin(e + f*x))/(f*(2*n + 3)*sqrt(a*sin(e + f*x) + a)*sin(e + f*x)**n) - 2*a**2*(d*sin(e + f*x))**(n + 1)*cos(e + f*x)/(d*f*(2*n + 3)*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_128():
    f = (d*sin(e + f*x))**n*sqrt(a*sin(e + f*x) + a)
    F = -2*a*(d*sin(e + f*x))**n*cos(e + f*x)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), 1 - sin(e + f*x))/(f*sqrt(a*sin(e + f*x) + a)*sin(e + f*x)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_129():
    f = (d*sin(e + f*x))**n/sqrt(a*sin(e + f*x) + a)
    F = -(d*sin(e + f*x))**n*cos(e + f*x)*appellf1(sympy.S.Half, 1, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, 1 - sin(e + f*x))/(f*sqrt(a*sin(e + f*x) + a)*sin(e + f*x)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_130():
    f = (d*sin(e + f*x))**n/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -(d*sin(e + f*x))**n*cos(e + f*x)*appellf1(sympy.S.Half, 2, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, 1 - sin(e + f*x))/(2*a*f*sqrt(a*sin(e + f*x) + a)*sin(e + f*x)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_131():
    f = (sin(e + f*x) + 1)**m*sin(e + f*x)**n
    F = -2**(m + sympy.S.Half)*cos(e + f*x)*appellf1(sympy.S.Half, -n, sympy.S.Half - m, sympy.S(3)/2, 1 - sin(e + f*x), sympy.S.Half - sin(e + f*x)/2)/(f*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_132():
    f = (-sin(e + f*x))**n*(1 - sin(e + f*x))**m
    F = 2**(m + sympy.S.Half)*cos(e + f*x)*appellf1(sympy.S.Half, -n, sympy.S.Half - m, sympy.S(3)/2, sin(e + f*x) + 1, sin(e + f*x)/2 + sympy.S.Half)/(f*sqrt(1 - sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_133():
    f = (d*sin(e + f*x))**n*(sin(e + f*x) + 1)**m
    F = -2**(m + sympy.S.Half)*(d*sin(e + f*x))**n*cos(e + f*x)*appellf1(sympy.S.Half, -n, sympy.S.Half - m, sympy.S(3)/2, 1 - sin(e + f*x), sympy.S.Half - sin(e + f*x)/2)/(f*sqrt(sin(e + f*x) + 1)*sin(e + f*x)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_134():
    f = (d*sin(e + f*x))**n*(1 - sin(e + f*x))**m
    F = 2**(m + sympy.S.Half)*(d*sin(e + f*x))**n*cos(e + f*x)*appellf1(sympy.S.Half, -n, sympy.S.Half - m, sympy.S(3)/2, sin(e + f*x) + 1, sin(e + f*x)/2 + sympy.S.Half)/(f*(-sin(e + f*x))**n*sqrt(1 - sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_135():
    f = (a*sin(e + f*x) + a)**m*sin(e + f*x)**n
    F = -2**(m + sympy.S.Half)*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(-m + sympy.S(-1)/2)*cos(e + f*x)*appellf1(sympy.S.Half, -n, sympy.S.Half - m, sympy.S(3)/2, 1 - sin(e + f*x), sympy.S.Half - sin(e + f*x)/2)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_136():
    f = (-sin(e + f*x))**n*(-a*sin(e + f*x) + a)**m
    F = 2**(m + sympy.S.Half)*(1 - sin(e + f*x))**(-m + sympy.S(-1)/2)*(-a*sin(e + f*x) + a)**m*cos(e + f*x)*appellf1(sympy.S.Half, -n, sympy.S.Half - m, sympy.S(3)/2, sin(e + f*x) + 1, sin(e + f*x)/2 + sympy.S.Half)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_137():
    f = (d*sin(e + f*x))**n*(a*sin(e + f*x) + a)**m
    F = -2**(m + sympy.S.Half)*(d*sin(e + f*x))**n*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(-m + sympy.S(-1)/2)*cos(e + f*x)*appellf1(sympy.S.Half, -n, sympy.S.Half - m, sympy.S(3)/2, 1 - sin(e + f*x), sympy.S.Half - sin(e + f*x)/2)/(f*sin(e + f*x)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_138():
    f = (d*sin(e + f*x))**n*(-a*sin(e + f*x) + a)**m
    F = 2**(m + sympy.S.Half)*(d*sin(e + f*x))**n*(1 - sin(e + f*x))**(-m + sympy.S(-1)/2)*(-a*sin(e + f*x) + a)**m*cos(e + f*x)*appellf1(sympy.S.Half, -n, sympy.S.Half - m, sympy.S(3)/2, sin(e + f*x) + 1, sin(e + f*x)/2 + sympy.S.Half)/(f*(-sin(e + f*x))**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_139():
    f = (a*sin(c + d*x) + a)**n*sin(c + d*x)**4
    F = -2**(n + sympy.S.Half)*(a*sin(c + d*x) + a)**n*(sin(c + d*x) + 1)**(-n + sympy.S(-1)/2)*(n**4 + 6*n**3 + 17*n**2 + 12*n + 9)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S.Half - n), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*(n + 1)*(n + 2)*(n + 3)*(n + 4)) - n*(a*sin(c + d*x) + a)**n*sin(c + d*x)**2*cos(c + d*x)/(d*(n + 3)*(n + 4)) - (a*sin(c + d*x) + a)**n*sin(c + d*x)**3*cos(c + d*x)/(d*(n + 4)) + (a*sin(c + d*x) + a)**n*(-n**2 - n + 9)*cos(c + d*x)/(d*(n + 1)*(n + 2)*(n + 3)*(n + 4)) - (a*sin(c + d*x) + a)**(n + 1)*(n**2 + 3*n + 9)*cos(c + d*x)/(a*d*(n + 2)*(n + 3)*(n + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_140():
    f = (a*sin(c + d*x) + a)**n*sin(c + d*x)**3
    F = -2**(n + sympy.S.Half)*n*(a*sin(c + d*x) + a)**n*(sin(c + d*x) + 1)**(-n + sympy.S(-1)/2)*(n**2 + 3*n + 5)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S.Half - n), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*(n + 1)*(n + 2)*(n + 3)) - (a*sin(c + d*x) + a)**n*sin(c + d*x)**2*cos(c + d*x)/(d*(n + 3)) - (n + 4)*(a*sin(c + d*x) + a)**n*cos(c + d*x)/(d*(n + 1)*(n + 2)*(n + 3)) - n*(a*sin(c + d*x) + a)**(n + 1)*cos(c + d*x)/(a*d*(n**2 + 5*n + 6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_141():
    f = (a*sin(c + d*x) + a)**n*sin(c + d*x)**2
    F = -2**(n + sympy.S.Half)*(a*sin(c + d*x) + a)**n*(sin(c + d*x) + 1)**(-n + sympy.S(-1)/2)*(n**2 + n + 1)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S.Half - n), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*(n + 1)*(n + 2)) + (a*sin(c + d*x) + a)**n*cos(c + d*x)/(d*(n**2 + 3*n + 2)) - (a*sin(c + d*x) + a)**(n + 1)*cos(c + d*x)/(a*d*(n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_142():
    f = (a*sin(c + d*x) + a)**n*sin(c + d*x)
    F = -2**(n + sympy.S.Half)*n*(a*sin(c + d*x) + a)**n*(sin(c + d*x) + 1)**(-n + sympy.S(-1)/2)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S.Half - n), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*(n + 1)) - (a*sin(c + d*x) + a)**n*cos(c + d*x)/(d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_143():
    f = (a*sin(c + d*x) + a)**n
    F = -2**(n + sympy.S.Half)*(a*sin(c + d*x) + a)**n*(sin(c + d*x) + 1)**(-n + sympy.S(-1)/2)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S.Half - n), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_144():
    f = (a*sin(c + d*x) + a)**n*csc(c + d*x)
    F = -2**(n + sympy.S.Half)*(a*sin(c + d*x) + a)**n*(sin(c + d*x) + 1)**(-n + sympy.S(-1)/2)*cos(c + d*x)*appellf1(sympy.S.Half, 1, sympy.S.Half - n, sympy.S(3)/2, 1 - sin(c + d*x), sympy.S.Half - sin(c + d*x)/2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_145():
    f = (a*sin(c + d*x) + a)**n*csc(c + d*x)**2
    F = -2**(n + sympy.S.Half)*(a*sin(c + d*x) + a)**n*(sin(c + d*x) + 1)**(-n + sympy.S(-1)/2)*cos(c + d*x)*appellf1(sympy.S.Half, 2, sympy.S.Half - n, sympy.S(3)/2, 1 - sin(c + d*x), sympy.S.Half - sin(c + d*x)/2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_146():
    f = (sin(c + d*x) + 1)**n
    F = -2**(n + sympy.S.Half)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S.Half - n), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*sqrt(sin(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_147():
    f = (1 - sin(c + d*x))**n
    F = 2**(n + sympy.S.Half)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S.Half - n), (sympy.S(3)/2,), sin(c + d*x)/2 + sympy.S.Half)/(d*sqrt(1 - sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_148():
    f = (a + b*sin(e + f*x))*sin(e + f*x)**3
    F = a*cos(e + f*x)**3/(3*f) - a*cos(e + f*x)/f + 3*b*x/8 - b*sin(e + f*x)**3*cos(e + f*x)/(4*f) - 3*b*sin(e + f*x)*cos(e + f*x)/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_149():
    f = (a + b*sin(e + f*x))*sin(e + f*x)**2
    F = a*x/2 - a*sin(e + f*x)*cos(e + f*x)/(2*f) + b*cos(e + f*x)**3/(3*f) - b*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_150():
    f = (a + b*sin(e + f*x))*sin(e + f*x)
    F = -a*cos(e + f*x)/f + b*x/2 - b*sin(e + f*x)*cos(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_151():
    f = a + b*sin(e + f*x)
    F = a*x - b*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_152():
    f = (a + b*sin(e + f*x))*csc(e + f*x)
    F = -a*atanh(cos(e + f*x))/f + b*x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_153():
    f = (a + b*sin(e + f*x))*csc(e + f*x)**2
    F = -a*cot(e + f*x)/f - b*atanh(cos(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_154():
    f = (a + b*sin(e + f*x))*csc(e + f*x)**3
    F = -a*cot(e + f*x)*csc(e + f*x)/(2*f) - a*atanh(cos(e + f*x))/(2*f) - b*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_155():
    f = (a + b*sin(e + f*x))*csc(e + f*x)**4
    F = -a*cot(e + f*x)**3/(3*f) - a*cot(e + f*x)/f - b*cot(e + f*x)*csc(e + f*x)/(2*f) - b*atanh(cos(e + f*x))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_156():
    f = (a + b*sin(e + f*x))**2*sin(e + f*x)**3
    F = 3*a*b*x/4 - a*b*sin(e + f*x)**3*cos(e + f*x)/(2*f) - 3*a*b*sin(e + f*x)*cos(e + f*x)/(4*f) - b**2*cos(e + f*x)**5/(5*f) - (a**2 + b**2)*cos(e + f*x)/f + (a**2 + 2*b**2)*cos(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_157():
    f = (a + b*sin(e + f*x))**2*sin(e + f*x)**2
    F = 2*a*b*cos(e + f*x)**3/(3*f) - 2*a*b*cos(e + f*x)/f - b**2*sin(e + f*x)**3*cos(e + f*x)/(4*f) + x*(a**2/2 + 3*b**2/8) - (4*a**2 + 3*b**2)*sin(e + f*x)*cos(e + f*x)/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_158():
    f = (a + b*sin(e + f*x))**2*sin(e + f*x)
    F = a*b*x - a*b*sin(e + f*x)*cos(e + f*x)/(3*f) - (a + b*sin(e + f*x))**2*cos(e + f*x)/(3*f) - (2*a**2 + 2*b**2)*cos(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_159():
    f = (a + b*sin(e + f*x))**2
    F = -2*a*b*cos(e + f*x)/f - b**2*sin(e + f*x)*cos(e + f*x)/(2*f) + x*(a**2 + b**2/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_160():
    f = (a + b*sin(e + f*x))**2*csc(e + f*x)
    F = -a**2*atanh(cos(e + f*x))/f + 2*a*b*x - b**2*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_161():
    f = (a + b*sin(e + f*x))**2*csc(e + f*x)**2
    F = -a**2*cot(e + f*x)/f - 2*a*b*atanh(cos(e + f*x))/f + b**2*x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_162():
    f = (a + b*sin(e + f*x))**2*csc(e + f*x)**3
    F = -a**2*cot(e + f*x)*csc(e + f*x)/(2*f) - 2*a*b*cot(e + f*x)/f - (a**2 + 2*b**2)*atanh(cos(e + f*x))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_163():
    f = (a + b*sin(e + f*x))**2*csc(e + f*x)**4
    F = -a**2*cot(e + f*x)*csc(e + f*x)**2/(3*f) - a*b*cot(e + f*x)*csc(e + f*x)/f - a*b*atanh(cos(e + f*x))/f - (2*a**2 + 3*b**2)*cot(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_164():
    f = (a + b*sin(e + f*x))**2*csc(e + f*x)**5
    F = -a**2*cot(e + f*x)*csc(e + f*x)**3/(4*f) - 2*a*b*cot(e + f*x)**3/(3*f) - 2*a*b*cot(e + f*x)/f - (3*a**2 + 4*b**2)*cot(e + f*x)*csc(e + f*x)/(8*f) - (3*a**2 + 4*b**2)*atanh(cos(e + f*x))/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_165():
    f = (a + b*sin(e + f*x))**3*sin(e + f*x)
    F = -a*(a + b*sin(e + f*x))**2*cos(e + f*x)/(4*f) - a*(a**2 + 4*b**2)*cos(e + f*x)/(2*f) + 3*b*x*(4*a**2 + b**2)/8 - b*(2*a**2 + 3*b**2)*sin(e + f*x)*cos(e + f*x)/(8*f) - (a + b*sin(e + f*x))**3*cos(e + f*x)/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_166():
    f = (a + b*sin(e + f*x))**3
    F = -5*a*b**2*sin(e + f*x)*cos(e + f*x)/(6*f) + a*x*(2*a**2 + 3*b**2)/2 - b*(a + b*sin(e + f*x))**2*cos(e + f*x)/(3*f) - 2*b*(4*a**2 + b**2)*cos(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_167():
    f = (a + b*sin(e + f*x))**3*csc(e + f*x)
    F = -a**3*atanh(cos(e + f*x))/f - 5*a*b**2*cos(e + f*x)/(2*f) - b**2*(a + b*sin(e + f*x))*cos(e + f*x)/(2*f) + b*x*(6*a**2 + b**2)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_168():
    f = (a + b*sin(e + f*x))**3*csc(e + f*x)**2
    F = -3*a**2*b*atanh(cos(e + f*x))/f - a**2*(a + b*sin(e + f*x))*cot(e + f*x)/f + 3*a*b**2*x + b*(a**2 - b**2)*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_169():
    f = (a + b*sin(e + f*x))**3*csc(e + f*x)**3
    F = -5*a**2*b*cot(e + f*x)/(2*f) - a**2*(a + b*sin(e + f*x))*cot(e + f*x)*csc(e + f*x)/(2*f) - a*(a**2 + 6*b**2)*atanh(cos(e + f*x))/(2*f) + b**3*x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_170():
    f = (a + b*sin(e + f*x))**3*csc(e + f*x)**4
    F = -7*a**2*b*cot(e + f*x)*csc(e + f*x)/(6*f) - a**2*(a + b*sin(e + f*x))*cot(e + f*x)*csc(e + f*x)**2/(3*f) - a*(2*a**2 + 9*b**2)*cot(e + f*x)/(3*f) - b*(3*a**2 + 2*b**2)*atanh(cos(e + f*x))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_171():
    f = (a + b*sin(e + f*x))**3*csc(e + f*x)**5
    F = -3*a**2*b*cot(e + f*x)*csc(e + f*x)**2/(4*f) - a**2*(a + b*sin(e + f*x))*cot(e + f*x)*csc(e + f*x)**3/(4*f) - 3*a*(a**2 + 4*b**2)*cot(e + f*x)*csc(e + f*x)/(8*f) - 3*a*(a**2 + 4*b**2)*atanh(cos(e + f*x))/(8*f) - b*(2*a**2 + b**2)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_172():
    f = (a + b*sin(e + f*x))**4
    F = -7*a*b*(a + b*sin(e + f*x))**2*cos(e + f*x)/(12*f) - a*b*(19*a**2 + 16*b**2)*cos(e + f*x)/(6*f) - b**2*(26*a**2 + 9*b**2)*sin(e + f*x)*cos(e + f*x)/(24*f) - b*(a + b*sin(e + f*x))**3*cos(e + f*x)/(4*f) + x*(a**4 + 3*a**2*b**2 + 3*b**4/8)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_173():
    f = sin(x)**4/(a + b*sin(x))
    F = 2*a**4*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(b**4*sqrt(a**2 - b**2)) + a*sin(x)*cos(x)/(2*b**2) - a*x*(2*a**2 + b**2)/(2*b**4) - sin(x)**2*cos(x)/(3*b) - (3*a**2 + 2*b**2)*cos(x)/(3*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_174():
    f = sin(x)**3/(a + b*sin(x))
    F = -2*a**3*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(b**3*sqrt(a**2 - b**2)) + a*cos(x)/b**2 - sin(x)*cos(x)/(2*b) + x*(2*a**2 + b**2)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_175():
    f = sin(x)**2/(a + b*sin(x))
    F = 2*a**2*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(b**2*sqrt(a**2 - b**2)) - a*x/b**2 - cos(x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_176():
    f = sin(x)/(a + b*sin(x))
    F = -2*a*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(b*sqrt(a**2 - b**2)) + x/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_177():
    f = 1/(a + b*sin(x))
    F = 2*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/sqrt(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_178():
    f = csc(x)/(a + b*sin(x))
    F = -2*b*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(a*sqrt(a**2 - b**2)) - atanh(cos(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_179():
    f = csc(x)**2/(a + b*sin(x))
    F = -cot(x)/a + 2*b**2*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(a**2*sqrt(a**2 - b**2)) + b*atanh(cos(x))/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_180():
    f = csc(x)**3/(a + b*sin(x))
    F = -cot(x)*csc(x)/(2*a) + b*cot(x)/a**2 - 2*b**3*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(a**3*sqrt(a**2 - b**2)) - (a**2 + 2*b**2)*atanh(cos(x))/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_181():
    f = csc(x)**4/(a + b*sin(x))
    F = -cot(x)*csc(x)**2/(3*a) + b*cot(x)*csc(x)/(2*a**2) - (2*a**2 + 3*b**2)*cot(x)/(3*a**3) + 2*b**4*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(a**4*sqrt(a**2 - b**2)) + b*(a**2 + 2*b**2)*atanh(cos(x))/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_182():
    f = sin(x)**4/(a + b*sin(x))**2
    F = -2*a**3*(3*a**2 - 4*b**2)*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(b**4*(a**2 - b**2)**(sympy.S(3)/2)) + a**2*sin(x)**2*cos(x)/(b*(a + b*sin(x))*(a**2 - b**2)) + a*(3*a**2 - 2*b**2)*cos(x)/(b**3*(a**2 - b**2)) - (3*a**2 - b**2)*sin(x)*cos(x)/(2*b**2*(a**2 - b**2)) + x*(6*a**2 + b**2)/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_183():
    f = sin(x)**3/(a + b*sin(x))**2
    F = a**2*sin(x)*cos(x)/(b*(a + b*sin(x))*(a**2 - b**2)) + 2*a**2*(2*a**2 - 3*b**2)*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(b**3*(a**2 - b**2)**(sympy.S(3)/2)) - 2*a*x/b**3 - (2*a**2 - b**2)*cos(x)/(b**2*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_184():
    f = sin(x)**2/(a + b*sin(x))**2
    F = a**2*cos(x)/(b*(a + b*sin(x))*(a**2 - b**2)) - 2*a*(a**2 - 2*b**2)*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(b**2*(a**2 - b**2)**(sympy.S(3)/2)) + x/b**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_185():
    f = sin(x)/(a + b*sin(x))**2
    F = -a*cos(x)/((a + b*sin(x))*(a**2 - b**2)) - 2*b*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_186():
    f = (a + b*sin(x))**(-2)
    F = 2*a*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(3)/2) + b*cos(x)/((a + b*sin(x))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_187():
    f = csc(x)/(a + b*sin(x))**2
    F = -b**2*cos(x)/(a*(a + b*sin(x))*(a**2 - b**2)) - 2*b*(2*a**2 - b**2)*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(a**2*(a**2 - b**2)**(sympy.S(3)/2)) - atanh(cos(x))/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_188():
    f = csc(x)**2/(a + b*sin(x))**2
    F = -b**2*cot(x)/(a*(a + b*sin(x))*(a**2 - b**2)) - (a**2 - 2*b**2)*cot(x)/(a**2*(a**2 - b**2)) + 2*b**2*(3*a**2 - 2*b**2)*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(a**3*(a**2 - b**2)**(sympy.S(3)/2)) + 2*b*atanh(cos(x))/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_189():
    f = csc(x)**3/(a + b*sin(x))**2
    F = -b**2*cot(x)*csc(x)/(a*(a + b*sin(x))*(a**2 - b**2)) - (a**2 - 3*b**2)*cot(x)*csc(x)/(2*a**2*(a**2 - b**2)) + b*(2*a**2 - 3*b**2)*cot(x)/(a**3*(a**2 - b**2)) - 2*b**3*(4*a**2 - 3*b**2)*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(a**4*(a**2 - b**2)**(sympy.S(3)/2)) - (a**2 + 6*b**2)*atanh(cos(x))/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_190():
    f = sin(x)**5/(a + b*sin(x))**3
    F = -a**3*(12*a**4 - 29*a**2*b**2 + 20*b**4)*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(b**5*(a**2 - b**2)**(sympy.S(5)/2)) + a**2*sin(x)**3*cos(x)/(2*b*(a + b*sin(x))**2*(a**2 - b**2)) + a**2*(4*a**2 - 7*b**2)*sin(x)**2*cos(x)/(2*b**2*(a + b*sin(x))*(a**2 - b**2)**2) + 3*a*(4*a**4 - 7*a**2*b**2 + 2*b**4)*cos(x)/(2*b**4*(a**2 - b**2)**2) - (6*a**4 - 10*a**2*b**2 + b**4)*sin(x)*cos(x)/(2*b**3*(a**2 - b**2)**2) + x*(12*a**2 + b**2)/(2*b**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_191():
    f = sin(x)**4/(a + b*sin(x))**3
    F = -3*a**3*(a**2 - 2*b**2)*cos(x)/(2*b**3*(a + b*sin(x))*(a**2 - b**2)**2) + a**2*sin(x)**2*cos(x)/(2*b*(a + b*sin(x))**2*(a**2 - b**2)) + 3*a**2*(2*a**4 - 5*a**2*b**2 + 4*b**4)*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(b**4*(a**2 - b**2)**(sympy.S(5)/2)) - 3*a*x/b**4 - (3*a**2 - 2*b**2)*cos(x)/(2*b**3*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_192():
    f = sin(x)**3/(a + b*sin(x))**3
    F = a**2*sin(x)*cos(x)/(2*b*(a + b*sin(x))**2*(a**2 - b**2)) + a**2*(2*a**2 - 5*b**2)*cos(x)/(2*b**2*(a + b*sin(x))*(a**2 - b**2)**2) - a*(2*a**4 - 5*a**2*b**2 + 6*b**4)*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(b**3*(a**2 - b**2)**(sympy.S(5)/2)) + x/b**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_193():
    f = sin(x)**2/(a + b*sin(x))**3
    F = a**2*cos(x)/(2*b*(a + b*sin(x))**2*(a**2 - b**2)) - a*(a**2 - 4*b**2)*cos(x)/(2*b*(a + b*sin(x))*(a**2 - b**2)**2) + (a**2 + 2*b**2)*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_194():
    f = sin(x)/(a + b*sin(x))**3
    F = -3*a*b*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(5)/2) - a*cos(x)/((a + b*sin(x))**2*(2*a**2 - 2*b**2)) - (a**2 + 2*b**2)*cos(x)/(2*(a + b*sin(x))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_195():
    f = (a + b*sin(x))**(-3)
    F = 3*a*b*cos(x)/(2*(a + b*sin(x))*(a**2 - b**2)**2) + b*cos(x)/((a + b*sin(x))**2*(2*a**2 - 2*b**2)) + (2*a**2 + b**2)*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_196():
    f = csc(x)/(a + b*sin(x))**3
    F = -b**2*cos(x)/(2*a*(a + b*sin(x))**2*(a**2 - b**2)) - b**2*(5*a**2 - 2*b**2)*cos(x)/(2*a**2*(a + b*sin(x))*(a**2 - b**2)**2) - b*(6*a**4 - 5*a**2*b**2 + 2*b**4)*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(a**3*(a**2 - b**2)**(sympy.S(5)/2)) - atanh(cos(x))/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_197():
    f = csc(x)**2/(a + b*sin(x))**3
    F = -b**2*cot(x)/(2*a*(a + b*sin(x))**2*(a**2 - b**2)) - 3*b**2*(2*a**2 - b**2)*cot(x)/(2*a**2*(a + b*sin(x))*(a**2 - b**2)**2) - (2*a**4 - 11*a**2*b**2 + 6*b**4)*cot(x)/(2*a**3*(a**2 - b**2)**2) + 3*b**2*(4*a**4 - 5*a**2*b**2 + 2*b**4)*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(a**4*(a**2 - b**2)**(sympy.S(5)/2)) + 3*b*atanh(cos(x))/a**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_198():
    f = csc(x)**3/(a + b*sin(x))**3
    F = -b**2*cot(x)*csc(x)/(2*a*(a + b*sin(x))**2*(a**2 - b**2)) - b**2*(7*a**2 - 4*b**2)*cot(x)*csc(x)/(2*a**2*(a + b*sin(x))*(a**2 - b**2)**2) - (a**4 - 10*a**2*b**2 + 6*b**4)*cot(x)*csc(x)/(2*a**3*(a**2 - b**2)**2) + 3*b*(2*a**4 - 7*a**2*b**2 + 4*b**4)*cot(x)/(2*a**4*(a**2 - b**2)**2) - b**3*(20*a**4 - 29*a**2*b**2 + 12*b**4)*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(a**5*(a**2 - b**2)**(sympy.S(5)/2)) - (a**2 + 12*b**2)*atanh(cos(x))/(2*a**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_199():
    f = (a + b*sin(c + d*x))**(-4)
    F = 5*a*b*cos(c + d*x)/(6*d*(a + b*sin(c + d*x))**2*(a**2 - b**2)**2) + a*(2*a**2 + 3*b**2)*atan((a*tan(c/2 + d*x/2) + b)/sqrt(a**2 - b**2))/(d*(a**2 - b**2)**(sympy.S(7)/2)) + b*(11*a**2 + 4*b**2)*cos(c + d*x)/(6*d*(a + b*sin(c + d*x))*(a**2 - b**2)**3) + b*cos(c + d*x)/(d*(a + b*sin(c + d*x))**3*(3*a**2 - 3*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_200():
    f = sqrt(a + b*sin(e + f*x))*sin(e + f*x)
    F = 2*a*sqrt(a + b*sin(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2*b/(a + b))/(3*b*f*sqrt((a + b*sin(e + f*x))/(a + b))) - 2*sqrt(a + b*sin(e + f*x))*cos(e + f*x)/(3*f) - sqrt((a + b*sin(e + f*x))/(a + b))*(2*a**2 - 2*b**2)*elliptic_f(e/2 + f*x/2 - pi/4, 2*b/(a + b))/(3*b*f*sqrt(a + b*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_201():
    f = sqrt(a + b*sin(e + f*x))
    F = 2*sqrt(a + b*sin(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2*b/(a + b))/(f*sqrt((a + b*sin(e + f*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_202():
    f = sqrt(a + b*sin(e + f*x))*csc(e + f*x)
    F = ((Integer(2) * Symbol('b') * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_203():
    f = sqrt(a + b*sin(e + f*x))*csc(e + f*x)**2
    F = (Integer(-1) * ((sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * (Symbol('f'))**(Integer(-1)))) + (Integer(-1) * ((sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('f') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('a') * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Symbol('b') * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_204():
    f = sin(e + f*x)/sqrt(a + b*sin(e + f*x))
    F = -2*a*sqrt((a + b*sin(e + f*x))/(a + b))*elliptic_f(e/2 + f*x/2 - pi/4, 2*b/(a + b))/(b*f*sqrt(a + b*sin(e + f*x))) + 2*sqrt(a + b*sin(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2*b/(a + b))/(b*f*sqrt((a + b*sin(e + f*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_205():
    f = 1/sqrt(a + b*sin(e + f*x))
    F = 2*sqrt((a + b*sin(e + f*x))/(a + b))*elliptic_f(e/2 + f*x/2 - pi/4, 2*b/(a + b))/(f*sqrt(a + b*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_206():
    f = csc(e + f*x)/sqrt(a + b*sin(e + f*x))
    F = (Integer(2) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_207():
    f = csc(e + f*x)**2/sqrt(a + b*sin(e + f*x))
    F = (Integer(-1) * ((sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('a') * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('a') * Symbol('f') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('a') * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_208():
    f = sqrt(a + b*sin(c + d*x))*sqrt(sin(c + d*x))
    F = (Integer(-1) * ((sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('d') * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.csc((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.csc((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.csc((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.csc((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))) + ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.csc((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.csc((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_209():
    f = 1/(sqrt(a + b*sin(c + d*x))*sqrt(sin(c + d*x)))
    F = -2*sqrt(a*(1 - csc(c + d*x))/(a + b))*sqrt(a*(csc(c + d*x) + 1)/(a - b))*sqrt(a + b)*tan(c + d*x)*elliptic_f(asin(sqrt(a + b*sin(c + d*x))/(sqrt(a + b)*sqrt(sin(c + d*x)))), -(a + b)/(a - b))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_210():
    f = (d*sin(e + f*x))**m*(a + b*sin(e + f*x))**3
    F = -a*b**2*(d*sin(e + f*x))**(m + 1)*(2*m + 7)*cos(e + f*x)/(d*f*(m + 2)*(m + 3)) + a*(d*sin(e + f*x))**(m + 1)*(a**2*(m + 2) + 3*b**2*(m + 1))*cos(e + f*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(d*f*(m + 1)*(m + 2)*sqrt(cos(e + f*x)**2)) - b**2*(d*sin(e + f*x))**(m + 1)*(a + b*sin(e + f*x))*cos(e + f*x)/(d*f*(m + 3)) + b*(d*sin(e + f*x))**(m + 2)*(3*a**2*(m + 3) + b**2*(m + 2))*cos(e + f*x)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), sin(e + f*x)**2)/(d**2*f*(m + 2)*(m + 3)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_211():
    f = (d*sin(e + f*x))**m*(a + b*sin(e + f*x))**2
    F = 2*a*b*(d*sin(e + f*x))**(m + 2)*cos(e + f*x)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), sin(e + f*x)**2)/(d**2*f*(m + 2)*sqrt(cos(e + f*x)**2)) - b**2*(d*sin(e + f*x))**(m + 1)*cos(e + f*x)/(d*f*(m + 2)) + (d*sin(e + f*x))**(m + 1)*(a**2*(m + 2) + b**2*(m + 1))*cos(e + f*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(d*f*(m + 1)*(m + 2)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_212():
    f = (d*sin(e + f*x))**m*(a + b*sin(e + f*x))
    F = a*(d*sin(e + f*x))**(m + 1)*cos(e + f*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(d*f*(m + 1)*sqrt(cos(e + f*x)**2)) + b*(d*sin(e + f*x))**(m + 2)*cos(e + f*x)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), sin(e + f*x)**2)/(d**2*f*(m + 2)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_213():
    f = (d*sin(e + f*x))**m/(a + b*sin(e + f*x))
    F = -a*d*(d*sin(e + f*x))**(m - 1)*(sin(e + f*x)**2)**(sympy.S.Half - m/2)*cos(e + f*x)*appellf1(sympy.S.Half, 1, sympy.S.Half - m/2, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(f*(a**2 - b**2)) + b*(d*sin(e + f*x))**m*cos(e + f*x)*appellf1(sympy.S.Half, 1, -m/2, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(f*(a**2 - b**2)*(sin(e + f*x)**2)**(m/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_214():
    f = (d*sin(e + f*x))**m/(a + b*sin(e + f*x))**2
    F = -a**2*d*(d*sin(e + f*x))**(m - 1)*(sin(e + f*x)**2)**(sympy.S.Half - m/2)*cos(e + f*x)*appellf1(sympy.S.Half, 2, sympy.S.Half - m/2, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(f*(a**2 - b**2)**2) + 2*a*b*(d*sin(e + f*x))**m*cos(e + f*x)*appellf1(sympy.S.Half, 2, -m/2, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(f*(a**2 - b**2)**2*(sin(e + f*x)**2)**(m/2)) - b**2*(d*sin(e + f*x))**(m + 1)*(sin(e + f*x)**2)**(-m/2 + sympy.S(-1)/2)*cos(e + f*x)*appellf1(sympy.S.Half, 2, -m/2 + sympy.S(-1)/2, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(d*f*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_215():
    f = (d*sin(e + f*x))**m/(a + b*sin(e + f*x))**3
    F = -a**3*d*(d*sin(e + f*x))**(m - 1)*(sin(e + f*x)**2)**(sympy.S.Half - m/2)*cos(e + f*x)*appellf1(sympy.S.Half, 3, sympy.S.Half - m/2, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(f*(a**2 - b**2)**3) + 3*a**2*b*(d*sin(e + f*x))**m*cos(e + f*x)*appellf1(sympy.S.Half, 3, -m/2, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(f*(a**2 - b**2)**3*(sin(e + f*x)**2)**(m/2)) - 3*a*b**2*(d*sin(e + f*x))**(m + 1)*(sin(e + f*x)**2)**(-m/2 + sympy.S(-1)/2)*cos(e + f*x)*appellf1(sympy.S.Half, 3, -m/2 + sympy.S(-1)/2, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(d*f*(a**2 - b**2)**3) + b**3*(d*sin(e + f*x))**m*cos(e + f*x)*appellf1(sympy.S.Half, 3, -m/2 - 1, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(f*(a**2 - b**2)**3*(sin(e + f*x)**2)**(m/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_216():
    f = (a + b*sin(c + d*x))**2*sin(c + d*x)**(-a**2/(a**2 + b**2) - 1)
    F = 2*a*(a**2 + b**2)*sin(c + d*x)**(b**2/(a**2 + b**2))*cos(c + d*x)*hyper((sympy.S.Half, b**2/(2*a**2 + 2*b**2)), (-a**2/(2*(a**2 + b**2)) + sympy.S(3)/2,), sin(c + d*x)**2)/(b*d*sqrt(cos(c + d*x)**2)) - (a**2 + b**2)*cos(c + d*x)/(d*sin(c + d*x)**(a**2/(a**2 + b**2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_217():
    f = (2*sin(c + d*x) + 1)**2/sin(c + d*x)**(sympy.S(6)/5)
    F = -5*cos(c + d*x)/(d*sin(c + d*x)**(sympy.S(1)/5)) + 5*sin(c + d*x)**(sympy.S(4)/5)*cos(c + d*x)*hyper((sympy.S(2)/5, sympy.S.Half), (sympy.S(7)/5,), sin(c + d*x)**2)/(d*sqrt(cos(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_218():
    f = (a + b*sin(c + d*x))**n*sin(c + d*x)**m
    F = sympy.Function('Unintegrable')(((sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_219():
    f = (a + b*sin(c + d*x))**n*sin(c + d*x)**3
    F = 2*a*(a + b*sin(c + d*x))**(n + 1)*cos(c + d*x)/(b**2*d*(n + 2)*(n + 3)) + sqrt(2)*a*(a + b*sin(c + d*x))**n*(2*a**2 + b**2*(n**2 + 5*n + 4))*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, -n, sympy.S(3)/2, sympy.S.Half - sin(c + d*x)/2, b*(1 - sin(c + d*x))/(a + b))/(b**3*d*((a + b*sin(c + d*x))/(a + b))**n*(n + 2)*(n + 3)*sqrt(sin(c + d*x) + 1)) - (a + b*sin(c + d*x))**(n + 1)*sin(c + d*x)*cos(c + d*x)/(b*d*(n + 3)) - sqrt(2)*(a + b)*(a + b*sin(c + d*x))**n*(2*a**2 + b**2*(n + 2)**2)*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, -n - 1, sympy.S(3)/2, sympy.S.Half - sin(c + d*x)/2, b*(1 - sin(c + d*x))/(a + b))/(b**3*d*((a + b*sin(c + d*x))/(a + b))**n*(n + 2)*(n + 3)*sqrt(sin(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_220():
    f = (a + b*sin(c + d*x))**n*sin(c + d*x)**2
    F = sqrt(2)*a*(a + b)*(a + b*sin(c + d*x))**n*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, -n - 1, sympy.S(3)/2, sympy.S.Half - sin(c + d*x)/2, b*(1 - sin(c + d*x))/(a + b))/(b**2*d*((a + b*sin(c + d*x))/(a + b))**n*(n + 2)*sqrt(sin(c + d*x) + 1)) - (a + b*sin(c + d*x))**(n + 1)*cos(c + d*x)/(b*d*(n + 2)) - sqrt(2)*(a + b*sin(c + d*x))**n*(a**2 + b**2*(n + 1))*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, -n, sympy.S(3)/2, sympy.S.Half - sin(c + d*x)/2, b*(1 - sin(c + d*x))/(a + b))/(b**2*d*((a + b*sin(c + d*x))/(a + b))**n*(n + 2)*sqrt(sin(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_221():
    f = (a + b*sin(c + d*x))**n*sin(c + d*x)
    F = sqrt(2)*a*(a + b*sin(c + d*x))**n*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, -n, sympy.S(3)/2, sympy.S.Half - sin(c + d*x)/2, b*(1 - sin(c + d*x))/(a + b))/(b*d*((a + b*sin(c + d*x))/(a + b))**n*sqrt(sin(c + d*x) + 1)) - sqrt(2)*(a + b)*(a + b*sin(c + d*x))**n*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, -n - 1, sympy.S(3)/2, sympy.S.Half - sin(c + d*x)/2, b*(1 - sin(c + d*x))/(a + b))/(b*d*((a + b*sin(c + d*x))/(a + b))**n*sqrt(sin(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_222():
    f = (a + b*sin(c + d*x))**n
    F = -sqrt(2)*(a + b*sin(c + d*x))**n*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, -n, sympy.S(3)/2, sympy.S.Half - sin(c + d*x)/2, b*(1 - sin(c + d*x))/(a + b))/(d*((a + b*sin(c + d*x))/(a + b))**n*sqrt(sin(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_223():
    f = (a + b*sin(c + d*x))**n*csc(c + d*x)
    F = sympy.Function('Unintegrable')((sympy.csc((Symbol('c') + (Symbol('d') * x))) * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_224():
    f = (a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**4
    F = 7*a*c**4*x/8 + 7*a*c**4*sin(e + f*x)*cos(e + f*x)/(8*f) + 7*a*c**4*cos(e + f*x)**3/(12*f) + a*(-c**2*sin(e + f*x) + c**2)**2*cos(e + f*x)**3/(5*f) + 7*a*(-c**4*sin(e + f*x) + c**4)*cos(e + f*x)**3/(20*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_225():
    f = (a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**3
    F = 5*a*c**3*x/8 + 5*a*c**3*sin(e + f*x)*cos(e + f*x)/(8*f) + 5*a*c**3*cos(e + f*x)**3/(12*f) + a*(-c**3*sin(e + f*x) + c**3)*cos(e + f*x)**3/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_226():
    f = (a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**2
    F = a*c**2*x/2 + a*c**2*sin(e + f*x)*cos(e + f*x)/(2*f) + a*c**2*cos(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_227():
    f = (a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)
    F = a*c*x/2 + a*c*sin(e + f*x)*cos(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_228():
    f = (a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)
    F = 2*a*cos(e + f*x)/(f*(-c*sin(e + f*x) + c)) - a*x/c
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_229():
    f = (a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)**2
    F = a*c*cos(e + f*x)**3/(3*f*(-c*sin(e + f*x) + c)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_230():
    f = (a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)**3
    F = a*c*cos(e + f*x)**3/(5*f*(-c*sin(e + f*x) + c)**4) + a*cos(e + f*x)**3/(15*f*(-c*sin(e + f*x) + c)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_231():
    f = (a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)**4
    F = a*c*cos(e + f*x)**3/(7*f*(-c*sin(e + f*x) + c)**5) + 2*a*cos(e + f*x)**3/(35*f*(-c*sin(e + f*x) + c)**4) + 2*a*cos(e + f*x)**3/(105*c*f*(-c*sin(e + f*x) + c)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_232():
    f = (a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)**5
    F = 2*a*c*cos(e + f*x)**3/(315*f*(-c**2*sin(e + f*x) + c**2)**3) + a*c*cos(e + f*x)**3/(9*f*(-c*sin(e + f*x) + c)**6) + a*cos(e + f*x)**3/(21*f*(-c*sin(e + f*x) + c)**5) + 2*a*cos(e + f*x)**3/(105*c*f*(-c*sin(e + f*x) + c)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_233():
    f = (a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**5
    F = 9*a**2*c**5*x/16 + 3*a**2*c**5*sin(e + f*x)*cos(e + f*x)**3/(8*f) + 9*a**2*c**5*sin(e + f*x)*cos(e + f*x)/(16*f) + 3*a**2*c**5*cos(e + f*x)**5/(10*f) + a**2*c**3*(-c*sin(e + f*x) + c)**2*cos(e + f*x)**5/(7*f) + 3*a**2*(-c**5*sin(e + f*x) + c**5)*cos(e + f*x)**5/(14*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_234():
    f = (a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**4
    F = 7*a**2*c**4*x/16 + 7*a**2*c**4*sin(e + f*x)*cos(e + f*x)**3/(24*f) + 7*a**2*c**4*sin(e + f*x)*cos(e + f*x)/(16*f) + 7*a**2*c**4*cos(e + f*x)**5/(30*f) + a**2*(-c**4*sin(e + f*x) + c**4)*cos(e + f*x)**5/(6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_235():
    f = (a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**3
    F = 3*a**2*c**3*x/8 + a**2*c**3*sin(e + f*x)*cos(e + f*x)**3/(4*f) + 3*a**2*c**3*sin(e + f*x)*cos(e + f*x)/(8*f) + a**2*c**3*cos(e + f*x)**5/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_236():
    f = (a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**2
    F = 3*a**2*c**2*x/8 + a**2*c**2*sin(e + f*x)*cos(e + f*x)**3/(4*f) + 3*a**2*c**2*sin(e + f*x)*cos(e + f*x)/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_237():
    f = (a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)
    F = a**2*c*x/2 + a**2*c*sin(e + f*x)*cos(e + f*x)/(2*f) - a**2*c*cos(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_238():
    f = (a*sin(e + f*x) + a)**2/(-c*sin(e + f*x) + c)
    F = 2*a**2*c*cos(e + f*x)**3/(f*(-c*sin(e + f*x) + c)**2) - 3*a**2*x/c + 3*a**2*cos(e + f*x)/(c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_239():
    f = (a*sin(e + f*x) + a)**2/(-c*sin(e + f*x) + c)**2
    F = 2*a**2*c*cos(e + f*x)**3/(3*f*(-c*sin(e + f*x) + c)**3) - 2*a**2*cos(e + f*x)/(f*(-c**2*sin(e + f*x) + c**2)) + a**2*x/c**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_240():
    f = (a*sin(e + f*x) + a)**2/(-c*sin(e + f*x) + c)**3
    F = a**2*c**2*cos(e + f*x)**5/(5*f*(-c*sin(e + f*x) + c)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_241():
    f = (a*sin(e + f*x) + a)**2/(-c*sin(e + f*x) + c)**4
    F = a**2*c**2*cos(e + f*x)**5/(7*f*(-c*sin(e + f*x) + c)**6) + a**2*c*cos(e + f*x)**5/(35*f*(-c*sin(e + f*x) + c)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_242():
    f = (a*sin(e + f*x) + a)**2/(-c*sin(e + f*x) + c)**5
    F = a**2*c**2*cos(e + f*x)**5/(9*f*(-c*sin(e + f*x) + c)**7) + 2*a**2*c*cos(e + f*x)**5/(63*f*(-c*sin(e + f*x) + c)**6) + 2*a**2*cos(e + f*x)**5/(315*f*(-c*sin(e + f*x) + c)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_243():
    f = (a*sin(e + f*x) + a)**2/(-c*sin(e + f*x) + c)**6
    F = a**2*c**2*cos(e + f*x)**5/(11*f*(-c*sin(e + f*x) + c)**8) + a**2*c*cos(e + f*x)**5/(33*f*(-c*sin(e + f*x) + c)**7) + 2*a**2*cos(e + f*x)**5/(231*f*(-c*sin(e + f*x) + c)**6) + 2*a**2*cos(e + f*x)**5/(1155*c*f*(-c*sin(e + f*x) + c)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_244():
    f = (a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**6
    F = 55*a**3*c**6*x/128 + 11*a**3*c**6*sin(e + f*x)*cos(e + f*x)**5/(48*f) + 55*a**3*c**6*sin(e + f*x)*cos(e + f*x)**3/(192*f) + 55*a**3*c**6*sin(e + f*x)*cos(e + f*x)/(128*f) + 11*a**3*c**6*cos(e + f*x)**7/(56*f) + a**3*(-c**3*sin(e + f*x) + c**3)**2*cos(e + f*x)**7/(9*f) + 11*a**3*(-c**6*sin(e + f*x) + c**6)*cos(e + f*x)**7/(72*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_245():
    f = (a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**5
    F = 45*a**3*c**5*x/128 + 3*a**3*c**5*sin(e + f*x)*cos(e + f*x)**5/(16*f) + 15*a**3*c**5*sin(e + f*x)*cos(e + f*x)**3/(64*f) + 45*a**3*c**5*sin(e + f*x)*cos(e + f*x)/(128*f) + 9*a**3*c**5*cos(e + f*x)**7/(56*f) + a**3*(-c**5*sin(e + f*x) + c**5)*cos(e + f*x)**7/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_246():
    f = (a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**4
    F = 5*a**3*c**4*x/16 + a**3*c**4*sin(e + f*x)*cos(e + f*x)**5/(6*f) + 5*a**3*c**4*sin(e + f*x)*cos(e + f*x)**3/(24*f) + 5*a**3*c**4*sin(e + f*x)*cos(e + f*x)/(16*f) + a**3*c**4*cos(e + f*x)**7/(7*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_247():
    f = (a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**3
    F = 5*a**3*c**3*x/16 + a**3*c**3*sin(e + f*x)*cos(e + f*x)**5/(6*f) + 5*a**3*c**3*sin(e + f*x)*cos(e + f*x)**3/(24*f) + 5*a**3*c**3*sin(e + f*x)*cos(e + f*x)/(16*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_248():
    f = (a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**2
    F = 3*a**3*c**2*x/8 + a**3*c**2*sin(e + f*x)*cos(e + f*x)**3/(4*f) + 3*a**3*c**2*sin(e + f*x)*cos(e + f*x)/(8*f) - a**3*c**2*cos(e + f*x)**5/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_249():
    f = (a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)
    F = 5*a**3*c*x/8 + 5*a**3*c*sin(e + f*x)*cos(e + f*x)/(8*f) - 5*a**3*c*cos(e + f*x)**3/(12*f) - c*(a**3*sin(e + f*x) + a**3)*cos(e + f*x)**3/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_250():
    f = (a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)
    F = 2*a**3*c**2*cos(e + f*x)**5/(f*(-c*sin(e + f*x) + c)**3) + 5*a**3*cos(e + f*x)**3/(2*f*(-c*sin(e + f*x) + c)) - 15*a**3*x/(2*c) + 15*a**3*cos(e + f*x)/(2*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_251():
    f = (a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**2
    F = 2*a**3*c**2*cos(e + f*x)**5/(3*f*(-c*sin(e + f*x) + c)**4) - 10*a**3*cos(e + f*x)**3/(3*f*(-c*sin(e + f*x) + c)**2) + 5*a**3*x/c**2 - 5*a**3*cos(e + f*x)/(c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_252():
    f = (a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**3
    F = 2*a**3*c**2*cos(e + f*x)**5/(5*f*(-c*sin(e + f*x) + c)**5) + 2*a**3*cos(e + f*x)/(f*(-c**3*sin(e + f*x) + c**3)) - 2*a**3*cos(e + f*x)**3/(3*f*(-c*sin(e + f*x) + c)**3) - a**3*x/c**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_253():
    f = (a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**4
    F = a**3*c**3*cos(e + f*x)**7/(7*f*(-c*sin(e + f*x) + c)**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_254():
    f = (a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**5
    F = a**3*c**3*cos(e + f*x)**7/(9*f*(-c*sin(e + f*x) + c)**8) + a**3*c**2*cos(e + f*x)**7/(63*f*(-c*sin(e + f*x) + c)**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_255():
    f = (a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**6
    F = a**3*c**3*cos(e + f*x)**7/(11*f*(-c*sin(e + f*x) + c)**9) + 2*a**3*c**2*cos(e + f*x)**7/(99*f*(-c*sin(e + f*x) + c)**8) + 2*a**3*c*cos(e + f*x)**7/(693*f*(-c*sin(e + f*x) + c)**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_256():
    f = (a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**7
    F = a**3*c**3*cos(e + f*x)**7/(13*f*(-c*sin(e + f*x) + c)**10) + 3*a**3*c**2*cos(e + f*x)**7/(143*f*(-c*sin(e + f*x) + c)**9) + 2*a**3*c*cos(e + f*x)**7/(429*f*(-c*sin(e + f*x) + c)**8) + 2*a**3*cos(e + f*x)**7/(3003*f*(-c*sin(e + f*x) + c)**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_257():
    f = (a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**8
    F = a**3*c**3*cos(e + f*x)**7/(15*f*(-c*sin(e + f*x) + c)**11) + 4*a**3*c**2*cos(e + f*x)**7/(195*f*(-c*sin(e + f*x) + c)**10) + 4*a**3*c*cos(e + f*x)**7/(715*f*(-c*sin(e + f*x) + c)**9) + 8*a**3*cos(e + f*x)**7/(6435*f*(-c*sin(e + f*x) + c)**8) + 8*a**3*cos(e + f*x)**7/(45045*c*f*(-c*sin(e + f*x) + c)**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_258():
    f = (-c*sin(e + f*x) + c)**4/(a*sin(e + f*x) + a)
    F = -2*a**3*c**4*cos(e + f*x)**7/(f*(a*sin(e + f*x) + a)**4) - 14*a*c**4*cos(e + f*x)**5/(f*(a*sin(e + f*x) + a)**2) - 35*c**4*x/(2*a) - 35*c**4*sin(e + f*x)*cos(e + f*x)/(2*a*f) - 35*c**4*cos(e + f*x)**3/(3*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_259():
    f = (-c*sin(e + f*x) + c)**3/(a*sin(e + f*x) + a)
    F = -2*a**2*c**3*cos(e + f*x)**5/(f*(a*sin(e + f*x) + a)**3) - 5*c**3*cos(e + f*x)**3/(2*f*(a*sin(e + f*x) + a)) - 15*c**3*x/(2*a) - 15*c**3*cos(e + f*x)/(2*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_260():
    f = (-c*sin(e + f*x) + c)**2/(a*sin(e + f*x) + a)
    F = -2*a*c**2*cos(e + f*x)**3/(f*(a*sin(e + f*x) + a)**2) - 3*c**2*x/a - 3*c**2*cos(e + f*x)/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_261():
    f = (-c*sin(e + f*x) + c)/(a*sin(e + f*x) + a)
    F = -2*c*cos(e + f*x)/(f*(a*sin(e + f*x) + a)) - c*x/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_262():
    f = 1/((a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c))
    F = tan(e + f*x)/(a*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_263():
    f = 1/((a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**2)
    F = sec(e + f*x)/(3*a*f*(-c**2*sin(e + f*x) + c**2)) + 2*tan(e + f*x)/(3*a*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_264():
    f = 1/((a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**3)
    F = sec(e + f*x)/(5*a*f*(-c**3*sin(e + f*x) + c**3)) + sec(e + f*x)/(5*a*c*f*(-c*sin(e + f*x) + c)**2) + 2*tan(e + f*x)/(5*a*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_265():
    f = 1/((a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**4)
    F = 4*sec(e + f*x)/(35*a*f*(-c**4*sin(e + f*x) + c**4)) + 4*sec(e + f*x)/(35*a*f*(-c**2*sin(e + f*x) + c**2)**2) + sec(e + f*x)/(7*a*c*f*(-c*sin(e + f*x) + c)**3) + 8*tan(e + f*x)/(35*a*c**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_266():
    f = (-c*sin(e + f*x) + c)**5/(a*sin(e + f*x) + a)**2
    F = -2*a**4*c**5*cos(e + f*x)**9/(3*f*(a*sin(e + f*x) + a)**6) + 6*a**2*c**5*cos(e + f*x)**7/(f*(a*sin(e + f*x) + a)**4) + 42*c**5*cos(e + f*x)**5/(f*(a*sin(e + f*x) + a)**2) + 105*c**5*x/(2*a**2) + 105*c**5*sin(e + f*x)*cos(e + f*x)/(2*a**2*f) + 35*c**5*cos(e + f*x)**3/(a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_267():
    f = (-c*sin(e + f*x) + c)**4/(a*sin(e + f*x) + a)**2
    F = 14*a**4*c**4*cos(e + f*x)**5/(3*f*(a**2*sin(e + f*x) + a**2)**3) - 2*a**3*c**4*cos(e + f*x)**7/(3*f*(a*sin(e + f*x) + a)**5) + 35*c**4*cos(e + f*x)**3/(6*f*(a**2*sin(e + f*x) + a**2)) + 35*c**4*x/(2*a**2) + 35*c**4*cos(e + f*x)/(2*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_268():
    f = (-c*sin(e + f*x) + c)**3/(a*sin(e + f*x) + a)**2
    F = -2*a**2*c**3*cos(e + f*x)**5/(3*f*(a*sin(e + f*x) + a)**4) + 10*c**3*cos(e + f*x)**3/(3*f*(a*sin(e + f*x) + a)**2) + 5*c**3*x/a**2 + 5*c**3*cos(e + f*x)/(a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_269():
    f = (-c*sin(e + f*x) + c)**2/(a*sin(e + f*x) + a)**2
    F = -2*a*c**2*cos(e + f*x)**3/(3*f*(a*sin(e + f*x) + a)**3) + 2*c**2*cos(e + f*x)/(f*(a**2*sin(e + f*x) + a**2)) + c**2*x/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_270():
    f = (-c*sin(e + f*x) + c)/(a*sin(e + f*x) + a)**2
    F = -a*c*cos(e + f*x)**3/(3*f*(a*sin(e + f*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_271():
    f = 1/((a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c))
    F = -sec(e + f*x)/(3*c*f*(a**2*sin(e + f*x) + a**2)) + 2*tan(e + f*x)/(3*a**2*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_272():
    f = 1/((a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**2)
    F = tan(e + f*x)**3/(3*a**2*c**2*f) + tan(e + f*x)/(a**2*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_273():
    f = 1/((a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**3)
    F = sec(e + f*x)**3/(5*a**2*f*(-c**3*sin(e + f*x) + c**3)) + 4*tan(e + f*x)**3/(15*a**2*c**3*f) + 4*tan(e + f*x)/(5*a**2*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_274():
    f = 1/((a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**4)
    F = sec(e + f*x)**3/(7*a**2*f*(-c**4*sin(e + f*x) + c**4)) + sec(e + f*x)**3/(7*a**2*f*(-c**2*sin(e + f*x) + c**2)**2) + 4*tan(e + f*x)**3/(21*a**2*c**4*f) + 4*tan(e + f*x)/(7*a**2*c**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_275():
    f = 1/((a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**5)
    F = 2*sec(e + f*x)**3/(21*a**2*f*(-c**5*sin(e + f*x) + c**5)) + sec(e + f*x)**3/(9*a**2*c**2*f*(-c*sin(e + f*x) + c)**3) + 2*sec(e + f*x)**3/(21*a**2*c**3*f*(-c*sin(e + f*x) + c)**2) + 8*tan(e + f*x)**3/(63*a**2*c**5*f) + 8*tan(e + f*x)/(21*a**2*c**5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_276():
    f = (-c*sin(e + f*x) + c)**5/(a*sin(e + f*x) + a)**3
    F = -2*a**4*c**5*cos(e + f*x)**9/(5*f*(a*sin(e + f*x) + a)**7) + 6*a**2*c**5*cos(e + f*x)**7/(5*f*(a*sin(e + f*x) + a)**5) - 21*c**5*cos(e + f*x)**3/(2*f*(a**3*sin(e + f*x) + a**3)) - 42*c**5*cos(e + f*x)**5/(5*f*(a*sin(e + f*x) + a)**3) - 63*c**5*x/(2*a**3) - 63*c**5*cos(e + f*x)/(2*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_277():
    f = (-c*sin(e + f*x) + c)**4/(a*sin(e + f*x) + a)**3
    F = -2*a**3*c**4*cos(e + f*x)**7/(5*f*(a*sin(e + f*x) + a)**6) + 14*a*c**4*cos(e + f*x)**5/(15*f*(a*sin(e + f*x) + a)**4) - 14*c**4*cos(e + f*x)**3/(3*a*f*(a*sin(e + f*x) + a)**2) - 7*c**4*x/a**3 - 7*c**4*cos(e + f*x)/(a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_278():
    f = (-c*sin(e + f*x) + c)**3/(a*sin(e + f*x) + a)**3
    F = -2*a**2*c**3*cos(e + f*x)**5/(5*f*(a*sin(e + f*x) + a)**5) - 2*c**3*cos(e + f*x)/(f*(a**3*sin(e + f*x) + a**3)) + 2*c**3*cos(e + f*x)**3/(3*f*(a*sin(e + f*x) + a)**3) - c**3*x/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_279():
    f = (-c*sin(e + f*x) + c)**2/(a*sin(e + f*x) + a)**3
    F = -a**2*c**2*cos(e + f*x)**5/(5*f*(a*sin(e + f*x) + a)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_280():
    f = (-c*sin(e + f*x) + c)/(a*sin(e + f*x) + a)**3
    F = -a*c*cos(e + f*x)**3/(5*f*(a*sin(e + f*x) + a)**4) - c*cos(e + f*x)**3/(15*f*(a*sin(e + f*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_281():
    f = 1/((a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c))
    F = -sec(e + f*x)/(5*c*f*(a**3*sin(e + f*x) + a**3)) - sec(e + f*x)/(5*a*c*f*(a*sin(e + f*x) + a)**2) + 2*tan(e + f*x)/(5*a**3*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_282():
    f = 1/((a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**2)
    F = -sec(e + f*x)**3/(5*c**2*f*(a**3*sin(e + f*x) + a**3)) + 4*tan(e + f*x)**3/(15*a**3*c**2*f) + 4*tan(e + f*x)/(5*a**3*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_283():
    f = 1/((a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**3)
    F = tan(e + f*x)**5/(5*a**3*c**3*f) + 2*tan(e + f*x)**3/(3*a**3*c**3*f) + tan(e + f*x)/(a**3*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_284():
    f = 1/((a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**4)
    F = sec(e + f*x)**5/(7*a**3*f*(-c**4*sin(e + f*x) + c**4)) + 6*tan(e + f*x)**5/(35*a**3*c**4*f) + 4*tan(e + f*x)**3/(7*a**3*c**4*f) + 6*tan(e + f*x)/(7*a**3*c**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_285():
    f = 1/((a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**5)
    F = sec(e + f*x)**5/(9*a**3*f*(-c**5*sin(e + f*x) + c**5)) + sec(e + f*x)**5/(9*a**3*c**3*f*(-c*sin(e + f*x) + c)**2) + 2*tan(e + f*x)**5/(15*a**3*c**5*f) + 4*tan(e + f*x)**3/(9*a**3*c**5*f) + 2*tan(e + f*x)/(3*a**3*c**5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_286():
    f = 1/((a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**6)
    F = 8*sec(e + f*x)**5/(99*a**3*f*(-c**6*sin(e + f*x) + c**6)) + 8*sec(e + f*x)**5/(99*a**3*f*(-c**3*sin(e + f*x) + c**3)**2) + sec(e + f*x)**5/(11*a**3*f*(-c**2*sin(e + f*x) + c**2)**3) + 16*tan(e + f*x)**5/(165*a**3*c**6*f) + 32*tan(e + f*x)**3/(99*a**3*c**6*f) + 16*tan(e + f*x)/(33*a**3*c**6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_287():
    f = (a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = 256*a*c**5*cos(e + f*x)**3/(315*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 64*a*c**4*cos(e + f*x)**3/(105*f*sqrt(-c*sin(e + f*x) + c)) + 8*a*c**3*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)**3/(21*f) + 2*a*c**2*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)**3/(9*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_288():
    f = (a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = 64*a*c**4*cos(e + f*x)**3/(105*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 16*a*c**3*cos(e + f*x)**3/(35*f*sqrt(-c*sin(e + f*x) + c)) + 2*a*c**2*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)**3/(7*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_289():
    f = (a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = 8*a*c**3*cos(e + f*x)**3/(15*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 2*a*c**2*cos(e + f*x)**3/(5*f*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_290():
    f = (a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)
    F = 2*a*c**2*cos(e + f*x)**3/(3*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_291():
    f = (a*sin(e + f*x) + a)/sqrt(-c*sin(e + f*x) + c)
    F = -2*a*cos(e + f*x)/(f*sqrt(-c*sin(e + f*x) + c)) + 2*sqrt(2)*a*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_292():
    f = (a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = a*cos(e + f*x)/(f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - sqrt(2)*a*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(2*c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_293():
    f = (a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = a*cos(e + f*x)/(2*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) - a*cos(e + f*x)/(8*c*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - sqrt(2)*a*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(16*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_294():
    f = (a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = a*cos(e + f*x)/(3*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) - a*cos(e + f*x)/(24*c*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) - a*cos(e + f*x)/(32*c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - sqrt(2)*a*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(64*c**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_295():
    f = (a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = 256*a**2*c**6*cos(e + f*x)**5/(1155*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + 64*a**2*c**5*cos(e + f*x)**5/(231*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 8*a**2*c**4*cos(e + f*x)**5/(33*f*sqrt(-c*sin(e + f*x) + c)) + 2*a**2*c**3*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)**5/(11*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_296():
    f = (a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = 64*a**2*c**5*cos(e + f*x)**5/(315*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + 16*a**2*c**4*cos(e + f*x)**5/(63*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 2*a**2*c**3*cos(e + f*x)**5/(9*f*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_297():
    f = (a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = 8*a**2*c**4*cos(e + f*x)**5/(35*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + 2*a**2*c**3*cos(e + f*x)**5/(7*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_298():
    f = (a*sin(e + f*x) + a)**2*sqrt(-c*sin(e + f*x) + c)
    F = 2*a**2*c**3*cos(e + f*x)**5/(5*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_299():
    f = (a*sin(e + f*x) + a)**2/sqrt(-c*sin(e + f*x) + c)
    F = -2*a**2*c*cos(e + f*x)**3/(3*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - 4*a**2*cos(e + f*x)/(f*sqrt(-c*sin(e + f*x) + c)) + 4*sqrt(2)*a**2*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_300():
    f = (a*sin(e + f*x) + a)**2/(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = a**2*c*cos(e + f*x)**3/(f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + 3*a**2*cos(e + f*x)/(c*f*sqrt(-c*sin(e + f*x) + c)) - 3*sqrt(2)*a**2*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_301():
    f = (a*sin(e + f*x) + a)**2/(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = a**2*c*cos(e + f*x)**3/(2*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) - 3*a**2*cos(e + f*x)/(4*c*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 3*sqrt(2)*a**2*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(8*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_302():
    f = (a*sin(e + f*x) + a)**2/(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = a**2*c*cos(e + f*x)**3/(3*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)) - a**2*cos(e + f*x)/(4*c*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + a**2*cos(e + f*x)/(16*c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + sqrt(2)*a**2*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(32*c**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_303():
    f = (a*sin(e + f*x) + a)**2/(-c*sin(e + f*x) + c)**(sympy.S(9)/2)
    F = a**2*c*cos(e + f*x)**3/(4*f*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)) - a**2*cos(e + f*x)/(8*c*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) + a**2*cos(e + f*x)/(64*c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + 3*a**2*cos(e + f*x)/(256*c**3*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 3*sqrt(2)*a**2*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(512*c**(sympy.S(9)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_304():
    f = (a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = 256*a**3*c**7*cos(e + f*x)**7/(3003*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) + 64*a**3*c**6*cos(e + f*x)**7/(429*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + 24*a**3*c**5*cos(e + f*x)**7/(143*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 2*a**3*c**4*cos(e + f*x)**7/(13*f*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_305():
    f = (a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = 64*a**3*c**6*cos(e + f*x)**7/(693*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) + 16*a**3*c**5*cos(e + f*x)**7/(99*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + 2*a**3*c**4*cos(e + f*x)**7/(11*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_306():
    f = (a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = 8*a**3*c**5*cos(e + f*x)**7/(63*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) + 2*a**3*c**4*cos(e + f*x)**7/(9*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_307():
    f = (a*sin(e + f*x) + a)**3*sqrt(-c*sin(e + f*x) + c)
    F = 2*a**3*c**4*cos(e + f*x)**7/(7*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_308():
    f = (a*sin(e + f*x) + a)**3/sqrt(-c*sin(e + f*x) + c)
    F = -2*a**3*c**2*cos(e + f*x)**5/(5*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) - 4*a**3*c*cos(e + f*x)**3/(3*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - 8*a**3*cos(e + f*x)/(f*sqrt(-c*sin(e + f*x) + c)) + 8*sqrt(2)*a**3*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_309():
    f = (a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = a**3*c**2*cos(e + f*x)**5/(f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) + 5*a**3*cos(e + f*x)**3/(3*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 10*a**3*cos(e + f*x)/(c*f*sqrt(-c*sin(e + f*x) + c)) - 10*sqrt(2)*a**3*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_310():
    f = (a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = a**3*c**2*cos(e + f*x)**5/(2*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)) - 5*a**3*cos(e + f*x)**3/(4*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) - 15*a**3*cos(e + f*x)/(4*c**2*f*sqrt(-c*sin(e + f*x) + c)) + 15*sqrt(2)*a**3*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(4*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_311():
    f = (a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = a**3*c**2*cos(e + f*x)**5/(3*f*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)) - 5*a**3*cos(e + f*x)**3/(12*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) + 5*a**3*cos(e + f*x)/(8*c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - 5*sqrt(2)*a**3*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(16*c**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_312():
    f = (a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**(sympy.S(9)/2)
    F = a**3*c**2*cos(e + f*x)**5/(4*f*(-c*sin(e + f*x) + c)**(sympy.S(13)/2)) - 5*a**3*cos(e + f*x)**3/(24*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)) + 5*a**3*cos(e + f*x)/(32*c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) - 5*a**3*cos(e + f*x)/(128*c**3*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - 5*sqrt(2)*a**3*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(256*c**(sympy.S(9)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_313():
    f = (a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**(sympy.S(11)/2)
    F = a**3*c**2*cos(e + f*x)**5/(5*f*(-c*sin(e + f*x) + c)**(sympy.S(15)/2)) - a**3*cos(e + f*x)**3/(8*f*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)) + a**3*cos(e + f*x)/(16*c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) - a**3*cos(e + f*x)/(128*c**3*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) - 3*a**3*cos(e + f*x)/(512*c**4*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - 3*sqrt(2)*a**3*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(1024*c**(sympy.S(11)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_314():
    f = (-c*sin(e + f*x) + c)**(sympy.S(7)/2)/(a*sin(e + f*x) + a)
    F = -256*c**3*sqrt(-c*sin(e + f*x) + c)*sec(e + f*x)/(5*a*f) + 64*c**2*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)/(5*a*f) + 8*c*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)/(5*a*f) + 2*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)/(5*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_315():
    f = (-c*sin(e + f*x) + c)**(sympy.S(5)/2)/(a*sin(e + f*x) + a)
    F = -64*c**2*sqrt(-c*sin(e + f*x) + c)*sec(e + f*x)/(3*a*f) + 16*c*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)/(3*a*f) + 2*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)/(3*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_316():
    f = (-c*sin(e + f*x) + c)**(sympy.S(3)/2)/(a*sin(e + f*x) + a)
    F = -8*c*sqrt(-c*sin(e + f*x) + c)*sec(e + f*x)/(a*f) + 2*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_317():
    f = sqrt(-c*sin(e + f*x) + c)/(a*sin(e + f*x) + a)
    F = -2*sqrt(-c*sin(e + f*x) + c)*sec(e + f*x)/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_318():
    f = 1/((a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    F = -sqrt(-c*sin(e + f*x) + c)*sec(e + f*x)/(a*c*f) + sqrt(2)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(2*a*sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_319():
    f = 1/((a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    F = 3*cos(e + f*x)/(4*a*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - sec(e + f*x)/(a*c*f*sqrt(-c*sin(e + f*x) + c)) + 3*sqrt(2)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(8*a*c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_320():
    f = 1/((a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    F = 15*cos(e + f*x)/(32*a*c*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + sec(e + f*x)/(4*a*c*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - 5*sec(e + f*x)/(8*a*c**2*f*sqrt(-c*sin(e + f*x) + c)) + 15*sqrt(2)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(64*a*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_321():
    f = (-c*sin(e + f*x) + c)**(sympy.S(9)/2)/(a*sin(e + f*x) + a)**2
    F = 4096*c**3*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)**3/(15*a**2*f) - 1024*c**2*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)**3/(5*a**2*f) + 128*c*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)**3/(5*a**2*f) + 32*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*sec(e + f*x)**3/(15*a**2*f) + 2*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)*sec(e + f*x)**3/(5*a**2*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_322():
    f = (-c*sin(e + f*x) + c)**(sympy.S(7)/2)/(a*sin(e + f*x) + a)**2
    F = 256*c**2*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)**3/(3*a**2*f) - 64*c*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)**3/(a**2*f) + 8*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)**3/(a**2*f) + 2*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*sec(e + f*x)**3/(3*a**2*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_323():
    f = (-c*sin(e + f*x) + c)**(sympy.S(5)/2)/(a*sin(e + f*x) + a)**2
    F = 64*c*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)**3/(3*a**2*f) - 16*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)**3/(a**2*f) + 2*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)**3/(a**2*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_324():
    f = (-c*sin(e + f*x) + c)**(sympy.S(3)/2)/(a*sin(e + f*x) + a)**2
    F = 8*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)**3/(3*a**2*f) - 2*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)**3/(a**2*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_325():
    f = sqrt(-c*sin(e + f*x) + c)/(a*sin(e + f*x) + a)**2
    F = -2*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)**3/(3*a**2*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_326():
    f = 1/((a*sin(e + f*x) + a)**2*sqrt(-c*sin(e + f*x) + c))
    F = -sqrt(-c*sin(e + f*x) + c)*sec(e + f*x)/(2*a**2*c*f) - (-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)**3/(3*a**2*c**2*f) + sqrt(2)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(4*a**2*sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_327():
    f = 1/((a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    F = 5*cos(e + f*x)/(8*a**2*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - 5*sec(e + f*x)/(6*a**2*c*f*sqrt(-c*sin(e + f*x) + c)) - sqrt(-c*sin(e + f*x) + c)*sec(e + f*x)**3/(3*a**2*c**2*f) + 5*sqrt(2)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(16*a**2*c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_328():
    f = 1/((a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    F = 35*cos(e + f*x)/(64*a**2*c*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 7*sec(e + f*x)/(24*a**2*c*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - sec(e + f*x)**3/(3*a**2*c**2*f*sqrt(-c*sin(e + f*x) + c)) - 35*sec(e + f*x)/(48*a**2*c**2*f*sqrt(-c*sin(e + f*x) + c)) + 35*sqrt(2)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(128*a**2*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_329():
    f = (-c*sin(e + f*x) + c)**(sympy.S(9)/2)/(a*sin(e + f*x) + a)**3
    F = -4096*c**2*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)**5/(15*a**3*f) + 1024*c*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)**5/(3*a**3*f) - 128*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*sec(e + f*x)**5/(a**3*f) + 32*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)*sec(e + f*x)**5/(3*a**3*c*f) + 2*(-c*sin(e + f*x) + c)**(sympy.S(13)/2)*sec(e + f*x)**5/(3*a**3*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_330():
    f = (-c*sin(e + f*x) + c)**(sympy.S(7)/2)/(a*sin(e + f*x) + a)**3
    F = -256*c*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)**5/(5*a**3*f) + 64*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)**5/(a**3*f) - 24*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*sec(e + f*x)**5/(a**3*c*f) + 2*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)*sec(e + f*x)**5/(a**3*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_331():
    f = (-c*sin(e + f*x) + c)**(sympy.S(5)/2)/(a*sin(e + f*x) + a)**3
    F = -64*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)**5/(15*a**3*f) + 16*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)**5/(3*a**3*c*f) - 2*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*sec(e + f*x)**5/(a**3*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_332():
    f = (-c*sin(e + f*x) + c)**(sympy.S(3)/2)/(a*sin(e + f*x) + a)**3
    F = 8*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)**5/(15*a**3*c*f) - 2*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)**5/(3*a**3*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_333():
    f = sqrt(-c*sin(e + f*x) + c)/(a*sin(e + f*x) + a)**3
    F = -2*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)**5/(5*a**3*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_334():
    f = 1/((a*sin(e + f*x) + a)**3*sqrt(-c*sin(e + f*x) + c))
    F = -sqrt(-c*sin(e + f*x) + c)*sec(e + f*x)/(4*a**3*c*f) - (-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)**3/(6*a**3*c**2*f) - (-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)**5/(5*a**3*c**3*f) + sqrt(2)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(8*a**3*sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_335():
    f = 1/((a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    F = 7*cos(e + f*x)/(16*a**3*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - 7*sec(e + f*x)/(12*a**3*c*f*sqrt(-c*sin(e + f*x) + c)) - 7*sqrt(-c*sin(e + f*x) + c)*sec(e + f*x)**3/(30*a**3*c**2*f) - (-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)**5/(5*a**3*c**3*f) + 7*sqrt(2)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(32*a**3*c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_336():
    f = 1/((a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    F = 63*cos(e + f*x)/(128*a**3*c*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 21*sec(e + f*x)/(80*a**3*c*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - 3*sec(e + f*x)**3/(10*a**3*c**2*f*sqrt(-c*sin(e + f*x) + c)) - 21*sec(e + f*x)/(32*a**3*c**2*f*sqrt(-c*sin(e + f*x) + c)) - sqrt(-c*sin(e + f*x) + c)*sec(e + f*x)**5/(5*a**3*c**3*f) + 63*sqrt(2)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(256*a**3*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_337():
    f = sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = -a*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(4*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_338():
    f = sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = -a*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(3*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_339():
    f = sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = -a*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(2*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_340():
    f = sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)
    F = -a*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_341():
    f = sqrt(a*sin(e + f*x) + a)/sqrt(-c*sin(e + f*x) + c)
    F = -a*log(1 - sin(e + f*x))*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_342():
    f = sqrt(a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = a*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_343():
    f = sqrt(a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = a*cos(e + f*x)/(2*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_344():
    f = sqrt(a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = a*cos(e + f*x)/(3*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_345():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = -a**2*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(10*f*sqrt(a*sin(e + f*x) + a)) - a*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_346():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = -a**2*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(6*f*sqrt(a*sin(e + f*x) + a)) - a*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_347():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = -a**2*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(3*f*sqrt(a*sin(e + f*x) + a)) - a*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_348():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sin(e + f*x) + c)
    F = c*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(2*f*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_349():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)/sqrt(-c*sin(e + f*x) + c)
    F = -2*a**2*log(1 - sin(e + f*x))*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - a*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(f*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_350():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)/(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = a**2*log(1 - sin(e + f*x))*cos(e + f*x)/(c*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) + a*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_351():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)/(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = (a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(4*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_352():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)/(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = (a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(6*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) + (a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(24*c*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_353():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)/(-c*sin(e + f*x) + c)**(sympy.S(9)/2)
    F = -a**2*cos(e + f*x)/(12*c*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) + a*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(4*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_354():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)/(-c*sin(e + f*x) + c)**(sympy.S(11)/2)
    F = -a**2*cos(e + f*x)/(20*c*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)) + a*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(5*f*(-c*sin(e + f*x) + c)**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_355():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = -a**3*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(15*f*sqrt(a*sin(e + f*x) + a)) - 2*a**2*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(15*f) - a*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_356():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = -2*a**3*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(15*f*sqrt(a*sin(e + f*x) + a)) - a**2*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(5*f) - a*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_357():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = c**2*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(6*f*sqrt(-c*sin(e + f*x) + c)) + c*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_358():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-c*sin(e + f*x) + c)
    F = c*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(3*f*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_359():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)/sqrt(-c*sin(e + f*x) + c)
    F = -4*a**3*log(1 - sin(e + f*x))*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - 2*a**2*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(f*sqrt(-c*sin(e + f*x) + c)) - a*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(2*f*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_360():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)/(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = 4*a**3*log(1 - sin(e + f*x))*cos(e + f*x)/(c*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) + 2*a**2*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(c*f*sqrt(-c*sin(e + f*x) + c)) + a*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_361():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)/(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = -a**3*log(1 - sin(e + f*x))*cos(e + f*x)/(c**2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - a**2*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(c*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + a*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(2*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_362():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)/(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = (a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(6*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_363():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)/(-c*sin(e + f*x) + c)**(sympy.S(9)/2)
    F = (a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(8*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)) + (a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(48*c*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_364():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)/(-c*sin(e + f*x) + c)**(sympy.S(11)/2)
    F = (a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(10*f*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)) + (a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(40*c*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)) + (a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(240*c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_365():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)/(-c*sin(e + f*x) + c)**(sympy.S(13)/2)
    F = a**3*cos(e + f*x)/(60*c**2*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)) - a**2*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(15*c*f*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)) + a*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(6*f*(-c*sin(e + f*x) + c)**(sympy.S(13)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_366():
    f = (a*sin(e + f*x) + a)**(sympy.S(7)/2)*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)
    F = -a**4*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*cos(e + f*x)/(35*f*sqrt(a*sin(e + f*x) + a)) - a**3*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*cos(e + f*x)/(14*f) - 3*a**2*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*cos(e + f*x)/(28*f) - a*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*cos(e + f*x)/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_367():
    f = (a*sin(e + f*x) + a)**(sympy.S(7)/2)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = -2*a**4*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(35*f*sqrt(a*sin(e + f*x) + a)) - 4*a**3*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(35*f) - a**2*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(7*f) - a*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(7*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_368():
    f = (a*sin(e + f*x) + a)**(sympy.S(7)/2)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = c**3*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(15*f*sqrt(-c*sin(e + f*x) + c)) + 2*c**2*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(15*f) + c*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_369():
    f = (a*sin(e + f*x) + a)**(sympy.S(7)/2)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = c**2*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(10*f*sqrt(-c*sin(e + f*x) + c)) + c*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_370():
    f = (a*sin(e + f*x) + a)**(sympy.S(7)/2)*sqrt(-c*sin(e + f*x) + c)
    F = c*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(4*f*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_371():
    f = (a*sin(e + f*x) + a)**(sympy.S(7)/2)/sqrt(-c*sin(e + f*x) + c)
    F = -8*a**4*log(1 - sin(e + f*x))*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - 4*a**3*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(f*sqrt(-c*sin(e + f*x) + c)) - a**2*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(f*sqrt(-c*sin(e + f*x) + c)) - a*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(3*f*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_372():
    f = (a*sin(e + f*x) + a)**(sympy.S(7)/2)/(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = 12*a**4*log(1 - sin(e + f*x))*cos(e + f*x)/(c*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) + 6*a**3*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(c*f*sqrt(-c*sin(e + f*x) + c)) + 3*a**2*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(2*c*f*sqrt(-c*sin(e + f*x) + c)) + a*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_373():
    f = (a*sin(e + f*x) + a)**(sympy.S(7)/2)/(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = -6*a**4*log(1 - sin(e + f*x))*cos(e + f*x)/(c**2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - 3*a**3*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(c**2*f*sqrt(-c*sin(e + f*x) + c)) - 3*a**2*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(2*c*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + a*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(2*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_374():
    f = (a*sin(e + f*x) + a)**(sympy.S(7)/2)/(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = a**4*log(1 - sin(e + f*x))*cos(e + f*x)/(c**3*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) + a**3*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - a**2*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(2*c*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + a*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(3*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_375():
    f = (a*sin(e + f*x) + a)**(sympy.S(7)/2)/(-c*sin(e + f*x) + c)**(sympy.S(9)/2)
    F = (a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(8*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_376():
    f = (a*sin(e + f*x) + a)**(sympy.S(7)/2)/(-c*sin(e + f*x) + c)**(sympy.S(11)/2)
    F = (a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(10*f*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)) + (a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(80*c*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_377():
    f = (a*sin(e + f*x) + a)**(sympy.S(7)/2)/(-c*sin(e + f*x) + c)**(sympy.S(13)/2)
    F = (a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(12*f*(-c*sin(e + f*x) + c)**(sympy.S(13)/2)) + (a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(60*c*f*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)) + (a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(480*c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_378():
    f = (a*sin(e + f*x) + a)**(sympy.S(7)/2)/(-c*sin(e + f*x) + c)**(sympy.S(15)/2)
    F = (a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(14*f*(-c*sin(e + f*x) + c)**(sympy.S(15)/2)) + (a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(56*c*f*(-c*sin(e + f*x) + c)**(sympy.S(13)/2)) + (a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(280*c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)) + (a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(2240*c**3*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_379():
    f = (a*sin(e + f*x) + a)**(sympy.S(7)/2)/(-c*sin(e + f*x) + c)**(sympy.S(17)/2)
    F = -a**4*cos(e + f*x)/(280*c**3*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)) + a**3*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(56*c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(13)/2)) - 3*a**2*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(56*c*f*(-c*sin(e + f*x) + c)**(sympy.S(15)/2)) + a*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(8*f*(-c*sin(e + f*x) + c)**(sympy.S(17)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_380():
    f = (-c*sin(e + f*x) + c)**(sympy.S(5)/2)/sqrt(a*sin(e + f*x) + a)
    F = 4*c**3*log(sin(e + f*x) + 1)*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) + 2*c**2*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)) + c*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(2*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_381():
    f = (-c*sin(e + f*x) + c)**(sympy.S(3)/2)/sqrt(a*sin(e + f*x) + a)
    F = 2*c**2*log(sin(e + f*x) + 1)*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) + c*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_382():
    f = sqrt(-c*sin(e + f*x) + c)/sqrt(a*sin(e + f*x) + a)
    F = c*log(sin(e + f*x) + 1)*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_383():
    f = 1/(sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    F = cos(e + f*x)*atanh(sin(e + f*x))/(f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_384():
    f = 1/(sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    F = cos(e + f*x)/(2*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + cos(e + f*x)*atanh(sin(e + f*x))/(2*c*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_385():
    f = 1/(sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    F = cos(e + f*x)/(4*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + cos(e + f*x)/(4*c*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + cos(e + f*x)*atanh(sin(e + f*x))/(4*c**2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_386():
    f = (-c*sin(e + f*x) + c)**(sympy.S(7)/2)/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -c*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - 12*c**4*log(sin(e + f*x) + 1)*cos(e + f*x)/(a*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - 6*c**3*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(a*f*sqrt(a*sin(e + f*x) + a)) - 3*c**2*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(2*a*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_387():
    f = (-c*sin(e + f*x) + c)**(sympy.S(5)/2)/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -c*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - 4*c**3*log(sin(e + f*x) + 1)*cos(e + f*x)/(a*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - 2*c**2*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(a*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_388():
    f = (-c*sin(e + f*x) + c)**(sympy.S(3)/2)/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -c*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - c**2*log(sin(e + f*x) + 1)*cos(e + f*x)/(a*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_389():
    f = sqrt(-c*sin(e + f*x) + c)/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -c*cos(e + f*x)/(f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_390():
    f = 1/((a*sin(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sin(e + f*x) + c))
    F = -cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sin(e + f*x) + c)) + cos(e + f*x)*atanh(sin(e + f*x))/(2*a*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_391():
    f = 1/((a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    F = -cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + cos(e + f*x)/(2*a*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + cos(e + f*x)*atanh(sin(e + f*x))/(2*a*c*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_392():
    f = 1/((a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    F = -cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + 3*cos(e + f*x)/(8*a*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + 3*cos(e + f*x)/(8*a*c*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 3*cos(e + f*x)*atanh(sin(e + f*x))/(8*a*c**2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_393():
    f = (-c*sin(e + f*x) + c)**(sympy.S(9)/2)/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -c*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) + 2*c**2*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + 24*c**5*log(sin(e + f*x) + 1)*cos(e + f*x)/(a**2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) + 12*c**4*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(a**2*f*sqrt(a*sin(e + f*x) + a)) + 3*c**3*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(a**2*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_394():
    f = (-c*sin(e + f*x) + c)**(sympy.S(7)/2)/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -c*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) + 3*c**2*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(2*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + 6*c**4*log(sin(e + f*x) + 1)*cos(e + f*x)/(a**2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) + 3*c**3*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(a**2*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_395():
    f = (-c*sin(e + f*x) + c)**(sympy.S(5)/2)/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -c*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) + c**2*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + c**3*log(sin(e + f*x) + 1)*cos(e + f*x)/(a**2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_396():
    f = (-c*sin(e + f*x) + c)**(sympy.S(3)/2)/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_397():
    f = sqrt(-c*sin(e + f*x) + c)/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -c*cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_398():
    f = 1/((a*sin(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-c*sin(e + f*x) + c))
    F = -cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-c*sin(e + f*x) + c)) - cos(e + f*x)/(4*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sin(e + f*x) + c)) + cos(e + f*x)*atanh(sin(e + f*x))/(4*a**2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_399():
    f = 1/((a*sin(e + f*x) + a)**(sympy.S(5)/2)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    F = -cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - 3*cos(e + f*x)/(8*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 3*cos(e + f*x)/(8*a**2*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 3*cos(e + f*x)*atanh(sin(e + f*x))/(8*a**2*c*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_400():
    f = 1/((a*sin(e + f*x) + a)**(sympy.S(5)/2)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    F = -cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) - cos(e + f*x)/(2*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + 3*cos(e + f*x)/(8*a**2*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + 3*cos(e + f*x)/(8*a**2*c*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 3*cos(e + f*x)*atanh(sin(e + f*x))/(8*a**2*c**2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_401():
    f = (a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**n
    F = 2**(n + sympy.S.Half)*c*(1 - sin(e + f*x))**(sympy.S.Half - n)*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(n - 1)*cos(e + f*x)*hyper((m + sympy.S.Half, sympy.S.Half - n), (m + sympy.S(3)/2,), sin(e + f*x)/2 + sympy.S.Half)/(f*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_402():
    f = (a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**3
    F = -2**(m + sympy.S.Half)*a**4*c**3*(a*sin(e + f*x) + a)**(m - 4)*(sin(e + f*x) + 1)**(sympy.S.Half - m)*cos(e + f*x)**7*hyper((sympy.S(7)/2, sympy.S.Half - m), (sympy.S(9)/2,), sympy.S.Half - sin(e + f*x)/2)/(7*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_403():
    f = (a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**2
    F = -2**(m + sympy.S.Half)*a**3*c**2*(a*sin(e + f*x) + a)**(m - 3)*(sin(e + f*x) + 1)**(sympy.S.Half - m)*cos(e + f*x)**5*hyper((sympy.S(5)/2, sympy.S.Half - m), (sympy.S(7)/2,), sympy.S.Half - sin(e + f*x)/2)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_404():
    f = (a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)
    F = -2**(m + sympy.S.Half)*a**2*c*(a*sin(e + f*x) + a)**(m - 2)*(sin(e + f*x) + 1)**(sympy.S.Half - m)*cos(e + f*x)**3*hyper((sympy.S(3)/2, sympy.S.Half - m), (sympy.S(5)/2,), sympy.S.Half - sin(e + f*x)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_405():
    f = (a*sin(e + f*x) + a)**m/(-c*sin(e + f*x) + c)
    F = 2**(m + sympy.S.Half)*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(sympy.S.Half - m)*hyper((sympy.S(-1)/2, sympy.S.Half - m), (sympy.S.Half,), sympy.S.Half - sin(e + f*x)/2)*sec(e + f*x)/(c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_406():
    f = (a*sin(e + f*x) + a)**m/(-c*sin(e + f*x) + c)**2
    F = 2**(m + sympy.S.Half)*(a*sin(e + f*x) + a)**(m + 1)*(sin(e + f*x) + 1)**(sympy.S.Half - m)*hyper((sympy.S(-3)/2, sympy.S.Half - m), (sympy.S(-1)/2,), sympy.S.Half - sin(e + f*x)/2)*sec(e + f*x)**3/(3*a*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_407():
    f = (a*sin(e + f*x) + a)**m/(-c*sin(e + f*x) + c)**3
    F = 2**(m + sympy.S.Half)*(a*sin(e + f*x) + a)**(m + 2)*(sin(e + f*x) + 1)**(sympy.S.Half - m)*hyper((sympy.S(-5)/2, sympy.S.Half - m), (sympy.S(-3)/2,), sympy.S.Half - sin(e + f*x)/2)*sec(e + f*x)**5/(5*a**2*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_408():
    f = (a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = 64*c**3*(a*sin(e + f*x) + a)**m*cos(e + f*x)/(f*(2*m + 5)*sqrt(-c*sin(e + f*x) + c)*(4*m**2 + 8*m + 3)) + 16*c**2*(a*sin(e + f*x) + a)**m*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(f*(4*m**2 + 16*m + 15)) + 2*c*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(f*(2*m + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_409():
    f = (a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = 8*c**2*(a*sin(e + f*x) + a)**m*cos(e + f*x)/(f*sqrt(-c*sin(e + f*x) + c)*(4*m**2 + 8*m + 3)) + 2*c*(a*sin(e + f*x) + a)**m*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(f*(2*m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_410():
    f = (a*sin(e + f*x) + a)**m*sqrt(-c*sin(e + f*x) + c)
    F = 2*c*(a*sin(e + f*x) + a)**m*cos(e + f*x)/(f*(2*m + 1)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_411():
    f = (a*sin(e + f*x) + a)**m/sqrt(-c*sin(e + f*x) + c)
    F = (a*sin(e + f*x) + a)**m*cos(e + f*x)*hyper((1, m + sympy.S.Half), (m + sympy.S(3)/2,), sin(e + f*x)/2 + sympy.S.Half)/(f*(2*m + 1)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_412():
    f = (a*sin(e + f*x) + a)**m/(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = (a*sin(e + f*x) + a)**m*cos(e + f*x)*hyper((2, m + sympy.S.Half), (m + sympy.S(3)/2,), sin(e + f*x)/2 + sympy.S.Half)/(2*c*f*(2*m + 1)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_413():
    f = (a*sin(e + f*x) + a)**m/(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = (a*sin(e + f*x) + a)**m*cos(e + f*x)*hyper((3, m + sympy.S.Half), (m + sympy.S(3)/2,), sin(e + f*x)/2 + sympy.S.Half)/(4*c**2*f*(2*m + 1)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_414():
    f = (a*sin(e + f*x) + a)**m/sqrt(-c*sin(e + f*x) + c)
    F = (a*sin(e + f*x) + a)**m*cos(e + f*x)*hyper((1, m + sympy.S.Half), (m + sympy.S(3)/2,), sin(e + f*x)/2 + sympy.S.Half)/(f*(2*m + 1)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_415():
    f = (c*sin(e + f*x) + c)**m/sqrt(-a*sin(e + f*x) + a)
    F = (c*sin(e + f*x) + c)**m*cos(e + f*x)*hyper((1, m + sympy.S.Half), (m + sympy.S(3)/2,), sin(e + f*x)/2 + sympy.S.Half)/(f*(2*m + 1)*sqrt(-a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_416():
    f = (a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 3)
    F = (a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 3)*cos(e + f*x)/(f*(2*m + 5)) + 2*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 2)*cos(e + f*x)/(c*f*(4*m**2 + 16*m + 15)) + 2*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 1)*cos(e + f*x)/(c**2*f*(2*m + 5)*(4*m**2 + 8*m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_417():
    f = (a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 2)
    F = (a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 2)*cos(e + f*x)/(f*(2*m + 3)) + (a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 1)*cos(e + f*x)/(c*f*(4*m**2 + 8*m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_418():
    f = (a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 1)
    F = (a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 1)*cos(e + f*x)/(f*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_419():
    f = (a*sin(e + f*x) + a)**m/(-c*sin(e + f*x) + c)**m
    F = 2**(sympy.S.Half - m)*c*(1 - sin(e + f*x))**(m + sympy.S.Half)*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 1)*cos(e + f*x)*hyper((m + sympy.S.Half, m + sympy.S.Half), (m + sympy.S(3)/2,), sin(e + f*x)/2 + sympy.S.Half)/(f*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_420():
    f = (a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(1 - m)
    F = 2**(sympy.S(3)/2 - m)*c**2*(1 - sin(e + f*x))**(m + sympy.S.Half)*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 1)*cos(e + f*x)*hyper((m + sympy.S(-1)/2, m + sympy.S.Half), (m + sympy.S(3)/2,), sin(e + f*x)/2 + sympy.S.Half)/(f*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_421():
    f = (a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(2 - m)
    F = 2**(sympy.S(5)/2 - m)*c**3*(1 - sin(e + f*x))**(m + sympy.S.Half)*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 1)*cos(e + f*x)*hyper((m + sympy.S(-3)/2, m + sympy.S.Half), (m + sympy.S(3)/2,), sin(e + f*x)/2 + sympy.S.Half)/(f*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_422():
    f = (c + d*sin(e + f*x))**4*(a*sin(e + f*x) + a)
    F = -a*d*(24*c**3 + 130*c**2*d + 116*c*d**2 + 45*d**3)*sin(e + f*x)*cos(e + f*x)/(120*f) + a*x*(8*c**4 + 16*c**3*d + 24*c**2*d**2 + 12*c*d**3 + 3*d**4)/8 - a*(c + d*sin(e + f*x))**4*cos(e + f*x)/(5*f) - a*(c + d*sin(e + f*x))**3*(4*c + 5*d)*cos(e + f*x)/(20*f) - a*(c + d*sin(e + f*x))**2*(12*c**2 + 35*c*d + 16*d**2)*cos(e + f*x)/(60*f) - a*(12*c**4 + 95*c**3*d + 112*c**2*d**2 + 80*c*d**3 + 16*d**4)*cos(e + f*x)/(30*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_423():
    f = (c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a)
    F = -a*d*(6*c**2 + 20*c*d + 9*d**2)*sin(e + f*x)*cos(e + f*x)/(24*f) + a*x*(8*c**3 + 12*c**2*d + 12*c*d**2 + 3*d**3)/8 - a*(c + d*sin(e + f*x))**3*cos(e + f*x)/(4*f) - a*(c + d*sin(e + f*x))**2*(3*c + 4*d)*cos(e + f*x)/(12*f) - a*(3*c**3 + 16*c**2*d + 12*c*d**2 + 4*d**3)*cos(e + f*x)/(6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_424():
    f = (c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)
    F = -a*d*(2*c + 3*d)*sin(e + f*x)*cos(e + f*x)/(6*f) + a*x*(2*c**2 + 2*c*d + d**2)/2 - a*(c + d*sin(e + f*x))**2*cos(e + f*x)/(3*f) - 2*a*(c**2 + 3*c*d + d**2)*cos(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_425():
    f = (c + d*sin(e + f*x))*(a*sin(e + f*x) + a)
    F = -a*d*sin(e + f*x)*cos(e + f*x)/(2*f) + a*x*(2*c + d)/2 - a*(c + d)*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_426():
    f = a*sin(e + f*x) + a
    F = a*x - a*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_427():
    f = (a*sin(e + f*x) + a)/(c + d*sin(e + f*x))
    F = a*x/d - 2*a*(c - d)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(d*f*sqrt(c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_428():
    f = (a*sin(e + f*x) + a)/(c + d*sin(e + f*x))**2
    F = 2*a*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*(c + d)*sqrt(c**2 - d**2)) - a*cos(e + f*x)/(f*(c + d)*(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_429():
    f = (a*sin(e + f*x) + a)/(c + d*sin(e + f*x))**3
    F = -a*(c - 2*d)*cos(e + f*x)/(f*(c + d)**2*(c + d*sin(e + f*x))*(2*c - 2*d)) - a*cos(e + f*x)/(f*(c + d*sin(e + f*x))**2*(2*c + 2*d)) + a*(2*c - d)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*(c + d)*(c**2 - d**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_430():
    f = (a*sin(e + f*x) + a)/(c + d*sin(e + f*x))**4
    F = -a*(c - 4*d)*(2*c - d)*cos(e + f*x)/(6*f*(c - d)**2*(c + d)**3*(c + d*sin(e + f*x))) - a*cos(e + f*x)/(f*(c + d*sin(e + f*x))**3*(3*c + 3*d)) + a*(2*c**2 - 2*c*d + d**2)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*(c + d)*(c**2 - d**2)**(sympy.S(5)/2)) - a*(2*c - 3*d)*cos(e + f*x)/(f*(c + d)**2*(c + d*sin(e + f*x))**2*(6*c - 6*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_431():
    f = (c + d*sin(e + f*x))**4*(a*sin(e + f*x) + a)**2
    F = a**2*x*(24*c**4 + 64*c**3*d + 84*c**2*d**2 + 48*c*d**3 + 11*d**4)/16 + a**2*(8*c**4 - 96*c**3*d - 438*c**2*d**2 - 464*c*d**3 - 165*d**4)*sin(e + f*x)*cos(e + f*x)/(240*f) + a**2*(c - 12*d)*(c + d*sin(e + f*x))**4*cos(e + f*x)/(30*d*f) - a**2*(c + d*sin(e + f*x))**5*cos(e + f*x)/(6*d*f) + a**2*(c + d*sin(e + f*x))**3*(4*c**2 - 48*c*d - 55*d**2)*cos(e + f*x)/(120*d*f) + a**2*(c + d*sin(e + f*x))**2*(4*c**3 - 48*c**2*d - 123*c*d**2 - 64*d**3)*cos(e + f*x)/(120*d*f) + a**2*(4*c**5 - 48*c**4*d - 311*c**3*d**2 - 448*c**2*d**3 - 288*c*d**4 - 64*d**5)*cos(e + f*x)/(60*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_432():
    f = (c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a)**2
    F = 3*a**2*x*(2*c + d)*(2*c**2 + 3*c*d + 2*d**2)/8 + a**2*(2*c**3 - 20*c**2*d - 57*c*d**2 - 30*d**3)*sin(e + f*x)*cos(e + f*x)/(40*f) + a**2*(c - 10*d)*(c + d*sin(e + f*x))**3*cos(e + f*x)/(20*d*f) - a**2*(c + d*sin(e + f*x))**4*cos(e + f*x)/(5*d*f) + a**2*(c + d*sin(e + f*x))**2*(c**2 - 10*c*d - 12*d**2)*cos(e + f*x)/(20*d*f) + a**2*(c**4 - 10*c**3*d - 44*c**2*d**2 - 40*c*d**3 - 12*d**4)*cos(e + f*x)/(10*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_433():
    f = (c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**2
    F = a**2*x*(12*c**2 + 16*c*d + 7*d**2)/8 - a**2*(12*c**2 + 16*c*d + 7*d**2)*sin(e + f*x)*cos(e + f*x)/(24*f) - a**2*(12*c**2 + 16*c*d + 7*d**2)*cos(e + f*x)/(6*f) - d*(8*c - d)*(a*sin(e + f*x) + a)**2*cos(e + f*x)/(12*f) - d**2*(a*sin(e + f*x) + a)**3*cos(e + f*x)/(4*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_434():
    f = (c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**2
    F = a**2*x*(3*c + 2*d)/2 - a**2*(3*c + 2*d)*sin(e + f*x)*cos(e + f*x)/(6*f) - 2*a**2*(3*c + 2*d)*cos(e + f*x)/(3*f) - d*(a*sin(e + f*x) + a)**2*cos(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_435():
    f = (a*sin(e + f*x) + a)**2
    F = 3*a**2*x/2 - a**2*sin(e + f*x)*cos(e + f*x)/(2*f) - 2*a**2*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_436():
    f = (a*sin(e + f*x) + a)**2/(c + d*sin(e + f*x))
    F = -a**2*cos(e + f*x)/(d*f) - a**2*x*(c - 2*d)/d**2 + 2*a**2*(c - d)**2*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(d**2*f*sqrt(c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_437():
    f = (a*sin(e + f*x) + a)**2/(c + d*sin(e + f*x))**3
    F = 3*a**2*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*(c + d)**2*sqrt(c**2 - d**2)) + a**2*(c - d)*cos(e + f*x)/(2*d*f*(c + d)*(c + d*sin(e + f*x))**2) - a**2*(c + 4*d)*cos(e + f*x)/(2*d*f*(c + d)**2*(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_438():
    f = (a*sin(e + f*x) + a)**2/(c + d*sin(e + f*x))**4
    F = a**2*(3*c - 2*d)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*(c - d)*(c + d)**3*sqrt(c**2 - d**2)) + a**2*(c - d)*cos(e + f*x)/(3*d*f*(c + d)*(c + d*sin(e + f*x))**3) - a**2*(c + 6*d)*cos(e + f*x)/(6*d*f*(c + d)**2*(c + d*sin(e + f*x))**2) - a**2*(c**2 + 6*c*d - 10*d**2)*cos(e + f*x)/(d*f*(c + d)**3*(c + d*sin(e + f*x))*(6*c - 6*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_439():
    f = (a*sin(e + f*x) + a)**2/(c + d*sin(e + f*x))**5
    F = a**2*(12*c**2 - 16*c*d + 7*d**2)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(4*f*(c - d)**2*(c + d)**4*sqrt(c**2 - d**2)) + a**2*(c - d)*cos(e + f*x)/(4*d*f*(c + d)*(c + d*sin(e + f*x))**4) - a**2*(c + 8*d)*cos(e + f*x)/(12*d*f*(c + d)**2*(c + d*sin(e + f*x))**3) - a**2*(2*c**2 + 16*c*d - 21*d**2)*cos(e + f*x)/(d*f*(c + d)**3*(c + d*sin(e + f*x))**2*(24*c - 24*d)) - a**2*(2*c**3 + 16*c**2*d - 59*c*d**2 + 32*d**3)*cos(e + f*x)/(24*d*f*(c - d)**2*(c + d)**4*(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_440():
    f = (a*sin(e + f*x) + a)**3
    F = 5*a**3*x/2 - 3*a**3*sin(e + f*x)*cos(e + f*x)/(2*f) + a**3*cos(e + f*x)**3/(3*f) - 4*a**3*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_441():
    f = (a*sin(e + f*x) + a)**3/(c + d*sin(e + f*x))
    F = a**3*(2*c - 5*d)*cos(e + f*x)/(2*d**2*f) + a**3*x*(2*c**2 - 6*c*d + 7*d**2)/(2*d**3) - 2*a**3*(c - d)**3*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(d**3*f*sqrt(c**2 - d**2)) - (a**3*sin(e + f*x) + a**3)*cos(e + f*x)/(2*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_442():
    f = (a*sin(e + f*x) + a)**3/(c + d*sin(e + f*x))**2
    F = -2*a**3*c*cos(e + f*x)/(d**2*f*(c + d)) - a**3*x*(2*c - 3*d)/d**3 + 2*a**3*(c - d)**2*(2*c + 3*d)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(d**3*f*(c + d)*sqrt(c**2 - d**2)) + (c - d)*(a**3*sin(e + f*x) + a**3)*cos(e + f*x)/(d*f*(c + d)*(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_443():
    f = (a*sin(e + f*x) + a)**3/(c + d*sin(e + f*x))**3
    F = a**3*(c - d)*(2*c + 5*d)*cos(e + f*x)/(2*d**2*f*(c + d)**2*(c + d*sin(e + f*x))) + a**3*x/d**3 - a**3*(c - d)*(2*c**2 + 6*c*d + 7*d**2)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(d**3*f*(c + d)**2*sqrt(c**2 - d**2)) + (c - d)*(a**3*sin(e + f*x) + a**3)*cos(e + f*x)/(2*d*f*(c + d)*(c + d*sin(e + f*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_444():
    f = (a*sin(e + f*x) + a)**3/(c + d*sin(e + f*x))**4
    F = 5*a**3*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*(c + d)**3*sqrt(c**2 - d**2)) + a**3*(c - d)*(2*c + 7*d)*cos(e + f*x)/(6*d**2*f*(c + d)**2*(c + d*sin(e + f*x))**2) - a**3*(2*c**2 + 9*c*d + 22*d**2)*cos(e + f*x)/(6*d**2*f*(c + d)**3*(c + d*sin(e + f*x))) + (c - d)*(a**3*sin(e + f*x) + a**3)*cos(e + f*x)/(3*d*f*(c + d)*(c + d*sin(e + f*x))**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_445():
    f = (a*sin(e + f*x) + a)**3/(c + d*sin(e + f*x))**5
    F = 5*a**3*(4*c - 3*d)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*(c + d)**4*(4*c - 4*d)*sqrt(c**2 - d**2)) + a**3*(c - d)*(2*c + 9*d)*cos(e + f*x)/(12*d**2*f*(c + d)**2*(c + d*sin(e + f*x))**3) - a**3*(2*c**2 + 12*c*d + 45*d**2)*cos(e + f*x)/(24*d**2*f*(c + d)**3*(c + d*sin(e + f*x))**2) - a**3*(2*c**3 + 12*c**2*d + 43*c*d**2 - 72*d**3)*cos(e + f*x)/(d**2*f*(c + d)**4*(c + d*sin(e + f*x))*(24*c - 24*d)) + (c - d)*(a**3*sin(e + f*x) + a**3)*cos(e + f*x)/(4*d*f*(c + d)*(c + d*sin(e + f*x))**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_446():
    f = (c + d*sin(e + f*x))**4/(a*sin(e + f*x) + a)
    F = -(c - d)*(c + d*sin(e + f*x))**3*cos(e + f*x)/(f*(a*sin(e + f*x) + a)) + d**2*(6*c**2 - 20*c*d + 9*d**2)*sin(e + f*x)*cos(e + f*x)/(6*a*f) + d*x*(8*c**3 - 12*c**2*d + 12*c*d**2 - 3*d**3)/(2*a) + d*(c + d*sin(e + f*x))**2*(3*c - 4*d)*cos(e + f*x)/(3*a*f) + 2*d*(3*c**3 - 16*c**2*d + 12*c*d**2 - 4*d**3)*cos(e + f*x)/(3*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_447():
    f = (c + d*sin(e + f*x))**3/(a*sin(e + f*x) + a)
    F = -(c - d)*(c + d*sin(e + f*x))**2*cos(e + f*x)/(f*(a*sin(e + f*x) + a)) + d**2*(2*c - 3*d)*sin(e + f*x)*cos(e + f*x)/(2*a*f) + 3*d*x*(2*c**2 - 2*c*d + d**2)/(2*a) + 2*d*(c**2 - 3*c*d + d**2)*cos(e + f*x)/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_448():
    f = (c + d*sin(e + f*x))**2/(a*sin(e + f*x) + a)
    F = -d**2*cos(e + f*x)/(a*f) + d*x*(2*c - d)/a - (c - d)**2*cos(e + f*x)/(a*f*(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_449():
    f = (c + d*sin(e + f*x))/(a*sin(e + f*x) + a)
    F = -(c - d)*cos(e + f*x)/(f*(a*sin(e + f*x) + a)) + d*x/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_450():
    f = 1/(a*sin(e + f*x) + a)
    F = -cos(e + f*x)/(f*(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_451():
    f = 1/((c + d*sin(e + f*x))*(a*sin(e + f*x) + a))
    F = -cos(e + f*x)/(f*(c - d)*(a*sin(e + f*x) + a)) - 2*d*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(a*f*(c - d)*sqrt(c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_452():
    f = 1/((c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a))
    F = -cos(e + f*x)/(f*(c - d)*(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)) - 2*d*(2*c + d)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(a*f*(c - d)*(c**2 - d**2)**(sympy.S(3)/2)) - d*(c + 2*d)*cos(e + f*x)/(a*f*(c - d)**2*(c + d)*(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_453():
    f = 1/((c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a))
    F = -cos(e + f*x)/(f*(c - d)*(c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)) - 3*d*(2*c**2 + 2*c*d + d**2)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(a*f*(c - d)*(c**2 - d**2)**(sympy.S(5)/2)) - d*(2*c + 3*d)*cos(e + f*x)/(2*a*f*(c - d)**2*(c + d)*(c + d*sin(e + f*x))**2) - d*(c + 4*d)*(2*c + d)*cos(e + f*x)/(2*a*f*(c - d)**3*(c + d)**2*(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_454():
    f = (c + d*sin(e + f*x))**5/(a*sin(e + f*x) + a)**2
    F = -(c - d)*(c + d*sin(e + f*x))**4*cos(e + f*x)/(3*f*(a*sin(e + f*x) + a)**2) + d**2*x*(10*c - 5*d)*(2*c**2 - 3*c*d + 2*d**2)/(2*a**2) + d**2*(2*c**3 + 20*c**2*d - 57*c*d**2 + 30*d**3)*sin(e + f*x)*cos(e + f*x)/(6*a**2*f) + d*(c + d*sin(e + f*x))**2*(c**2 + 10*c*d - 12*d**2)*cos(e + f*x)/(3*a**2*f) + 2*d*(c**4 + 10*c**3*d - 44*c**2*d**2 + 40*c*d**3 - 12*d**4)*cos(e + f*x)/(3*a**2*f) - (c - d)*(c + 10*d)*(c + d*sin(e + f*x))**3*cos(e + f*x)/(3*a**2*f*(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_455():
    f = (c + d*sin(e + f*x))**4/(a*sin(e + f*x) + a)**2
    F = -(c - d)*(c + d*sin(e + f*x))**3*cos(e + f*x)/(3*f*(a*sin(e + f*x) + a)**2) + d**2*x*(12*c**2 - 16*c*d + 7*d**2)/(2*a**2) + d**2*(2*c**2 + 16*c*d - 21*d**2)*sin(e + f*x)*cos(e + f*x)/(6*a**2*f) + 2*d*(c**3 + 8*c**2*d - 20*c*d**2 + 8*d**3)*cos(e + f*x)/(3*a**2*f) - (c - d)*(c + 8*d)*(c + d*sin(e + f*x))**2*cos(e + f*x)/(3*a**2*f*(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_456():
    f = (c + d*sin(e + f*x))**3/(a*sin(e + f*x) + a)**2
    F = -(c - d)*(c + d*sin(e + f*x))**2*cos(e + f*x)/(3*f*(a*sin(e + f*x) + a)**2) + d**2*x*(3*c - 2*d)/a**2 + d**2*(c - 4*d)*cos(e + f*x)/(3*a**2*f) - (c - d)**2*(c + 6*d)*cos(e + f*x)/(3*a**2*f*(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_457():
    f = (c + d*sin(e + f*x))**2/(a*sin(e + f*x) + a)**2
    F = -(c - d)*(c + d*sin(e + f*x))*cos(e + f*x)/(3*f*(a*sin(e + f*x) + a)**2) + d**2*x/a**2 - (c - d)*(c + 4*d)*cos(e + f*x)/(3*a**2*f*(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_458():
    f = (c + d*sin(e + f*x))/(a*sin(e + f*x) + a)**2
    F = -(c - d)*cos(e + f*x)/(3*f*(a*sin(e + f*x) + a)**2) - (c + 2*d)*cos(e + f*x)/(3*f*(a**2*sin(e + f*x) + a**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_459():
    f = (a*sin(e + f*x) + a)**(-2)
    F = -cos(e + f*x)/(3*f*(a**2*sin(e + f*x) + a**2)) - cos(e + f*x)/(3*f*(a*sin(e + f*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_460():
    f = 1/((c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**2)
    F = -cos(e + f*x)/(f*(3*c - 3*d)*(a*sin(e + f*x) + a)**2) + 2*d**2*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(a**2*f*(c - d)**2*sqrt(c**2 - d**2)) - (c - 4*d)*cos(e + f*x)/(3*a**2*f*(c - d)**2*(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_461():
    f = 1/((c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**2)
    F = -cos(e + f*x)/(f*(c + d*sin(e + f*x))*(3*c - 3*d)*(a*sin(e + f*x) + a)**2) + 2*d**2*(3*c + 2*d)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(a**2*f*(c - d)**3*(c + d)*sqrt(c**2 - d**2)) - d*(c**2 - 6*c*d - 10*d**2)*cos(e + f*x)/(3*a**2*f*(c - d)**3*(c + d)*(c + d*sin(e + f*x))) - (c - 6*d)*cos(e + f*x)/(3*a**2*f*(c - d)**2*(c + d*sin(e + f*x))*(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_462():
    f = 1/((c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a)**2)
    F = -cos(e + f*x)/(f*(c + d*sin(e + f*x))**2*(3*c - 3*d)*(a*sin(e + f*x) + a)**2) + d**2*(12*c**2 + 16*c*d + 7*d**2)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(a**2*f*(c - d)**4*(c + d)**2*sqrt(c**2 - d**2)) - d*(2*c**2 - 16*c*d - 21*d**2)*cos(e + f*x)/(6*a**2*f*(c - d)**3*(c + d)*(c + d*sin(e + f*x))**2) - d*(2*c**3 - 16*c**2*d - 59*c*d**2 - 32*d**3)*cos(e + f*x)/(6*a**2*f*(c - d)**4*(c + d)**2*(c + d*sin(e + f*x))) - (c - 8*d)*cos(e + f*x)/(3*a**2*f*(c - d)**2*(c + d*sin(e + f*x))**2*(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_463():
    f = (c + d*sin(e + f*x))**6/(a*sin(e + f*x) + a)**3
    F = -(c - d)*(c + d*sin(e + f*x))**5*cos(e + f*x)/(5*f*(a*sin(e + f*x) + a)**3) - (c - d)*(c + d*sin(e + f*x))**3*(2*c**2 + 18*c*d + 115*d**2)*cos(e + f*x)/(15*f*(a**3*sin(e + f*x) + a**3)) - (c - d)*(c + d*sin(e + f*x))**4*(2*c + 13*d)*cos(e + f*x)/(15*a*f*(a*sin(e + f*x) + a)**2) + d**3*x*(40*c**3 - 90*c**2*d + 78*c*d**2 - 23*d**3)/(2*a**3) + d**2*(4*c**4 + 36*c**3*d + 216*c**2*d**2 - 626*c*d**3 + 345*d**4)*sin(e + f*x)*cos(e + f*x)/(30*a**3*f) + d*(c + d*sin(e + f*x))**2*(2*c**3 + 18*c**2*d + 111*c*d**2 - 136*d**3)*cos(e + f*x)/(15*a**3*f) + 2*d*(2*c**5 + 18*c**4*d + 107*c**3*d**2 - 472*c**2*d**3 + 456*c*d**4 - 136*d**5)*cos(e + f*x)/(15*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_464():
    f = (c + d*sin(e + f*x))**5/(a*sin(e + f*x) + a)**3
    F = -(c - d)*(c + d*sin(e + f*x))**4*cos(e + f*x)/(5*f*(a*sin(e + f*x) + a)**3) - (c - d)*(c + d*sin(e + f*x))**2*(2*c**2 + 15*c*d + 76*d**2)*cos(e + f*x)/(15*f*(a**3*sin(e + f*x) + a**3)) - (c - d)*(c + d*sin(e + f*x))**3*(2*c + 11*d)*cos(e + f*x)/(15*a*f*(a*sin(e + f*x) + a)**2) + d**3*x*(20*c**2 - 30*c*d + 13*d**2)/(2*a**3) + d**2*(4*c**3 + 30*c**2*d + 146*c*d**2 - 195*d**3)*sin(e + f*x)*cos(e + f*x)/(30*a**3*f) + 2*d*(2*c**4 + 15*c**3*d + 72*c**2*d**2 - 180*c*d**3 + 76*d**4)*cos(e + f*x)/(15*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_465():
    f = (c + d*sin(e + f*x))**4/(a*sin(e + f*x) + a)**3
    F = -(c - d)**2*(2*c**2 + 12*c*d + 45*d**2)*cos(e + f*x)/(15*f*(a**3*sin(e + f*x) + a**3)) - (c - d)*(c + d*sin(e + f*x))**3*cos(e + f*x)/(5*f*(a*sin(e + f*x) + a)**3) - (c - d)*(c + d*sin(e + f*x))**2*(2*c + 9*d)*cos(e + f*x)/(15*a*f*(a*sin(e + f*x) + a)**2) + d**3*x*(4*c - 3*d)/a**3 + d**2*(2*c**2 + 10*c*d - 27*d**2)*cos(e + f*x)/(15*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_466():
    f = (c + d*sin(e + f*x))**3/(a*sin(e + f*x) + a)**3
    F = -(c - d)*(c + d*sin(e + f*x))**2*cos(e + f*x)/(5*f*(a*sin(e + f*x) + a)**3) - (c - d)*(2*c**2 + 11*c*d + 29*d**2)*cos(e + f*x)/(15*f*(a**3*sin(e + f*x) + a**3)) - (c - d)**2*(2*c + 7*d)*cos(e + f*x)/(15*a*f*(a*sin(e + f*x) + a)**2) + d**3*x/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_467():
    f = (c + d*sin(e + f*x))**2/(a*sin(e + f*x) + a)**3
    F = -(c - d)*(c + d*sin(e + f*x))*cos(e + f*x)/(5*f*(a*sin(e + f*x) + a)**3) - (2*c**2 + 6*c*d + 7*d**2)*cos(e + f*x)/(15*f*(a**3*sin(e + f*x) + a**3)) - (c - d)*(2*c + 5*d)*cos(e + f*x)/(15*a*f*(a*sin(e + f*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_468():
    f = (c + d*sin(e + f*x))/(a*sin(e + f*x) + a)**3
    F = -(c - d)*cos(e + f*x)/(5*f*(a*sin(e + f*x) + a)**3) - (2*c + 3*d)*cos(e + f*x)/(15*f*(a**3*sin(e + f*x) + a**3)) - (2*c + 3*d)*cos(e + f*x)/(15*a*f*(a*sin(e + f*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_469():
    f = (a*sin(e + f*x) + a)**(-3)
    F = -2*cos(e + f*x)/(15*f*(a**3*sin(e + f*x) + a**3)) - cos(e + f*x)/(5*f*(a*sin(e + f*x) + a)**3) - 2*cos(e + f*x)/(15*a*f*(a*sin(e + f*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_470():
    f = 1/((c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**3)
    F = -cos(e + f*x)/(f*(5*c - 5*d)*(a*sin(e + f*x) + a)**3) - (2*c**2 - 9*c*d + 22*d**2)*cos(e + f*x)/(15*f*(c - d)**3*(a**3*sin(e + f*x) + a**3)) - (2*c - 7*d)*cos(e + f*x)/(15*a*f*(c - d)**2*(a*sin(e + f*x) + a)**2) - 2*d**3*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(a**3*f*(c - d)**3*sqrt(c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_471():
    f = 1/((c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**3)
    F = -cos(e + f*x)/(f*(c + d*sin(e + f*x))*(5*c - 5*d)*(a*sin(e + f*x) + a)**3) - (2*c**2 - 12*c*d + 45*d**2)*cos(e + f*x)/(15*f*(c - d)**3*(c + d*sin(e + f*x))*(a**3*sin(e + f*x) + a**3)) - (2*c - 9*d)*cos(e + f*x)/(15*a*f*(c - d)**2*(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**2) - 2*d**3*(4*c + 3*d)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(a**3*f*(c - d)**4*(c + d)*sqrt(c**2 - d**2)) - d*(2*c**3 - 12*c**2*d + 43*c*d**2 + 72*d**3)*cos(e + f*x)/(15*a**3*f*(c - d)**4*(c + d)*(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_472():
    f = 1/((c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a)**3)
    F = -cos(e + f*x)/(f*(c + d*sin(e + f*x))**2*(5*c - 5*d)*(a*sin(e + f*x) + a)**3) - (2*c**2 - 15*c*d + 76*d**2)*cos(e + f*x)/(15*f*(c - d)**3*(c + d*sin(e + f*x))**2*(a**3*sin(e + f*x) + a**3)) - (2*c - 11*d)*cos(e + f*x)/(15*a*f*(c - d)**2*(c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**2) - d**3*(20*c**2 + 30*c*d + 13*d**2)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(a**3*f*(c - d)**5*(c + d)**2*sqrt(c**2 - d**2)) - d*(4*c**3 - 30*c**2*d + 146*c*d**2 + 195*d**3)*cos(e + f*x)/(30*a**3*f*(c - d)**4*(c + d)*(c + d*sin(e + f*x))**2) - d*(4*c**4 - 30*c**3*d + 142*c**2*d**2 + 525*c*d**3 + 304*d**4)*cos(e + f*x)/(30*a**3*f*(c - d)**5*(c + d)**2*(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_473():
    f = (A + B*sin(x))/(sin(x) + 1)**4
    F = -(A - B)*cos(x)/(7*(sin(x) + 1)**4) - (3*A + 4*B)*cos(x)/(35*(sin(x) + 1)**3) - (6*A + 8*B)*cos(x)/(105*sin(x) + 105) - (6*A + 8*B)*cos(x)/(105*(sin(x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_474():
    f = (A + B*sin(x))/(1 - sin(x))**4
    F = (6*A - 8*B)*cos(x)/(105 - 105*sin(x)) + (6*A - 8*B)*cos(x)/(105*(1 - sin(x))**2) + (3*A - 4*B)*cos(x)/(35*(1 - sin(x))**3) + (A + B)*cos(x)/(7*(1 - sin(x))**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_475():
    f = (c + d*sin(e + f*x))**(sympy.S(5)/2)*(a*sin(e + f*x) + a)
    F = -2*a*(c + d*sin(e + f*x))**(sympy.S(5)/2)*cos(e + f*x)/(7*f) - 2*a*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(5*c + 7*d)*cos(e + f*x)/(35*f) - 2*a*sqrt(c + d*sin(e + f*x))*(15*c**2 + 56*c*d + 25*d**2)*cos(e + f*x)/(105*f) - 2*a*sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2)*(15*c**2 + 56*c*d + 25*d**2)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(105*d*f*sqrt(c + d*sin(e + f*x))) + 2*a*sqrt(c + d*sin(e + f*x))*(15*c**3 + 161*c**2*d + 145*c*d**2 + 63*d**3)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(105*d*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_476():
    f = (c + d*sin(e + f*x))**(sympy.S(3)/2)*(a*sin(e + f*x) + a)
    F = -2*a*(c + d*sin(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(5*f) - 2*a*sqrt(c + d*sin(e + f*x))*(3*c + 5*d)*cos(e + f*x)/(15*f) - 2*a*sqrt((c + d*sin(e + f*x))/(c + d))*(3*c + 5*d)*(c**2 - d**2)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d*f*sqrt(c + d*sin(e + f*x))) + 2*a*sqrt(c + d*sin(e + f*x))*(3*c**2 + 20*c*d + 9*d**2)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_477():
    f = sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)
    F = -2*a*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(3*f) - 2*a*sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d*f*sqrt(c + d*sin(e + f*x))) + 2*a*(c + 3*d)*sqrt(c + d*sin(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_478():
    f = (a*sin(e + f*x) + a)/sqrt(c + d*sin(e + f*x))
    F = -2*a*sqrt((c + d*sin(e + f*x))/(c + d))*(c - d)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(d*f*sqrt(c + d*sin(e + f*x))) + 2*a*sqrt(c + d*sin(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(d*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_479():
    f = (a*sin(e + f*x) + a)/(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = -2*a*cos(e + f*x)/(f*(c + d)*sqrt(c + d*sin(e + f*x))) + 2*a*sqrt((c + d*sin(e + f*x))/(c + d))*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(d*f*sqrt(c + d*sin(e + f*x))) - 2*a*sqrt(c + d*sin(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(d*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_480():
    f = (a*sin(e + f*x) + a)/(c + d*sin(e + f*x))**(sympy.S(5)/2)
    F = -2*a*(c - 3*d)*cos(e + f*x)/(f*(c + d)**2*sqrt(c + d*sin(e + f*x))*(3*c - 3*d)) - 2*a*cos(e + f*x)/(f*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(3*c + 3*d)) + 2*a*sqrt((c + d*sin(e + f*x))/(c + d))*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d*f*(c + d)*sqrt(c + d*sin(e + f*x))) - 2*a*(c - 3*d)*sqrt(c + d*sin(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(d*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c + d)**2*(3*c - 3*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_481():
    f = (a*sin(e + f*x) + a)/(c + d*sin(e + f*x))**(sympy.S(7)/2)
    F = -2*a*cos(e + f*x)/(f*(c + d*sin(e + f*x))**(sympy.S(5)/2)*(5*c + 5*d)) - 2*a*(3*c - 5*d)*cos(e + f*x)/(f*(c + d)**2*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(15*c - 15*d)) - 2*a*(3*c**2 - 20*c*d + 9*d**2)*cos(e + f*x)/(15*f*(c - d)**2*(c + d)**3*sqrt(c + d*sin(e + f*x))) + 2*a*sqrt((c + d*sin(e + f*x))/(c + d))*(3*c - 5*d)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(d*f*(c + d)**2*sqrt(c + d*sin(e + f*x))*(15*c - 15*d)) - 2*a*sqrt(c + d*sin(e + f*x))*(3*c**2 - 20*c*d + 9*d**2)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c - d)**2*(c + d)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_482():
    f = (c + d*sin(e + f*x))**(sympy.S(5)/2)*(a*sin(e + f*x) + a)**2
    F = 4*a**2*(c - 9*d)*(c + d*sin(e + f*x))**(sympy.S(5)/2)*cos(e + f*x)/(63*d*f) - 2*a**2*(c + d*sin(e + f*x))**(sympy.S(7)/2)*cos(e + f*x)/(9*d*f) + 4*a**2*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(5*c*(c - 9*d) - 56*d**2)*cos(e + f*x)/(315*d*f) + 4*a**2*sqrt(c + d*sin(e + f*x))*(5*c**3 - 45*c**2*d - 141*c*d**2 - 75*d**3)*cos(e + f*x)/(315*d*f) + 4*a**2*sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2)*(5*c**3 - 45*c**2*d - 141*c*d**2 - 75*d**3)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(315*d**2*f*sqrt(c + d*sin(e + f*x))) - 4*a**2*sqrt(c + d*sin(e + f*x))*(5*c**4 - 45*c**3*d - 381*c**2*d**2 - 435*c*d**3 - 168*d**4)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(315*d**2*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_483():
    f = (c + d*sin(e + f*x))**(sympy.S(3)/2)*(a*sin(e + f*x) + a)**2
    F = 4*a**2*(c - 7*d)*(c + d*sin(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(35*d*f) - 2*a**2*(c + d*sin(e + f*x))**(sympy.S(5)/2)*cos(e + f*x)/(7*d*f) + 4*a**2*sqrt(c + d*sin(e + f*x))*(c**2 - 7*c*d - 10*d**2)*cos(e + f*x)/(35*d*f) + 4*a**2*sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2)*(c**2 - 7*c*d - 10*d**2)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(35*d**2*f*sqrt(c + d*sin(e + f*x))) - 4*a**2*(c + 3*d)*sqrt(c + d*sin(e + f*x))*(c**2 - 10*c*d - 7*d**2)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(35*d**2*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_484():
    f = sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**2
    F = 4*a**2*(c - 5*d)*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(15*d*f) - 2*a**2*(c + d*sin(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(5*d*f) + 4*a**2*sqrt((c + d*sin(e + f*x))/(c + d))*(c - 5*d)*(c**2 - d**2)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d**2*f*sqrt(c + d*sin(e + f*x))) - 4*a**2*sqrt(c + d*sin(e + f*x))*(c**2 - 5*c*d - 12*d**2)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d**2*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_485():
    f = (a*sin(e + f*x) + a)**2/sqrt(c + d*sin(e + f*x))
    F = -2*a**2*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(3*d*f) + 4*a**2*sqrt((c + d*sin(e + f*x))/(c + d))*(c - 2*d)*(c - d)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d**2*f*sqrt(c + d*sin(e + f*x))) - 4*a**2*(c - 3*d)*sqrt(c + d*sin(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d**2*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_486():
    f = (a*sin(e + f*x) + a)**2/(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = 4*a**2*c*sqrt(c + d*sin(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(d**2*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c + d)) + 2*a**2*(c - d)*cos(e + f*x)/(d*f*(c + d)*sqrt(c + d*sin(e + f*x))) - 4*a**2*sqrt((c + d*sin(e + f*x))/(c + d))*(c - d)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(d**2*f*sqrt(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_487():
    f = (a*sin(e + f*x) + a)**2/(c + d*sin(e + f*x))**(sympy.S(5)/2)
    F = 2*a**2*(c - d)*cos(e + f*x)/(3*d*f*(c + d)*(c + d*sin(e + f*x))**(sympy.S(3)/2)) - 4*a**2*(c + 3*d)*cos(e + f*x)/(3*d*f*(c + d)**2*sqrt(c + d*sin(e + f*x))) + 4*a**2*sqrt((c + d*sin(e + f*x))/(c + d))*(c + 2*d)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d**2*f*(c + d)*sqrt(c + d*sin(e + f*x))) - 4*a**2*(c + 3*d)*sqrt(c + d*sin(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d**2*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c + d)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_488():
    f = (a*sin(e + f*x) + a)**2/(c + d*sin(e + f*x))**(sympy.S(7)/2)
    F = 2*a**2*(c - d)*cos(e + f*x)/(5*d*f*(c + d)*(c + d*sin(e + f*x))**(sympy.S(5)/2)) - 4*a**2*(c + 5*d)*cos(e + f*x)/(15*d*f*(c + d)**2*(c + d*sin(e + f*x))**(sympy.S(3)/2)) - 4*a**2*(c**2 + 5*c*d - 12*d**2)*cos(e + f*x)/(d*f*(c + d)**3*sqrt(c + d*sin(e + f*x))*(15*c - 15*d)) + 4*a**2*sqrt((c + d*sin(e + f*x))/(c + d))*(c + 5*d)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d**2*f*(c + d)**2*sqrt(c + d*sin(e + f*x))) - 4*a**2*sqrt(c + d*sin(e + f*x))*(c**2 + 5*c*d - 12*d**2)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(d**2*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c + d)**3*(15*c - 15*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_489():
    f = (c + d*sin(e + f*x))**(sympy.S(5)/2)*(a*sin(e + f*x) + a)**3
    F = 8*a**3*(c - 6*d)*(c + d*sin(e + f*x))**(sympy.S(7)/2)*cos(e + f*x)/(99*d**2*f) - 4*a**3*(c + d*sin(e + f*x))**(sympy.S(5)/2)*(4*c**2 - 33*c*d + 189*d**2)*cos(e + f*x)/(693*d**2*f) - 4*a**3*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(4*c**3 - 33*c**2*d + 182*c*d**2 + 231*d**3)*cos(e + f*x)/(693*d**2*f) - 4*a**3*sqrt(c + d*sin(e + f*x))*(4*c**4 - 33*c**3*d + 177*c**2*d**2 + 561*c*d**3 + 315*d**4)*cos(e + f*x)/(693*d**2*f) - 4*a**3*sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2)*(4*c**4 - 33*c**3*d + 177*c**2*d**2 + 561*c*d**3 + 315*d**4)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(693*d**3*f*sqrt(c + d*sin(e + f*x))) + 4*a**3*(c + 3*d)*sqrt(c + d*sin(e + f*x))*(4*c**4 - 45*c**3*d + 309*c**2*d**2 + 525*c*d**3 + 231*d**4)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(693*d**3*f*sqrt((c + d*sin(e + f*x))/(c + d))) - 2*(c + d*sin(e + f*x))**(sympy.S(7)/2)*(a**3*sin(e + f*x) + a**3)*cos(e + f*x)/(11*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_490():
    f = (c + d*sin(e + f*x))**(sympy.S(3)/2)*(a*sin(e + f*x) + a)**3
    F = 8*a**3*(c - 5*d)*(c + d*sin(e + f*x))**(sympy.S(5)/2)*cos(e + f*x)/(63*d**2*f) - 4*a**3*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(4*c**2 - 27*c*d + 119*d**2)*cos(e + f*x)/(315*d**2*f) - 4*a**3*sqrt(c + d*sin(e + f*x))*(4*c**3 - 27*c**2*d + 114*c*d**2 + 165*d**3)*cos(e + f*x)/(315*d**2*f) - 4*a**3*sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2)*(4*c**3 - 27*c**2*d + 114*c*d**2 + 165*d**3)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(315*d**3*f*sqrt(c + d*sin(e + f*x))) + 4*a**3*sqrt(c + d*sin(e + f*x))*(4*c**4 - 27*c**3*d + 111*c**2*d**2 + 579*c*d**3 + 357*d**4)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(315*d**3*f*sqrt((c + d*sin(e + f*x))/(c + d))) - 2*(c + d*sin(e + f*x))**(sympy.S(5)/2)*(a**3*sin(e + f*x) + a**3)*cos(e + f*x)/(9*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_491():
    f = sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**3
    F = 8*a**3*(c - 4*d)*(c + d*sin(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(35*d**2*f) - 4*a**3*sqrt(c + d*sin(e + f*x))*(4*c**2 - 21*c*d + 65*d**2)*cos(e + f*x)/(105*d**2*f) - 4*a**3*sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2)*(4*c**2 - 21*c*d + 65*d**2)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(105*d**3*f*sqrt(c + d*sin(e + f*x))) + 4*a**3*sqrt(c + d*sin(e + f*x))*(4*c**3 - 21*c**2*d + 62*c*d**2 + 147*d**3)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(105*d**3*f*sqrt((c + d*sin(e + f*x))/(c + d))) - 2*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(a**3*sin(e + f*x) + a**3)*cos(e + f*x)/(7*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_492():
    f = (a*sin(e + f*x) + a)**3/sqrt(c + d*sin(e + f*x))
    F = 8*a**3*(c - 3*d)*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(15*d**2*f) - 4*a**3*sqrt((c + d*sin(e + f*x))/(c + d))*(c - d)*(4*c**2 - 11*c*d + 15*d**2)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d**3*f*sqrt(c + d*sin(e + f*x))) + 4*a**3*sqrt(c + d*sin(e + f*x))*(4*c**2 - 15*c*d + 27*d**2)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d**3*f*sqrt((c + d*sin(e + f*x))/(c + d))) - 2*sqrt(c + d*sin(e + f*x))*(a**3*sin(e + f*x) + a**3)*cos(e + f*x)/(5*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_493():
    f = (a*sin(e + f*x) + a)**3/(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = -4*a**3*sqrt(c + d*sin(e + f*x))*(2*c - d)*cos(e + f*x)/(3*d**2*f*(c + d)) + 4*a**3*sqrt((c + d*sin(e + f*x))/(c + d))*(c - d)*(4*c - 5*d)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d**3*f*sqrt(c + d*sin(e + f*x))) - 4*a**3*sqrt(c + d*sin(e + f*x))*(4*c**2 - 5*c*d - 3*d**2)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d**3*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c + d)) + (2*c - 2*d)*(a**3*sin(e + f*x) + a**3)*cos(e + f*x)/(d*f*(c + d)*sqrt(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_494():
    f = (a*sin(e + f*x) + a)**3/(c + d*sin(e + f*x))**(sympy.S(5)/2)
    F = 8*a**3*(c - d)*(c + 2*d)*cos(e + f*x)/(3*d**2*f*(c + d)**2*sqrt(c + d*sin(e + f*x))) - 4*a**3*sqrt((c + d*sin(e + f*x))/(c + d))*(c - d)*(4*c + 5*d)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d**3*f*(c + d)*sqrt(c + d*sin(e + f*x))) + 4*a**3*sqrt(c + d*sin(e + f*x))*(4*c**2 + 5*c*d - 3*d**2)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d**3*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c + d)**2) + (2*c - 2*d)*(a**3*sin(e + f*x) + a**3)*cos(e + f*x)/(3*d*f*(c + d)*(c + d*sin(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_495():
    f = (a*sin(e + f*x) + a)**3/(c + d*sin(e + f*x))**(sympy.S(7)/2)
    F = 8*a**3*(c - d)*(c + 3*d)*cos(e + f*x)/(15*d**2*f*(c + d)**2*(c + d*sin(e + f*x))**(sympy.S(3)/2)) - 4*a**3*(4*c**2 + 15*c*d + 27*d**2)*cos(e + f*x)/(15*d**2*f*(c + d)**3*sqrt(c + d*sin(e + f*x))) + 4*a**3*sqrt((c + d*sin(e + f*x))/(c + d))*(4*c**2 + 11*c*d + 15*d**2)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d**3*f*(c + d)**2*sqrt(c + d*sin(e + f*x))) - 4*a**3*sqrt(c + d*sin(e + f*x))*(4*c**2 + 15*c*d + 27*d**2)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d**3*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c + d)**3) + (2*c - 2*d)*(a**3*sin(e + f*x) + a**3)*cos(e + f*x)/(5*d*f*(c + d)*(c + d*sin(e + f*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_496():
    f = (a*sin(e + f*x) + a)**3/(c + d*sin(e + f*x))**(sympy.S(9)/2)
    F = 8*a**3*(c - d)*(c + 4*d)*cos(e + f*x)/(35*d**2*f*(c + d)**2*(c + d*sin(e + f*x))**(sympy.S(5)/2)) - 4*a**3*(4*c**2 + 21*c*d + 65*d**2)*cos(e + f*x)/(105*d**2*f*(c + d)**3*(c + d*sin(e + f*x))**(sympy.S(3)/2)) - 4*a**3*(4*c**3 + 21*c**2*d + 62*c*d**2 - 147*d**3)*cos(e + f*x)/(d**2*f*(c + d)**4*sqrt(c + d*sin(e + f*x))*(105*c - 105*d)) + 4*a**3*sqrt((c + d*sin(e + f*x))/(c + d))*(4*c**2 + 21*c*d + 65*d**2)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(105*d**3*f*(c + d)**3*sqrt(c + d*sin(e + f*x))) - 4*a**3*sqrt(c + d*sin(e + f*x))*(4*c**3 + 21*c**2*d + 62*c*d**2 - 147*d**3)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(d**3*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c + d)**4*(105*c - 105*d)) + (2*c - 2*d)*(a**3*sin(e + f*x) + a**3)*cos(e + f*x)/(7*d*f*(c + d)*(c + d*sin(e + f*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_497():
    f = (c + d*sin(e + f*x))**(sympy.S(5)/2)/(a*sin(e + f*x) + a)
    F = -(c - d)*(c + d*sin(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(f*(a*sin(e + f*x) + a)) + d*sqrt(c + d*sin(e + f*x))*(3*c - 5*d)*cos(e + f*x)/(3*a*f) + sqrt((c + d*sin(e + f*x))/(c + d))*(3*c - 5*d)*(c**2 - d**2)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*a*f*sqrt(c + d*sin(e + f*x))) - sqrt(c + d*sin(e + f*x))*(3*c**2 - 20*c*d + 9*d**2)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*a*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_498():
    f = (c + d*sin(e + f*x))**(sympy.S(3)/2)/(a*sin(e + f*x) + a)
    F = -(c - d)*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(f*(a*sin(e + f*x) + a)) + sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(a*f*sqrt(c + d*sin(e + f*x))) - (c - 3*d)*sqrt(c + d*sin(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(a*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_499():
    f = sqrt(c + d*sin(e + f*x))/(a*sin(e + f*x) + a)
    F = -sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(f*(a*sin(e + f*x) + a)) + sqrt((c + d*sin(e + f*x))/(c + d))*(c + d)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(a*f*sqrt(c + d*sin(e + f*x))) - sqrt(c + d*sin(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(a*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_500():
    f = 1/(sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a))
    F = -sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(f*(c - d)*(a*sin(e + f*x) + a)) + sqrt((c + d*sin(e + f*x))/(c + d))*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(a*f*sqrt(c + d*sin(e + f*x))) - sqrt(c + d*sin(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(a*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c - d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_501():
    f = 1/((c + d*sin(e + f*x))**(sympy.S(3)/2)*(a*sin(e + f*x) + a))
    F = -cos(e + f*x)/(f*(c - d)*sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)) - d*(c + 3*d)*cos(e + f*x)/(a*f*(c - d)**2*(c + d)*sqrt(c + d*sin(e + f*x))) + sqrt((c + d*sin(e + f*x))/(c + d))*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(a*f*(c - d)*sqrt(c + d*sin(e + f*x))) - (c + 3*d)*sqrt(c + d*sin(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(a*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c - d)**2*(c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_502():
    f = 1/((c + d*sin(e + f*x))**(sympy.S(5)/2)*(a*sin(e + f*x) + a))
    F = -cos(e + f*x)/(f*(c - d)*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(a*sin(e + f*x) + a)) - d*(3*c + 5*d)*cos(e + f*x)/(3*a*f*(c - d)**2*(c + d)*(c + d*sin(e + f*x))**(sympy.S(3)/2)) - d*(3*c**2 + 20*c*d + 9*d**2)*cos(e + f*x)/(3*a*f*(c - d)**3*(c + d)**2*sqrt(c + d*sin(e + f*x))) + sqrt((c + d*sin(e + f*x))/(c + d))*(3*c + 5*d)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*a*f*(c - d)**2*(c + d)*sqrt(c + d*sin(e + f*x))) - sqrt(c + d*sin(e + f*x))*(3*c**2 + 20*c*d + 9*d**2)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*a*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c - d)**3*(c + d)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_503():
    f = (c + d*sin(e + f*x))**(sympy.S(5)/2)/(a*sin(e + f*x) + a)**2
    F = -(c - d)*(c + d*sin(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(3*f*(a*sin(e + f*x) + a)**2) + sqrt((c + d*sin(e + f*x))/(c + d))*(c + 5*d)*(c**2 - d**2)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*a**2*f*sqrt(c + d*sin(e + f*x))) - (c - d)*(c + 5*d)*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(3*a**2*f*(sin(e + f*x) + 1)) - sqrt(c + d*sin(e + f*x))*(c**2 + 5*c*d - 12*d**2)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*a**2*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_504():
    f = (c + d*sin(e + f*x))**(sympy.S(3)/2)/(a*sin(e + f*x) + a)**2
    F = -(c - d)*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(3*f*(a*sin(e + f*x) + a)**2) + sqrt((c + d*sin(e + f*x))/(c + d))*(c + d)*(c + 2*d)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*a**2*f*sqrt(c + d*sin(e + f*x))) - (c + 3*d)*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(3*a**2*f*(sin(e + f*x) + 1)) - (c + 3*d)*sqrt(c + d*sin(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*a**2*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_505():
    f = sqrt(c + d*sin(e + f*x))/(a*sin(e + f*x) + a)**2
    F = -sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(3*f*(a*sin(e + f*x) + a)**2) - c*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(3*a**2*f*(c - d)*(sin(e + f*x) + 1)) - c*sqrt(c + d*sin(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*a**2*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c - d)) + sqrt((c + d*sin(e + f*x))/(c + d))*(c + d)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*a**2*f*sqrt(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_506():
    f = 1/(sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**2)
    F = -sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(f*(3*c - 3*d)*(a*sin(e + f*x) + a)**2) + sqrt((c + d*sin(e + f*x))/(c + d))*(c - 2*d)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*a**2*f*(c - d)*sqrt(c + d*sin(e + f*x))) - (c - 3*d)*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(3*a**2*f*(c - d)**2*(sin(e + f*x) + 1)) - (c - 3*d)*sqrt(c + d*sin(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*a**2*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c - d)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_507():
    f = 1/((c + d*sin(e + f*x))**(sympy.S(3)/2)*(a*sin(e + f*x) + a)**2)
    F = -cos(e + f*x)/(f*sqrt(c + d*sin(e + f*x))*(3*c - 3*d)*(a*sin(e + f*x) + a)**2) - d*(c**2 - 5*c*d - 12*d**2)*cos(e + f*x)/(3*a**2*f*(c - d)**3*(c + d)*sqrt(c + d*sin(e + f*x))) + sqrt((c + d*sin(e + f*x))/(c + d))*(c - 5*d)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*a**2*f*(c - d)**2*sqrt(c + d*sin(e + f*x))) - (c - 5*d)*cos(e + f*x)/(3*a**2*f*(c - d)**2*sqrt(c + d*sin(e + f*x))*(sin(e + f*x) + 1)) - sqrt(c + d*sin(e + f*x))*(c**2 - 5*c*d - 12*d**2)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*a**2*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c - d)**3*(c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_508():
    f = 1/((c + d*sin(e + f*x))**(sympy.S(5)/2)*(a*sin(e + f*x) + a)**2)
    F = -cos(e + f*x)/(f*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(3*c - 3*d)*(a*sin(e + f*x) + a)**2) - d*(c**2 - 7*c*d - 10*d**2)*cos(e + f*x)/(3*a**2*f*(c - d)**3*(c + d)*(c + d*sin(e + f*x))**(sympy.S(3)/2)) - d*(c + 3*d)*(c**2 - 10*c*d - 7*d**2)*cos(e + f*x)/(3*a**2*f*(c - d)**4*(c + d)**2*sqrt(c + d*sin(e + f*x))) + sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - 7*c*d - 10*d**2)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*a**2*f*(c - d)**3*(c + d)*sqrt(c + d*sin(e + f*x))) - (c - 7*d)*cos(e + f*x)/(3*a**2*f*(c - d)**2*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(sin(e + f*x) + 1)) - (c + 3*d)*sqrt(c + d*sin(e + f*x))*(c**2 - 10*c*d - 7*d**2)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*a**2*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c - d)**4*(c + d)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_509():
    f = (c + d*sin(e + f*x))**(sympy.S(5)/2)/(a*sin(e + f*x) + a)**3
    F = -(c - d)*(c + d*sin(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(5*f*(a*sin(e + f*x) + a)**3) - sqrt(c + d*sin(e + f*x))*(4*c**2 + 15*c*d + 27*d**2)*cos(e + f*x)/(30*f*(a**3*sin(e + f*x) + a**3)) + (-2*c + 2*d)*(c + 3*d)*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(15*a*f*(a*sin(e + f*x) + a)**2) + sqrt((c + d*sin(e + f*x))/(c + d))*(c + d)*(4*c**2 + 11*c*d + 15*d**2)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(30*a**3*f*sqrt(c + d*sin(e + f*x))) - sqrt(c + d*sin(e + f*x))*(4*c**2 + 15*c*d + 27*d**2)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(30*a**3*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_510():
    f = (c + d*sin(e + f*x))**(sympy.S(3)/2)/(a*sin(e + f*x) + a)**3
    F = -(c - d)*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(5*f*(a*sin(e + f*x) + a)**3) - sqrt(c + d*sin(e + f*x))*(4*c**2 + 5*c*d - 3*d**2)*cos(e + f*x)/(f*(30*c - 30*d)*(a**3*sin(e + f*x) + a**3)) - sqrt(c + d*sin(e + f*x))*(2*c + 4*d)*cos(e + f*x)/(15*a*f*(a*sin(e + f*x) + a)**2) + sqrt((c + d*sin(e + f*x))/(c + d))*(c + d)*(4*c + 5*d)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(30*a**3*f*sqrt(c + d*sin(e + f*x))) - sqrt(c + d*sin(e + f*x))*(4*c**2 + 5*c*d - 3*d**2)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(30*a**3*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c - d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_511():
    f = sqrt(c + d*sin(e + f*x))/(a*sin(e + f*x) + a)**3
    F = -sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(5*f*(a*sin(e + f*x) + a)**3) - sqrt(c + d*sin(e + f*x))*(4*c**2 - 5*c*d - 3*d**2)*cos(e + f*x)/(30*f*(c - d)**2*(a**3*sin(e + f*x) + a**3)) - sqrt(c + d*sin(e + f*x))*(2*c - d)*cos(e + f*x)/(15*a*f*(c - d)*(a*sin(e + f*x) + a)**2) + sqrt((c + d*sin(e + f*x))/(c + d))*(c + d)*(4*c - 5*d)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(30*a**3*f*(c - d)*sqrt(c + d*sin(e + f*x))) - sqrt(c + d*sin(e + f*x))*(4*c**2 - 5*c*d - 3*d**2)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(30*a**3*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c - d)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_512():
    f = 1/(sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**3)
    F = -sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(f*(5*c - 5*d)*(a*sin(e + f*x) + a)**3) - sqrt(c + d*sin(e + f*x))*(4*c**2 - 15*c*d + 27*d**2)*cos(e + f*x)/(30*f*(c - d)**3*(a**3*sin(e + f*x) + a**3)) - sqrt(c + d*sin(e + f*x))*(2*c - 6*d)*cos(e + f*x)/(15*a*f*(c - d)**2*(a*sin(e + f*x) + a)**2) + sqrt((c + d*sin(e + f*x))/(c + d))*(4*c**2 - 11*c*d + 15*d**2)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(30*a**3*f*(c - d)**2*sqrt(c + d*sin(e + f*x))) - sqrt(c + d*sin(e + f*x))*(4*c**2 - 15*c*d + 27*d**2)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(30*a**3*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c - d)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_513():
    f = 1/((c + d*sin(e + f*x))**(sympy.S(3)/2)*(a*sin(e + f*x) + a)**3)
    F = -cos(e + f*x)/(f*sqrt(c + d*sin(e + f*x))*(5*c - 5*d)*(a*sin(e + f*x) + a)**3) - (4*c**2 - 21*c*d + 65*d**2)*cos(e + f*x)/(30*f*(c - d)**3*sqrt(c + d*sin(e + f*x))*(a**3*sin(e + f*x) + a**3)) - (2*c - 8*d)*cos(e + f*x)/(15*a*f*(c - d)**2*sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**2) - d*(4*c**3 - 21*c**2*d + 62*c*d**2 + 147*d**3)*cos(e + f*x)/(30*a**3*f*(c - d)**4*(c + d)*sqrt(c + d*sin(e + f*x))) + sqrt((c + d*sin(e + f*x))/(c + d))*(4*c**2 - 21*c*d + 65*d**2)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(30*a**3*f*(c - d)**3*sqrt(c + d*sin(e + f*x))) - sqrt(c + d*sin(e + f*x))*(4*c**3 - 21*c**2*d + 62*c*d**2 + 147*d**3)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(30*a**3*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c - d)**4*(c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_514():
    f = 1/((c + d*sin(e + f*x))**(sympy.S(5)/2)*(a*sin(e + f*x) + a)**3)
    F = -cos(e + f*x)/(f*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(5*c - 5*d)*(a*sin(e + f*x) + a)**3) - (4*c**2 - 27*c*d + 119*d**2)*cos(e + f*x)/(30*f*(c - d)**3*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(a**3*sin(e + f*x) + a**3)) - (2*c - 10*d)*cos(e + f*x)/(15*a*f*(c - d)**2*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(a*sin(e + f*x) + a)**2) - d*(4*c**3 - 27*c**2*d + 114*c*d**2 + 165*d**3)*cos(e + f*x)/(30*a**3*f*(c - d)**4*(c + d)*(c + d*sin(e + f*x))**(sympy.S(3)/2)) - d*(4*c**4 - 27*c**3*d + 111*c**2*d**2 + 579*c*d**3 + 357*d**4)*cos(e + f*x)/(30*a**3*f*(c - d)**5*(c + d)**2*sqrt(c + d*sin(e + f*x))) + sqrt((c + d*sin(e + f*x))/(c + d))*(4*c**3 - 27*c**2*d + 114*c*d**2 + 165*d**3)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(30*a**3*f*(c - d)**4*(c + d)*sqrt(c + d*sin(e + f*x))) - sqrt(c + d*sin(e + f*x))*(4*c**4 - 27*c**3*d + 111*c**2*d**2 + 579*c*d**3 + 357*d**4)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(30*a**3*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c - d)**5*(c + d)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_515():
    f = (c + d*sin(e + f*x))**3*sqrt(a*sin(e + f*x) + a)
    F = -4*a*(c + d)*(15*c**2 + 10*c*d + 7*d**2)*cos(e + f*x)/(35*f*sqrt(a*sin(e + f*x) + a)) - 2*a*(c + d*sin(e + f*x))**3*cos(e + f*x)/(7*f*sqrt(a*sin(e + f*x) + a)) - d*(c + d)*(40*c - 8*d)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(35*f) - 12*d**2*(c + d)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(35*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_516():
    f = (c + d*sin(e + f*x))**2*sqrt(a*sin(e + f*x) + a)
    F = -2*a*(15*c**2 + 10*c*d + 7*d**2)*cos(e + f*x)/(15*f*sqrt(a*sin(e + f*x) + a)) - d*(20*c - 4*d)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(15*f) - 2*d**2*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(5*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_517():
    f = (c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)
    F = -2*a*(3*c + d)*cos(e + f*x)/(3*f*sqrt(a*sin(e + f*x) + a)) - 2*d*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_518():
    f = sqrt(a*sin(e + f*x) + a)
    F = -2*a*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_519():
    f = sqrt(a*sin(e + f*x) + a)/(c + d*sin(e + f*x))
    F = -2*sqrt(a)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(sqrt(d)*f*sqrt(c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_520():
    f = sqrt(a*sin(e + f*x) + a)/(c + d*sin(e + f*x))**2
    F = -sqrt(a)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(sqrt(d)*f*(c + d)**(sympy.S(3)/2)) - a*cos(e + f*x)/(f*(c + d)*(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_521():
    f = sqrt(a*sin(e + f*x) + a)/(c + d*sin(e + f*x))**3
    F = -3*sqrt(a)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(4*sqrt(d)*f*(c + d)**(sympy.S(5)/2)) - a*cos(e + f*x)/(f*(c + d*sin(e + f*x))**2*(2*c + 2*d)*sqrt(a*sin(e + f*x) + a)) - 3*a*cos(e + f*x)/(4*f*(c + d)**2*(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_522():
    f = (c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = 4*a**2*(c - 17*d)*(c + d)*(15*c**2 + 10*c*d + 7*d**2)*cos(e + f*x)/(315*d*f*sqrt(a*sin(e + f*x) + a)) + 2*a**2*(c - 17*d)*(c + d*sin(e + f*x))**3*cos(e + f*x)/(63*d*f*sqrt(a*sin(e + f*x) + a)) - 2*a**2*(c + d*sin(e + f*x))**4*cos(e + f*x)/(9*d*f*sqrt(a*sin(e + f*x) + a)) + 8*a*(c - 17*d)*(c + d)*(5*c - d)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(315*f) + d*(c + d)*(4*c - 68*d)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(105*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_523():
    f = (c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -8*a**2*(35*c**2 + 42*c*d + 19*d**2)*cos(e + f*x)/(105*f*sqrt(a*sin(e + f*x) + a)) - 2*a*sqrt(a*sin(e + f*x) + a)*(35*c**2 + 42*c*d + 19*d**2)*cos(e + f*x)/(105*f) - d*(28*c - 4*d)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(35*f) - 2*d**2*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(7*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_524():
    f = (c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -8*a**2*(5*c + 3*d)*cos(e + f*x)/(15*f*sqrt(a*sin(e + f*x) + a)) - 2*a*(5*c + 3*d)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(15*f) - 2*d*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_525():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -8*a**2*cos(e + f*x)/(3*f*sqrt(a*sin(e + f*x) + a)) - 2*a*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_526():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)/(c + d*sin(e + f*x))
    F = 2*a**(sympy.S(3)/2)*(c - d)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(d**(sympy.S(3)/2)*f*sqrt(c + d)) - 2*a**2*cos(e + f*x)/(d*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_527():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)/(c + d*sin(e + f*x))**2
    F = -a**(sympy.S(3)/2)*(c + 3*d)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(d**(sympy.S(3)/2)*f*(c + d)**(sympy.S(3)/2)) + a**2*(c - d)*cos(e + f*x)/(d*f*(c + d)*(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_528():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)/(c + d*sin(e + f*x))**3
    F = -a**(sympy.S(3)/2)*(c + 7*d)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(4*d**(sympy.S(3)/2)*f*(c + d)**(sympy.S(5)/2)) + a**2*(c - d)*cos(e + f*x)/(2*d*f*(c + d)*(c + d*sin(e + f*x))**2*sqrt(a*sin(e + f*x) + a)) - a**2*(c + 7*d)*cos(e + f*x)/(4*d*f*(c + d)**2*(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_529():
    f = (c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -4*a**3*(c + d)*(3*c**2 - 38*c*d + 355*d**2)*(15*c**2 + 10*c*d + 7*d**2)*cos(e + f*x)/(3465*d**2*f*sqrt(a*sin(e + f*x) + a)) + 2*a**3*(c + d*sin(e + f*x))**4*(3*c - 23*d)*cos(e + f*x)/(99*d**2*f*sqrt(a*sin(e + f*x) + a)) - 2*a**3*(c + d*sin(e + f*x))**3*(3*c**2 - 38*c*d + 355*d**2)*cos(e + f*x)/(693*d**2*f*sqrt(a*sin(e + f*x) + a)) - 8*a**2*(c + d)*(5*c - d)*sqrt(a*sin(e + f*x) + a)*(3*c**2 - 38*c*d + 355*d**2)*cos(e + f*x)/(3465*d*f) - 2*a**2*(c + d*sin(e + f*x))**4*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(11*d*f) - 4*a*(c + d)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(3*c**2 - 38*c*d + 355*d**2)*cos(e + f*x)/(1155*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_530():
    f = (c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -64*a**3*(21*c**2 + 30*c*d + 13*d**2)*cos(e + f*x)/(315*f*sqrt(a*sin(e + f*x) + a)) - 16*a**2*sqrt(a*sin(e + f*x) + a)*(21*c**2 + 30*c*d + 13*d**2)*cos(e + f*x)/(315*f) - 2*a*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(21*c**2 + 30*c*d + 13*d**2)*cos(e + f*x)/(105*f) - d*(36*c - 4*d)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(63*f) - 2*d**2*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(9*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_531():
    f = (c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -64*a**3*(7*c + 5*d)*cos(e + f*x)/(105*f*sqrt(a*sin(e + f*x) + a)) - 16*a**2*(7*c + 5*d)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(105*f) - 2*a*(7*c + 5*d)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(35*f) - 2*d*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(7*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_532():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -64*a**3*cos(e + f*x)/(15*f*sqrt(a*sin(e + f*x) + a)) - 16*a**2*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(15*f) - 2*a*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_533():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)/(c + d*sin(e + f*x))
    F = -2*a**(sympy.S(5)/2)*(c - d)**2*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(d**(sympy.S(5)/2)*f*sqrt(c + d)) + 2*a**3*(3*c - 7*d)*cos(e + f*x)/(3*d**2*f*sqrt(a*sin(e + f*x) + a)) - 2*a**2*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(3*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_534():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)/(c + d*sin(e + f*x))**2
    F = a**(sympy.S(5)/2)*(c - d)*(3*c + 5*d)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(d**(sympy.S(5)/2)*f*(c + d)**(sympy.S(3)/2)) - a**3*(3*c + d)*cos(e + f*x)/(d**2*f*(c + d)*sqrt(a*sin(e + f*x) + a)) + a**2*(c - d)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(d*f*(c + d)*(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_535():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)/(c + d*sin(e + f*x))**3
    F = -a**(sympy.S(5)/2)*(3*c**2 + 10*c*d + 19*d**2)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(4*d**(sympy.S(5)/2)*f*(c + d)**(sympy.S(5)/2)) + 3*a**3*(c - d)*(c + 3*d)*cos(e + f*x)/(4*d**2*f*(c + d)**2*(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)) + a**2*(c - d)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(2*d*f*(c + d)*(c + d*sin(e + f*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_536():
    f = (c + d*sin(e + f*x))**3/sqrt(a*sin(e + f*x) + a)
    F = -2*d*(c + d*sin(e + f*x))**2*cos(e + f*x)/(5*f*sqrt(a*sin(e + f*x) + a)) - 4*d*(21*c**2 - 12*c*d + 7*d**2)*cos(e + f*x)/(15*f*sqrt(a*sin(e + f*x) + a)) - d**2*(18*c - 2*d)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(15*a*f) - sqrt(2)*(c - d)**3*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_537():
    f = (c + d*sin(e + f*x))**2/sqrt(a*sin(e + f*x) + a)
    F = -d*(12*c - 4*d)*cos(e + f*x)/(3*f*sqrt(a*sin(e + f*x) + a)) - 2*d**2*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(3*a*f) - sqrt(2)*(c - d)**2*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_538():
    f = (c + d*sin(e + f*x))/sqrt(a*sin(e + f*x) + a)
    F = -2*d*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)) - sqrt(2)*(c - d)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_539():
    f = 1/sqrt(a*sin(e + f*x) + a)
    F = -sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_540():
    f = 1/((c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a))
    F = 2*sqrt(d)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f*(c - d)*sqrt(c + d)) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f*(c - d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_541():
    f = 1/((c + d*sin(e + f*x))**2*sqrt(a*sin(e + f*x) + a))
    F = d*cos(e + f*x)/(f*(c + d*sin(e + f*x))*(c**2 - d**2)*sqrt(a*sin(e + f*x) + a)) + sqrt(d)*(3*c + d)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f*(c - d)**2*(c + d)**(sympy.S(3)/2)) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f*(c - d)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_542():
    f = 1/((c + d*sin(e + f*x))**3*sqrt(a*sin(e + f*x) + a))
    F = d*(7*c + d)*cos(e + f*x)/(4*f*(c + d*sin(e + f*x))*(c**2 - d**2)**2*sqrt(a*sin(e + f*x) + a)) + d*cos(e + f*x)/(f*(c + d*sin(e + f*x))**2*(2*c**2 - 2*d**2)*sqrt(a*sin(e + f*x) + a)) + sqrt(d)*(15*c**2 + 10*c*d + 7*d**2)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(4*sqrt(a)*f*(c - d)**3*(c + d)**(sympy.S(5)/2)) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f*(c - d)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_543():
    f = (c + d*sin(e + f*x))**3/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -(c - d)*(c + d*sin(e + f*x))**2*cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + d*(3*c**2 - 24*c*d + 13*d**2)*cos(e + f*x)/(3*a*f*sqrt(a*sin(e + f*x) + a)) + d**2*(3*c - 7*d)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(6*a**2*f) - sqrt(2)*(c - d)**2*(c + 11*d)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_544():
    f = (c + d*sin(e + f*x))**2/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -(c - d)*(c + d*sin(e + f*x))*cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + d*(c - 5*d)*cos(e + f*x)/(2*a*f*sqrt(a*sin(e + f*x) + a)) - sqrt(2)*(c - d)*(c + 7*d)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_545():
    f = (c + d*sin(e + f*x))/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -(c - d)*cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - sqrt(2)*(c + 3*d)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_546():
    f = (a*sin(e + f*x) + a)**(sympy.S(-3)/2)
    F = -cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_547():
    f = 1/((c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2))
    F = -cos(e + f*x)/(f*(2*c - 2*d)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - 2*d**(sympy.S(3)/2)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(a**(sympy.S(3)/2)*f*(c - d)**2*sqrt(c + d)) - sqrt(2)*(c - 5*d)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f*(c - d)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_548():
    f = 1/((c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**(sympy.S(3)/2))
    F = -cos(e + f*x)/(f*(c + d*sin(e + f*x))*(2*c - 2*d)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - d*(c + 3*d)*cos(e + f*x)/(2*a*f*(c - d)**2*(c + d)*(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)) - d**(sympy.S(3)/2)*(5*c + 3*d)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(a**(sympy.S(3)/2)*f*(c - d)**3*(c + d)**(sympy.S(3)/2)) - sqrt(2)*(c - 9*d)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f*(c - d)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_549():
    f = 1/((c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a)**(sympy.S(3)/2))
    F = -cos(e + f*x)/(f*(c + d*sin(e + f*x))**2*(2*c - 2*d)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - d*(c + 2*d)*cos(e + f*x)/(2*a*f*(c - d)**2*(c + d)*(c + d*sin(e + f*x))**2*sqrt(a*sin(e + f*x) + a)) - d*(c + 7*d)*(2*c + d)*cos(e + f*x)/(4*a*f*(c - d)**3*(c + d)**2*(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)) - d**(sympy.S(3)/2)*(35*c**2 + 42*c*d + 19*d**2)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f*(c - d)**4*(c + d)**(sympy.S(5)/2)) - sqrt(2)*(c - 13*d)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f*(c - d)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_550():
    f = (c + d*sin(e + f*x))**3/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -(c - d)*(c + d*sin(e + f*x))**2*cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) - (c - d)**2*(3*c + 13*d)*cos(e + f*x)/(16*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + d**2*(c - 9*d)*cos(e + f*x)/(4*a**2*f*sqrt(a*sin(e + f*x) + a)) - sqrt(2)*(3*c - 3*d)*(c**2 + 6*c*d + 25*d**2)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(32*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_551():
    f = (c + d*sin(e + f*x))**2/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -(c - d)*(c + d*sin(e + f*x))*cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) - (c + 3*d)*(3*c - 3*d)*cos(e + f*x)/(16*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - sqrt(2)*(3*c**2 + 10*c*d + 19*d**2)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(32*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_552():
    f = (c + d*sin(e + f*x))/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -(c - d)*cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) - (3*c + 5*d)*cos(e + f*x)/(16*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - sqrt(2)*(3*c + 5*d)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(32*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_553():
    f = (a*sin(e + f*x) + a)**(sympy.S(-5)/2)
    F = -cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) - 3*cos(e + f*x)/(16*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - 3*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(32*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_554():
    f = 1/((c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2))
    F = -cos(e + f*x)/(f*(4*c - 4*d)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) - (3*c - 11*d)*cos(e + f*x)/(16*a*f*(c - d)**2*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + 2*d**(sympy.S(5)/2)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(a**(sympy.S(5)/2)*f*(c - d)**3*sqrt(c + d)) - sqrt(2)*(3*c**2 - 14*c*d + 43*d**2)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(32*a**(sympy.S(5)/2)*f*(c - d)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_555():
    f = 1/((c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**(sympy.S(5)/2))
    F = -cos(e + f*x)/(f*(c + d*sin(e + f*x))*(4*c - 4*d)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) - (3*c - 15*d)*cos(e + f*x)/(16*a*f*(c - d)**2*(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - d*(c - 7*d)*(3*c + 5*d)*cos(e + f*x)/(16*a**2*f*(c - d)**3*(c + d)*(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)) + d**(sympy.S(5)/2)*(7*c + 5*d)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(a**(sympy.S(5)/2)*f*(c - d)**4*(c + d)**(sympy.S(3)/2)) - sqrt(2)*(3*c**2 - 22*c*d + 115*d**2)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(32*a**(sympy.S(5)/2)*f*(c - d)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_556():
    f = 1/((c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a)**(sympy.S(5)/2))
    F = -cos(e + f*x)/(f*(c + d*sin(e + f*x))**2*(4*c - 4*d)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) - (3*c - 19*d)*cos(e + f*x)/(16*a*f*(c - d)**2*(c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - d*(3*c**2 - 20*c*d - 31*d**2)*cos(e + f*x)/(16*a**2*f*(c - d)**3*(c + d)*(c + d*sin(e + f*x))**2*sqrt(a*sin(e + f*x) + a)) - 3*d*(c + 3*d)*(c**2 - 10*c*d - 7*d**2)*cos(e + f*x)/(16*a**2*f*(c - d)**4*(c + d)**2*(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)) + 3*d**(sympy.S(5)/2)*(21*c**2 + 30*c*d + 13*d**2)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(5)/2)*f*(c - d)**5*(c + d)**(sympy.S(5)/2)) - sqrt(2)*(3*c**2 - 30*c*d + 219*d**2)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(32*a**(sympy.S(5)/2)*f*(c - d)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_557():
    f = (c + d*sin(e + f*x))**(sympy.S(5)/2)*sqrt(a*sin(e + f*x) + a)
    F = -5*sqrt(a)*(c + d)**3*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(8*sqrt(d)*f) - 5*a*(c + d)**2*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(8*f*sqrt(a*sin(e + f*x) + a)) - 5*a*(c + d)*(c + d*sin(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(12*f*sqrt(a*sin(e + f*x) + a)) - a*(c + d*sin(e + f*x))**(sympy.S(5)/2)*cos(e + f*x)/(3*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_558():
    f = (c + d*sin(e + f*x))**(sympy.S(3)/2)*sqrt(a*sin(e + f*x) + a)
    F = -3*sqrt(a)*(c + d)**2*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(4*sqrt(d)*f) - 3*a*(c + d)*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(4*f*sqrt(a*sin(e + f*x) + a)) - a*(c + d*sin(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(2*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_559():
    f = sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)
    F = -sqrt(a)*(c + d)*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(sqrt(d)*f) - a*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_560():
    f = sqrt(a*sin(e + f*x) + a)/sqrt(c + d*sin(e + f*x))
    F = -2*sqrt(a)*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_561():
    f = sqrt(a*sin(e + f*x) + a)/(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = -2*a*cos(e + f*x)/(f*(c + d)*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_562():
    f = sqrt(a*sin(e + f*x) + a)/(c + d*sin(e + f*x))**(sympy.S(5)/2)
    F = -2*a*cos(e + f*x)/(f*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(3*c + 3*d)*sqrt(a*sin(e + f*x) + a)) - 4*a*cos(e + f*x)/(3*f*(c + d)**2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_563():
    f = sqrt(a*sin(e + f*x) + a)/(c + d*sin(e + f*x))**(sympy.S(7)/2)
    F = -2*a*cos(e + f*x)/(f*(c + d*sin(e + f*x))**(sympy.S(5)/2)*(5*c + 5*d)*sqrt(a*sin(e + f*x) + a)) - 8*a*cos(e + f*x)/(15*f*(c + d)**2*(c + d*sin(e + f*x))**(sympy.S(3)/2)*sqrt(a*sin(e + f*x) + a)) - 16*a*cos(e + f*x)/(15*f*(c + d)**3*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_564():
    f = (c + d*sin(e + f*x))**(sympy.S(5)/2)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = 5*a**(sympy.S(3)/2)*(c - 15*d)*(c + d)**3*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(64*d**(sympy.S(3)/2)*f) + 5*a**2*(c - 15*d)*(c + d)**2*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(64*d*f*sqrt(a*sin(e + f*x) + a)) + 5*a**2*(c - 15*d)*(c + d)*(c + d*sin(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(96*d*f*sqrt(a*sin(e + f*x) + a)) + a**2*(c - 15*d)*(c + d*sin(e + f*x))**(sympy.S(5)/2)*cos(e + f*x)/(24*d*f*sqrt(a*sin(e + f*x) + a)) - a**2*(c + d*sin(e + f*x))**(sympy.S(7)/2)*cos(e + f*x)/(4*d*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_565():
    f = (c + d*sin(e + f*x))**(sympy.S(3)/2)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = a**(sympy.S(3)/2)*(c - 11*d)*(c + d)**2*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(8*d**(sympy.S(3)/2)*f) + a**2*(c - 11*d)*(c + d)*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(8*d*f*sqrt(a*sin(e + f*x) + a)) + a**2*(c - 11*d)*(c + d*sin(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(12*d*f*sqrt(a*sin(e + f*x) + a)) - a**2*(c + d*sin(e + f*x))**(sympy.S(5)/2)*cos(e + f*x)/(3*d*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_566():
    f = sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = a**(sympy.S(3)/2)*(c - 7*d)*(c + d)*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(4*d**(sympy.S(3)/2)*f) + a**2*(c - 7*d)*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(4*d*f*sqrt(a*sin(e + f*x) + a)) - a**2*(c + d*sin(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(2*d*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_567():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)/sqrt(c + d*sin(e + f*x))
    F = a**(sympy.S(3)/2)*(c - 3*d)*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(d**(sympy.S(3)/2)*f) - a**2*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(d*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_568():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)/(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = -2*a**(sympy.S(3)/2)*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(d**(sympy.S(3)/2)*f) + 2*a**2*(c - d)*cos(e + f*x)/(d*f*(c + d)*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_569():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)/(c + d*sin(e + f*x))**(sympy.S(5)/2)
    F = 2*a**2*(c - d)*cos(e + f*x)/(3*d*f*(c + d)*(c + d*sin(e + f*x))**(sympy.S(3)/2)*sqrt(a*sin(e + f*x) + a)) - 2*a**2*(c + 5*d)*cos(e + f*x)/(3*d*f*(c + d)**2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_570():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)/(c + d*sin(e + f*x))**(sympy.S(7)/2)
    F = 2*a**2*(c - d)*cos(e + f*x)/(5*d*f*(c + d)*(c + d*sin(e + f*x))**(sympy.S(5)/2)*sqrt(a*sin(e + f*x) + a)) - 2*a**2*(c + 9*d)*cos(e + f*x)/(15*d*f*(c + d)**2*(c + d*sin(e + f*x))**(sympy.S(3)/2)*sqrt(a*sin(e + f*x) + a)) - 4*a**2*(c + 9*d)*cos(e + f*x)/(15*d*f*(c + d)**3*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_571():
    f = (a*sin(e + f*x) + a)**(sympy.S(3)/2)/(c + d*sin(e + f*x))**(sympy.S(9)/2)
    F = 2*a**2*(c - d)*cos(e + f*x)/(7*d*f*(c + d)*(c + d*sin(e + f*x))**(sympy.S(7)/2)*sqrt(a*sin(e + f*x) + a)) - 2*a**2*(c + 13*d)*cos(e + f*x)/(35*d*f*(c + d)**2*(c + d*sin(e + f*x))**(sympy.S(5)/2)*sqrt(a*sin(e + f*x) + a)) - 8*a**2*(c + 13*d)*cos(e + f*x)/(105*d*f*(c + d)**3*(c + d*sin(e + f*x))**(sympy.S(3)/2)*sqrt(a*sin(e + f*x) + a)) - 16*a**2*(c + 13*d)*cos(e + f*x)/(105*d*f*(c + d)**4*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_572():
    f = (c + d*sin(e + f*x))**(sympy.S(5)/2)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -a**(sympy.S(5)/2)*(c + d)**3*(3*c**2 - 34*c*d + 283*d**2)*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(128*d**(sympy.S(5)/2)*f) + 3*a**3*(c - 7*d)*(c + d*sin(e + f*x))**(sympy.S(7)/2)*cos(e + f*x)/(40*d**2*f*sqrt(a*sin(e + f*x) + a)) - a**3*(c + d)**2*sqrt(c + d*sin(e + f*x))*(3*c**2 - 34*c*d + 283*d**2)*cos(e + f*x)/(128*d**2*f*sqrt(a*sin(e + f*x) + a)) - a**3*(c + d)*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(3*c**2 - 34*c*d + 283*d**2)*cos(e + f*x)/(192*d**2*f*sqrt(a*sin(e + f*x) + a)) - a**3*(c + d*sin(e + f*x))**(sympy.S(5)/2)*(3*c**2 - 34*c*d + 283*d**2)*cos(e + f*x)/(240*d**2*f*sqrt(a*sin(e + f*x) + a)) - a**2*(c + d*sin(e + f*x))**(sympy.S(7)/2)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(5*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_573():
    f = (c + d*sin(e + f*x))**(sympy.S(3)/2)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -a**(sympy.S(5)/2)*(c + d)**2*(3*c**2 - 26*c*d + 163*d**2)*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(64*d**(sympy.S(5)/2)*f) - a**3*(c + d)*sqrt(c + d*sin(e + f*x))*(3*c**2 - 26*c*d + 163*d**2)*cos(e + f*x)/(64*d**2*f*sqrt(a*sin(e + f*x) + a)) + a**3*(c + d*sin(e + f*x))**(sympy.S(5)/2)*(3*c - 17*d)*cos(e + f*x)/(24*d**2*f*sqrt(a*sin(e + f*x) + a)) - a**3*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(3*c**2 - 26*c*d + 163*d**2)*cos(e + f*x)/(96*d**2*f*sqrt(a*sin(e + f*x) + a)) - a**2*(c + d*sin(e + f*x))**(sympy.S(5)/2)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(4*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_574():
    f = sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -a**(sympy.S(5)/2)*(c + d)*(c**2 - 6*c*d + 25*d**2)*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(8*d**(sympy.S(5)/2)*f) + a**3*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(3*c - 13*d)*cos(e + f*x)/(12*d**2*f*sqrt(a*sin(e + f*x) + a)) - a**3*sqrt(c + d*sin(e + f*x))*(c**2 - 6*c*d + 25*d**2)*cos(e + f*x)/(8*d**2*f*sqrt(a*sin(e + f*x) + a)) - a**2*(c + d*sin(e + f*x))**(sympy.S(3)/2)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(3*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_575():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)/sqrt(c + d*sin(e + f*x))
    F = -a**(sympy.S(5)/2)*(3*c**2 - 10*c*d + 19*d**2)*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(4*d**(sympy.S(5)/2)*f) + 3*a**3*(c - 3*d)*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(4*d**2*f*sqrt(a*sin(e + f*x) + a)) - a**2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(2*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_576():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)/(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = a**(sympy.S(5)/2)*(3*c - 5*d)*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(d**(sympy.S(5)/2)*f) - a**3*sqrt(c + d*sin(e + f*x))*(3*c - d)*cos(e + f*x)/(d**2*f*(c + d)*sqrt(a*sin(e + f*x) + a)) + 2*a**2*(c - d)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(d*f*(c + d)*sqrt(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_577():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)/(c + d*sin(e + f*x))**(sympy.S(5)/2)
    F = -2*a**(sympy.S(5)/2)*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(d**(sympy.S(5)/2)*f) + 2*a**3*(c - d)*(3*c + 7*d)*cos(e + f*x)/(3*d**2*f*(c + d)**2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)) + 2*a**2*(c - d)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(3*d*f*(c + d)*(c + d*sin(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_578():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)/(c + d*sin(e + f*x))**(sympy.S(7)/2)
    F = 2*a**3*(c - d)*(3*c + 11*d)*cos(e + f*x)/(15*d**2*f*(c + d)**2*(c + d*sin(e + f*x))**(sympy.S(3)/2)*sqrt(a*sin(e + f*x) + a)) - 2*a**3*(3*c**2 + 14*c*d + 43*d**2)*cos(e + f*x)/(15*d**2*f*(c + d)**3*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)) + 2*a**2*(c - d)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(5*d*f*(c + d)*(c + d*sin(e + f*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_579():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)/(c + d*sin(e + f*x))**(sympy.S(9)/2)
    F = 6*a**3*(c - d)*(c + 5*d)*cos(e + f*x)/(35*d**2*f*(c + d)**2*(c + d*sin(e + f*x))**(sympy.S(5)/2)*sqrt(a*sin(e + f*x) + a)) - 2*a**3*(3*c**2 + 22*c*d + 115*d**2)*cos(e + f*x)/(105*d**2*f*(c + d)**3*(c + d*sin(e + f*x))**(sympy.S(3)/2)*sqrt(a*sin(e + f*x) + a)) - 4*a**3*(3*c**2 + 22*c*d + 115*d**2)*cos(e + f*x)/(105*d**2*f*(c + d)**4*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)) + 2*a**2*(c - d)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(7*d*f*(c + d)*(c + d*sin(e + f*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_580():
    f = (a*sin(e + f*x) + a)**(sympy.S(5)/2)/(c + d*sin(e + f*x))**(sympy.S(11)/2)
    F = 2*a**3*(c - d)*(3*c + 19*d)*cos(e + f*x)/(63*d**2*f*(c + d)**2*(c + d*sin(e + f*x))**(sympy.S(7)/2)*sqrt(a*sin(e + f*x) + a)) - 2*a**3*(c**2 + 10*c*d + 73*d**2)*cos(e + f*x)/(105*d**2*f*(c + d)**3*(c + d*sin(e + f*x))**(sympy.S(5)/2)*sqrt(a*sin(e + f*x) + a)) - 8*a**3*(c**2 + 10*c*d + 73*d**2)*cos(e + f*x)/(315*d**2*f*(c + d)**4*(c + d*sin(e + f*x))**(sympy.S(3)/2)*sqrt(a*sin(e + f*x) + a)) - 16*a**3*(c**2 + 10*c*d + 73*d**2)*cos(e + f*x)/(315*d**2*f*(c + d)**5*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)) + 2*a**2*(c - d)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(9*d*f*(c + d)*(c + d*sin(e + f*x))**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_581():
    f = (c + d*sin(e + f*x))**(sympy.S(5)/2)/sqrt(a*sin(e + f*x) + a)
    F = -d*(c + d*sin(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(2*f*sqrt(a*sin(e + f*x) + a)) - d*sqrt(c + d*sin(e + f*x))*(7*c - d)*cos(e + f*x)/(4*f*sqrt(a*sin(e + f*x) + a)) - sqrt(d)*(15*c**2 - 10*c*d + 7*d**2)*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(4*sqrt(a)*f) - sqrt(2)*(c - d)**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(c - d)*cos(e + f*x)/(2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_582():
    f = (c + d*sin(e + f*x))**(sympy.S(3)/2)/sqrt(a*sin(e + f*x) + a)
    F = -d*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)) - sqrt(d)*(3*c - d)*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f) - sqrt(2)*(c - d)**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(c - d)*cos(e + f*x)/(2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_583():
    f = sqrt(c + d*sin(e + f*x))/sqrt(a*sin(e + f*x) + a)
    F = -2*sqrt(d)*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f) - sqrt(2)*sqrt(c - d)*atanh(sqrt(2)*sqrt(a)*sqrt(c - d)*cos(e + f*x)/(2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_584():
    f = 1/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a))
    F = -sqrt(2)*atanh(sqrt(2)*sqrt(a)*sqrt(c - d)*cos(e + f*x)/(2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f*sqrt(c - d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_585():
    f = 1/((c + d*sin(e + f*x))**(sympy.S(3)/2)*sqrt(a*sin(e + f*x) + a))
    F = 2*d*cos(e + f*x)/(f*sqrt(c + d*sin(e + f*x))*(c**2 - d**2)*sqrt(a*sin(e + f*x) + a)) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*sqrt(c - d)*cos(e + f*x)/(2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f*(c - d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_586():
    f = 1/((c + d*sin(e + f*x))**(sympy.S(5)/2)*sqrt(a*sin(e + f*x) + a))
    F = 2*d*(5*c + d)*cos(e + f*x)/(3*f*sqrt(c + d*sin(e + f*x))*(c**2 - d**2)**2*sqrt(a*sin(e + f*x) + a)) + 2*d*cos(e + f*x)/(f*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(3*c**2 - 3*d**2)*sqrt(a*sin(e + f*x) + a)) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*sqrt(c - d)*cos(e + f*x)/(2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f*(c - d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_587():
    f = (c + d*sin(e + f*x))**(sympy.S(5)/2)/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -(c - d)*(c + d*sin(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + d*(c - 3*d)*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(2*a*f*sqrt(a*sin(e + f*x) + a)) - d**(sympy.S(3)/2)*(5*c - 3*d)*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(a**(sympy.S(3)/2)*f) - sqrt(2)*(c - d)**(sympy.S(3)/2)*(c + 9*d)*atanh(sqrt(2)*sqrt(a)*sqrt(c - d)*cos(e + f*x)/(2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_588():
    f = (c + d*sin(e + f*x))**(sympy.S(3)/2)/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -(c - d)*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - 2*d**(sympy.S(3)/2)*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(a**(sympy.S(3)/2)*f) - sqrt(2)*sqrt(c - d)*(c + 5*d)*atanh(sqrt(2)*sqrt(a)*sqrt(c - d)*cos(e + f*x)/(2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_589():
    f = sqrt(c + d*sin(e + f*x))/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - sqrt(2)*(c + d)*atanh(sqrt(2)*sqrt(a)*sqrt(c - d)*cos(e + f*x)/(2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f*sqrt(c - d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_590():
    f = 1/(sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2))
    F = -sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(f*(2*c - 2*d)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - sqrt(2)*(c - 3*d)*atanh(sqrt(2)*sqrt(a)*sqrt(c - d)*cos(e + f*x)/(2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f*(c - d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_591():
    f = 1/((c + d*sin(e + f*x))**(sympy.S(3)/2)*(a*sin(e + f*x) + a)**(sympy.S(3)/2))
    F = -cos(e + f*x)/(f*sqrt(c + d*sin(e + f*x))*(2*c - 2*d)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - d*(c + 5*d)*cos(e + f*x)/(2*a*f*(c - d)**2*(c + d)*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)) - sqrt(2)*(c - 7*d)*atanh(sqrt(2)*sqrt(a)*sqrt(c - d)*cos(e + f*x)/(2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f*(c - d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_592():
    f = 1/((c + d*sin(e + f*x))**(sympy.S(5)/2)*(a*sin(e + f*x) + a)**(sympy.S(3)/2))
    F = -cos(e + f*x)/(f*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(2*c - 2*d)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - d*(3*c + 7*d)*cos(e + f*x)/(6*a*f*(c - d)**2*(c + d)*(c + d*sin(e + f*x))**(sympy.S(3)/2)*sqrt(a*sin(e + f*x) + a)) - d*(3*c**2 + 38*c*d + 19*d**2)*cos(e + f*x)/(6*a*f*(c - d)**3*(c + d)**2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)) - sqrt(2)*(c - 11*d)*atanh(sqrt(2)*sqrt(a)*sqrt(c - d)*cos(e + f*x)/(2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f*(c - d)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_593():
    f = (c + d*sin(e + f*x))**(sympy.S(5)/2)/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -(c - d)*(c + d*sin(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) - (c - d)*sqrt(c + d*sin(e + f*x))*(3*c + 11*d)*cos(e + f*x)/(16*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - 2*d**(sympy.S(5)/2)*atan(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(a**(sympy.S(5)/2)*f) - sqrt(2)*sqrt(c - d)*(3*c**2 + 14*c*d + 43*d**2)*atanh(sqrt(2)*sqrt(a)*sqrt(c - d)*cos(e + f*x)/(2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(32*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_594():
    f = (c + d*sin(e + f*x))**(sympy.S(3)/2)/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -(c - d)*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) - sqrt(c + d*sin(e + f*x))*(3*c + 7*d)*cos(e + f*x)/(16*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - 3*sqrt(2)*(c + d)**2*atanh(sqrt(2)*sqrt(a)*sqrt(c - d)*cos(e + f*x)/(2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(32*a**(sympy.S(5)/2)*f*sqrt(c - d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_595():
    f = sqrt(c + d*sin(e + f*x))/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) - sqrt(c + d*sin(e + f*x))*(3*c - d)*cos(e + f*x)/(16*a*f*(c - d)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - sqrt(2)*(c + d)*(3*c - 5*d)*atanh(sqrt(2)*sqrt(a)*sqrt(c - d)*cos(e + f*x)/(2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(32*a**(sympy.S(5)/2)*f*(c - d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_596():
    f = 1/(sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2))
    F = -sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(f*(4*c - 4*d)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) - sqrt(c + d*sin(e + f*x))*(3*c - 9*d)*cos(e + f*x)/(16*a*f*(c - d)**2*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - sqrt(2)*(3*c**2 - 10*c*d + 19*d**2)*atanh(sqrt(2)*sqrt(a)*sqrt(c - d)*cos(e + f*x)/(2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(32*a**(sympy.S(5)/2)*f*(c - d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_597():
    f = 1/((c + d*sin(e + f*x))**(sympy.S(3)/2)*(a*sin(e + f*x) + a)**(sympy.S(5)/2))
    F = -cos(e + f*x)/(f*sqrt(c + d*sin(e + f*x))*(4*c - 4*d)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) - (3*c - 13*d)*cos(e + f*x)/(16*a*f*(c - d)**2*sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - d*(c - 7*d)*(3*c + 7*d)*cos(e + f*x)/(16*a**2*f*(c - d)**3*(c + d)*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)) - sqrt(2)*(3*c**2 - 18*c*d + 75*d**2)*atanh(sqrt(2)*sqrt(a)*sqrt(c - d)*cos(e + f*x)/(2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(32*a**(sympy.S(5)/2)*f*(c - d)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_598():
    f = 1/((c + d*sin(e + f*x))**(sympy.S(5)/2)*(a*sin(e + f*x) + a)**(sympy.S(5)/2))
    F = -cos(e + f*x)/(f*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(4*c - 4*d)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) - (3*c - 17*d)*cos(e + f*x)/(16*a*f*(c - d)**2*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - d*(9*c**2 - 54*c*d - 95*d**2)*cos(e + f*x)/(48*a**2*f*(c - d)**3*(c + d)*(c + d*sin(e + f*x))**(sympy.S(3)/2)*sqrt(a*sin(e + f*x) + a)) - d*(9*c**3 - 57*c**2*d - 493*c*d**2 - 299*d**3)*cos(e + f*x)/(48*a**2*f*(c - d)**4*(c + d)**2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)) - sqrt(2)*(3*c**2 - 26*c*d + 163*d**2)*atanh(sqrt(2)*sqrt(a)*sqrt(c - d)*cos(e + f*x)/(2*sqrt(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)))/(32*a**(sympy.S(5)/2)*f*(c - d)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_599():
    f = (c + d*sin(e + f*x))**n*(a*sin(e + f*x) + a)**m
    F = sqrt(2)*(c + d*sin(e + f*x))**n*(a*sin(e + f*x) + a)**m*cos(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, -n, m + sympy.S(3)/2, sin(e + f*x)/2 + sympy.S.Half, -d*(sin(e + f*x) + 1)/(c - d))/(f*((c + d*sin(e + f*x))/(c - d))**n*sqrt(1 - sin(e + f*x))*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_600():
    f = (c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a)**m
    F = -2**(m + sympy.S.Half)*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(-m + sympy.S(-1)/2)*(c**3*(m**3 + 6*m**2 + 11*m + 6) + 3*c**2*d*m*(m**2 + 5*m + 6) + 3*c*d**2*(m**3 + 4*m**2 + 4*m + 3) + d**3*m*(m**2 + 3*m + 5))*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - sin(e + f*x)/2)/(f*(m + 1)*(m + 2)*(m + 3)) - d*(c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**m*cos(e + f*x)/(f*(m + 3)) - d*(a*sin(e + f*x) + a)**m*(2*c**2*(m**2 + 6*m + 8) - c*d*(-2*m**2 - 3*m + 5) + d**2*(m + 4))*cos(e + f*x)/(f*(m + 1)*(m + 2)*(m + 3)) - d**2*(a*sin(e + f*x) + a)**(m + 1)*(c*(m + 5) + d*m)*cos(e + f*x)/(a*f*(m + 2)*(m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_601():
    f = (c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**m
    F = -2**(m + sympy.S.Half)*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(-m + sympy.S(-1)/2)*(c**2*(m**2 + 3*m + 2) + 2*c*d*m*(m + 2) + d**2*(m**2 + m + 1))*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - sin(e + f*x)/2)/(f*(m + 1)*(m + 2)) + d*(a*sin(e + f*x) + a)**m*(-2*c*(m + 2) + d)*cos(e + f*x)/(f*(m + 1)*(m + 2)) - d**2*(a*sin(e + f*x) + a)**(m + 1)*cos(e + f*x)/(a*f*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_602():
    f = (c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**m
    F = -2**(m + sympy.S.Half)*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(-m + sympy.S(-1)/2)*(c*m + c + d*m)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - sin(e + f*x)/2)/(f*(m + 1)) - d*(a*sin(e + f*x) + a)**m*cos(e + f*x)/(f*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_603():
    f = (a*sin(e + f*x) + a)**m
    F = -2**(m + sympy.S.Half)*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(-m + sympy.S(-1)/2)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - sin(e + f*x)/2)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_604():
    f = (a*sin(e + f*x) + a)**m/(c + d*sin(e + f*x))
    F = sqrt(2)*(a*sin(e + f*x) + a)**m*cos(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, 1, m + sympy.S(3)/2, sin(e + f*x)/2 + sympy.S.Half, -d*(sin(e + f*x) + 1)/(c - d))/(f*sqrt(1 - sin(e + f*x))*(c - d)*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_605():
    f = (a*sin(e + f*x) + a)**m/(c + d*sin(e + f*x))**2
    F = sqrt(2)*(a*sin(e + f*x) + a)**m*cos(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, 2, m + sympy.S(3)/2, sin(e + f*x)/2 + sympy.S.Half, -d*(sin(e + f*x) + 1)/(c - d))/(f*sqrt(1 - sin(e + f*x))*(c - d)**2*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_606():
    f = (a*sin(e + f*x) + a)**m/(c + d*sin(e + f*x))**3
    F = sqrt(2)*(a*sin(e + f*x) + a)**m*cos(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, 3, m + sympy.S(3)/2, sin(e + f*x)/2 + sympy.S.Half, -d*(sin(e + f*x) + 1)/(c - d))/(f*sqrt(1 - sin(e + f*x))*(c - d)**3*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_607():
    f = (c + d*sin(e + f*x))**(sympy.S(5)/2)*(a*sin(e + f*x) + a)**m
    F = sqrt(2)*(c - d)**2*sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**m*cos(e + f*x)*appellf1(m + sympy.S.Half, sympy.S(-5)/2, sympy.S.Half, m + sympy.S(3)/2, -d*(sin(e + f*x) + 1)/(c - d), sin(e + f*x)/2 + sympy.S.Half)/(f*sqrt((c + d*sin(e + f*x))/(c - d))*sqrt(1 - sin(e + f*x))*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_608():
    f = (c + d*sin(e + f*x))**(sympy.S(3)/2)*(a*sin(e + f*x) + a)**m
    F = sqrt(2)*(c - d)*sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**m*cos(e + f*x)*appellf1(m + sympy.S.Half, sympy.S(-3)/2, sympy.S.Half, m + sympy.S(3)/2, -d*(sin(e + f*x) + 1)/(c - d), sin(e + f*x)/2 + sympy.S.Half)/(f*sqrt((c + d*sin(e + f*x))/(c - d))*sqrt(1 - sin(e + f*x))*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_609():
    f = sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**m
    F = sqrt(2)*sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**m*cos(e + f*x)*appellf1(m + sympy.S.Half, sympy.S(-1)/2, sympy.S.Half, m + sympy.S(3)/2, -d*(sin(e + f*x) + 1)/(c - d), sin(e + f*x)/2 + sympy.S.Half)/(f*sqrt((c + d*sin(e + f*x))/(c - d))*sqrt(1 - sin(e + f*x))*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_610():
    f = (a*sin(e + f*x) + a)**m/sqrt(c + d*sin(e + f*x))
    F = sqrt(2)*sqrt((c + d*sin(e + f*x))/(c - d))*(a*sin(e + f*x) + a)**m*cos(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, sympy.S.Half, m + sympy.S(3)/2, -d*(sin(e + f*x) + 1)/(c - d), sin(e + f*x)/2 + sympy.S.Half)/(f*sqrt(1 - sin(e + f*x))*sqrt(c + d*sin(e + f*x))*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_611():
    f = (a*sin(e + f*x) + a)**m/(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = sqrt(2)*sqrt((c + d*sin(e + f*x))/(c - d))*(a*sin(e + f*x) + a)**m*cos(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, sympy.S(3)/2, m + sympy.S(3)/2, sin(e + f*x)/2 + sympy.S.Half, -d*(sin(e + f*x) + 1)/(c - d))/(f*sqrt(1 - sin(e + f*x))*(c - d)*sqrt(c + d*sin(e + f*x))*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_612():
    f = (a*sin(e + f*x) + a)**m/(c + d*sin(e + f*x))**(sympy.S(5)/2)
    F = sqrt(2)*sqrt((c + d*sin(e + f*x))/(c - d))*(a*sin(e + f*x) + a)**m*cos(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, sympy.S(5)/2, m + sympy.S(3)/2, sin(e + f*x)/2 + sympy.S.Half, -d*(sin(e + f*x) + 1)/(c - d))/(f*sqrt(1 - sin(e + f*x))*(c - d)**2*sqrt(c + d*sin(e + f*x))*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_613():
    f = (sin(e + f*x) + 1)**m*(3*sin(e + f*x) + 3)**(-m - 1)
    F = -3**(-m - 1)*cos(e + f*x)/(f*(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_614():
    f = (sin(e + f*x) + 1)**m*(2*sin(e + f*x) + 3)**(-m - 1)
    F = -2**(m + sympy.S.Half)*5**(-m + sympy.S(-1)/2)*((sin(e + f*x) + 1)/(2*sin(e + f*x) + 3))**(sympy.S.Half - m)*(sin(e + f*x) + 1)**(m - 1)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), (1 - sin(e + f*x))/(4*sin(e + f*x) + 6))/(f*(2*sin(e + f*x) + 3)**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_615():
    f = (sin(e + f*x) + 1)**m*(sin(e + f*x) + 3)**(-m - 1)
    F = -2**(-m + sympy.S(-1)/2)*((sin(e + f*x) + 1)/(sin(e + f*x) + 3))**(sympy.S.Half - m)*(sin(e + f*x) + 1)**(m - 1)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), (1 - sin(e + f*x))/(sin(e + f*x) + 3))/(f*(sin(e + f*x) + 3)**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_616():
    f = 3**(-m - 1)*(sin(e + f*x) + 1)**m
    F = -2**(m + sympy.S.Half)*3**(-m - 1)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - sin(e + f*x)/2)/(f*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_617():
    f = (3 - sin(e + f*x))**(-m - 1)*(sin(e + f*x) + 1)**m
    F = -((3 - sin(e + f*x))/(sin(e + f*x) + 1))**(m + 1)*(3 - sin(e + f*x))**(-m - 1)*(sin(e + f*x) + 1)**m*cos(e + f*x)*hyper((sympy.S.Half, m + 1), (sympy.S(3)/2,), -(2 - 2*sin(e + f*x))/(sin(e + f*x) + 1))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_618():
    f = (3 - 2*sin(e + f*x))**(-m - 1)*(sin(e + f*x) + 1)**m
    F = sqrt(5)*sqrt(-(1 - sin(e + f*x))/(sin(e + f*x) + 1))*(sin(e + f*x) + 1)**m*cos(e + f*x)*hyper((sympy.S.Half, -m), (1 - m,), (6 - 4*sin(e + f*x))/(sin(e + f*x) + 1))/(5*f*m*(1 - sin(e + f*x))*(3 - 2*sin(e + f*x))**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_619():
    f = (3 - 3*sin(e + f*x))**(-m - 1)*(sin(e + f*x) + 1)**m
    F = (3 - 3*sin(e + f*x))**(-m - 1)*(sin(e + f*x) + 1)**m*cos(e + f*x)/(f*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_620():
    f = (a*sin(e + f*x) + a)**m*(3*sin(e + f*x) + 3)**(-m - 1)
    F = -(a*sin(e + f*x) + a)**m*(3*sin(e + f*x) + 3)**(-m - 1)*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_621():
    f = 3**(-m - 1)*(a*sin(e + f*x) + a)**m
    F = -2**(m + sympy.S.Half)*3**(-m - 1)*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(-m + sympy.S(-1)/2)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - sin(e + f*x)/2)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_622():
    f = (3 - 3*sin(e + f*x))**(-m - 1)*(a*sin(e + f*x) + a)**m
    F = (3 - 3*sin(e + f*x))**(-m - 1)*(a*sin(e + f*x) + a)**m*cos(e + f*x)/(f*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_623():
    f = (3 - 4*sin(e + f*x))**(-m - 1)*(a*sin(e + f*x) + a)**m
    F = sqrt(7)*sqrt((1 - sin(e + f*x))/(sin(e + f*x) + 1))*(a*sin(e + f*x) + a)**m*cos(e + f*x)*hyper((sympy.S.Half, -m), (1 - m,), -(6 - 8*sin(e + f*x))/(sin(e + f*x) + 1))/(7*f*m*(1 - sin(e + f*x))*(3 - 4*sin(e + f*x))**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_624():
    f = (3 - 5*sin(e + f*x))**(-m - 1)*(a*sin(e + f*x) + a)**m
    F = sqrt((1 - sin(e + f*x))/(sin(e + f*x) + 1))*(a*sin(e + f*x) + a)**m*cos(e + f*x)*hyper((sympy.S.Half, -m), (1 - m,), -(3 - 5*sin(e + f*x))/(sin(e + f*x) + 1))/(4*f*m*(1 - sin(e + f*x))*(3 - 5*sin(e + f*x))**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_625():
    f = (a*sin(e + f*x) + a)**m*(3*sin(e + f*x) - 3)**(-m - 1)
    F = (a*sin(e + f*x) + a)**m*(3*sin(e + f*x) - 3)**(-m - 1)*cos(e + f*x)/(f*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_626():
    f = (a*sin(e + f*x) + a)**m*(2*sin(e + f*x) - 3)**(-m - 1)
    F = -sqrt(5)*sqrt(-(1 - sin(e + f*x))/(sin(e + f*x) + 1))*(a*sin(e + f*x) + a)**m*cos(e + f*x)*hyper((sympy.S.Half, -m), (1 - m,), (6 - 4*sin(e + f*x))/(sin(e + f*x) + 1))/(5*f*m*(1 - sin(e + f*x))*(2*sin(e + f*x) - 3)**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_627():
    f = (a*sin(e + f*x) + a)**m*(sin(e + f*x) - 3)**(-m - 1)
    F = -sqrt(2)*sqrt(-(1 - sin(e + f*x))/(sin(e + f*x) + 1))*(a*sin(e + f*x) + a)**m*cos(e + f*x)*hyper((sympy.S.Half, -m), (1 - m,), (3 - sin(e + f*x))/(sin(e + f*x) + 1))/(4*f*m*(1 - sin(e + f*x))*(sin(e + f*x) - 3)**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_628():
    f = (-3)**(-m - 1)*(a*sin(e + f*x) + a)**m
    F = -(-3)**(-m - 1)*2**(m + sympy.S.Half)*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(-m + sympy.S(-1)/2)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - sin(e + f*x)/2)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_629():
    f = (a*sin(e + f*x) + a)**m*(-sin(e + f*x) - 3)**(-m - 1)
    F = -sqrt(2)*sqrt(-(1 - sin(e + f*x))/(sin(e + f*x) + 1))*(a*sin(e + f*x) + a)**m*cos(e + f*x)*hyper((sympy.S.Half, -m), (1 - m,), (sin(e + f*x) + 3)/(2*sin(e + f*x) + 2))/(4*f*m*(1 - sin(e + f*x))*(-sin(e + f*x) - 3)**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_630():
    f = (a*sin(e + f*x) + a)**m*(-2*sin(e + f*x) - 3)**(-m - 1)
    F = -sqrt(5)*sqrt(-(1 - sin(e + f*x))/(sin(e + f*x) + 1))*(a*sin(e + f*x) + a)**m*cos(e + f*x)*hyper((sympy.S.Half, -m), (1 - m,), (4*sin(e + f*x) + 6)/(5*sin(e + f*x) + 5))/(5*f*m*(1 - sin(e + f*x))*(-2*sin(e + f*x) - 3)**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_631():
    f = (a*sin(e + f*x) + a)**m*(-3*sin(e + f*x) - 3)**(-m - 1)
    F = -(a*sin(e + f*x) + a)**m*(-3*sin(e + f*x) - 3)**(-m - 1)*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_632():
    f = (a*sin(e + f*x) + a)**m*(-4*sin(e + f*x) - 3)**(-m - 1)
    F = sqrt(7)*sqrt((1 - sin(e + f*x))/(sin(e + f*x) + 1))*(a*sin(e + f*x) + a)**m*cos(e + f*x)*hyper((sympy.S.Half, -m), (1 - m,), (8*sin(e + f*x) + 6)/(7*sin(e + f*x) + 7))/(7*f*m*(1 - sin(e + f*x))*(-4*sin(e + f*x) - 3)**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_633():
    f = (a*sin(e + f*x) + a)**m*(-5*sin(e + f*x) - 3)**(-m - 1)
    F = sqrt((1 - sin(e + f*x))/(sin(e + f*x) + 1))*(a*sin(e + f*x) + a)**m*cos(e + f*x)*hyper((sympy.S.Half, -m), (1 - m,), (5*sin(e + f*x) + 3)/(4*sin(e + f*x) + 4))/(4*f*m*(1 - sin(e + f*x))*(-5*sin(e + f*x) - 3)**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_634():
    f = (d*sin(e + f*x))**(-m - 1)*(a*sin(e + f*x) + a)**m
    F = -((sin(e + f*x) + 1)/(1 - sin(e + f*x)))**(sympy.S.Half - m)*(a*sin(e + f*x) + a)**m*cos(e + f*x)*hyper((-m, sympy.S.Half - m), (1 - m,), -2*sin(e + f*x)/(1 - sin(e + f*x)))/(d*f*m*(d*sin(e + f*x))**m*(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_635():
    f = (c + d*sin(e + f*x))**(-m - 1)*(a*sin(e + f*x) + a)**m
    F = -2**(m + sympy.S.Half)*a*((c + d)*(sin(e + f*x) + 1)/(c + d*sin(e + f*x)))**(sympy.S.Half - m)*(a*sin(e + f*x) + a)**(m - 1)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), (1 - sin(e + f*x))*(c - d)/(2*c + 2*d*sin(e + f*x)))/(f*(c + d)*(c + d*sin(e + f*x))**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_636():
    f = (c + d*sin(e + f*x))**n*(a*sin(e + f*x) + a)**3
    F = -8*sqrt(2)*a**3*(c + d*sin(e + f*x))**n*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S(-5)/2, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, d*(1 - sin(e + f*x))/(c + d))/(f*((c + d*sin(e + f*x))/(c + d))**n*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_637():
    f = (c + d*sin(e + f*x))**n*(a*sin(e + f*x) + a)**2
    F = -4*sqrt(2)*a**2*(c + d*sin(e + f*x))**n*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S(-3)/2, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, d*(1 - sin(e + f*x))/(c + d))/(f*((c + d*sin(e + f*x))/(c + d))**n*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_638():
    f = (c + d*sin(e + f*x))**n*(a*sin(e + f*x) + a)
    F = -2*sqrt(2)*a*(c + d*sin(e + f*x))**n*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S(-1)/2, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, d*(1 - sin(e + f*x))/(c + d))/(f*((c + d*sin(e + f*x))/(c + d))**n*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_639():
    f = (c + d*sin(e + f*x))**n
    F = -sqrt(2)*(c + d*sin(e + f*x))**n*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, d*(1 - sin(e + f*x))/(c + d))/(f*((c + d*sin(e + f*x))/(c + d))**n*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_640():
    f = (c + d*sin(e + f*x))**n/(a*sin(e + f*x) + a)
    F = -sqrt(2)*(c + d*sin(e + f*x))**n*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S(3)/2, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, d*(1 - sin(e + f*x))/(c + d))/(2*a*f*((c + d*sin(e + f*x))/(c + d))**n*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_641():
    f = (c + d*sin(e + f*x))**n/(a*sin(e + f*x) + a)**2
    F = -sqrt(2)*(c + d*sin(e + f*x))**n*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S(5)/2, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, d*(1 - sin(e + f*x))/(c + d))/(4*a**2*f*((c + d*sin(e + f*x))/(c + d))**n*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_642():
    f = (c + d*sin(e + f*x))**n/(a*sin(e + f*x) + a)**3
    F = -sqrt(2)*(c + d*sin(e + f*x))**n*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S(7)/2, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, d*(1 - sin(e + f*x))/(c + d))/(8*a**3*f*((c + d*sin(e + f*x))/(c + d))**n*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_643():
    f = (c + d*sin(e + f*x))**n*(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = 2*a**3*(c + d*sin(e + f*x))**(n + 1)*(3*c - d*(4*n + 11))*cos(e + f*x)/(d**2*f*(2*n + 3)*(2*n + 5)*sqrt(a*sin(e + f*x) + a)) - 2*a**3*(c + d*sin(e + f*x))**n*(3*c**2 - 2*c*d*(4*n + 7) + d**2*(16*n**2 + 56*n + 43))*cos(e + f*x)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), d*(1 - sin(e + f*x))/(c + d))/(d**2*f*((c + d*sin(e + f*x))/(c + d))**n*(2*n + 3)*(2*n + 5)*sqrt(a*sin(e + f*x) + a)) - 2*a**2*(c + d*sin(e + f*x))**(n + 1)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(d*f*(2*n + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_644():
    f = (c + d*sin(e + f*x))**n*(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -2*a**2*(c + d*sin(e + f*x))**(n + 1)*cos(e + f*x)/(d*f*(2*n + 3)*sqrt(a*sin(e + f*x) + a)) + 2*a**2*(c - d*(4*n + 5))*(c + d*sin(e + f*x))**n*cos(e + f*x)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), d*(1 - sin(e + f*x))/(c + d))/(d*f*((c + d*sin(e + f*x))/(c + d))**n*(2*n + 3)*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_645():
    f = (c + d*sin(e + f*x))**n*sqrt(a*sin(e + f*x) + a)
    F = -2*a*(c + d*sin(e + f*x))**n*cos(e + f*x)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), d*(1 - sin(e + f*x))/(c + d))/(f*((c + d*sin(e + f*x))/(c + d))**n*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_646():
    f = (c + d*sin(e + f*x))**(sympy.S(1)/3)*(a*sin(e + f*x) + a)
    F = -2*sqrt(2)*a*(c + d*sin(e + f*x))**(sympy.S(1)/3)*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S(-1)/2, sympy.S(-1)/3, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, d*(1 - sin(e + f*x))/(c + d))/(f*((c + d*sin(e + f*x))/(c + d))**(sympy.S(1)/3)*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_647():
    f = (a*sin(e + f*x) + a)/(c + d*sin(e + f*x))**(sympy.S(1)/3)
    F = -2*sqrt(2)*a*((c + d*sin(e + f*x))/(c + d))**(sympy.S(1)/3)*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S(-1)/2, sympy.S(1)/3, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, d*(1 - sin(e + f*x))/(c + d))/(f*(c + d*sin(e + f*x))**(sympy.S(1)/3)*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_648():
    f = (a*sin(e + f*x) + a)/(c + d*sin(e + f*x))**(sympy.S(4)/3)
    F = -2*sqrt(2)*a*((c + d*sin(e + f*x))/(c + d))**(sympy.S(1)/3)*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S(-1)/2, sympy.S(4)/3, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, d*(1 - sin(e + f*x))/(c + d))/(f*(c + d)*(c + d*sin(e + f*x))**(sympy.S(1)/3)*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_649():
    f = (a + b*sin(e + f*x))*(c + d*sin(e + f*x))**3
    F = -b*(c + d*sin(e + f*x))**3*cos(e + f*x)/(4*f) - d*(20*a*c*d + 6*b*c**2 + 9*b*d**2)*sin(e + f*x)*cos(e + f*x)/(24*f) + x*(a*c**3 + 3*a*c*d**2/2 + 3*b*c**2*d/2 + 3*b*d**3/8) - (c + d*sin(e + f*x))**2*(4*a*d + 3*b*c)*cos(e + f*x)/(12*f) - (4*a*d*(4*c**2 + d**2) + 3*b*(c**3 + 4*c*d**2))*cos(e + f*x)/(6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_650():
    f = (a + b*sin(e + f*x))*(c + d*sin(e + f*x))**2
    F = -b*(c + d*sin(e + f*x))**2*cos(e + f*x)/(3*f) - d*(3*a*d + 2*b*c)*sin(e + f*x)*cos(e + f*x)/(6*f) + x*(a*(2*c**2 + d**2)/2 + b*c*d) - (6*a*c*d + 2*b*(c**2 + d**2))*cos(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_651():
    f = (a + b*sin(e + f*x))*(c + d*sin(e + f*x))
    F = -b*d*sin(e + f*x)*cos(e + f*x)/(2*f) + x*(a*c + b*d/2) - (a*d + b*c)*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_652():
    f = a + b*sin(e + f*x)
    F = a*x - b*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_653():
    f = (a + b*sin(e + f*x))/(c + d*sin(e + f*x))
    F = b*x/d - (-2*a*d + 2*b*c)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(d*f*sqrt(c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_654():
    f = (a + b*sin(e + f*x))/(c + d*sin(e + f*x))**2
    F = (2*a*c - 2*b*d)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*(c**2 - d**2)**(sympy.S(3)/2)) - (-a*d + b*c)*cos(e + f*x)/(f*(c + d*sin(e + f*x))*(c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_655():
    f = (a + b*sin(e + f*x))/(c + d*sin(e + f*x))**3
    F = -(-a*(2*c**2 + d**2) + 3*b*c*d)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*(c**2 - d**2)**(sympy.S(5)/2)) + (3*a*c*d - b*(c**2 + 2*d**2))*cos(e + f*x)/(2*f*(c + d*sin(e + f*x))*(c**2 - d**2)**2) - (-a*d + b*c)*cos(e + f*x)/(f*(c + d*sin(e + f*x))**2*(2*c**2 - 2*d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_656():
    f = (a + b*sin(e + f*x))**2*(c + d*sin(e + f*x))**3
    F = -b**2*(c + d*sin(e + f*x))**4*cos(e + f*x)/(5*d*f) + b*(c + d*sin(e + f*x))**3*(-10*a*d + b*c)*cos(e + f*x)/(20*d*f) + x*(a**2*(2*c**3 + 3*c*d**2)/2 + 3*a*b*d*(4*c**2 + d**2)/4 + b**2*c*(4*c**2 + 9*d**2)/8) - (100*a**2*c*d**2 + 30*a*b*d*(2*c**2 + 3*d**2) - b**2*(6*c**3 - 71*c*d**2))*sin(e + f*x)*cos(e + f*x)/(120*f) - (c + d*sin(e + f*x))**2*(-3*b*c*(-10*a*d + b*c) + d**2*(20*a**2 + 16*b**2))*cos(e + f*x)/(60*d*f) - (20*a**2*d**2*(4*c**2 + d**2) + 30*a*b*c*d*(c**2 + 4*d**2) - b**2*(3*c**4 - 52*c**2*d**2 - 16*d**4))*cos(e + f*x)/(30*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_657():
    f = (a + b*sin(e + f*x))**2*(c + d*sin(e + f*x))**2
    F = x*(a**2*(2*c**2 + d**2)/2 + 2*a*b*c*d + b**2*(4*c**2 + 3*d**2)/8) - (2*a*d*(-a*d + 8*b*c) + 3*b**2*(4*c**2 + 3*d**2))*sin(e + f*x)*cos(e + f*x)/(24*f) - d**2*(a + b*sin(e + f*x))**3*cos(e + f*x)/(4*b*f) - d*(a + b*sin(e + f*x))**2*(-a*d + 8*b*c)*cos(e + f*x)/(12*b*f) - (-a**3*d**2 + 8*a**2*b*c*d + 4*a*b**2*(3*c**2 + 2*d**2) + 8*b**3*c*d)*cos(e + f*x)/(6*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_658():
    f = (a + b*sin(e + f*x))**2*(c + d*sin(e + f*x))
    F = -b*(2*a*d + 3*b*c)*sin(e + f*x)*cos(e + f*x)/(6*f) - d*(a + b*sin(e + f*x))**2*cos(e + f*x)/(3*f) + x*(a**2*c + a*b*d + b**2*c/2) - (2*a**2*d + 6*a*b*c + 2*b**2*d)*cos(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_659():
    f = (a + b*sin(e + f*x))**2
    F = -2*a*b*cos(e + f*x)/f - b**2*sin(e + f*x)*cos(e + f*x)/(2*f) + x*(a**2 + b**2/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_660():
    f = (a + b*sin(e + f*x))**2/(c + d*sin(e + f*x))
    F = -b**2*cos(e + f*x)/(d*f) - b*x*(-2*a*d + b*c)/d**2 + 2*(-a*d + b*c)**2*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(d**2*f*sqrt(c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_661():
    f = (a + b*sin(e + f*x))**2/(c + d*sin(e + f*x))**2
    F = b**2*x/d**2 + (-a*d + b*c)**2*cos(e + f*x)/(d*f*(c + d*sin(e + f*x))*(c**2 - d**2)) - (-2*a*d + 2*b*c)*(a*c*d + b*(c**2 - 2*d**2))*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(d**2*f*(c**2 - d**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_662():
    f = (a + b*sin(e + f*x))**2/(c + d*sin(e + f*x))**3
    F = -(-a**2*(2*c**2 + d**2) + 6*a*b*c*d - b**2*(c**2 + 2*d**2))*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*(c**2 - d**2)**(sympy.S(5)/2)) - (-a*d + b*c)*(3*a*c*d + b*(c**2 - 4*d**2))*cos(e + f*x)/(2*d*f*(c + d*sin(e + f*x))*(c**2 - d**2)**2) + (-a*d + b*c)**2*cos(e + f*x)/(2*d*f*(c + d*sin(e + f*x))**2*(c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_663():
    f = (a + b*sin(e + f*x))**2/(c + d*sin(e + f*x))**4
    F = -(-a**2*(2*c**3 + 3*c*d**2) + 2*a*b*d*(4*c**2 + d**2) - b**2*c*(c**2 + 4*d**2))*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*(c**2 - d**2)**(sympy.S(7)/2)) + (a**2*d**2*(11*c**2 + 4*d**2) - a*b*(4*c**3*d + 26*c*d**3) - b**2*(c**4 - 10*c**2*d**2 - 6*d**4))*cos(e + f*x)/(6*d*f*(c + d*sin(e + f*x))*(c**2 - d**2)**3) - (-a*d + b*c)*(5*a*c*d + b*(c**2 - 6*d**2))*cos(e + f*x)/(6*d*f*(c + d*sin(e + f*x))**2*(c**2 - d**2)**2) + (-a*d + b*c)**2*cos(e + f*x)/(3*d*f*(c + d*sin(e + f*x))**3*(c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_664():
    f = (a + b*sin(e + f*x))**3*(c + d*sin(e + f*x))**2
    F = x*(a**3*(2*c**2 + d**2)/2 + 3*a**2*b*c*d + 3*a*b**2*(4*c**2 + 3*d**2)/8 + 3*b**3*c*d/4) - (-6*a**3*d**2 + 60*a**2*b*c*d + a*b**2*(100*c**2 + 71*d**2) + 90*b**3*c*d)*sin(e + f*x)*cos(e + f*x)/(120*f) - d**2*(a + b*sin(e + f*x))**4*cos(e + f*x)/(5*b*f) - d*(a + b*sin(e + f*x))**3*(-a*d + 10*b*c)*cos(e + f*x)/(20*b*f) - (a + b*sin(e + f*x))**2*(3*a*d*(-a*d + 10*b*c) + 4*b**2*(5*c**2 + 4*d**2))*cos(e + f*x)/(60*b*f) - (-3*a**4*d**2 + 30*a**3*b*c*d + 4*a**2*b**2*(20*c**2 + 13*d**2) + 120*a*b**3*c*d + 4*b**4*(5*c**2 + 4*d**2))*cos(e + f*x)/(30*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_665():
    f = (a + b*sin(e + f*x))**3*(c + d*sin(e + f*x))
    F = -b*(6*a**2*d + 20*a*b*c + 9*b**2*d)*sin(e + f*x)*cos(e + f*x)/(24*f) - d*(a + b*sin(e + f*x))**3*cos(e + f*x)/(4*f) + x*(a**3*c + 3*a**2*b*d/2 + 3*a*b**2*c/2 + 3*b**3*d/8) - (a + b*sin(e + f*x))**2*(3*a*d + 4*b*c)*cos(e + f*x)/(12*f) - (3*a**3*d + 16*a**2*b*c + 12*a*b**2*d + 4*b**3*c)*cos(e + f*x)/(6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_666():
    f = (a + b*sin(e + f*x))**3
    F = -5*a*b**2*sin(e + f*x)*cos(e + f*x)/(6*f) + a*x*(2*a**2 + 3*b**2)/2 - b*(a + b*sin(e + f*x))**2*cos(e + f*x)/(3*f) - 2*b*(4*a**2 + b**2)*cos(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_667():
    f = (a + b*sin(e + f*x))**3/(c + d*sin(e + f*x))
    F = -b**2*(a + b*sin(e + f*x))*cos(e + f*x)/(2*d*f) + b**2*(-5*a*d + 2*b*c)*cos(e + f*x)/(2*d**2*f) - b*x*(-6*a**2*d**2 + 6*a*b*c*d - b**2*(2*c**2 + d**2))/(2*d**3) - 2*(-a*d + b*c)**3*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(d**3*f*sqrt(c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_668():
    f = (a + b*sin(e + f*x))**3/(c + d*sin(e + f*x))**2
    F = -b**2*x*(-3*a*d + 2*b*c)/d**3 + b*(-a**2*d**2 + 2*a*b*c*d - b**2*(2*c**2 - d**2))*cos(e + f*x)/(d**2*f*(c**2 - d**2)) + (a + b*sin(e + f*x))*(-a*d + b*c)**2*cos(e + f*x)/(d*f*(c + d*sin(e + f*x))*(c**2 - d**2)) + 2*(-a*d + b*c)**2*(a*c*d + 2*b*c**2 - 3*b*d**2)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(d**3*f*(c**2 - d**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_669():
    f = (a + b*sin(e + f*x))**3/(c + d*sin(e + f*x))**3
    F = b**3*x/d**3 + (a + b*sin(e + f*x))*(-a*d + b*c)**2*cos(e + f*x)/(2*d*f*(c + d*sin(e + f*x))**2*(c**2 - d**2)) + (-a*d + b*c)**2*(3*a*c*d + 2*b*c**2 - 5*b*d**2)*cos(e + f*x)/(2*d**2*f*(c + d*sin(e + f*x))*(c**2 - d**2)**2) - (-a**3*d**3*(2*c**2 + d**2) + 9*a**2*b*c*d**4 - 3*a*b**2*d**3*(c**2 + 2*d**2) + b**3*(2*c**5 - 5*c**3*d**2 + 6*c*d**4))*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(d**3*f*(c**2 - d**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_670():
    f = (a + b*sin(e + f*x))**3/(c + d*sin(e + f*x))**4
    F = -(a*c - b*d)*(-a**2*(2*c**2 + 3*d**2) + 10*a*b*c*d - b**2*(3*c**2 + 2*d**2))*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*(c**2 - d**2)**(sympy.S(7)/2)) + (a + b*sin(e + f*x))*(-a*d + b*c)**2*cos(e + f*x)/(3*d*f*(c + d*sin(e + f*x))**3*(c**2 - d**2)) - (-a*d + b*c)*(a**2*d**2*(11*c**2 + 4*d**2) + 5*a*b*c*d*(c**2 - 7*d**2) + b**2*(2*c**4 - 5*c**2*d**2 + 18*d**4))*cos(e + f*x)/(6*d**2*f*(c + d*sin(e + f*x))*(c**2 - d**2)**3) + (-a*d + b*c)**2*(5*a*c*d + 2*b*c**2 - 7*b*d**2)*cos(e + f*x)/(6*d**2*f*(c + d*sin(e + f*x))**2*(c**2 - d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_671():
    f = (B*sin(x) + B*b/a)/(a + b*sin(x))
    F = B*x/b - 2*B*sqrt(a**2 - b**2)*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(a*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_672():
    f = (B*a/b + B*sin(x))/(a + b*sin(x))
    F = B*x/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_673():
    f = (a + b*sin(x))/(a*sin(x) + b)**2
    F = -cos(x)/(a*sin(x) + b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_674():
    f = (2 - sin(x))/(sin(x) + 2)
    F = -x + 4*sqrt(3)*x/3 + 8*sqrt(3)*atan(cos(x)/(sin(x) + sqrt(3) + 2))/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_675():
    f = (c + d*sin(e + f*x))**4/(a + b*sin(e + f*x))
    F = -d**2*(c + d*sin(e + f*x))**2*cos(e + f*x)/(3*b*f) - d**3*(-3*a*d + 8*b*c)*sin(e + f*x)*cos(e + f*x)/(6*b**2*f) + d**2*(-3*a**2*d**2 + 12*a*b*c*d - b**2*(17*c**2 + 2*d**2))*cos(e + f*x)/(3*b**3*f) + d*x*(-2*a**3*d**3 + 8*a**2*b*c*d**2 - a*b**2*d*(12*c**2 + d**2) + 4*b**3*c*(2*c**2 + d**2))/(2*b**4) + 2*(-a*d + b*c)**4*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(b**4*f*sqrt(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_676():
    f = (c + d*sin(e + f*x))**3/(a + b*sin(e + f*x))
    F = -d**2*(c + d*sin(e + f*x))*cos(e + f*x)/(2*b*f) - d**2*(-2*a*d + 5*b*c)*cos(e + f*x)/(2*b**2*f) - d*x*(-2*a**2*d**2 + 6*a*b*c*d - b**2*(6*c**2 + d**2))/(2*b**3) + 2*(-a*d + b*c)**3*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(b**3*f*sqrt(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_677():
    f = (c + d*sin(e + f*x))**2/(a + b*sin(e + f*x))
    F = -d**2*cos(e + f*x)/(b*f) + d*x*(-a*d + 2*b*c)/b**2 + 2*(-a*d + b*c)**2*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(b**2*f*sqrt(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_678():
    f = (c + d*sin(e + f*x))/(a + b*sin(e + f*x))
    F = d*x/b + (-2*a*d + 2*b*c)*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(b*f*sqrt(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_679():
    f = 1/(a + b*sin(e + f*x))
    F = 2*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(f*sqrt(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_680():
    f = 1/((a + b*sin(e + f*x))*(c + d*sin(e + f*x)))
    F = 2*b*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(f*sqrt(a**2 - b**2)*(-a*d + b*c)) - 2*d*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*sqrt(c**2 - d**2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_681():
    f = 1/((a + b*sin(e + f*x))*(c + d*sin(e + f*x))**2)
    F = 2*b**2*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(f*sqrt(a**2 - b**2)*(-a*d + b*c)**2) - d**2*cos(e + f*x)/(f*(c + d*sin(e + f*x))*(c**2 - d**2)*(-a*d + b*c)) + 2*d*(a*c*d - b*(2*c**2 - d**2))*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*(c**2 - d**2)**(sympy.S(3)/2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_682():
    f = 1/((a + b*sin(e + f*x))*(c + d*sin(e + f*x))**3)
    F = 2*b**3*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(f*sqrt(a**2 - b**2)*(-a*d + b*c)**3) - d**2*(-3*a*c*d + 5*b*c**2 - 2*b*d**2)*cos(e + f*x)/(2*f*(c + d*sin(e + f*x))*(c**2 - d**2)**2*(-a*d + b*c)**2) - d**2*cos(e + f*x)/(f*(c + d*sin(e + f*x))**2*(c**2 - d**2)*(-2*a*d + 2*b*c)) + d*(-a**2*d**2*(2*c**2 + d**2) + 6*a*b*c**3*d - b**2*(6*c**4 - 5*c**2*d**2 + 2*d**4))*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*(c**2 - d**2)**(sympy.S(5)/2)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_683():
    f = (c + d*sin(e + f*x))**4/(a + b*sin(e + f*x))**2
    F = (c + d*sin(e + f*x))**2*(-a*d + b*c)**2*cos(e + f*x)/(b*f*(a + b*sin(e + f*x))*(a**2 - b**2)) + d**2*(-3*a**2*d**2 + 4*a*b*c*d - b**2*(2*c**2 - d**2))*sin(e + f*x)*cos(e + f*x)/(2*b**2*f*(a**2 - b**2)) + d*(-a*d + 2*b*c)*(-3*a**2*d**2 + 2*a*b*c*d - b**2*(c**2 - 2*d**2))*cos(e + f*x)/(b**3*f*(a**2 - b**2)) - d**2*x*(-6*a**2*d**2 + 16*a*b*c*d - b**2*(12*c**2 + d**2))/(2*b**4) + 2*(-a*d + b*c)**3*(3*a**2*d + a*b*c - 4*b**2*d)*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(b**4*f*(a**2 - b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_684():
    f = (c + d*sin(e + f*x))**3/(a + b*sin(e + f*x))**2
    F = (c + d*sin(e + f*x))*(-a*d + b*c)**2*cos(e + f*x)/(b*f*(a + b*sin(e + f*x))*(a**2 - b**2)) + d*(-2*a**2*d**2 + 2*a*b*c*d - b**2*(c**2 - d**2))*cos(e + f*x)/(b**2*f*(a**2 - b**2)) + d**2*x*(-2*a*d + 3*b*c)/b**3 + 2*(-a*d + b*c)**2*(2*a**2*d + a*b*c - 3*b**2*d)*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(b**3*f*(a**2 - b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_685():
    f = (c + d*sin(e + f*x))**2/(a + b*sin(e + f*x))**2
    F = (-a*d + b*c)**2*cos(e + f*x)/(b*f*(a + b*sin(e + f*x))*(a**2 - b**2)) + d**2*x/b**2 + (-2*a*d + 2*b*c)*(a**2*d + a*b*c - 2*b**2*d)*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(b**2*f*(a**2 - b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_686():
    f = (c + d*sin(e + f*x))/(a + b*sin(e + f*x))**2
    F = (2*a*c - 2*b*d)*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(f*(a**2 - b**2)**(sympy.S(3)/2)) + (-a*d + b*c)*cos(e + f*x)/(f*(a + b*sin(e + f*x))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_687():
    f = (a + b*sin(e + f*x))**(-2)
    F = 2*a*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(f*(a**2 - b**2)**(sympy.S(3)/2)) + b*cos(e + f*x)/(f*(a + b*sin(e + f*x))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_688():
    f = 1/((a + b*sin(e + f*x))**2*(c + d*sin(e + f*x)))
    F = b**2*cos(e + f*x)/(f*(a + b*sin(e + f*x))*(a**2 - b**2)*(-a*d + b*c)) + 2*b*(-2*a**2*d + a*b*c + b**2*d)*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(f*(a**2 - b**2)**(sympy.S(3)/2)*(-a*d + b*c)**2) + 2*d**2*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*sqrt(c**2 - d**2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_689():
    f = 1/((a + b*sin(e + f*x))**2*(c + d*sin(e + f*x))**2)
    F = 2*b**2*(-3*a**2*d + a*b*c + 2*b**2*d)*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(f*(a**2 - b**2)**(sympy.S(3)/2)*(-a*d + b*c)**3) + b**2*cos(e + f*x)/(f*(a + b*sin(e + f*x))*(a**2 - b**2)*(c + d*sin(e + f*x))*(-a*d + b*c)) + 2*d**2*(-a*c*d + 3*b*c**2 - 2*b*d**2)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*(c**2 - d**2)**(sympy.S(3)/2)*(-a*d + b*c)**3) + d*(a**2*d**2 + b**2*(c**2 - 2*d**2))*cos(e + f*x)/(f*(a**2 - b**2)*(c + d*sin(e + f*x))*(c**2 - d**2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_690():
    f = 1/((a + b*sin(e + f*x))**2*(c + d*sin(e + f*x))**3)
    F = 2*b**3*(-4*a**2*d + a*b*c + 3*b**2*d)*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(f*(a**2 - b**2)**(sympy.S(3)/2)*(-a*d + b*c)**4) + b**2*cos(e + f*x)/(f*(a + b*sin(e + f*x))*(a**2 - b**2)*(c + d*sin(e + f*x))**2*(-a*d + b*c)) - d**2*(-a**2*d**2*(2*c**2 + d**2) + 2*a*b*c*d*(4*c**2 - d**2) - 3*b**2*(4*c**4 - 5*c**2*d**2 + 2*d**4))*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*(c**2 - d**2)**(sympy.S(5)/2)*(-a*d + b*c)**4) + d*(a**2*d**2 + b**2*(2*c**2 - 3*d**2))*cos(e + f*x)/(f*(2*a**2 - 2*b**2)*(c + d*sin(e + f*x))**2*(c**2 - d**2)*(-a*d + b*c)**2) - (3*a**3*c*d**4 - a**2*b*d**3*(7*c**2 - 4*d**2) - 3*a*b**2*c*d**4 - b**3*(2*c**4*d - 11*c**2*d**3 + 6*d**5))*cos(e + f*x)/(f*(2*a**2 - 2*b**2)*(c + d*sin(e + f*x))*(c**2 - d**2)**2*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_691():
    f = (c + d*sin(e + f*x))**5/(a + b*sin(e + f*x))**3
    F = (c + d*sin(e + f*x))**3*(-a*d + b*c)**2*cos(e + f*x)/(2*b*f*(a + b*sin(e + f*x))**2*(a**2 - b**2)) + (c + d*sin(e + f*x))**2*(-a*d + b*c)**2*(4*a**2*d + 3*a*b*c - 7*b**2*d)*cos(e + f*x)/(2*b**2*f*(a + b*sin(e + f*x))*(a**2 - b**2)**2) + d**2*(-6*a**4*d**3 + 7*a**3*b*c*d**2 + a**2*b**2*d*(c**2 + 10*d**2) - a*b**3*c*(3*c**2 + 16*d**2) + b**4*d*(8*c**2 - d**2))*sin(e + f*x)*cos(e + f*x)/(2*b**3*f*(a**2 - b**2)**2) - d*(-12*a**5*d**4 + 30*a**4*b*c*d**3 - a**3*b**2*d**2*(16*c**2 - 21*d**2) - a**2*b**3*c*d*(4*c**2 + 55*d**2) + a*b**4*(6*c**4 + 43*c**2*d**2 - 6*d**4) - b**5*c*d*(17*c**2 - 10*d**2))*cos(e + f*x)/(2*b**4*f*(a**2 - b**2)**2) - d**3*x*(-12*a**2*d**2 + 30*a*b*c*d - b**2*(20*c**2 + d**2))/(2*b**5) + (-a*d + b*c)**3*(12*a**4*d**2 + 6*a**3*b*c*d + a**2*b**2*(2*c**2 - 29*d**2) - 12*a*b**3*c*d + b**4*(c**2 + 20*d**2))*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(b**5*f*(a**2 - b**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_692():
    f = (c + d*sin(e + f*x))**4/(a + b*sin(e + f*x))**3
    F = (c + d*sin(e + f*x))**2*(-a*d + b*c)**2*cos(e + f*x)/(2*b*f*(a + b*sin(e + f*x))**2*(a**2 - b**2)) + d**2*(-3*a**2*d**2 + 2*a*b*c*d - b**2*(c**2 - 2*d**2))*cos(e + f*x)/(2*b**3*f*(a**2 - b**2)) + 3*(-a*d + b*c)**3*(a**2*d + a*b*c - 2*b**2*d)*cos(e + f*x)/(2*b**3*f*(a + b*sin(e + f*x))*(a**2 - b**2)**2) + d**3*x*(-3*a*d + 4*b*c)/b**4 + (-a*d + b*c)**2*(6*a**4*d**2 + 4*a**3*b*c*d + a**2*b**2*(2*c**2 - 15*d**2) - 10*a*b**3*c*d + b**4*(c**2 + 12*d**2))*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(b**4*f*(a**2 - b**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_693():
    f = (c + d*sin(e + f*x))**3/(a + b*sin(e + f*x))**3
    F = (c + d*sin(e + f*x))*(-a*d + b*c)**2*cos(e + f*x)/(2*b*f*(a + b*sin(e + f*x))**2*(a**2 - b**2)) + (-a*d + b*c)**2*(2*a**2*d + 3*a*b*c - 5*b**2*d)*cos(e + f*x)/(2*b**2*f*(a + b*sin(e + f*x))*(a**2 - b**2)**2) + d**3*x/b**3 + (-a*d + b*c)*(2*a**4*d**2 + 2*a**3*b*c*d + a**2*b**2*(2*c**2 - 5*d**2) - 8*a*b**3*c*d + b**4*(c**2 + 6*d**2))*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(b**3*f*(a**2 - b**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_694():
    f = (c + d*sin(e + f*x))**2/(a + b*sin(e + f*x))**3
    F = -(-a**2*(2*c**2 + d**2) + 6*a*b*c*d - b**2*(c**2 + 2*d**2))*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(f*(a**2 - b**2)**(sympy.S(5)/2)) + (-a*d + b*c)*(a**2*d + 3*a*b*c - 4*b**2*d)*cos(e + f*x)/(2*b*f*(a + b*sin(e + f*x))*(a**2 - b**2)**2) + (-a*d + b*c)**2*cos(e + f*x)/(2*b*f*(a + b*sin(e + f*x))**2*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_695():
    f = (c + d*sin(e + f*x))/(a + b*sin(e + f*x))**3
    F = (2*a**2*c - 3*a*b*d + b**2*c)*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(f*(a**2 - b**2)**(sympy.S(5)/2)) + (-a**2*d + 3*a*b*c - 2*b**2*d)*cos(e + f*x)/(2*f*(a + b*sin(e + f*x))*(a**2 - b**2)**2) + (-a*d + b*c)*cos(e + f*x)/(f*(a + b*sin(e + f*x))**2*(2*a**2 - 2*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_696():
    f = (a + b*sin(e + f*x))**(-3)
    F = 3*a*b*cos(e + f*x)/(2*f*(a + b*sin(e + f*x))*(a**2 - b**2)**2) + b*cos(e + f*x)/(f*(a + b*sin(e + f*x))**2*(2*a**2 - 2*b**2)) + (2*a**2 + b**2)*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(f*(a**2 - b**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_697():
    f = 1/((a + b*sin(e + f*x))**3*(c + d*sin(e + f*x)))
    F = b**2*(-5*a**2*d + 3*a*b*c + 2*b**2*d)*cos(e + f*x)/(2*f*(a + b*sin(e + f*x))*(a**2 - b**2)**2*(-a*d + b*c)**2) + b**2*cos(e + f*x)/(f*(a + b*sin(e + f*x))**2*(2*a**2 - 2*b**2)*(-a*d + b*c)) - b*(-6*a**4*d**2 + 6*a**3*b*c*d - a**2*b**2*(2*c**2 - 5*d**2) - b**4*(c**2 + 2*d**2))*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(f*(a**2 - b**2)**(sympy.S(5)/2)*(-a*d + b*c)**3) - 2*d**3*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*sqrt(c**2 - d**2)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_698():
    f = 1/((a + b*sin(e + f*x))**3*(c + d*sin(e + f*x))**2)
    F = -b**2*(-12*a**4*d**2 + 8*a**3*b*c*d - a**2*b**2*(2*c**2 - 15*d**2) - 2*a*b**3*c*d - b**4*(c**2 + 6*d**2))*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(f*(a**2 - b**2)**(sympy.S(5)/2)*(-a*d + b*c)**4) + 3*b**2*(-2*a**2*d + a*b*c + b**2*d)*cos(e + f*x)/(2*f*(a + b*sin(e + f*x))*(a**2 - b**2)**2*(c + d*sin(e + f*x))*(-a*d + b*c)**2) + b**2*cos(e + f*x)/(f*(a + b*sin(e + f*x))**2*(2*a**2 - 2*b**2)*(c + d*sin(e + f*x))*(-a*d + b*c)) - 2*d**3*(-a*c*d + 4*b*c**2 - 3*b*d**2)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*(c**2 - d**2)**(sympy.S(3)/2)*(-a*d + b*c)**4) - d*(2*a**4*d**3 + a**2*b**2*d*(7*c**2 - 11*d**2) - 3*a*b**3*c*(c**2 - d**2) - 2*b**4*d*(2*c**2 - 3*d**2))*cos(e + f*x)/(2*f*(a**2 - b**2)**2*(c + d*sin(e + f*x))*(c**2 - d**2)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_699():
    f = 1/((a + b*sin(e + f*x))**3*(c + d*sin(e + f*x))**3)
    F = -b**3*(-20*a**4*d**2 + 10*a**3*b*c*d - a**2*b**2*(2*c**2 - 29*d**2) - 4*a*b**3*c*d - b**4*(c**2 + 12*d**2))*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(f*(a**2 - b**2)**(sympy.S(5)/2)*(-a*d + b*c)**5) + b**2*(-7*a**2*d + 3*a*b*c + 4*b**2*d)*cos(e + f*x)/(2*f*(a + b*sin(e + f*x))*(a**2 - b**2)**2*(c + d*sin(e + f*x))**2*(-a*d + b*c)**2) + b**2*cos(e + f*x)/(f*(a + b*sin(e + f*x))**2*(2*a**2 - 2*b**2)*(c + d*sin(e + f*x))**2*(-a*d + b*c)) - d**3*(a**2*d**2*(2*c**2 + d**2) - a*b*(10*c**3*d - 4*c*d**3) + b**2*(20*c**4 - 29*c**2*d**2 + 12*d**4))*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*(c**2 - d**2)**(sympy.S(5)/2)*(-a*d + b*c)**5) + 3*d*(a**5*c*d**4 - a**4*b*(3*c**2*d**3 - 2*d**5) - 2*a**3*b**2*c*d**4 - a**2*b**3*d*(3*c**4 - 12*c**2*d**2 + 7*d**4) + a*b**4*c*(c**4 - 2*c**2*d**2 + 2*d**4) + b**5*d*(2*c**4 - 7*c**2*d**2 + 4*d**4))*cos(e + f*x)/(2*f*(a**2 - b**2)**2*(c + d*sin(e + f*x))*(c**2 - d**2)**2*(-a*d + b*c)**4) - d*(a**4*d**3 + 2*a**2*b**2*d*(4*c**2 - 5*d**2) - 3*a*b**3*c*(c**2 - d**2) - b**4*d*(5*c**2 - 6*d**2))*cos(e + f*x)/(2*f*(a**2 - b**2)**2*(c + d*sin(e + f*x))**2*(c**2 - d**2)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_700():
    f = (a + b*sin(e + f*x))*(c + d*sin(e + f*x))**(sympy.S(5)/2)
    F = -2*b*(c + d*sin(e + f*x))**(sympy.S(5)/2)*cos(e + f*x)/(7*f) - (c + d*sin(e + f*x))**(sympy.S(3)/2)*(14*a*d + 10*b*c)*cos(e + f*x)/(35*f) + sqrt(c + d*sin(e + f*x))*(-112*a*c*d - 30*b*c**2 - 50*b*d**2)*cos(e + f*x)/(105*f) - sqrt((c + d*sin(e + f*x))/(c + d))*(2*c**2 - 2*d**2)*(56*a*c*d + 15*b*c**2 + 25*b*d**2)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(105*d*f*sqrt(c + d*sin(e + f*x))) + sqrt(c + d*sin(e + f*x))*(322*a*c**2*d + 126*a*d**3 + 30*b*c**3 + 290*b*c*d**2)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(105*d*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_701():
    f = (a + b*sin(e + f*x))*(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = -2*b*(c + d*sin(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(5*f) + sqrt(c + d*sin(e + f*x))*(-10*a*d - 6*b*c)*cos(e + f*x)/(15*f) - sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2)*(10*a*d + 6*b*c)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d*f*sqrt(c + d*sin(e + f*x))) + sqrt(c + d*sin(e + f*x))*(40*a*c*d + 6*b*(c**2 + 3*d**2))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_702():
    f = (a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))
    F = -2*b*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(3*f) - 2*b*sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d*f*sqrt(c + d*sin(e + f*x))) + sqrt(c + d*sin(e + f*x))*(6*a*d + 2*b*c)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_703():
    f = (a + b*sin(e + f*x))/sqrt(c + d*sin(e + f*x))
    F = 2*b*sqrt(c + d*sin(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(d*f*sqrt((c + d*sin(e + f*x))/(c + d))) - sqrt((c + d*sin(e + f*x))/(c + d))*(-2*a*d + 2*b*c)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(d*f*sqrt(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_704():
    f = (a + b*sin(e + f*x))/(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = 2*b*sqrt((c + d*sin(e + f*x))/(c + d))*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(d*f*sqrt(c + d*sin(e + f*x))) + (2*a*d - 2*b*c)*cos(e + f*x)/(f*sqrt(c + d*sin(e + f*x))*(c**2 - d**2)) - sqrt(c + d*sin(e + f*x))*(-2*a*d + 2*b*c)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(d*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_705():
    f = (a + b*sin(e + f*x))/(c + d*sin(e + f*x))**(sympy.S(5)/2)
    F = (8*a*c*d - 2*b*(c**2 + 3*d**2))*cos(e + f*x)/(3*f*sqrt(c + d*sin(e + f*x))*(c**2 - d**2)**2) + (2*a*d - 2*b*c)*cos(e + f*x)/(f*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(3*c**2 - 3*d**2)) + sqrt((c + d*sin(e + f*x))/(c + d))*(-2*a*d + 2*b*c)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d*f*sqrt(c + d*sin(e + f*x))*(c**2 - d**2)) + sqrt(c + d*sin(e + f*x))*(8*a*c*d - 2*b*(c**2 + 3*d**2))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_706():
    f = (a + b*sin(e + f*x))/(c + d*sin(e + f*x))**(sympy.S(7)/2)
    F = -(-46*a*c**2*d - 18*a*d**3 + 6*b*c**3 + 58*b*c*d**2)*cos(e + f*x)/(15*f*sqrt(c + d*sin(e + f*x))*(c**2 - d**2)**3) - (-16*a*c*d + 6*b*c**2 + 10*b*d**2)*cos(e + f*x)/(15*f*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(c**2 - d**2)**2) + (2*a*d - 2*b*c)*cos(e + f*x)/(f*(c + d*sin(e + f*x))**(sympy.S(5)/2)*(5*c**2 - 5*d**2)) + sqrt((c + d*sin(e + f*x))/(c + d))*(-16*a*c*d + 6*b*c**2 + 10*b*d**2)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d*f*sqrt(c + d*sin(e + f*x))*(c**2 - d**2)**2) - sqrt(c + d*sin(e + f*x))*(-46*a*c**2*d - 18*a*d**3 + 6*b*c**3 + 58*b*c*d**2)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_707():
    f = (a + b*sin(e + f*x))**2*(c + d*sin(e + f*x))**(sympy.S(5)/2)
    F = -2*b**2*(c + d*sin(e + f*x))**(sympy.S(7)/2)*cos(e + f*x)/(9*d*f) + 4*b*(c + d*sin(e + f*x))**(sympy.S(5)/2)*(-9*a*d + b*c)*cos(e + f*x)/(63*d*f) - (c + d*sin(e + f*x))**(sympy.S(3)/2)*(-20*b*c*(-9*a*d + b*c) + 2*d**2*(63*a**2 + 49*b**2))*cos(e + f*x)/(315*d*f) + sqrt(c + d*sin(e + f*x))*(-336*a**2*c*d**2 - 60*a*b*d*(3*c**2 + 5*d**2) + 4*b**2*(5*c**3 - 57*c*d**2))*cos(e + f*x)/(315*d*f) + sqrt((c + d*sin(e + f*x))/(c + d))*(4*c**2 - 4*d**2)*(-84*a**2*c*d**2 - 45*a*b*c**2*d - 75*a*b*d**3 + 5*b**2*c**3 - 57*b**2*c*d**2)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(315*d**2*f*sqrt(c + d*sin(e + f*x))) + sqrt(c + d*sin(e + f*x))*(42*a**2*d**2*(23*c**2 + 9*d**2) + 60*a*b*d*(3*c**3 + 29*c*d**2) - 2*b**2*(10*c**4 - 279*c**2*d**2 - 147*d**4))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(315*d**2*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_708():
    f = (a + b*sin(e + f*x))**2*(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = -2*b**2*(c + d*sin(e + f*x))**(sympy.S(5)/2)*cos(e + f*x)/(7*d*f) + 4*b*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(-7*a*d + b*c)*cos(e + f*x)/(35*d*f) + sqrt(c + d*sin(e + f*x))*(12*b*c*(-7*a*d + b*c) - 2*d**2*(35*a**2 + 25*b**2))*cos(e + f*x)/(105*d*f) - sqrt((c + d*sin(e + f*x))/(c + d))*(2*c**2 - 2*d**2)*(35*a**2*d**2 + 42*a*b*c*d - b**2*(6*c**2 - 25*d**2))*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(105*d**2*f*sqrt(c + d*sin(e + f*x))) + sqrt(c + d*sin(e + f*x))*(280*a**2*c*d**2 + 84*a*b*d*(c**2 + 3*d**2) - 4*b**2*(3*c**3 - 41*c*d**2))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(105*d**2*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_709():
    f = (a + b*sin(e + f*x))**2*sqrt(c + d*sin(e + f*x))
    F = -2*b**2*(c + d*sin(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(5*d*f) + 4*b*sqrt(c + d*sin(e + f*x))*(-5*a*d + b*c)*cos(e + f*x)/(15*d*f) + 4*b*sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2)*(-5*a*d + b*c)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d**2*f*sqrt(c + d*sin(e + f*x))) + sqrt(c + d*sin(e + f*x))*(-4*b*c*(-5*a*d + b*c) + 2*d**2*(15*a**2 + 9*b**2))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d**2*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_710():
    f = (a + b*sin(e + f*x))**2/sqrt(c + d*sin(e + f*x))
    F = -2*b**2*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(3*d*f) - 4*b*sqrt(c + d*sin(e + f*x))*(-3*a*d + b*c)*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d**2*f*sqrt((c + d*sin(e + f*x))/(c + d))) + sqrt((c + d*sin(e + f*x))/(c + d))*(4*b*c*(-3*a*d + b*c) + 2*d**2*(3*a**2 + b**2))*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d**2*f*sqrt(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_711():
    f = (a + b*sin(e + f*x))**2/(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = -4*b*sqrt((c + d*sin(e + f*x))/(c + d))*(-a*d + b*c)*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(d**2*f*sqrt(c + d*sin(e + f*x))) + 2*(-a*d + b*c)**2*cos(e + f*x)/(d*f*sqrt(c + d*sin(e + f*x))*(c**2 - d**2)) + sqrt(c + d*sin(e + f*x))*(-4*a*b*c*d + 4*b**2*c**2 + 2*d**2*(a**2 - b**2))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(d**2*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_712():
    f = (a + b*sin(e + f*x))**2/(c + d*sin(e + f*x))**(sympy.S(5)/2)
    F = -(-4*a*d + 4*b*c)*(2*a*c*d + b*(c**2 - 3*d**2))*cos(e + f*x)/(3*d*f*sqrt(c + d*sin(e + f*x))*(c**2 - d**2)**2) + 2*(-a*d + b*c)**2*cos(e + f*x)/(3*d*f*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(c**2 - d**2)) + sqrt((c + d*sin(e + f*x))/(c + d))*(-2*a**2*d**2 + 4*a*b*c*d + 2*b**2*(2*c**2 - 3*d**2))*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d**2*f*sqrt(c + d*sin(e + f*x))*(c**2 - d**2)) - sqrt(c + d*sin(e + f*x))*(-4*a*d + 4*b*c)*(2*a*c*d + b*(c**2 - 3*d**2))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d**2*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_713():
    f = (a + b*sin(e + f*x))**2/(c + d*sin(e + f*x))**(sympy.S(7)/2)
    F = (2*a**2*d**2*(23*c**2 + 9*d**2) - 2*a*b*(6*c**3*d + 58*c*d**3) - 2*b**2*(2*c**4 - 19*c**2*d**2 - 15*d**4))*cos(e + f*x)/(15*d*f*sqrt(c + d*sin(e + f*x))*(c**2 - d**2)**3) - (-4*a*d + 4*b*c)*(4*a*c*d + b*(c**2 - 5*d**2))*cos(e + f*x)/(15*d*f*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(c**2 - d**2)**2) + 2*(-a*d + b*c)**2*cos(e + f*x)/(5*d*f*(c + d*sin(e + f*x))**(sympy.S(5)/2)*(c**2 - d**2)) + sqrt((c + d*sin(e + f*x))/(c + d))*(-4*a*d + 4*b*c)*(4*a*c*d + b*(c**2 - 5*d**2))*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d**2*f*sqrt(c + d*sin(e + f*x))*(c**2 - d**2)**2) + sqrt(c + d*sin(e + f*x))*(2*a**2*d**2*(23*c**2 + 9*d**2) - 2*a*b*(6*c**3*d + 58*c*d**3) - 2*b**2*(2*c**4 - 19*c**2*d**2 - 15*d**4))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d**2*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_714():
    f = (a + b*sin(e + f*x))**3*(c + d*sin(e + f*x))**(sympy.S(5)/2)
    F = -2*b**2*(a + b*sin(e + f*x))*(c + d*sin(e + f*x))**(sympy.S(7)/2)*cos(e + f*x)/(11*d*f) + 8*b**2*(c + d*sin(e + f*x))**(sympy.S(7)/2)*(-6*a*d + b*c)*cos(e + f*x)/(99*d**2*f) + 2*b*(c + d*sin(e + f*x))**(sympy.S(5)/2)*(-297*a**2*d**2 + 66*a*b*c*d - b**2*(8*c**2 + 81*d**2))*cos(e + f*x)/(693*d**2*f) - (c + d*sin(e + f*x))**(sympy.S(3)/2)*(1386*a**3*d**3 + 2970*a**2*b*c*d**2 - 66*a*b**2*d*(10*c**2 - 49*d**2) + 10*b**3*(8*c**3 + 67*c*d**2))*cos(e + f*x)/(3465*d**2*f) + sqrt(c + d*sin(e + f*x))*(-3696*a**3*c*d**3 - 990*a**2*b*d**2*(3*c**2 + 5*d**2) + 132*a*b**2*d*(5*c**3 - 57*c*d**2) - 10*b**3*(8*c**4 + 57*c**2*d**2 + 135*d**4))*cos(e + f*x)/(3465*d**2*f) - sqrt((c + d*sin(e + f*x))/(c + d))*(2*c**2 - 2*d**2)*(1848*a**3*c*d**3 + 495*a**2*b*d**2*(3*c**2 + 5*d**2) - 66*a*b**2*d*(5*c**3 - 57*c*d**2) + 5*b**3*(8*c**4 + 57*c**2*d**2 + 135*d**4))*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3465*d**3*f*sqrt(c + d*sin(e + f*x))) + sqrt(c + d*sin(e + f*x))*(462*a**3*d**3*(23*c**2 + 9*d**2) + 990*a**2*b*c*d**2*(3*c**2 + 29*d**2) - 66*a*b**2*d*(10*c**4 - 279*c**2*d**2 - 147*d**4) + 10*b**3*(8*c**5 + 51*c**3*d**2 + 741*c*d**4))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3465*d**3*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_715():
    f = (a + b*sin(e + f*x))**3*(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = -2*b**2*(a + b*sin(e + f*x))*(c + d*sin(e + f*x))**(sympy.S(5)/2)*cos(e + f*x)/(9*d*f) + 8*b**2*(c + d*sin(e + f*x))**(sympy.S(5)/2)*(-5*a*d + b*c)*cos(e + f*x)/(63*d**2*f) + 2*b*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(-189*a**2*d**2 + 54*a*b*c*d - b**2*(8*c**2 + 49*d**2))*cos(e + f*x)/(315*d**2*f) + sqrt(c + d*sin(e + f*x))*(-210*a**3*d**3 - 378*a**2*b*c*d**2 + 18*a*b**2*d*(6*c**2 - 25*d**2) - 2*b**3*(8*c**3 + 39*c*d**2))*cos(e + f*x)/(315*d**2*f) - sqrt((c + d*sin(e + f*x))/(c + d))*(2*c**2 - 2*d**2)*(105*a**3*d**3 + 189*a**2*b*c*d**2 - 9*a*b**2*d*(6*c**2 - 25*d**2) + b**3*(8*c**3 + 39*c*d**2))*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(315*d**3*f*sqrt(c + d*sin(e + f*x))) + sqrt(c + d*sin(e + f*x))*(840*a**3*c*d**3 + 378*a**2*b*d**2*(c**2 + 3*d**2) - 2*a*b**2*(54*c**3*d - 738*c*d**3) + 2*b**3*(8*c**4 + 33*c**2*d**2 + 147*d**4))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(315*d**3*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_716():
    f = (a + b*sin(e + f*x))**3*sqrt(c + d*sin(e + f*x))
    F = -2*b**2*(a + b*sin(e + f*x))*(c + d*sin(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(7*d*f) + 8*b**2*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(-4*a*d + b*c)*cos(e + f*x)/(35*d**2*f) + 2*b*sqrt(c + d*sin(e + f*x))*(-105*a**2*d**2 + 42*a*b*c*d - b**2*(8*c**2 + 25*d**2))*cos(e + f*x)/(105*d**2*f) + 2*b*sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2)*(-105*a**2*d**2 + 42*a*b*c*d - b**2*(8*c**2 + 25*d**2))*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(105*d**3*f*sqrt(c + d*sin(e + f*x))) + sqrt(c + d*sin(e + f*x))*(210*a**3*d**3 + 210*a**2*b*c*d**2 - 42*a*b**2*d*(2*c**2 - 9*d**2) + 2*b**3*(8*c**3 + 19*c*d**2))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(105*d**3*f*sqrt((c + d*sin(e + f*x))/(c + d)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_717():
    f = (a + b*sin(e + f*x))**3/sqrt(c + d*sin(e + f*x))
    F = -2*b**2*(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(5*d*f) + 8*b**2*sqrt(c + d*sin(e + f*x))*(-3*a*d + b*c)*cos(e + f*x)/(15*d**2*f) - 2*b*sqrt(c + d*sin(e + f*x))*(-45*a**2*d**2 + 30*a*b*c*d - b**2*(8*c**2 + 9*d**2))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d**3*f*sqrt((c + d*sin(e + f*x))/(c + d))) - sqrt((c + d*sin(e + f*x))/(c + d))*(-30*a**3*d**3 + 90*a**2*b*c*d**2 - 30*a*b**2*d*(2*c**2 + d**2) + 2*b**3*(8*c**3 + 7*c*d**2))*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d**3*f*sqrt(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_718():
    f = (a + b*sin(e + f*x))**3/(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = 2*b*sqrt(c + d*sin(e + f*x))*(-3*a**2*d**2 + 6*a*b*c*d - b**2*(4*c**2 - d**2))*cos(e + f*x)/(3*d**2*f*(c**2 - d**2)) - 2*b*sqrt((c + d*sin(e + f*x))/(c + d))*(-9*a**2*d**2 + 18*a*b*c*d - b**2*(8*c**2 + d**2))*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d**3*f*sqrt(c + d*sin(e + f*x))) + 2*(a + b*sin(e + f*x))*(-a*d + b*c)**2*cos(e + f*x)/(d*f*sqrt(c + d*sin(e + f*x))*(c**2 - d**2)) - sqrt(c + d*sin(e + f*x))*(-6*a**3*d**3 + 18*a**2*b*c*d**2 - 18*a*b**2*d*(2*c**2 - d**2) + 2*b**3*(8*c**3 - 5*c*d**2))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d**3*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_719():
    f = (a + b*sin(e + f*x))**3/(c + d*sin(e + f*x))**(sympy.S(5)/2)
    F = 2*(a + b*sin(e + f*x))*(-a*d + b*c)**2*cos(e + f*x)/(3*d*f*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(c**2 - d**2)) + 8*(-a*d + b*c)**2*(a*c*d + b*(c**2 - 2*d**2))*cos(e + f*x)/(3*d**2*f*sqrt(c + d*sin(e + f*x))*(c**2 - d**2)**2) - sqrt((c + d*sin(e + f*x))/(c + d))*(-2*a*d + 2*b*c)*(-a**2*d**2 + 2*a*b*c*d + b**2*(8*c**2 - 9*d**2))*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d**3*f*sqrt(c + d*sin(e + f*x))*(c**2 - d**2)) + sqrt(c + d*sin(e + f*x))*(8*a**3*c*d**3 - 6*a**2*b*d**2*(c**2 + 3*d**2) - 12*a*b**2*c*d*(c**2 - 3*d**2) + 2*b**3*(8*c**4 - 15*c**2*d**2 + 3*d**4))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(3*d**3*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_720():
    f = (a + b*sin(e + f*x))**3/(c + d*sin(e + f*x))**(sympy.S(7)/2)
    F = 2*(a + b*sin(e + f*x))*(-a*d + b*c)**2*cos(e + f*x)/(5*d*f*(c + d*sin(e + f*x))**(sympy.S(5)/2)*(c**2 - d**2)) - (-2*a*d + 2*b*c)*(a**2*d**2*(23*c**2 + 9*d**2) + 2*a*b*d*(7*c**3 - 39*c*d**2) + b**2*(8*c**4 - 21*c**2*d**2 + 45*d**4))*cos(e + f*x)/(15*d**2*f*sqrt(c + d*sin(e + f*x))*(c**2 - d**2)**3) + 8*(-a*d + b*c)**2*(2*a*c*d + b*(c**2 - 3*d**2))*cos(e + f*x)/(15*d**2*f*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(c**2 - d**2)**2) - sqrt((c + d*sin(e + f*x))/(c + d))*(16*a**3*c*d**3 - 6*a**2*b*d**2*(3*c**2 + 5*d**2) - 12*a*b**2*c*d*(c**2 - 5*d**2) - 2*b**3*(8*c**4 - 15*c**2*d**2 + 15*d**4))*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d**3*f*sqrt(c + d*sin(e + f*x))*(c**2 - d**2)**2) - sqrt(c + d*sin(e + f*x))*(-2*a*d + 2*b*c)*(a**2*d**2*(23*c**2 + 9*d**2) + 2*a*b*d*(7*c**3 - 39*c*d**2) + b**2*(8*c**4 - 21*c**2*d**2 + 45*d**4))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(15*d**3*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_721():
    f = (a + b*sin(e + f*x))**3/(c + d*sin(e + f*x))**(sympy.S(9)/2)
    F = 2*(a + b*sin(e + f*x))*(-a*d + b*c)**2*cos(e + f*x)/(7*d*f*(c + d*sin(e + f*x))**(sympy.S(7)/2)*(c**2 - d**2)) + (32*a**3*c*d**3*(11*c**2 + 13*d**2) - 18*a**2*b*d**2*(5*c**4 + 102*c**2*d**2 + 21*d**4) - 12*a*b**2*c*d*(3*c**4 - 62*c**2*d**2 - 133*d**4) - 2*b**3*(8*c**6 - 23*c**4*d**2 + 294*c**2*d**4 + 105*d**6))*cos(e + f*x)/(105*d**2*f*sqrt(c + d*sin(e + f*x))*(c**2 - d**2)**4) - (-2*a*d + 2*b*c)*(a**2*d**2*(71*c**2 + 25*d**2) + a*b*(26*c**3*d - 218*c*d**3) + b**2*(8*c**4 - 17*c**2*d**2 + 105*d**4))*cos(e + f*x)/(105*d**2*f*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(c**2 - d**2)**3) + 8*(-a*d + b*c)**2*(3*a*c*d + b*(c**2 - 4*d**2))*cos(e + f*x)/(35*d**2*f*(c + d*sin(e + f*x))**(sympy.S(5)/2)*(c**2 - d**2)**2) + sqrt((c + d*sin(e + f*x))/(c + d))*(-2*a*d + 2*b*c)*(a**2*d**2*(71*c**2 + 25*d**2) + a*b*(26*c**3*d - 218*c*d**3) + b**2*(8*c**4 - 17*c**2*d**2 + 105*d**4))*elliptic_f(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(105*d**3*f*sqrt(c + d*sin(e + f*x))*(c**2 - d**2)**3) + sqrt(c + d*sin(e + f*x))*(32*a**3*c*d**3*(11*c**2 + 13*d**2) - 18*a**2*b*d**2*(5*c**4 + 102*c**2*d**2 + 21*d**4) - 12*a*b**2*c*d*(3*c**4 - 62*c**2*d**2 - 133*d**4) - 2*b**3*(8*c**6 - 23*c**4*d**2 + 294*c**2*d**4 + 105*d**6))*elliptic_e(e/2 + f*x/2 - pi/4, 2*d/(c + d))/(105*d**3*f*sqrt((c + d*sin(e + f*x))/(c + d))*(c**2 - d**2)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_722():
    f = (c + d*sin(e + f*x))**(sympy.S(5)/2)/(a + b*sin(e + f*x))
    F = ((Integer(-2) * (Symbol('d'))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(3) * Symbol('b') * Symbol('f')))**(Integer(-1))) + ((Integer(2) * Symbol('d') * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.elliptic_e(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('d') * ((Integer(6) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * (Symbol('c'))**(Integer(2))) + (Symbol('d'))**(Integer(2)))))) * sympy.elliptic_f(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('a') + Symbol('b')) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_723():
    f = (c + d*sin(e + f*x))**(sympy.S(3)/2)/(a + b*sin(e + f*x))
    F = ((Integer(2) * Symbol('d') * sympy.elliptic_e(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('b') * Symbol('f') * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.elliptic_f(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('a') + Symbol('b')) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_724():
    f = sqrt(c + d*sin(e + f*x))/(a + b*sin(e + f*x))
    F = ((Integer(2) * Symbol('d') * sympy.elliptic_f(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Symbol('b') * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Symbol('b') * (Symbol('a') + Symbol('b')) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_725():
    f = 1/((a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x)))
    F = (Integer(2) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * (((Symbol('a') + Symbol('b')) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_726():
    f = 1/((a + b*sin(e + f*x))*(c + d*sin(e + f*x))**(sympy.S(3)/2))
    F = ((Integer(-2) * (Symbol('d'))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x)))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('d') * sympy.elliptic_e(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))) * Symbol('f') * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * (((Symbol('a') + Symbol('b')) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_727():
    f = 1/((a + b*sin(e + f*x))*(c + d*sin(e + f*x))**(sympy.S(5)/2))
    F = ((Integer(-2) * (Symbol('d'))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x)))) * ((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))) * Symbol('f') * ((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * ((Integer(7) * Symbol('b') * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c') * Symbol('d'))) + (Integer(-1) * (Integer(3) * Symbol('b') * (Symbol('d'))**(Integer(2))))) * sympy.cos((Symbol('e') + (Symbol('f') * x)))) * ((Integer(3) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * ((Integer(7) * Symbol('b') * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c') * Symbol('d'))) + (Integer(-1) * (Integer(3) * Symbol('b') * (Symbol('d'))**(Integer(2))))) * sympy.elliptic_e(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(3) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))))**(Integer(2)) * Symbol('f') * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * Symbol('d') * sympy.elliptic_f(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * (((Symbol('a') + Symbol('b')) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_728():
    f = (c + d*sin(e + f*x))**(sympy.S(7)/2)/(a + b*sin(e + f*x))**2
    F = ((Symbol('d') * ((Integer(6) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('d'))**(Integer(2)))))))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f')))**(Integer(-1))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * ((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + ((((Integer(29) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3)))) + ((Symbol('b'))**(Integer(3)) * ((Integer(3) * (Symbol('c'))**(Integer(3))) + (Integer(-1) * (Integer(20) * Symbol('c') * (Symbol('d'))**(Integer(2)))))) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * ((Integer(9) * (Symbol('c'))**(Integer(2)) * Symbol('d')) + (Integer(-1) * (Integer(12) * (Symbol('d'))**(Integer(3)))))))) * sympy.elliptic_e(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(24) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(3))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(4)) * (Symbol('d'))**(Integer(4)))) + (Integer(-1) * (Integer(12) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('c') * Symbol('d') * ((Symbol('c'))**(Integer(2)) + (Integer(3) * (Symbol('d'))**(Integer(2)))))) + (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * ((Symbol('c'))**(Integer(2)) + (Integer(8) * (Symbol('d'))**(Integer(2))))) + ((Symbol('b'))**(Integer(4)) * ((Integer(3) * (Symbol('c'))**(Integer(4))) + (Integer(16) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(2) * (Symbol('d'))**(Integer(4)))))) * sympy.elliptic_f(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('c')) + (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('d')) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2)) * Symbol('d')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(4)) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_729():
    f = (c + d*sin(e + f*x))**(sympy.S(5)/2)/(a + b*sin(e + f*x))**2
    F = (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(2) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Integer(2) * (Symbol('d'))**(Integer(2)))))))) * sympy.elliptic_e(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))))**(Integer(-1)))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Symbol('c'))**(Integer(2)) + (Integer(4) * (Symbol('d'))**(Integer(2))))))) * sympy.elliptic_f(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * (((Symbol('b'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('c')) + (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('d')) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2)) * Symbol('d')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(3)) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_730():
    f = (c + d*sin(e + f*x))**(sympy.S(3)/2)/(a + b*sin(e + f*x))**2
    F = ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.elliptic_e(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(2) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + ((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Symbol('c'))**(Integer(2)) + (Integer(2) * (Symbol('d'))**(Integer(2))))))) * sympy.elliptic_f(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('c')) + ((Symbol('a'))**(Integer(2)) * Symbol('d')) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('d')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_731():
    f = sqrt(c + d*sin(e + f*x))/(a + b*sin(e + f*x))**2
    F = ((Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + ((sympy.elliptic_e(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.elliptic_f(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((((Integer(2) * Symbol('a') * Symbol('b') * Symbol('c')) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * Symbol('d'))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * Symbol('d')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_732():
    f = 1/((a + b*sin(e + f*x))**2*sqrt(c + d*sin(e + f*x)))
    F = (((Symbol('b'))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + ((Symbol('b') * sympy.elliptic_e(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.elliptic_f(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((((Integer(2) * Symbol('a') * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('d'))) + ((Symbol('b'))**(Integer(2)) * Symbol('d'))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_733():
    f = 1/((a + b*sin(e + f*x))**2*(c + d*sin(e + f*x))**(sympy.S(3)/2))
    F = ((Symbol('d') * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + ((Symbol('b'))**(Integer(2)) * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Integer(3) * (Symbol('d'))**(Integer(2))))))) * sympy.cos((Symbol('e') + (Symbol('f') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + ((Symbol('b'))**(Integer(2)) * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Integer(3) * (Symbol('d'))**(Integer(2))))))) * sympy.elliptic_e(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))) * Symbol('f') * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.elliptic_f(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((Symbol('b') * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('d'))) + (Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('d'))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_734():
    f = 1/((a + b*sin(e + f*x))**2*(c + d*sin(e + f*x))**(sympy.S(5)/2))
    F = ((Symbol('d') * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + ((Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('d'))**(Integer(2))))))) * sympy.cos((Symbol('e') + (Symbol('f') * x)))) * ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))) * Symbol('f') * ((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(8) * (Symbol('a'))**(Integer(3)) * Symbol('c') * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Integer(8) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('c') * (Symbol('d'))**(Integer(4)))) + (Integer(-1) * (Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('d'))**(Integer(3)) * ((Integer(5) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('d'))**(Integer(2))))))) + (Integer(-1) * ((Symbol('b'))**(Integer(3)) * ((Integer(3) * (Symbol('c'))**(Integer(4)) * Symbol('d')) + (Integer(-1) * (Integer(26) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(3)))) + (Integer(15) * (Symbol('d'))**(Integer(5))))))) * sympy.cos((Symbol('e') + (Symbol('f') * x)))) * ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(8) * (Symbol('a'))**(Integer(3)) * Symbol('c') * (Symbol('d'))**(Integer(3))) + (Integer(-1) * (Integer(8) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('c') * (Symbol('d'))**(Integer(3)))) + (Integer(-1) * (Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('d'))**(Integer(2)) * ((Integer(5) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('d'))**(Integer(2))))))) + (Integer(-1) * ((Symbol('b'))**(Integer(3)) * ((Integer(3) * (Symbol('c'))**(Integer(4))) + (Integer(-1) * (Integer(26) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + (Integer(15) * (Symbol('d'))**(Integer(4))))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))))**(Integer(2)) * Symbol('f') * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + ((Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('d'))**(Integer(2))))))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('d'))) + (Integer(5) * (Symbol('b'))**(Integer(2)) * Symbol('d'))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_735():
    f = (c + d*sin(e + f*x))**(sympy.S(9)/2)/(a + b*sin(e + f*x))**3
    F = ((Symbol('d') * ((Integer(36) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(35) * (Symbol('a'))**(Integer(4)) * (Symbol('d'))**(Integer(3)))) + ((Symbol('b'))**(Integer(4)) * Symbol('d') * ((Integer(45) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(8) * (Symbol('d'))**(Integer(2)))))) + (Integer(-1) * (Integer(18) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('c') * ((Symbol('c'))**(Integer(2)) + (Integer(5) * (Symbol('d'))**(Integer(2)))))) + ((Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('d') * ((Integer(9) * (Symbol('c'))**(Integer(2))) + (Integer(61) * (Symbol('d'))**(Integer(2)))))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(12) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('f')))**(Integer(-1))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Integer(6) * Symbol('a') * Symbol('b') * Symbol('c')) + (Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('d')) + (Integer(-1) * (Integer(13) * (Symbol('b'))**(Integer(2)) * Symbol('d')))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * ((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * ((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1))) + ((((Integer(185) * (Symbol('a'))**(Integer(4)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(3))) + (Integer(-1) * (Integer(105) * (Symbol('a'))**(Integer(5)) * (Symbol('d'))**(Integer(4)))) + (Integer(-1) * ((Symbol('b'))**(Integer(5)) * Symbol('c') * Symbol('d') * ((Integer(51) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(104) * (Symbol('d'))**(Integer(2))))))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * ((Integer(3) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(13) * (Symbol('d'))**(Integer(2))))))) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * Symbol('c') * Symbol('d') * ((Integer(21) * (Symbol('c'))**(Integer(2))) + (Integer(361) * (Symbol('d'))**(Integer(2)))))) + (Integer(9) * Symbol('a') * (Symbol('b'))**(Integer(4)) * ((Integer(2) * (Symbol('c'))**(Integer(4))) + (Integer(17) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(8) * (Symbol('d'))**(Integer(4))))))) * sympy.elliptic_e(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(12) * (Symbol('b'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('f') * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(150) * (Symbol('a'))**(Integer(5)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Integer(105) * (Symbol('a'))**(Integer(6)) * (Symbol('d'))**(Integer(5)))) + (Integer(-1) * (Integer(12) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(3)) * Symbol('c') * (Symbol('d'))**(Integer(2)) * ((Integer(4) * (Symbol('c'))**(Integer(2))) + (Integer(29) * (Symbol('d'))**(Integer(2)))))) + ((Symbol('a'))**(Integer(4)) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(3)) * ((Integer(26) * (Symbol('c'))**(Integer(2))) + (Integer(223) * (Symbol('d'))**(Integer(2))))) + (Integer(-1) * ((Symbol('b'))**(Integer(6)) * Symbol('d') * ((Integer(57) * (Symbol('c'))**(Integer(4))) + (Integer(136) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(8) * (Symbol('d'))**(Integer(4)))))) + (Integer(6) * Symbol('a') * (Symbol('b'))**(Integer(5)) * Symbol('c') * ((Integer(3) * (Symbol('c'))**(Integer(4))) + (Integer(38) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(48) * (Symbol('d'))**(Integer(4))))) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(4)) * Symbol('d') * ((Integer(33) * (Symbol('c'))**(Integer(4))) + (Integer(70) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(128) * (Symbol('d'))**(Integer(4))))))) * sympy.elliptic_f(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Integer(12) * (Symbol('b'))**(Integer(5)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * ((Integer(20) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(44) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('c') * Symbol('d'))) + (Integer(35) * (Symbol('a'))**(Integer(4)) * (Symbol('d'))**(Integer(2))) + (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * ((Integer(4) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(43) * (Symbol('d'))**(Integer(2)))))) + ((Symbol('b'))**(Integer(4)) * ((Integer(4) * (Symbol('c'))**(Integer(2))) + (Integer(63) * (Symbol('d'))**(Integer(2)))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(5)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_736():
    f = (c + d*sin(e + f*x))**(sympy.S(7)/2)/(a + b*sin(e + f*x))**3
    F = (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Integer(6) * Symbol('a') * Symbol('b') * Symbol('c')) + (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('d')) + (Integer(-1) * (Integer(11) * (Symbol('b'))**(Integer(2)) * Symbol('d')))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * ((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((((Integer(8) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(4)) * (Symbol('d'))**(Integer(3)))) + ((Symbol('b'))**(Integer(4)) * Symbol('d') * ((Integer(13) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(8) * (Symbol('d'))**(Integer(2)))))) + (Integer(-1) * (Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('c') * ((Integer(3) * (Symbol('c'))**(Integer(2))) + (Integer(13) * (Symbol('d'))**(Integer(2)))))) + ((Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('d') * ((Integer(5) * (Symbol('c'))**(Integer(2))) + (Integer(29) * (Symbol('d'))**(Integer(2)))))) * sympy.elliptic_e(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('f') * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(5) * (Symbol('a'))**(Integer(4)) * (Symbol('d'))**(Integer(3))) + ((Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('d') * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Integer(11) * (Symbol('d'))**(Integer(2)))))) + (Integer(-1) * (Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('c') * ((Symbol('c'))**(Integer(2)) + (Integer(5) * (Symbol('d'))**(Integer(2)))))) + ((Symbol('b'))**(Integer(4)) * Symbol('d') * ((Integer(5) * (Symbol('c'))**(Integer(2))) + (Integer(8) * (Symbol('d'))**(Integer(2)))))) * sympy.elliptic_f(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Integer(12) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(36) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('c') * Symbol('d'))) + (Integer(15) * (Symbol('a'))**(Integer(4)) * (Symbol('d'))**(Integer(2))) + (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * ((Integer(4) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(19) * (Symbol('d'))**(Integer(2)))))) + ((Symbol('b'))**(Integer(4)) * ((Integer(4) * (Symbol('c'))**(Integer(2))) + (Integer(35) * (Symbol('d'))**(Integer(2)))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(4)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_737():
    f = (c + d*sin(e + f*x))**(sympy.S(5)/2)/(a + b*sin(e + f*x))**3
    F = (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('c')) + ((Symbol('a'))**(Integer(2)) * Symbol('d')) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('d')))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(4) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + ((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('c')) + ((Symbol('a'))**(Integer(2)) * Symbol('d')) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('d')))) * sympy.elliptic_e(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('f') * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(3) * (Symbol('a'))**(Integer(4)) * (Symbol('d'))**(Integer(3))) + ((Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('d') * ((Integer(7) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('d'))**(Integer(2)))))) + ((Symbol('b'))**(Integer(4)) * Symbol('d') * ((Integer(11) * (Symbol('c'))**(Integer(2))) + (Integer(8) * (Symbol('d'))**(Integer(2))))) + (Integer(-1) * (Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('c') * ((Integer(3) * (Symbol('c'))**(Integer(2))) + (Integer(11) * (Symbol('d'))**(Integer(2))))))) * sympy.elliptic_f(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(28) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('c') * Symbol('d'))) + (Integer(3) * (Symbol('a'))**(Integer(4)) * (Symbol('d'))**(Integer(2))) + (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * ((Integer(4) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('d'))**(Integer(2)))))) + ((Symbol('b'))**(Integer(4)) * ((Integer(4) * (Symbol('c'))**(Integer(2))) + (Integer(15) * (Symbol('d'))**(Integer(2)))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_738():
    f = (c + d*sin(e + f*x))**(sympy.S(3)/2)/(a + b*sin(e + f*x))**3
    F = ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1))) + ((((Integer(6) * Symbol('a') * Symbol('b') * Symbol('c')) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * Symbol('d'))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2)) * Symbol('d')))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + ((((Integer(6) * Symbol('a') * Symbol('b') * Symbol('c')) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * Symbol('d'))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2)) * Symbol('d')))) * sympy.elliptic_e(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(4) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('f') * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Integer(6) * Symbol('a') * Symbol('b') * Symbol('c')) + ((Symbol('a'))**(Integer(2)) * Symbol('d')) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2)) * Symbol('d')))) * sympy.elliptic_f(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(20) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('c') * Symbol('d')) + ((Symbol('a'))**(Integer(4)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * ((Symbol('b'))**(Integer(4)) * ((Integer(4) * (Symbol('c'))**(Integer(2))) + (Integer(3) * (Symbol('d'))**(Integer(2)))))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * ((Integer(4) * (Symbol('c'))**(Integer(2))) + (Integer(5) * (Symbol('d'))**(Integer(2))))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_739():
    f = sqrt(c + d*sin(e + f*x))/(a + b*sin(e + f*x))**3
    F = ((Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * ((Integer(6) * Symbol('a') * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('d'))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * Symbol('d')))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + ((((Integer(6) * Symbol('a') * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('d'))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * Symbol('d')))) * sympy.elliptic_e(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('c')) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * Symbol('d'))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * Symbol('d')))) * sympy.elliptic_f(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(12) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(12) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(4)) * (Symbol('d'))**(Integer(2)))) + (Integer(-1) * ((Symbol('b'))**(Integer(4)) * ((Integer(4) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Symbol('d'))**(Integer(2)))))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * ((Integer(4) * (Symbol('c'))**(Integer(2))) + (Integer(5) * (Symbol('d'))**(Integer(2))))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(3)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_740():
    f = 1/((a + b*sin(e + f*x))**3*sqrt(c + d*sin(e + f*x)))
    F = (((Symbol('b'))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('d'))) + ((Symbol('b'))**(Integer(2)) * Symbol('d'))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('d'))) + ((Symbol('b'))**(Integer(2)) * Symbol('d'))) * sympy.elliptic_e(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('f') * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(6) * Symbol('a') * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('d'))) + ((Symbol('b'))**(Integer(2)) * Symbol('d'))) * sympy.elliptic_f(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(20) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(4)) * (Symbol('d'))**(Integer(2)))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * ((Integer(4) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('d'))**(Integer(2))))))) + (Integer(-1) * ((Symbol('b'))**(Integer(4)) * ((Integer(4) * (Symbol('c'))**(Integer(2))) + (Integer(3) * (Symbol('d'))**(Integer(2))))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_741():
    f = 1/((a + b*sin(e + f*x))**3*(c + d*sin(e + f*x))**(sympy.S(3)/2))
    F = ((Integer(-1) * (Symbol('d') * ((Integer(8) * (Symbol('a'))**(Integer(4)) * (Symbol('d'))**(Integer(3))) + ((Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('d') * ((Integer(13) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(29) * (Symbol('d'))**(Integer(2)))))) + (Integer(-1) * ((Symbol('b'))**(Integer(4)) * Symbol('d') * ((Integer(7) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(15) * (Symbol('d'))**(Integer(2))))))) + (Integer(-1) * (Integer(6) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('c') * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2))))))) * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x)))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * ((Integer(6) * Symbol('a') * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(11) * (Symbol('a'))**(Integer(2)) * Symbol('d'))) + (Integer(5) * (Symbol('b'))**(Integer(2)) * Symbol('d'))) * sympy.cos((Symbol('e') + (Symbol('f') * x)))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(8) * (Symbol('a'))**(Integer(4)) * (Symbol('d'))**(Integer(3))) + ((Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('d') * ((Integer(13) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(29) * (Symbol('d'))**(Integer(2)))))) + (Integer(-1) * ((Symbol('b'))**(Integer(4)) * Symbol('d') * ((Integer(7) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(15) * (Symbol('d'))**(Integer(2))))))) + (Integer(-1) * (Integer(6) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('c') * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2))))))) * sympy.elliptic_e(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))) * Symbol('f') * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * ((Integer(6) * Symbol('a') * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(11) * (Symbol('a'))**(Integer(2)) * Symbol('d'))) + (Integer(5) * (Symbol('b'))**(Integer(2)) * Symbol('d'))) * sympy.elliptic_f(((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * ((Integer(28) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('c') * Symbol('d'))) + (Integer(-1) * (Integer(35) * (Symbol('a'))**(Integer(4)) * (Symbol('d'))**(Integer(2)))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * ((Integer(4) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(19) * (Symbol('d'))**(Integer(2))))))) + (Integer(-1) * ((Symbol('b'))**(Integer(4)) * ((Integer(4) * (Symbol('c'))**(Integer(2))) + (Integer(15) * (Symbol('d'))**(Integer(2))))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('e') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('f') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_742():
    f = sqrt(a + b*sin(e + f*x))*(c + d*sin(e + f*x))**(sympy.S(5)/2)
    F = ((sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('c') + (Integer(-1) * Symbol('d'))) * sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(14) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + ((Symbol('b'))**(Integer(2)) * ((Integer(33) * (Symbol('c'))**(Integer(2))) + (Integer(16) * (Symbol('d'))**(Integer(2)))))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(24) * (Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * (((Integer(8) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * Symbol('f')))**(Integer(-1)) * (sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3)))) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * ((Integer(15) * (Symbol('c'))**(Integer(2))) + (Integer(4) * (Symbol('d'))**(Integer(2)))))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(3)) * ((Symbol('c'))**(Integer(3)) + (Integer(4) * Symbol('c') * (Symbol('d'))**(Integer(2))))))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))) + (Integer(-1) * ((((Integer(14) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + ((Symbol('b'))**(Integer(2)) * ((Integer(33) * (Symbol('c'))**(Integer(2))) + (Integer(16) * (Symbol('d'))**(Integer(2)))))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(24) * Symbol('b') * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * ((Integer(13) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(12) * Symbol('b') * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(3) * Symbol('b') * Symbol('f')))**(Integer(-1)))) + ((((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(6) * Symbol('a') * Symbol('b') * Symbol('d') * ((Integer(2) * Symbol('c')) + Symbol('d')))) + ((Symbol('b'))**(Integer(2)) * ((Integer(33) * (Symbol('c'))**(Integer(2))) + (Integer(26) * Symbol('c') * Symbol('d')) + (Integer(16) * (Symbol('d'))**(Integer(2)))))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(24) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_743():
    f = sqrt(a + b*sin(e + f*x))*(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = ((sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('c') + (Integer(-1) * Symbol('d'))) * sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(4) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1))) + ((sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(6) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + ((Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('c'))**(Integer(2))) + (Integer(4) * (Symbol('d'))**(Integer(2)))))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * Symbol('f')))**(Integer(-1))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(2) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(5) * Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(4) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))) + (Integer(2) * Symbol('b') * Symbol('d'))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x))) * ((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_744():
    f = sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))
    F = ((sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('c') + (Integer(-1) * Symbol('d'))) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1))) + ((sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('b') * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_745():
    f = sqrt(a + b*sin(e + f*x))/sqrt(c + d*sin(e + f*x))
    F = (Integer(2) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * Symbol('f')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_746():
    f = sqrt(a + b*sin(e + f*x))/(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = -sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*sqrt(a + b)*(2*a - 2*b)*(c + d*sin(e + f*x))*elliptic_e(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(f*(c - d)*sqrt(c + d)*(-a*d + b*c)) + sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*sqrt(a + b)*(2*a - 2*b)*(c + d*sin(e + f*x))*elliptic_f(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(f*(c - d)*sqrt(c + d)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_747():
    f = sqrt(a + b*sin(e + f*x))/(c + d*sin(e + f*x))**(sympy.S(5)/2)
    F = 2*d*sqrt(a + b*sin(e + f*x))*cos(e + f*x)/(f*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(3*c**2 - 3*d**2)) + sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*sqrt(a + b)*(2*a - 2*b)*(c + d*sin(e + f*x))*(3*c + d)*elliptic_f(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(3*f*(c - d)**2*(c + d)**(sympy.S(3)/2)*(-a*d + b*c)) + sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*sqrt(a + b)*(2*a - 2*b)*(c + d*sin(e + f*x))*(4*a*c*d - b*(3*c**2 + d**2))*elliptic_e(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(3*f*(c - d)**2*(c + d)**(sympy.S(3)/2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_748():
    f = (a + b*sin(e + f*x))**(sympy.S(3)/2)*(c + d*sin(e + f*x))**(sympy.S(5)/2)
    F = (((Integer(192) * (Symbol('b'))**(Integer(2)) * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1)) * (sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('c') + (Integer(-1) * Symbol('d'))) * sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(57) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(9) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3)))) + (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * ((Integer(337) * (Symbol('c'))**(Integer(2))) + (Integer(156) * (Symbol('d'))**(Integer(2))))) + ((Symbol('b'))**(Integer(3)) * ((Integer(15) * (Symbol('c'))**(Integer(3))) + (Integer(284) * Symbol('c') * (Symbol('d'))**(Integer(2)))))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) + (Integer(-1) * (((Integer(64) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('d'))**(Integer(2)) * Symbol('f')))**(Integer(-1)) * (sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(20) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(3))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(4)) * (Symbol('d'))**(Integer(4)))) + (Integer(-1) * (Integer(60) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('c') * Symbol('d') * ((Symbol('c'))**(Integer(2)) + (Integer(4) * (Symbol('d'))**(Integer(2)))))) + (Integer(-1) * (Integer(6) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * ((Integer(15) * (Symbol('c'))**(Integer(2))) + (Integer(4) * (Symbol('d'))**(Integer(2)))))) + ((Symbol('b'))**(Integer(4)) * ((Integer(5) * (Symbol('c'))**(Integer(4))) + (Integer(-1) * (Integer(120) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + (Integer(-1) * (Integer(48) * (Symbol('d'))**(Integer(4))))))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))) + (Integer(-1) * ((((Integer(57) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(9) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3)))) + (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * ((Integer(337) * (Symbol('c'))**(Integer(2))) + (Integer(156) * (Symbol('d'))**(Integer(2))))) + ((Symbol('b'))**(Integer(3)) * ((Integer(15) * (Symbol('c'))**(Integer(3))) + (Integer(284) * Symbol('c') * (Symbol('d'))**(Integer(2)))))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(192) * Symbol('b') * Symbol('d') * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(54) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(9) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + ((Symbol('b'))**(Integer(2)) * ((Integer(59) * (Symbol('c'))**(Integer(2))) + (Integer(36) * (Symbol('d'))**(Integer(2)))))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(96) * Symbol('b') * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * ((Integer(17) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(24) * Symbol('b') * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(4) * Symbol('b') * Symbol('f')))**(Integer(-1)))) + (((Integer(192) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1)) * (((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(9) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('d'))**(Integer(2)) * ((Integer(17) * Symbol('c')) + (Integer(6) * Symbol('d'))))) + (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * ((Integer(73) * (Symbol('c'))**(Integer(2))) + (Integer(36) * Symbol('c') * Symbol('d')) + (Integer(28) * (Symbol('d'))**(Integer(2))))) + ((Symbol('b'))**(Integer(3)) * ((Integer(15) * (Symbol('c'))**(Integer(3))) + (Integer(118) * (Symbol('c'))**(Integer(2)) * Symbol('d')) + (Integer(284) * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(72) * (Symbol('d'))**(Integer(3)))))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_749():
    f = (a + b*sin(e + f*x))**(sympy.S(3)/2)*(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = ((sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('c') + (Integer(-1) * Symbol('d'))) * sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(38) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + ((Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('c'))**(Integer(2))) + (Integer(16) * (Symbol('d'))**(Integer(2)))))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(24) * Symbol('b') * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1))) + ((sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * ((Integer(10) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Integer(12) * (Symbol('d'))**(Integer(2)))))))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('d'))**(Integer(2)) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((((Integer(38) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + ((Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('c'))**(Integer(2))) + (Integer(16) * (Symbol('d'))**(Integer(2)))))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(24) * Symbol('d') * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(7) * Symbol('a') * Symbol('d'))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(12) * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(6) * Symbol('a') * Symbol('b') * Symbol('d') * ((Integer(4) * Symbol('c')) + Symbol('d')))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('c'))**(Integer(2))) + (Integer(14) * Symbol('c') * Symbol('d')) + (Integer(16) * (Symbol('d'))**(Integer(2))))))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(24) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * Symbol('f')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_750():
    f = (a + b*sin(e + f*x))**(sympy.S(3)/2)*sqrt(c + d*sin(e + f*x))
    F = ((sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('c') + (Integer(-1) * Symbol('d'))) * sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Symbol('b') * Symbol('c')) + (Integer(5) * Symbol('a') * Symbol('d'))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(4) * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1))) + ((sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(6) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Integer(4) * (Symbol('d'))**(Integer(2)))))))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(4) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('d'))**(Integer(2)) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(5) * Symbol('a') * Symbol('d'))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(4) * Symbol('d') * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(2) * Symbol('f')))**(Integer(-1)))) + ((((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(3) * Symbol('a') * Symbol('d')) + (Symbol('b') * (Symbol('c') + (Integer(2) * Symbol('d'))))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_751():
    f = (a + b*sin(e + f*x))**(sympy.S(3)/2)/sqrt(c + d*sin(e + f*x))
    F = (Integer(-1) * ((Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('b') * (Symbol('c') + (Integer(-1) * Symbol('d')))) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.Function('EllipticPi')((((Symbol('a') + Symbol('b')) * Symbol('d')) * ((Symbol('b') * (Symbol('c') + Symbol('d'))))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_752():
    f = (a + b*sin(e + f*x))**(sympy.S(3)/2)/(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * Symbol('d') * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1)) * Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('b') * (Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('d'))))) + (Symbol('a') * Symbol('d'))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) + ((Integer(2) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.Function('EllipticPi')((((Symbol('a') + Symbol('b')) * Symbol('d')) * ((Symbol('b') * (Symbol('c') + Symbol('d'))))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_753():
    f = (a + b*sin(e + f*x))**(sympy.S(3)/2)/(c + d*sin(e + f*x))**(sympy.S(5)/2)
    F = sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*sqrt(a + b)*(2*a - 2*b)*(c + d*sin(e + f*x))*(a*(3*c + d) - b*(c + 3*d))*elliptic_f(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(3*f*(c - d)**2*(c + d)**(sympy.S(3)/2)*(-a*d + b*c)) - sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*sqrt(a + b)*(8*a - 8*b)*(c + d*sin(e + f*x))*(a*c - b*d)*elliptic_e(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(3*f*(c - d)**2*(c + d)**(sympy.S(3)/2)*(-a*d + b*c)) - sqrt(a + b*sin(e + f*x))*(-2*a*d + 2*b*c)*cos(e + f*x)/(f*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(3*c**2 - 3*d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_754():
    f = (a + b*sin(e + f*x))**(sympy.S(5)/2)*(c + d*sin(e + f*x))**(sympy.S(5)/2)
    F = (((Integer(1920) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1)) * (sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('c') + (Integer(-1) * Symbol('d'))) * sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(360) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(3))) + (Integer(-1) * (Integer(45) * (Symbol('a'))**(Integer(4)) * (Symbol('d'))**(Integer(4)))) + (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * ((Integer(1877) * (Symbol('c'))**(Integer(2))) + (Integer(846) * (Symbol('d'))**(Integer(2))))) + (Integer(8) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('d') * ((Integer(45) * (Symbol('c'))**(Integer(3))) + (Integer(791) * Symbol('c') * (Symbol('d'))**(Integer(2))))) + (Integer(-1) * ((Symbol('b'))**(Integer(4)) * ((Integer(45) * (Symbol('c'))**(Integer(4))) + (Integer(-1) * (Integer(1692) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + (Integer(-1) * (Integer(1024) * (Symbol('d'))**(Integer(4)))))))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) + (Integer(-1) * (((Integer(128) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('d'))**(Integer(3)) * Symbol('f')))**(Integer(-1)) * (sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * ((Integer(28) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(3))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(4)) * (Symbol('d'))**(Integer(4)))) + (Integer(28) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('c') * Symbol('d') * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Integer(20) * (Symbol('d'))**(Integer(2)))))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * ((Integer(89) * (Symbol('c'))**(Integer(2))) + (Integer(20) * (Symbol('d'))**(Integer(2)))))) + (Integer(-1) * ((Symbol('b'))**(Integer(4)) * ((Integer(3) * (Symbol('c'))**(Integer(4))) + (Integer(40) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(240) * (Symbol('d'))**(Integer(4))))))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))) + (Integer(-1) * ((((Integer(360) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(3))) + (Integer(-1) * (Integer(45) * (Symbol('a'))**(Integer(4)) * (Symbol('d'))**(Integer(4)))) + (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * ((Integer(1877) * (Symbol('c'))**(Integer(2))) + (Integer(846) * (Symbol('d'))**(Integer(2))))) + (Integer(8) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('d') * ((Integer(45) * (Symbol('c'))**(Integer(3))) + (Integer(791) * Symbol('c') * (Symbol('d'))**(Integer(2))))) + (Integer(-1) * ((Symbol('b'))**(Integer(4)) * ((Integer(45) * (Symbol('c'))**(Integer(4))) + (Integer(-1) * (Integer(1692) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + (Integer(-1) * (Integer(1024) * (Symbol('d'))**(Integer(4)))))))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(1920) * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(917) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(15) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3))) + (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * ((Integer(345) * (Symbol('c'))**(Integer(2))) + (Integer(772) * (Symbol('d'))**(Integer(2))))) + (Integer(-1) * ((Symbol('b'))**(Integer(3)) * ((Integer(45) * (Symbol('c'))**(Integer(3))) + (Integer(-1) * (Integer(516) * Symbol('c') * (Symbol('d'))**(Integer(2)))))))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(960) * Symbol('b') * Symbol('d') * Symbol('f')))**(Integer(-1)))) + (((Integer(1920) * (Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1)) * (((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(45) * (Symbol('a'))**(Integer(4)) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Integer(30) * (Symbol('a'))**(Integer(3)) * Symbol('b') * (Symbol('d'))**(Integer(3)) * ((Integer(11) * Symbol('c')) + (Integer(3) * Symbol('d'))))) + (Integer(30) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * ((Integer(64) * (Symbol('c'))**(Integer(2))) + (Integer(23) * Symbol('c') * Symbol('d')) + (Integer(22) * (Symbol('d'))**(Integer(2))))) + (Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('d') * ((Integer(165) * (Symbol('c'))**(Integer(3))) + (Integer(917) * (Symbol('c'))**(Integer(2)) * Symbol('d')) + (Integer(2392) * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(516) * (Symbol('d'))**(Integer(3))))) + (Integer(-1) * ((Symbol('b'))**(Integer(4)) * ((Integer(45) * (Symbol('c'))**(Integer(4))) + (Integer(-1) * (Integer(30) * (Symbol('c'))**(Integer(3)) * Symbol('d'))) + (Integer(-1) * (Integer(1692) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + (Integer(-1) * (Integer(1544) * Symbol('c') * (Symbol('d'))**(Integer(3)))) + (Integer(-1) * (Integer(1024) * (Symbol('d'))**(Integer(4)))))))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) + (Integer(-1) * ((((Integer(110) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(93) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(15) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(64) * (Symbol('d'))**(Integer(2)))))))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(240) * Symbol('d') * Symbol('f')))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(40) * Symbol('d') * Symbol('f')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(7) * (Integer(2))**(Integer(-1))))) * ((Integer(5) * Symbol('d') * Symbol('f')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_755():
    f = (a + b*sin(e + f*x))**(sympy.S(5)/2)*(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = (((Integer(192) * Symbol('b') * (Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1)) * (sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('c') + (Integer(-1) * Symbol('d'))) * sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(337) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(15) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3))) + (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * ((Integer(57) * (Symbol('c'))**(Integer(2))) + (Integer(284) * (Symbol('d'))**(Integer(2))))) + (Integer(-1) * ((Symbol('b'))**(Integer(3)) * ((Integer(9) * (Symbol('c'))**(Integer(3))) + (Integer(-1) * (Integer(156) * Symbol('c') * (Symbol('d'))**(Integer(2)))))))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) + (((Integer(64) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('d'))**(Integer(3)) * Symbol('f')))**(Integer(-1)) * (sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(60) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(3))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(4)) * (Symbol('d'))**(Integer(4)))) + (Integer(-1) * (Integer(20) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('c') * Symbol('d') * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Integer(12) * (Symbol('d'))**(Integer(2))))))) + (Integer(3) * (Symbol('b'))**(Integer(4)) * (((Symbol('c'))**(Integer(2)) + (Integer(4) * (Symbol('d'))**(Integer(2)))))**(Integer(2))) + (Integer(30) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * ((Integer(3) * (Symbol('c'))**(Integer(2))) + (Integer(4) * (Symbol('d'))**(Integer(2)))))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) + (Integer(-1) * ((((Integer(337) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(15) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3))) + (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * ((Integer(57) * (Symbol('c'))**(Integer(2))) + (Integer(284) * (Symbol('d'))**(Integer(2))))) + (Integer(-1) * ((Symbol('b'))**(Integer(3)) * ((Integer(9) * (Symbol('c'))**(Integer(3))) + (Integer(-1) * (Integer(156) * Symbol('c') * (Symbol('d'))**(Integer(2)))))))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(192) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(54) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(59) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(9) * (Symbol('b'))**(Integer(2)) * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Integer(4) * (Symbol('d'))**(Integer(2)))))))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(96) * Symbol('d') * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * (((Integer(192) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1)) * (((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(15) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('d'))**(Integer(2)) * ((Integer(11) * Symbol('c')) + (Integer(2) * Symbol('d'))))) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * ((Integer(51) * (Symbol('c'))**(Integer(2))) + (Integer(172) * Symbol('c') * Symbol('d')) + (Integer(212) * (Symbol('d'))**(Integer(2)))))) + ((Symbol('b'))**(Integer(3)) * ((Integer(9) * (Symbol('c'))**(Integer(3))) + (Integer(-1) * (Integer(6) * (Symbol('c'))**(Integer(2)) * Symbol('d'))) + (Integer(-1) * (Integer(156) * Symbol('c') * (Symbol('d'))**(Integer(2)))) + (Integer(-1) * (Integer(72) * (Symbol('d'))**(Integer(3))))))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))) + ((Symbol('b') * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(17) * Symbol('a') * Symbol('d')))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(24) * Symbol('d') * Symbol('f')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * Symbol('f')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_756():
    f = (a + b*sin(e + f*x))**(sympy.S(5)/2)*sqrt(c + d*sin(e + f*x))
    F = ((sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('c') + (Integer(-1) * Symbol('d'))) * sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(14) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(33) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(16) * (Symbol('d'))**(Integer(2)))))))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(24) * (Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1))) + ((sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(15) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(5) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3))) + (Integer(-1) * (Integer(5) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Integer(4) * (Symbol('d'))**(Integer(2))))))) + ((Symbol('b'))**(Integer(3)) * ((Symbol('c'))**(Integer(3)) + (Integer(4) * Symbol('c') * (Symbol('d'))**(Integer(2)))))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(8) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('d'))**(Integer(3)) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(14) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(33) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(16) * (Symbol('d'))**(Integer(2)))))))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(24) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((Symbol('b') * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(13) * Symbol('a') * Symbol('d')))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(12) * Symbol('d') * Symbol('f')))**(Integer(-1))) + ((((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(15) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(6) * Symbol('a') * Symbol('b') * Symbol('d') * ((Integer(2) * Symbol('c')) + (Integer(3) * Symbol('d')))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d'))) + (Integer(-1) * (Integer(16) * (Symbol('d'))**(Integer(2)))))))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(24) * Symbol('b') * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * Symbol('f')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_757():
    f = (a + b*sin(e + f*x))**(sympy.S(5)/2)/sqrt(c + d*sin(e + f*x))
    F = (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('c') + (Integer(-1) * Symbol('d'))) * sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(10) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('c'))**(Integer(2))) + (Integer(4) * (Symbol('d'))**(Integer(2))))))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(4) * sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('d'))**(Integer(3)) * Symbol('f')))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(2) * Symbol('d') * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d'))) + (Integer(-1) * (Integer(2) * Symbol('b') * Symbol('d')))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_758():
    f = (a + b*sin(e + f*x))**(sympy.S(5)/2)/(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Symbol('d'))**(Integer(2))))))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('d'))**(Integer(3)) * Symbol('f')))**(Integer(-1)))) + ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('d') * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(4) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Symbol('d'))**(Integer(2))))))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('d'))**(Integer(2)) * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * (Symbol('b') * ((Integer(3) * Symbol('c')) + Symbol('d'))))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('d'))**(Integer(2)) * ((Symbol('c') + Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_759():
    f = (a + b*sin(e + f*x))**(sympy.S(5)/2)/(c + d*sin(e + f*x))**(sympy.S(5)/2)
    F = ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(3) * Symbol('d') * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))) * Symbol('f') * ((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(3) * Symbol('b') * (Symbol('c'))**(Integer(2))) + (Integer(4) * Symbol('a') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(7) * Symbol('b') * (Symbol('d'))**(Integer(2))))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(3) * ((Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * ((Symbol('c') + Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * (((Integer(3) * ((Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(2)) * (Symbol('d'))**(Integer(3)) * ((Symbol('c') + Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f')))**(Integer(-1)) * (Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * ((Integer(3) * Symbol('c')) + Symbol('d'))) + (Symbol('a') * Symbol('b') * Symbol('d') * ((Integer(3) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('d'))) + (Integer(-1) * (Integer(7) * (Symbol('d'))**(Integer(2)))))) + ((Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('c'))**(Integer(3))) + (Integer(-1) * (Integer(6) * (Symbol('c'))**(Integer(2)) * Symbol('d'))) + (Integer(-1) * (Integer(2) * Symbol('c') * (Symbol('d'))**(Integer(2)))) + (Integer(9) * (Symbol('d'))**(Integer(3)))))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.Function('EllipticPi')((((Symbol('a') + Symbol('b')) * Symbol('d')) * ((Symbol('b') * (Symbol('c') + Symbol('d'))))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_760():
    f = (c + d*sin(e + f*x))**(sympy.S(5)/2)/sqrt(a + b*sin(e + f*x))
    F = ((Integer(3) * sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('c') + (Integer(-1) * Symbol('d'))) * Symbol('d') * sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(10) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(15) * (Symbol('c'))**(Integer(2))) + (Integer(4) * (Symbol('d'))**(Integer(2))))))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('d') * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(4) * Symbol('b') * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(2) * Symbol('b') * Symbol('f')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * Symbol('b') * Symbol('d') * ((Integer(7) * Symbol('c')) + (Integer(3) * Symbol('d'))))) + ((Symbol('b'))**(Integer(2)) * ((Integer(8) * (Symbol('c'))**(Integer(2))) + (Integer(9) * Symbol('c') * Symbol('d')) + (Integer(2) * (Symbol('d'))**(Integer(2)))))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_761():
    f = (c + d*sin(e + f*x))**(sympy.S(3)/2)/sqrt(a + b*sin(e + f*x))
    F = ((sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('c') + (Integer(-1) * Symbol('d'))) * Symbol('d') * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1))) + ((sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('a') * Symbol('d')) + (Integer(-1) * (Symbol('b') * ((Integer(2) * Symbol('c')) + Symbol('d'))))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_762():
    f = sqrt(c + d*sin(e + f*x))/sqrt(a + b*sin(e + f*x))
    F = (Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.Function('EllipticPi')((((Symbol('a') + Symbol('b')) * Symbol('d')) * ((Symbol('b') * (Symbol('c') + Symbol('d'))))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('b') * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_763():
    f = 1/(sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x)))
    F = 2*sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*sqrt(a + b)*(c + d*sin(e + f*x))*elliptic_f(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(f*sqrt(c + d)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_764():
    f = 1/(sqrt(a + b*sin(e + f*x))*(c + d*sin(e + f*x))**(sympy.S(3)/2))
    F = d*sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*sqrt(a + b)*(2*a - 2*b)*(c + d*sin(e + f*x))*elliptic_e(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(f*(c - d)*sqrt(c + d)*(-a*d + b*c)**2) + 2*sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*sqrt(a + b)*(c + d*sin(e + f*x))*elliptic_f(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(f*(c - d)*sqrt(c + d)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_765():
    f = 1/(sqrt(a + b*sin(e + f*x))*(c + d*sin(e + f*x))**(sympy.S(5)/2))
    F = -2*d**2*sqrt(a + b*sin(e + f*x))*cos(e + f*x)/(f*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(c**2 - d**2)*(-3*a*d + 3*b*c)) - d*sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*sqrt(a + b)*(4*a - 4*b)*(c + d*sin(e + f*x))*(2*a*c*d - b*(3*c**2 - d**2))*elliptic_e(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(3*f*(c - d)**2*(c + d)**(sympy.S(3)/2)*(-a*d + b*c)**3) - 2*sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*sqrt(a + b)*(c + d*sin(e + f*x))*(a*d*(3*c + d) - b*(3*c**2 + 3*c*d - 2*d**2))*elliptic_f(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(3*f*(c - d)**2*(c + d)**(sympy.S(3)/2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_766():
    f = (c + d*sin(e + f*x))**(sympy.S(5)/2)/(a + b*sin(e + f*x))**(sympy.S(3)/2)
    F = (((Symbol('c') + (Integer(-1) * Symbol('d'))) * sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1))) + ((Symbol('d') * sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('f')))**(Integer(-1))) + ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((((Integer(4) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Symbol('d'))**(Integer(2))))))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * Symbol('d') * (Symbol('c') + (Integer(3) * Symbol('d'))))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(2) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(6) * Symbol('c') * Symbol('d'))) + (Integer(-1) * (Symbol('d'))**(Integer(2))))))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_767():
    f = (c + d*sin(e + f*x))**(sympy.S(3)/2)/(a + b*sin(e + f*x))**(sympy.S(3)/2)
    F = ((((Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('f')))**(Integer(-1)) * Integer(2) * (Symbol('c') + (Integer(-1) * Symbol('d'))) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) + ((((Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('f')))**(Integer(-1)) * Integer(2) * Symbol('d') * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) + ((((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1)) * Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('b') * (Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('d'))))) + (Symbol('a') * Symbol('d'))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_768():
    f = sqrt(c + d*sin(e + f*x))/(a + b*sin(e + f*x))**(sympy.S(3)/2)
    F = 2*sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*sqrt(a + b)*(c - d)*(c + d*sin(e + f*x))*elliptic_f(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(f*(a - b)*sqrt(c + d)*(-a*d + b*c)) + sqrt((-a*d + b*c)*(sin(e + f*x) + 1)/((a + b*sin(e + f*x))*(c - d)))*sqrt(-(1 - sin(e + f*x))*(-a*d + b*c)/((a + b*sin(e + f*x))*(c + d)))*(a + b*sin(e + f*x))*sqrt(c + d)*(2*c - 2*d)*elliptic_e(asin(sqrt(a + b)*sqrt(c + d*sin(e + f*x))/(sqrt(a + b*sin(e + f*x))*sqrt(c + d))), (a - b)*(c + d)/((a + b)*(c - d)))*sec(e + f*x)/(f*(a - b)*sqrt(a + b)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_769():
    f = 1/((a + b*sin(e + f*x))**(sympy.S(3)/2)*sqrt(c + d*sin(e + f*x)))
    F = 2*b*sqrt((-a*d + b*c)*(sin(e + f*x) + 1)/((a + b*sin(e + f*x))*(c - d)))*sqrt(-(1 - sin(e + f*x))*(-a*d + b*c)/((a + b*sin(e + f*x))*(c + d)))*(a + b*sin(e + f*x))*(c - d)*sqrt(c + d)*elliptic_e(asin(sqrt(a + b)*sqrt(c + d*sin(e + f*x))/(sqrt(a + b*sin(e + f*x))*sqrt(c + d))), (a - b)*(c + d)/((a + b)*(c - d)))*sec(e + f*x)/(f*(a - b)*sqrt(a + b)*(-a*d + b*c)**2) + 2*sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*sqrt(a + b)*(c + d*sin(e + f*x))*elliptic_f(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(f*(a - b)*sqrt(c + d)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_770():
    f = 1/((a + b*sin(e + f*x))**(sympy.S(3)/2)*(c + d*sin(e + f*x))**(sympy.S(3)/2))
    F = 2*b**2*cos(e + f*x)/(f*sqrt(a + b*sin(e + f*x))*(a**2 - b**2)*sqrt(c + d*sin(e + f*x))*(-a*d + b*c)) + sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*(c + d*sin(e + f*x))*(-2*a*d + 2*b*(c - 2*d))*elliptic_f(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(f*sqrt(a + b)*(c - d)*sqrt(c + d)*(-a*d + b*c)**2) - sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*(c + d*sin(e + f*x))*(2*a**2*d**2 + 2*b**2*(c**2 - 2*d**2))*elliptic_e(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(f*sqrt(a + b)*(c - d)*sqrt(c + d)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_771():
    f = 1/((a + b*sin(e + f*x))**(sympy.S(3)/2)*(c + d*sin(e + f*x))**(sympy.S(5)/2))
    F = 2*b**2*cos(e + f*x)/(f*sqrt(a + b*sin(e + f*x))*(a**2 - b**2)*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(-a*d + b*c)) + 2*d*sqrt(a + b*sin(e + f*x))*(a**2*d**2 + b**2*(3*c**2 - 4*d**2))*cos(e + f*x)/(f*(3*a**2 - 3*b**2)*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(c**2 - d**2)*(-a*d + b*c)**2) + sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*(c + d*sin(e + f*x))*(2*a**2*d**2*(3*c + d) - 12*a*b*d*(c**2 - d**2) + 2*b**2*(3*c**3 - 9*c**2*d - 6*c*d**2 + 8*d**3))*elliptic_f(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(3*f*sqrt(a + b)*(c - d)**2*(c + d)**(sympy.S(3)/2)*(-a*d + b*c)**3) + sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*(c + d*sin(e + f*x))*(8*a**3*c*d**3 - 2*a**2*b*d**2*(9*c**2 - 5*d**2) - 8*a*b**2*c*d**3 - 2*b**3*(3*c**4 - 15*c**2*d**2 + 8*d**4))*elliptic_e(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(3*f*sqrt(a + b)*(c - d)**2*(c + d)**(sympy.S(3)/2)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_772():
    f = (c + d*sin(e + f*x))**(sympy.S(5)/2)/(a + b*sin(e + f*x))**(sympy.S(5)/2)
    F = ((Integer(2) * (Symbol('c') + (Integer(-1) * Symbol('d'))) * sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(4) * Symbol('a') * Symbol('b') * Symbol('c')) + (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('d')) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2)) * Symbol('d')))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(3) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f')))**(Integer(-1))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('f')))**(Integer(-1))) + ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(3) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('d')))) * Symbol('d')) + (Integer(3) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('d'))) + (Integer(-1) * (Integer(2) * (Symbol('d'))**(Integer(2)))))) + ((Symbol('b'))**(Integer(3)) * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Integer(7) * Symbol('c') * Symbol('d'))) + (Integer(9) * (Symbol('d'))**(Integer(2)))))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(3) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_773():
    f = (c + d*sin(e + f*x))**(sympy.S(3)/2)/(a + b*sin(e + f*x))**(sympy.S(5)/2)
    F = sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*(c + d*sin(e + f*x))*(2*c - 2*d)*(3*a*c - a*d + b*c - 3*b*d)*elliptic_f(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(3*f*(a - b)**2*sqrt(a + b)*sqrt(c + d)*(-a*d + b*c)) + sqrt((-a*d + b*c)*(sin(e + f*x) + 1)/((a + b*sin(e + f*x))*(c - d)))*sqrt(-(1 - sin(e + f*x))*(-a*d + b*c)/((a + b*sin(e + f*x))*(c + d)))*(a + b*sin(e + f*x))*sqrt(c + d)*(8*c - 8*d)*(a*c - b*d)*elliptic_e(asin(sqrt(a + b)*sqrt(c + d*sin(e + f*x))/(sqrt(a + b*sin(e + f*x))*sqrt(c + d))), (a - b)*(c + d)/((a + b)*(c - d)))*sec(e + f*x)/(3*f*(a - b)**2*(a + b)**(sympy.S(3)/2)*(-a*d + b*c)) + sqrt(c + d*sin(e + f*x))*(-2*a*d + 2*b*c)*cos(e + f*x)/(f*(a + b*sin(e + f*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_774():
    f = sqrt(c + d*sin(e + f*x))/(a + b*sin(e + f*x))**(sympy.S(5)/2)
    F = 2*b*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(f*(a + b*sin(e + f*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)) + sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*(6*a + 2*b)*(c - d)*(c + d*sin(e + f*x))*elliptic_f(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(3*f*(a - b)**2*sqrt(a + b)*sqrt(c + d)*(-a*d + b*c)) + sqrt((-a*d + b*c)*(sin(e + f*x) + 1)/((a + b*sin(e + f*x))*(c - d)))*sqrt(-(1 - sin(e + f*x))*(-a*d + b*c)/((a + b*sin(e + f*x))*(c + d)))*(a + b*sin(e + f*x))*sqrt(c + d)*(2*c - 2*d)*(-3*a**2*d + 4*a*b*c - b**2*d)*elliptic_e(asin(sqrt(a + b)*sqrt(c + d*sin(e + f*x))/(sqrt(a + b*sin(e + f*x))*sqrt(c + d))), (a - b)*(c + d)/((a + b)*(c - d)))*sec(e + f*x)/(3*f*(a - b)**2*(a + b)**(sympy.S(3)/2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_775():
    f = 1/((a + b*sin(e + f*x))**(sympy.S(5)/2)*sqrt(c + d*sin(e + f*x)))
    F = 2*b**2*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(f*(a + b*sin(e + f*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)*(-a*d + b*c)) + 4*b*sqrt((-a*d + b*c)*(sin(e + f*x) + 1)/((a + b*sin(e + f*x))*(c - d)))*sqrt(-(1 - sin(e + f*x))*(-a*d + b*c)/((a + b*sin(e + f*x))*(c + d)))*(a + b*sin(e + f*x))*(c - d)*sqrt(c + d)*(-3*a**2*d + 2*a*b*c + b**2*d)*elliptic_e(asin(sqrt(a + b)*sqrt(c + d*sin(e + f*x))/(sqrt(a + b*sin(e + f*x))*sqrt(c + d))), (a - b)*(c + d)/((a + b)*(c - d)))*sec(e + f*x)/(3*f*(a - b)**2*(a + b)**(sympy.S(3)/2)*(-a*d + b*c)**3) + sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*(c + d*sin(e + f*x))*(-6*a**2*d + 6*a*b*(c - d) + 2*b**2*(c + 2*d))*elliptic_f(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(3*f*(a - b)**2*sqrt(a + b)*sqrt(c + d)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_776():
    f = 1/((a + b*sin(e + f*x))**(sympy.S(5)/2)*(c + d*sin(e + f*x))**(sympy.S(3)/2))
    F = 8*b**2*(-2*a**2*d + a*b*c + b**2*d)*cos(e + f*x)/(3*f*sqrt(a + b*sin(e + f*x))*(a**2 - b**2)**2*sqrt(c + d*sin(e + f*x))*(-a*d + b*c)**2) + 2*b**2*cos(e + f*x)/(f*(a + b*sin(e + f*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)*sqrt(c + d*sin(e + f*x))*(-a*d + b*c)) - sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*(c + d*sin(e + f*x))*(-6*a**3*d**2 + 6*a**2*b*d*(2*c - 3*d) - 6*a*b**2*(c**2 - 2*d**2) + 2*b**3*(c**2 - 6*c*d + 8*d**2))*elliptic_f(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(3*f*sqrt(a + b)*(a**2 - b**2)*(c - d)*sqrt(c + d)*(-a*d + b*c)**3) + sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*(c + d*sin(e + f*x))*(6*a**4*d**3 + 6*a**2*b**2*d*(3*c**2 - 5*d**2) - 8*a*b**3*c*(c**2 - d**2) - 2*b**4*d*(5*c**2 - 8*d**2))*elliptic_e(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(3*f*sqrt(a + b)*(a**2 - b**2)*(c - d)*sqrt(c + d)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_777():
    f = 1/((a + b*sin(e + f*x))**(sympy.S(5)/2)*(c + d*sin(e + f*x))**(sympy.S(5)/2))
    F = 4*b**2*(-5*a**2*d + 2*a*b*c + 3*b**2*d)*cos(e + f*x)/(3*f*sqrt(a + b*sin(e + f*x))*(a**2 - b**2)**2*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(-a*d + b*c)**2) + 2*b**2*cos(e + f*x)/(f*(a + b*sin(e + f*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(-a*d + b*c)) - 2*d*sqrt(a + b*sin(e + f*x))*(a**4*d**3 + a**2*b**2*d*(11*c**2 - 13*d**2) - 4*a*b**3*c*(c**2 - d**2) - b**4*d*(7*c**2 - 8*d**2))*cos(e + f*x)/(3*f*(a**2 - b**2)**2*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(c**2 - d**2)*(-a*d + b*c)**3) - sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*(c + d*sin(e + f*x))*(2*a**4*d**3*(3*c + d) - 18*a**3*b*d**2*(c**2 - d**2) + 2*a**2*b**2*d*(9*c**3 - 18*c**2*d - 15*c*d**2 + 16*d**3) - 6*a*b**3*(c**4 - 5*c**2*d**2 + 4*d**4) + 2*b**4*(c**4 - 9*c**3*d + 16*c**2*d**2 + 12*c*d**3 - 16*d**4))*elliptic_f(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(3*f*sqrt(a + b)*(a**2 - b**2)*(c - d)**2*(c + d)**(sympy.S(3)/2)*(-a*d + b*c)**4) - sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*(c + d*sin(e + f*x))*(8*a**5*c*d**4 - 8*a**4*b*(3*c**2*d**3 - 2*d**5) - 16*a**3*b**2*c*d**4 - 8*a**2*b**3*d*(3*c**4 - 12*c**2*d**2 + 7*d**4) + 8*a*b**4*c*(c**4 - 2*c**2*d**2 + 2*d**4) + 8*b**5*d*(2*c**4 - 7*c**2*d**2 + 4*d**4))*elliptic_e(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(3*f*sqrt(a + b)*(a**2 - b**2)*(c - d)**2*(c + d)**(sympy.S(3)/2)*(-a*d + b*c)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_778():
    f = (a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m')) * ((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Symbol('n'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_779():
    f = (a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**2
    F = -d**2*(a + b*sin(e + f*x))**(m + 1)*cos(e + f*x)/(b*f*(m + 2)) + sqrt(2)*d*(a + b)*(a + b*sin(e + f*x))**m*(a*d - 2*b*c*(m + 2))*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m - 1, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, b*(1 - sin(e + f*x))/(a + b))/(b**2*f*((a + b*sin(e + f*x))/(a + b))**m*(m + 2)*sqrt(sin(e + f*x) + 1)) - sqrt(2)*(a + b*sin(e + f*x))**m*(a*d*(a*d - 2*b*c*(m + 2)) + b**2*(c**2*(m + 2) + d**2*(m + 1)))*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, b*(1 - sin(e + f*x))/(a + b))/(b**2*f*((a + b*sin(e + f*x))/(a + b))**m*(m + 2)*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_780():
    f = (a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))
    F = -sqrt(2)*d*(a + b)*(a + b*sin(e + f*x))**m*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m - 1, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, b*(1 - sin(e + f*x))/(a + b))/(b*f*((a + b*sin(e + f*x))/(a + b))**m*sqrt(sin(e + f*x) + 1)) - sqrt(2)*(a + b*sin(e + f*x))**m*(-a*d + b*c)*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, b*(1 - sin(e + f*x))/(a + b))/(b*f*((a + b*sin(e + f*x))/(a + b))**m*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_781():
    f = (a + b*sin(e + f*x))**m
    F = -sqrt(2)*(a + b*sin(e + f*x))**m*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, b*(1 - sin(e + f*x))/(a + b))/(f*((a + b*sin(e + f*x))/(a + b))**m*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_782():
    f = (a + b*sin(e + f*x))**m/(c + d*sin(e + f*x))
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m')) * ((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_783():
    f = (a + b*sin(e + f*x))**m/(c + d*sin(e + f*x))**2
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m')) * (((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_784():
    f = (a + b*sin(e + f*x))**m/(c + d*sin(e + f*x))**3
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m')) * (((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_785():
    f = (a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(sympy.S(5)/2)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m')) * ((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_786():
    f = (a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m')) * ((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_787():
    f = (a + b*sin(e + f*x))**m*sqrt(c + d*sin(e + f*x))
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m')) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_788():
    f = (a + b*sin(e + f*x))**m/sqrt(c + d*sin(e + f*x))
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m')) * (sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_789():
    f = (a + b*sin(e + f*x))**m/(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m')) * (((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_790():
    f = (a + b*sin(e + f*x))**m/(c + d*sin(e + f*x))**(sympy.S(5)/2)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m')) * (((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_791():
    f = (d*csc(e + f*x))**n*(a*sin(e + f*x) + a)**3
    F = a**3*d**4*(d*csc(e + f*x))**(n - 4)*(11 - 4*n)*cos(e + f*x)*hyper((sympy.S.Half, 2 - n/2), (3 - n/2,), sin(e + f*x)**2)/(f*(2 - n)*(4 - n)*sqrt(cos(e + f*x)**2)) + a**3*d**3*(d*csc(e + f*x))**(n - 3)*(1 - 2*n)*cot(e + f*x)/(f*(1 - n)*(2 - n)) + a**3*d**3*(d*csc(e + f*x))**(n - 3)*(5 - 4*n)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S(3)/2 - n/2), (sympy.S(5)/2 - n/2,), sin(e + f*x)**2)/(f*(1 - n)*(3 - n)*sqrt(cos(e + f*x)**2)) + d**3*(d*csc(e + f*x))**(n - 3)*(a**3*csc(e + f*x) + a**3)*cot(e + f*x)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_792():
    f = (d*csc(e + f*x))**n*(a*sin(e + f*x) + a)**2
    F = a**2*d**3*(d*csc(e + f*x))**(n - 3)*(3 - 2*n)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S(3)/2 - n/2), (sympy.S(5)/2 - n/2,), sin(e + f*x)**2)/(f*(1 - n)*(3 - n)*sqrt(cos(e + f*x)**2)) + 2*a**2*d**2*(d*csc(e + f*x))**(n - 2)*cos(e + f*x)*hyper((sympy.S.Half, 1 - n/2), (2 - n/2,), sin(e + f*x)**2)/(f*(2 - n)*sqrt(cos(e + f*x)**2)) + a**2*d**2*(d*csc(e + f*x))**(n - 2)*cot(e + f*x)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_793():
    f = (d*csc(e + f*x))**n*(a*sin(e + f*x) + a)
    F = a*d**2*(d*csc(e + f*x))**(n - 2)*cos(e + f*x)*hyper((sympy.S.Half, 1 - n/2), (2 - n/2,), sin(e + f*x)**2)/(f*(2 - n)*sqrt(cos(e + f*x)**2)) + a*d*(d*csc(e + f*x))**(n - 1)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), sin(e + f*x)**2)/(f*(1 - n)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_794():
    f = (d*csc(e + f*x))**n/(a*sin(e + f*x) + a)
    F = -(d*csc(e + f*x))**n*cot(e + f*x)/(f*(a*csc(e + f*x) + a)) + d*n*(d*csc(e + f*x))**(n - 1)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), sin(e + f*x)**2)/(a*f*(1 - n)*sqrt(cos(e + f*x)**2)) + (d*csc(e + f*x))**n*cos(e + f*x)*hyper((sympy.S.Half, -n/2), (1 - n/2,), sin(e + f*x)**2)/(a*f*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_795():
    f = (d*csc(e + f*x))**n/(a*sin(e + f*x) + a)**2
    F = (d*csc(e + f*x))**(n + 2)*cot(e + f*x)/(3*d**2*f*(a*csc(e + f*x) + a)**2) - (d*csc(e + f*x))**(n + 1)*(2*n + 1)*cos(e + f*x)*hyper((sympy.S.Half, -n/2 + sympy.S(-1)/2), (sympy.S.Half - n/2,), sin(e + f*x)**2)/(3*a**2*d*f*sqrt(cos(e + f*x)**2)) + 2*n*(d*csc(e + f*x))**(n + 2)*cos(e + f*x)*hyper((sympy.S.Half, -n/2 - 1), (-n/2,), sin(e + f*x)**2)/(3*a**2*d**2*f*sqrt(cos(e + f*x)**2)) - 2*n*(d*csc(e + f*x))**(n + 2)*cot(e + f*x)/(3*a**2*d**2*f*(csc(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_796():
    f = (c*(d*sin(e + f*x))**p)**n*(a*sin(e + f*x) + a)**m
    F = -2**(m + sympy.S.Half)*(c*(d*sin(e + f*x))**p)**n*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(-m + sympy.S(-1)/2)*cos(e + f*x)*appellf1(sympy.S.Half, -n*p, sympy.S.Half - m, sympy.S(3)/2, 1 - sin(e + f*x), sympy.S.Half - sin(e + f*x)/2)/(f*sin(e + f*x)**(n*p))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_797():
    f = (c*(d*sin(e + f*x))**p)**n*(a*sin(e + f*x) + a)**3
    F = -a**3*(c*(d*sin(e + f*x))**p)**n*(2*n*p + 7)*sin(e + f*x)*cos(e + f*x)/(f*(n*p + 2)*(n*p + 3)) + a**3*(c*(d*sin(e + f*x))**p)**n*(4*n*p + 11)*sin(e + f*x)**2*cos(e + f*x)*hyper((sympy.S.Half, n*p/2 + 1), (n*p/2 + 2,), sin(e + f*x)**2)/(f*(n*p + 2)*(n*p + 3)*sqrt(cos(e + f*x)**2)) + a**3*(c*(d*sin(e + f*x))**p)**n*(4*n*p + 5)*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(f*(n*p + 1)*(n*p + 2)*sqrt(cos(e + f*x)**2)) - (c*(d*sin(e + f*x))**p)**n*(a**3*sin(e + f*x) + a**3)*sin(e + f*x)*cos(e + f*x)/(f*(n*p + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_798():
    f = (c*(d*sin(e + f*x))**p)**n*(a*sin(e + f*x) + a)**2
    F = -a**2*(c*(d*sin(e + f*x))**p)**n*sin(e + f*x)*cos(e + f*x)/(f*(n*p + 2)) + 2*a**2*(c*(d*sin(e + f*x))**p)**n*sin(e + f*x)**2*cos(e + f*x)*hyper((sympy.S.Half, n*p/2 + 1), (n*p/2 + 2,), sin(e + f*x)**2)/(f*(n*p + 2)*sqrt(cos(e + f*x)**2)) + a**2*(c*(d*sin(e + f*x))**p)**n*(2*n*p + 3)*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(f*(n*p + 1)*(n*p + 2)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_799():
    f = (c*(d*sin(e + f*x))**p)**n*(a*sin(e + f*x) + a)
    F = a*(c*(d*sin(e + f*x))**p)**n*sin(e + f*x)**2*cos(e + f*x)*hyper((sympy.S.Half, n*p/2 + 1), (n*p/2 + 2,), sin(e + f*x)**2)/(f*(n*p + 2)*sqrt(cos(e + f*x)**2)) + a*(c*(d*sin(e + f*x))**p)**n*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(f*(n*p + 1)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_800():
    f = (c*(d*sin(e + f*x))**p)**n/(a*sin(e + f*x) + a)
    F = -(c*(d*sin(e + f*x))**p)**n*cos(e + f*x)/(f*(a*sin(e + f*x) + a)) - n*p*(c*(d*sin(e + f*x))**p)**n*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(a*f*(n*p + 1)*sqrt(cos(e + f*x)**2)) + (c*(d*sin(e + f*x))**p)**n*cos(e + f*x)*hyper((sympy.S.Half, n*p/2), (n*p/2 + 1,), sin(e + f*x)**2)/(a*f*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_801():
    f = (c*(d*sin(e + f*x))**p)**n/(a*sin(e + f*x) + a)**2
    F = (c*(d*sin(e + f*x))**p)**n*sin(e + f*x)*cos(e + f*x)/(3*f*(a*sin(e + f*x) + a)**2) - n*p*(c*(d*sin(e + f*x))**p)**n*(-2*n*p + 1)*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(3*a**2*f*(n*p + 1)*sqrt(cos(e + f*x)**2)) + (c*(d*sin(e + f*x))**p)**n*(-2*n*p + 2)*sin(e + f*x)*cos(e + f*x)/(3*a**2*f*(sin(e + f*x) + 1)) + (c*(d*sin(e + f*x))**p)**n*(-2*n**2*p**2 + 2)*sin(e + f*x)**2*cos(e + f*x)*hyper((sympy.S.Half, n*p/2 + 1), (n*p/2 + 2,), sin(e + f*x)**2)/(3*a**2*f*(n*p + 2)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_802():
    f = (d*csc(e + f*x))**n*(a + b*sin(e + f*x))**3
    F = a**2*b*d**3*(d*csc(e + f*x))**(n - 3)*(1 - 2*n)*cot(e + f*x)/(f*(1 - n)*(2 - n)) + a**2*d**3*(d*csc(e + f*x))**(n - 3)*(a*csc(e + f*x) + b)*cot(e + f*x)/(f*(1 - n)) + a*d**3*(d*csc(e + f*x))**(n - 3)*(a**2*(2 - n) + 3*b**2*(1 - n))*cos(e + f*x)*hyper((sympy.S.Half, sympy.S(3)/2 - n/2), (sympy.S(5)/2 - n/2,), sin(e + f*x)**2)/(f*(1 - n)*(3 - n)*sqrt(cos(e + f*x)**2)) + b*d**4*(d*csc(e + f*x))**(n - 4)*(3*a**2*(3 - n) + b**2*(2 - n))*cos(e + f*x)*hyper((sympy.S.Half, 2 - n/2), (3 - n/2,), sin(e + f*x)**2)/(f*(2 - n)*(4 - n)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_803():
    f = (d*csc(e + f*x))**n*(a + b*sin(e + f*x))**2
    F = a**2*d**2*(d*csc(e + f*x))**(n - 2)*cot(e + f*x)/(f*(1 - n)) + 2*a*b*d**2*(d*csc(e + f*x))**(n - 2)*cos(e + f*x)*hyper((sympy.S.Half, 1 - n/2), (2 - n/2,), sin(e + f*x)**2)/(f*(2 - n)*sqrt(cos(e + f*x)**2)) + d**3*(d*csc(e + f*x))**(n - 3)*(a**2*(2 - n) + b**2*(1 - n))*cos(e + f*x)*hyper((sympy.S.Half, sympy.S(3)/2 - n/2), (sympy.S(5)/2 - n/2,), sin(e + f*x)**2)/(f*(1 - n)*(3 - n)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_804():
    f = (d*csc(e + f*x))**n*(a + b*sin(e + f*x))
    F = a*d*(d*csc(e + f*x))**(n - 1)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), sin(e + f*x)**2)/(f*(1 - n)*sqrt(cos(e + f*x)**2)) + b*d**2*(d*csc(e + f*x))**(n - 2)*cos(e + f*x)*hyper((sympy.S.Half, 1 - n/2), (2 - n/2,), sin(e + f*x)**2)/(f*(2 - n)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_805():
    f = (d*csc(e + f*x))**n/(a + b*sin(e + f*x))
    F = -a*(d*csc(e + f*x))**(n + 1)*(sin(e + f*x)**2)**(n/2 + sympy.S.Half)*cos(e + f*x)*appellf1(sympy.S.Half, 1, n/2 + sympy.S.Half, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(d*f*(a**2 - b**2)) + b*(d*csc(e + f*x))**(n + 1)*(sin(e + f*x)**2)**(n/2)*sin(e + f*x)*cos(e + f*x)*appellf1(sympy.S.Half, 1, n/2, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(d*f*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_806():
    f = (d*csc(e + f*x))**n/(a + b*sin(e + f*x))**2
    F = -a**2*(d*csc(e + f*x))**(n + 2)*(sin(e + f*x)**2)**(n/2 + sympy.S.Half)*sin(e + f*x)*cos(e + f*x)*appellf1(sympy.S.Half, 2, n/2 + sympy.S.Half, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(d**2*f*(a**2 - b**2)**2) + 2*a*b*(d*csc(e + f*x))**(n + 2)*(sin(e + f*x)**2)**(n/2 + 1)*cos(e + f*x)*appellf1(sympy.S.Half, 2, n/2, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(d**2*f*(a**2 - b**2)**2) - b**2*(d*csc(e + f*x))**(n + 2)*(sin(e + f*x)**2)**(n/2 + sympy.S(-1)/2)*sin(e + f*x)**3*cos(e + f*x)*appellf1(sympy.S.Half, 2, n/2 + sympy.S(-1)/2, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(d**2*f*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_807():
    f = (d*csc(e + f*x))**n/(a + b*sin(e + f*x))**3
    F = -a**3*(d*csc(e + f*x))**(n + 3)*(sin(e + f*x)**2)**(n/2 + sympy.S(3)/2)*cos(e + f*x)*appellf1(sympy.S.Half, 3, n/2 + sympy.S.Half, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(d**3*f*(a**2 - b**2)**3) + 3*a**2*b*(d*csc(e + f*x))**(n + 3)*(sin(e + f*x)**2)**(n/2)*sin(e + f*x)**3*cos(e + f*x)*appellf1(sympy.S.Half, 3, n/2, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(d**3*f*(a**2 - b**2)**3) - 3*a*b**2*(d*csc(e + f*x))**(n + 3)*(sin(e + f*x)**2)**(n/2 + sympy.S(-1)/2)*sin(e + f*x)**4*cos(e + f*x)*appellf1(sympy.S.Half, 3, n/2 + sympy.S(-1)/2, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(d**3*f*(a**2 - b**2)**3) + b**3*(d*csc(e + f*x))**(n + 3)*(sin(e + f*x)**2)**(n/2)*sin(e + f*x)**3*cos(e + f*x)*appellf1(sympy.S.Half, 3, n/2 - 1, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(d**3*f*(a**2 - b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_808():
    f = (c*(d*sin(e + f*x))**p)**n*(a + b*sin(e + f*x))**m
    F = (((Symbol('c') * ((Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))**(Symbol('p'))))**(Symbol('n')) * sympy.Function('Unintegrable')((((Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))**((Symbol('n') * Symbol('p'))) * ((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m'))), x)) * (((Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))**((Symbol('n') * Symbol('p'))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_809():
    f = (c*(d*sin(e + f*x))**p)**n*(a + b*sin(e + f*x))**3
    F = -a*b**2*(c*(d*sin(e + f*x))**p)**n*(2*n*p + 7)*sin(e + f*x)*cos(e + f*x)/(f*(n*p + 2)*(n*p + 3)) + a*(c*(d*sin(e + f*x))**p)**n*(a**2*(n*p + 2) + 3*b**2*(n*p + 1))*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(f*(n*p + 1)*(n*p + 2)*sqrt(cos(e + f*x)**2)) - b**2*(c*(d*sin(e + f*x))**p)**n*(a + b*sin(e + f*x))*sin(e + f*x)*cos(e + f*x)/(f*(n*p + 3)) + b*(c*(d*sin(e + f*x))**p)**n*(3*a**2*(n*p + 3) + b**2*(n*p + 2))*sin(e + f*x)**2*cos(e + f*x)*hyper((sympy.S.Half, n*p/2 + 1), (n*p/2 + 2,), sin(e + f*x)**2)/(f*(n*p + 2)*(n*p + 3)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_810():
    f = (c*(d*sin(e + f*x))**p)**n*(a + b*sin(e + f*x))**2
    F = 2*a*b*(c*(d*sin(e + f*x))**p)**n*sin(e + f*x)**2*cos(e + f*x)*hyper((sympy.S.Half, n*p/2 + 1), (n*p/2 + 2,), sin(e + f*x)**2)/(f*(n*p + 2)*sqrt(cos(e + f*x)**2)) - b**2*(c*(d*sin(e + f*x))**p)**n*sin(e + f*x)*cos(e + f*x)/(f*(n*p + 2)) + (c*(d*sin(e + f*x))**p)**n*(a**2*(n*p + 2) + b**2*(n*p + 1))*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(f*(n*p + 1)*(n*p + 2)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_811():
    f = (c*(d*sin(e + f*x))**p)**n*(a + b*sin(e + f*x))
    F = a*(c*(d*sin(e + f*x))**p)**n*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(f*(n*p + 1)*sqrt(cos(e + f*x)**2)) + b*(c*(d*sin(e + f*x))**p)**n*sin(e + f*x)**2*cos(e + f*x)*hyper((sympy.S.Half, n*p/2 + 1), (n*p/2 + 2,), sin(e + f*x)**2)/(f*(n*p + 2)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_812():
    f = (c*(d*sin(e + f*x))**p)**n/(a + b*sin(e + f*x))
    F = -a*(c*(d*sin(e + f*x))**p)**n*(sin(e + f*x)**2)**(-n*p/2 + sympy.S.Half)*cot(e + f*x)*appellf1(sympy.S.Half, 1, -n*p/2 + sympy.S.Half, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(f*(a**2 - b**2)) + b*(c*(d*sin(e + f*x))**p)**n*cos(e + f*x)*appellf1(sympy.S.Half, 1, -n*p/2, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(f*(a**2 - b**2)*(sin(e + f*x)**2)**(n*p/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_813():
    f = (c*(d*sin(e + f*x))**p)**n/(a + b*sin(e + f*x))**2
    F = -a**2*(c*(d*sin(e + f*x))**p)**n*(sin(e + f*x)**2)**(-n*p/2 + sympy.S.Half)*cot(e + f*x)*appellf1(sympy.S.Half, 2, -n*p/2 + sympy.S.Half, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(f*(a**2 - b**2)**2) + 2*a*b*(c*(d*sin(e + f*x))**p)**n*cos(e + f*x)*appellf1(sympy.S.Half, 2, -n*p/2, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(f*(a**2 - b**2)**2*(sin(e + f*x)**2)**(n*p/2)) - b**2*(c*(d*sin(e + f*x))**p)**n*(sin(e + f*x)**2)**(-n*p/2 + sympy.S(-1)/2)*sin(e + f*x)*cos(e + f*x)*appellf1(sympy.S.Half, 2, -n*p/2 + sympy.S(-1)/2, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(f*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_2_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_814():
    f = (c*(d*sin(e + f*x))**p)**n/(a + b*sin(e + f*x))**3
    F = -a**3*(c*(d*sin(e + f*x))**p)**n*(sin(e + f*x)**2)**(-n*p/2 + sympy.S.Half)*cot(e + f*x)*appellf1(sympy.S.Half, 3, -n*p/2 + sympy.S.Half, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(f*(a**2 - b**2)**3) + 3*a**2*b*(c*(d*sin(e + f*x))**p)**n*cos(e + f*x)*appellf1(sympy.S.Half, 3, -n*p/2, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(f*(a**2 - b**2)**3*(sin(e + f*x)**2)**(n*p/2)) - 3*a*b**2*(c*(d*sin(e + f*x))**p)**n*(sin(e + f*x)**2)**(-n*p/2 + sympy.S(-1)/2)*sin(e + f*x)*cos(e + f*x)*appellf1(sympy.S.Half, 3, -n*p/2 + sympy.S(-1)/2, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(f*(a**2 - b**2)**3) + b**3*(c*(d*sin(e + f*x))**p)**n*cos(e + f*x)*appellf1(sympy.S.Half, 3, -n*p/2 - 1, sympy.S(3)/2, -b**2*cos(e + f*x)**2/(a**2 - b**2), cos(e + f*x)**2)/(f*(a**2 - b**2)**3*(sin(e + f*x)**2)**(n*p/2))
    assert integrate(f, x) == F

