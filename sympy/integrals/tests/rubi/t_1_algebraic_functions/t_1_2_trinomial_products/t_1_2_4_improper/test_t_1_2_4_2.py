"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.2 Trinomial products/1.2.4 Improper/1.2.4.2 (d x)^m (a x^q+b x^n+c x^(2 n-q))^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, m, n, q = symbols('a b c d e m n q')

def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_1():
    f = x**2*(a*x**2 + b*x**3 + c*x**4)
    F = a*x**5/5 + b*x**6/6 + c*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_2():
    f = x*(a*x**2 + b*x**3 + c*x**4)
    F = a*x**4/4 + b*x**5/5 + c*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_3():
    f = a*x**2 + b*x**3 + c*x**4
    F = a*x**3/3 + b*x**4/4 + c*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_4():
    f = (a*x**2 + b*x**3 + c*x**4)/x
    F = a*x**2/2 + b*x**3/3 + c*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_5():
    f = (a*x**2 + b*x**3 + c*x**4)/x**2
    F = a*x + b*x**2/2 + c*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_6():
    f = x**2*(a*x**2 + b*x**3 + c*x**4)**2
    F = a**2*x**7/7 + a*b*x**8/4 + b*c*x**10/5 + c**2*x**11/11 + x**9*(2*a*c + b**2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_7():
    f = x*(a*x**2 + b*x**3 + c*x**4)**2
    F = a**2*x**6/6 + 2*a*b*x**7/7 + 2*b*c*x**9/9 + c**2*x**10/10 + x**8*(2*a*c + b**2)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_8():
    f = (a*x**2 + b*x**3 + c*x**4)**2
    F = a**2*x**5/5 + a*b*x**6/3 + b*c*x**8/4 + c**2*x**9/9 + x**7*(2*a*c + b**2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_9():
    f = (a*x**2 + b*x**3 + c*x**4)**2/x
    F = a**2*x**4/4 + 2*a*b*x**5/5 + 2*b*c*x**7/7 + c**2*x**8/8 + x**6*(2*a*c + b**2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_10():
    f = (a*x**2 + b*x**3 + c*x**4)**2/x**2
    F = a**2*x**3/3 + a*b*x**4/2 + b*c*x**6/3 + c**2*x**7/7 + x**5*(2*a*c + b**2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_11():
    f = x**5/(a*x**2 + b*x**3 + c*x**4)
    F = -b*x/c**2 + b*(-3*a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**3*sqrt(-4*a*c + b**2)) + x**2/(2*c) + (-a*c + b**2)*log(a + b*x + c*x**2)/(2*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_12():
    f = x**4/(a*x**2 + b*x**3 + c*x**4)
    F = -b*log(a + b*x + c*x**2)/(2*c**2) + x/c - (-2*a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_13():
    f = x**3/(a*x**2 + b*x**3 + c*x**4)
    F = b*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c*sqrt(-4*a*c + b**2)) + log(a + b*x + c*x**2)/(2*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_14():
    f = x**2/(a*x**2 + b*x**3 + c*x**4)
    F = -2*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_15():
    f = x/(a*x**2 + b*x**3 + c*x**4)
    F = b*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(a*sqrt(-4*a*c + b**2)) + log(x)/a - log(a + b*x + c*x**2)/(2*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_16():
    f = 1/(a*x**2 + b*x**3 + c*x**4)
    F = -1/(a*x) - b*log(x)/a**2 + b*log(a + b*x + c*x**2)/(2*a**2) - (-2*a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(a**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_17():
    f = 1/(x*(a*x**2 + b*x**3 + c*x**4))
    F = -1/(2*a*x**2) + b/(a**2*x) + b*(-3*a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(a**3*sqrt(-4*a*c + b**2)) + (-a*c + b**2)*log(x)/a**3 - (-a*c + b**2)*log(a + b*x + c*x**2)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_18():
    f = 1/(x**2*(a*x**2 + b*x**3 + c*x**4))
    F = -1/(3*a*x**3) + b/(2*a**2*x**2) - (-a*c + b**2)/(a**3*x) - b*(-2*a*c + b**2)*log(x)/a**4 + b*(-2*a*c + b**2)*log(a + b*x + c*x**2)/(2*a**4) - (2*a**2*c**2 - 4*a*b**2*c + b**4)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(a**4*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_19():
    f = x**8/(a*x**2 + b*x**3 + c*x**4)**2
    F = -b*x**2/(c*(-4*a*c + b**2)) - b*log(a + b*x + c*x**2)/c**3 + x**3*(2*a + b*x)/((-4*a*c + b**2)*(a + b*x + c*x**2)) + x*(-6*a*c + 2*b**2)/(c**2*(-4*a*c + b**2)) - (12*a**2*c**2 - 12*a*b**2*c + 2*b**4)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**3*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_20():
    f = x**7/(a*x**2 + b*x**3 + c*x**4)**2
    F = -b*x/(c*(-4*a*c + b**2)) + b*(-6*a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**2*(-4*a*c + b**2)**(sympy.S(3)/2)) + x**2*(2*a + b*x)/((-4*a*c + b**2)*(a + b*x + c*x**2)) + log(a + b*x + c*x**2)/(2*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_21():
    f = x**6/(a*x**2 + b*x**3 + c*x**4)**2
    F = 4*a*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) + x*(2*a + b*x)/((-4*a*c + b**2)*(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_22():
    f = x**5/(a*x**2 + b*x**3 + c*x**4)**2
    F = -2*b*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) + (2*a + b*x)/((-4*a*c + b**2)*(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_23():
    f = x**4/(a*x**2 + b*x**3 + c*x**4)**2
    F = 4*c*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) - (b + 2*c*x)/((-4*a*c + b**2)*(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_24():
    f = x**3/(a*x**2 + b*x**3 + c*x**4)**2
    F = (-2*a*c + b**2 + b*c*x)/(a*(-4*a*c + b**2)*(a + b*x + c*x**2)) + b*(-6*a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(a**2*(-4*a*c + b**2)**(sympy.S(3)/2)) + log(x)/a**2 - log(a + b*x + c*x**2)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_25():
    f = x**2/(a*x**2 + b*x**3 + c*x**4)**2
    F = (-2*a*c + b**2 + b*c*x)/(a*x*(-4*a*c + b**2)*(a + b*x + c*x**2)) + (6*a*c - 2*b**2)/(a**2*x*(-4*a*c + b**2)) - 2*b*log(x)/a**3 + b*log(a + b*x + c*x**2)/a**3 - (12*a**2*c**2 - 12*a*b**2*c + 2*b**4)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(a**3*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_26():
    f = x/(a*x**2 + b*x**3 + c*x**4)**2
    F = (-2*a*c + b**2 + b*c*x)/(a*x**2*(-4*a*c + b**2)*(a + b*x + c*x**2)) - (-8*a*c + 3*b**2)/(2*a**2*x**2*(-4*a*c + b**2)) + b*(-11*a*c + 3*b**2)/(a**3*x*(-4*a*c + b**2)) + b*(30*a**2*c**2 - 20*a*b**2*c + 3*b**4)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(a**4*(-4*a*c + b**2)**(sympy.S(3)/2)) + (-2*a*c + 3*b**2)*log(x)/a**4 - (-2*a*c + 3*b**2)*log(a + b*x + c*x**2)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_27():
    f = (a*x**2 + b*x**3 + c*x**4)**(-2)
    F = (-2*a*c + b**2 + b*c*x)/(a*x**3*(-4*a*c + b**2)*(a + b*x + c*x**2)) - (-10*a*c + 4*b**2)/(3*a**2*x**3*(-4*a*c + b**2)) + b*(-7*a*c + 2*b**2)/(a**3*x**2*(-4*a*c + b**2)) - (10*a**2*c**2 - 18*a*b**2*c + 4*b**4)/(a**4*x*(-4*a*c + b**2)) - 2*b*(-3*a*c + 2*b**2)*log(x)/a**5 + b*(-3*a*c + 2*b**2)*log(a + b*x + c*x**2)/a**5 - (-20*a**3*c**3 + 60*a**2*b**2*c**2 - 30*a*b**4*c + 4*b**6)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(a**5*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_28():
    f = 1/(x*(a*x**2 + b*x**3 + c*x**4)**2)
    F = (-2*a*c + b**2 + b*c*x)/(a*x**4*(-4*a*c + b**2)*(a + b*x + c*x**2)) - (-12*a*c + 5*b**2)/(4*a**2*x**4*(-4*a*c + b**2)) + b*(-17*a*c + 5*b**2)/(3*a**3*x**3*(-4*a*c + b**2)) - (12*a**2*c**2 - 22*a*b**2*c + 5*b**4)/(2*a**4*x**2*(-4*a*c + b**2)) + b*(29*a**2*c**2 - 27*a*b**2*c + 5*b**4)/(a**5*x*(-4*a*c + b**2)) + b*(-70*a**3*c**3 + 105*a**2*b**2*c**2 - 42*a*b**4*c + 5*b**6)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(a**6*(-4*a*c + b**2)**(sympy.S(3)/2)) + (3*a**2*c**2 - 12*a*b**2*c + 5*b**4)*log(x)/a**6 - (3*a**2*c**2 - 12*a*b**2*c + 5*b**4)*log(a + b*x + c*x**2)/(2*a**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_29():
    f = x**2*sqrt(a*x**2 + b*x**3 + c*x**4)
    F = b*(-116*a*c + 35*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(960*c**3) + b*x*(-12*a*c + 7*b**2)*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(256*c**(sympy.S(9)/2)*sqrt(a*x**2 + b*x**3 + c*x**4)) + x**2*(b + 8*c*x)*sqrt(a*x**2 + b*x**3 + c*x**4)/(40*c) - x*(-16*a*c + 7*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(240*c**2) - sqrt(a*x**2 + b*x**3 + c*x**4)*(256*a**2*c**2 - 460*a*b**2*c + 105*b**4)/(1920*c**4*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_30():
    f = x*sqrt(a*x**2 + b*x**3 + c*x**4)
    F = b*(-52*a*c + 15*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(192*c**3*x) + x*(b + 6*c*x)*sqrt(a*x**2 + b*x**3 + c*x**4)/(24*c) - (-12*a*c + 5*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(96*c**2) - x*(-4*a*c + b**2)*(-4*a*c + 5*b**2)*sqrt(a + b*x + c*x**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(128*c**(sympy.S(7)/2)*sqrt(a*x**2 + b*x**3 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_31():
    f = sqrt(a*x**2 + b*x**3 + c*x**4)
    F = -b*(b + 2*c*x)*sqrt(a*x**2 + b*x**3 + c*x**4)/(8*c**2*x) + b*(-4*a*c + b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(5)/2)*x*sqrt(a + b*x + c*x**2)) + (a + b*x + c*x**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(3*c*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_32():
    f = sqrt(a*x**2 + b*x**3 + c*x**4)/x
    F = (b + 2*c*x)*sqrt(a*x**2 + b*x**3 + c*x**4)/(4*c*x) - x*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*c**(sympy.S(3)/2)*sqrt(a*x**2 + b*x**3 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_33():
    f = sqrt(a*x**2 + b*x**3 + c*x**4)/x**2
    F = -sqrt(a)*x*sqrt(a + b*x + c*x**2)*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/sqrt(a*x**2 + b*x**3 + c*x**4) + b*x*sqrt(a + b*x + c*x**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*sqrt(c)*sqrt(a*x**2 + b*x**3 + c*x**4)) + sqrt(a*x**2 + b*x**3 + c*x**4)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_34():
    f = sqrt(a*x**2 + b*x**3 + c*x**4)/x**3
    F = sqrt(c)*x*sqrt(a + b*x + c*x**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/sqrt(a*x**2 + b*x**3 + c*x**4) - sqrt(a*x**2 + b*x**3 + c*x**4)/x**2 - b*x*sqrt(a + b*x + c*x**2)*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(2*sqrt(a)*sqrt(a*x**2 + b*x**3 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_35():
    f = sqrt(a*x**2 + b*x**3 + c*x**4)/x**4
    F = -sqrt(a*x**2 + b*x**3 + c*x**4)/(2*x**3) - b*sqrt(a*x**2 + b*x**3 + c*x**4)/(4*a*x**2) + (-4*a*c + b**2)*atanh(x*(2*a + b*x)/(2*sqrt(a)*sqrt(a*x**2 + b*x**3 + c*x**4)))/(8*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_36():
    f = sqrt(a*x**2 + b*x**3 + c*x**4)/x**5
    F = -sqrt(a*x**2 + b*x**3 + c*x**4)/(3*x**4) - b*sqrt(a*x**2 + b*x**3 + c*x**4)/(12*a*x**3) + (-8*a*c + 3*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(24*a**2*x**2) - b*(-4*a*c + b**2)*atanh(x*(2*a + b*x)/(2*sqrt(a)*sqrt(a*x**2 + b*x**3 + c*x**4)))/(16*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_37():
    f = sqrt(a*x**2 + b*x**3 + c*x**4)/x**6
    F = -sqrt(a*x**2 + b*x**3 + c*x**4)/(4*x**5) - b*sqrt(a*x**2 + b*x**3 + c*x**4)/(24*a*x**4) + (-12*a*c + 5*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(96*a**2*x**3) - b*(-52*a*c + 15*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(192*a**3*x**2) + (-4*a*c + b**2)*(-4*a*c + 5*b**2)*atanh(x*(2*a + b*x)/(2*sqrt(a)*sqrt(a*x**2 + b*x**3 + c*x**4)))/(128*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_38():
    f = x*(a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)
    F = -b*x*sqrt(a*x**2 + b*x**3 + c*x**4)*(2416*a**2*c**2 - 1560*a*b**2*c + 231*b**4)/(71680*c**4) - b*sqrt(a*x**2 + b*x**3 + c*x**4)*(-58816*a**3*c**3 + 81648*a**2*b**2*c**2 - 30660*a*b**4*c + 3465*b**6)/(573440*c**6*x) + x*(3*b + 14*c*x)*(a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)/(112*c) - x**3*(b*(68*a*c + 11*b**2) + 10*c*x*(-28*a*c + 11*b**2))*sqrt(a*x**2 + b*x**3 + c*x**4)/(4480*c**2) + x**2*sqrt(a*x**2 + b*x**3 + c*x**4)*(560*a**2*c**2 - 568*a*b**2*c + 99*b**4)/(35840*c**3) + sqrt(a*x**2 + b*x**3 + c*x**4)*(-6720*a**3*c**3 + 18896*a**2*b**2*c**2 - 8988*a*b**4*c + 1155*b**6)/(286720*c**5) + 3*x*(-4*a*c + b**2)**2*sqrt(a + b*x + c*x**2)*(16*a**2*c**2 - 72*a*b**2*c + 33*b**4)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(32768*c**(sympy.S(13)/2)*sqrt(a*x**2 + b*x**3 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_39():
    f = (a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)
    F = -b*x**2*(-44*a*c + 9*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(2240*c**2) - b*sqrt(a*x**2 + b*x**3 + c*x**4)*(1168*a**2*c**2 - 728*a*b**2*c + 105*b**4)/(17920*c**4) - 3*b*x*(-4*a*c + b**2)**2*(-4*a*c + 3*b**2)*sqrt(a + b*x + c*x**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2048*c**(sympy.S(11)/2)*sqrt(a*x**2 + b*x**3 + c*x**4)) + x*(a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)/7 + x**3*(24*a*c + b**2 + 10*b*c*x)*sqrt(a*x**2 + b*x**3 + c*x**4)/(280*c) + x*(-32*a*c + 7*b**2)*(-4*a*c + 3*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(4480*c**3) + sqrt(a*x**2 + b*x**3 + c*x**4)*(-2048*a**3*c**3 + 5488*a**2*b**2*c**2 - 2520*a*b**4*c + 315*b**6)/(35840*c**5*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_40():
    f = (a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)/x
    F = -b*sqrt(a*x**2 + b*x**3 + c*x**4)*(1296*a**2*c**2 - 760*a*b**2*c + 105*b**4)/(7680*c**4*x) + (3*b + 10*c*x)*(a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)/(60*c*x) - x*(b*(12*a*c + 7*b**2) + 6*c*x*(-20*a*c + 7*b**2))*sqrt(a*x**2 + b*x**3 + c*x**4)/(960*c**2) + sqrt(a*x**2 + b*x**3 + c*x**4)*(240*a**2*c**2 - 216*a*b**2*c + 35*b**4)/(3840*c**3) + x*(-4*a*c + b**2)**2*(-4*a*c + 7*b**2)*sqrt(a + b*x + c*x**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(1024*c**(sympy.S(9)/2)*sqrt(a*x**2 + b*x**3 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_41():
    f = (a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)/x**2
    F = -b*(b + 2*c*x)*(a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)/(16*c**2*x**3) + 3*b*(b + 2*c*x)*(-4*a*c + b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(128*c**3*x) - 3*b*x*(-4*a*c + b**2)**2*sqrt(a + b*x + c*x**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(256*c**(sympy.S(7)/2)*sqrt(a*x**2 + b*x**3 + c*x**4)) + (a*x**2 + b*x**3 + c*x**4)**(sympy.S(5)/2)/(5*c*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_42():
    f = (a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)/x**3
    F = (b + 2*c*x)*(a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)/(8*c*x**3) - (b + 2*c*x)*(-12*a*c + 3*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(64*c**2*x) + 3*x*(-4*a*c + b**2)**2*sqrt(a + b*x + c*x**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(128*c**(sympy.S(5)/2)*sqrt(a*x**2 + b*x**3 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_43():
    f = (a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)/x**4
    F = -a**(sympy.S(3)/2)*x*sqrt(a + b*x + c*x**2)*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/sqrt(a*x**2 + b*x**3 + c*x**4) - b*x*(-12*a*c + b**2)*sqrt(a + b*x + c*x**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(3)/2)*sqrt(a*x**2 + b*x**3 + c*x**4)) + (a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)/(3*x**3) + (8*a*c + b**2 + 2*b*c*x)*sqrt(a*x**2 + b*x**3 + c*x**4)/(8*c*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_44():
    f = (a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)/x**5
    F = -3*sqrt(a)*b*x*sqrt(a + b*x + c*x**2)*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(2*sqrt(a*x**2 + b*x**3 + c*x**4)) + (9*b + 6*c*x)*sqrt(a*x**2 + b*x**3 + c*x**4)/(4*x) - (a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)/x**4 + x*(12*a*c + 3*b**2)*sqrt(a + b*x + c*x**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*sqrt(c)*sqrt(a*x**2 + b*x**3 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_45():
    f = (a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)/x**6
    F = 3*b*sqrt(c)*x*sqrt(a + b*x + c*x**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*sqrt(a*x**2 + b*x**3 + c*x**4)) - (3*b - 6*c*x)*sqrt(a*x**2 + b*x**3 + c*x**4)/(4*x**2) - (a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)/(2*x**5) - x*(12*a*c + 3*b**2)*sqrt(a + b*x + c*x**2)*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(8*sqrt(a)*sqrt(a*x**2 + b*x**3 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_46():
    f = (a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)/x**7
    F = c**(sympy.S(3)/2)*x*sqrt(a + b*x + c*x**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/sqrt(a*x**2 + b*x**3 + c*x**4) - (a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)/(3*x**6) - b*(a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)/(4*a*x**5) + (-8*a*c + b**2 + 2*b*c*x)*sqrt(a*x**2 + b*x**3 + c*x**4)/(8*a*x**2) + b*x*(-12*a*c + b**2)*sqrt(a + b*x + c*x**2)*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(16*a**(sympy.S(3)/2)*sqrt(a*x**2 + b*x**3 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_47():
    f = (a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)/x**8
    F = -(b + 6*c*x)*sqrt(a*x**2 + b*x**3 + c*x**4)/(8*x**4) - (a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)/(4*x**7) - (-12*a*c + b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(32*a*x**3) + b*(-20*a*c + 3*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(64*a**2*x**2) - 3*(-4*a*c + b**2)**2*atanh(x*(2*a + b*x)/(2*sqrt(a)*sqrt(a*x**2 + b*x**3 + c*x**4)))/(128*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_48():
    f = (a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)/x**9
    F = -(3*b + 12*c*x)*sqrt(a*x**2 + b*x**3 + c*x**4)/(40*x**5) - (a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)/(5*x**8) - (-8*a*c + b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(80*a*x**4) + b*(-28*a*c + 5*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(320*a**2*x**3) - sqrt(a*x**2 + b*x**3 + c*x**4)*(128*a**2*c**2 - 100*a*b**2*c + 15*b**4)/(640*a**3*x**2) + 3*b*(-4*a*c + b**2)**2*atanh(x*(2*a + b*x)/(2*sqrt(a)*sqrt(a*x**2 + b*x**3 + c*x**4)))/(256*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_49():
    f = x**3/sqrt(a*x**2 + b*x**3 + c*x**4)
    F = -3*b*sqrt(a*x**2 + b*x**3 + c*x**4)/(4*c**2*x) + sqrt(a*x**2 + b*x**3 + c*x**4)/(2*c) + x*(-4*a*c + 3*b**2)*sqrt(a + b*x + c*x**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*c**(sympy.S(5)/2)*sqrt(a*x**2 + b*x**3 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_50():
    f = x**2/sqrt(a*x**2 + b*x**3 + c*x**4)
    F = -b*x*sqrt(a + b*x + c*x**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*c**(sympy.S(3)/2)*sqrt(a*x**2 + b*x**3 + c*x**4)) + sqrt(a*x**2 + b*x**3 + c*x**4)/(c*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_51():
    f = x/sqrt(a*x**2 + b*x**3 + c*x**4)
    F = x*sqrt(a + b*x + c*x**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(sqrt(c)*sqrt(a*x**2 + b*x**3 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_52():
    f = 1/sqrt(a*x**2 + b*x**3 + c*x**4)
    F = -atanh(x*(2*a + b*x)/(2*sqrt(a)*sqrt(a*x**2 + b*x**3 + c*x**4)))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_53():
    f = 1/(x*sqrt(a*x**2 + b*x**3 + c*x**4))
    F = -sqrt(a*x**2 + b*x**3 + c*x**4)/(a*x**2) + b*atanh(x*(2*a + b*x)/(2*sqrt(a)*sqrt(a*x**2 + b*x**3 + c*x**4)))/(2*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_54():
    f = 1/(x**2*sqrt(a*x**2 + b*x**3 + c*x**4))
    F = -sqrt(a*x**2 + b*x**3 + c*x**4)/(2*a*x**3) + 3*b*sqrt(a*x**2 + b*x**3 + c*x**4)/(4*a**2*x**2) - (-4*a*c + 3*b**2)*atanh(x*(2*a + b*x)/(2*sqrt(a)*sqrt(a*x**2 + b*x**3 + c*x**4)))/(8*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_55():
    f = x**7/(a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)
    F = -2*b*x*sqrt(a*x**2 + b*x**3 + c*x**4)/(c*(-4*a*c + b**2)) - b*(-52*a*c + 15*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(4*c**3*x*(-4*a*c + b**2)) + 2*x**4*(2*a + b*x)/((-4*a*c + b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)) + (-12*a*c + 5*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(2*c**2*(-4*a*c + b**2)) + x*(-12*a*c + 15*b**2)*sqrt(a + b*x + c*x**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*c**(sympy.S(7)/2)*sqrt(a*x**2 + b*x**3 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_56():
    f = x**6/(a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)
    F = -2*b*sqrt(a*x**2 + b*x**3 + c*x**4)/(c*(-4*a*c + b**2)) - 3*b*x*sqrt(a + b*x + c*x**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*c**(sympy.S(5)/2)*sqrt(a*x**2 + b*x**3 + c*x**4)) + 2*x**3*(2*a + b*x)/((-4*a*c + b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)) + (-8*a*c + 3*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(c**2*x*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_57():
    f = x**5/(a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)
    F = -2*b*sqrt(a*x**2 + b*x**3 + c*x**4)/(c*x*(-4*a*c + b**2)) + 2*x**2*(2*a + b*x)/((-4*a*c + b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)) + x*sqrt(a + b*x + c*x**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(c**(sympy.S(3)/2)*sqrt(a*x**2 + b*x**3 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_58():
    f = x**4/(a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)
    F = 2*x*(2*a + b*x)/((-4*a*c + b**2)*sqrt(a*x**2 + b*x**3 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_59():
    f = x**3/(a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)
    F = -2*x*(b + 2*c*x)/((-4*a*c + b**2)*sqrt(a*x**2 + b*x**3 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_60():
    f = x**2/(a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)
    F = 2*x*(-2*a*c + b**2 + b*c*x)/(a*(-4*a*c + b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)) - atanh(x*(2*a + b*x)/(2*sqrt(a)*sqrt(a*x**2 + b*x**3 + c*x**4)))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_61():
    f = x/(a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2)
    F = (-4*a*c + 2*b**2 + 2*b*c*x)/(a*(-4*a*c + b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)) - (-8*a*c + 3*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(a**2*x**2*(-4*a*c + b**2)) + 3*b*atanh(x*(2*a + b*x)/(2*sqrt(a)*sqrt(a*x**2 + b*x**3 + c*x**4)))/(2*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_62():
    f = (a*x**2 + b*x**3 + c*x**4)**(sympy.S(-3)/2)
    F = (-4*a*c + 2*b**2 + 2*b*c*x)/(a*x*(-4*a*c + b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)) - (-12*a*c + 5*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(2*a**2*x**3*(-4*a*c + b**2)) + b*(-52*a*c + 15*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(4*a**3*x**2*(-4*a*c + b**2)) - (-12*a*c + 15*b**2)*atanh(x*(2*a + b*x)/(2*sqrt(a)*sqrt(a*x**2 + b*x**3 + c*x**4)))/(8*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_63():
    f = 1/(x*(a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2))
    F = (-4*a*c + 2*b**2 + 2*b*c*x)/(a*x**2*(-4*a*c + b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)) - (-16*a*c + 7*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(3*a**2*x**4*(-4*a*c + b**2)) + b*(-116*a*c + 35*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(12*a**3*x**3*(-4*a*c + b**2)) - sqrt(a*x**2 + b*x**3 + c*x**4)*(256*a**2*c**2 - 460*a*b**2*c + 105*b**4)/(24*a**4*x**2*(-4*a*c + b**2)) + 5*b*(-12*a*c + 7*b**2)*atanh(x*(2*a + b*x)/(2*sqrt(a)*sqrt(a*x**2 + b*x**3 + c*x**4)))/(16*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_64():
    f = 1/(x**2*(a*x**2 + b*x**3 + c*x**4)**(sympy.S(3)/2))
    F = (-4*a*c + 2*b**2 + 2*b*c*x)/(a*x**3*(-4*a*c + b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)) - (-20*a*c + 9*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(4*a**2*x**5*(-4*a*c + b**2)) + b*(-68*a*c + 21*b**2)*sqrt(a*x**2 + b*x**3 + c*x**4)/(8*a**3*x**4*(-4*a*c + b**2)) - sqrt(a*x**2 + b*x**3 + c*x**4)*(240*a**2*c**2 - 448*a*b**2*c + 105*b**4)/(32*a**4*x**3*(-4*a*c + b**2)) + b*sqrt(a*x**2 + b*x**3 + c*x**4)*(1808*a**2*c**2 - 1680*a*b**2*c + 315*b**4)/(64*a**5*x**2*(-4*a*c + b**2)) - (240*a**2*c**2 - 840*a*b**2*c + 315*b**4)*atanh(x*(2*a + b*x)/(2*sqrt(a)*sqrt(a*x**2 + b*x**3 + c*x**4)))/(128*a**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_65():
    f = x**m*(a*x + b*x**3 + c*x**5)
    F = a*x**(m + 2)/(m + 2) + b*x**(m + 4)/(m + 4) + c*x**(m + 6)/(m + 6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_66():
    f = x**2*(a*x + b*x**3 + c*x**5)
    F = a*x**4/4 + b*x**6/6 + c*x**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_67():
    f = x*(a*x + b*x**3 + c*x**5)
    F = a*x**3/3 + b*x**5/5 + c*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_68():
    f = a*x + b*x**3 + c*x**5
    F = a*x**2/2 + b*x**4/4 + c*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_69():
    f = (a*x + b*x**3 + c*x**5)/x
    F = a*x + b*x**3/3 + c*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_70():
    f = (a*x + b*x**3 + c*x**5)/x**2
    F = a*log(x) + b*x**2/2 + c*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_71():
    f = (a*x + b*x**3 + c*x**5)/x**3
    F = -a/x + b*x + c*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_72():
    f = x**m*(a*x + b*x**3 + c*x**5)**2
    F = a**2*x**(m + 3)/(m + 3) + 2*a*b*x**(m + 5)/(m + 5) + 2*b*c*x**(m + 9)/(m + 9) + c**2*x**(m + 11)/(m + 11) + x**(m + 7)*(2*a*c + b**2)/(m + 7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_73():
    f = x**2*(a*x + b*x**3 + c*x**5)**2
    F = a**2*x**5/5 + 2*a*b*x**7/7 + 2*b*c*x**11/11 + c**2*x**13/13 + x**9*(2*a*c + b**2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_74():
    f = x*(a*x + b*x**3 + c*x**5)**2
    F = a**2*x**4/4 + a*b*x**6/3 + b*c*x**10/5 + c**2*x**12/12 + x**8*(2*a*c + b**2)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_75():
    f = (a*x + b*x**3 + c*x**5)**2
    F = a**2*x**3/3 + 2*a*b*x**5/5 + 2*b*c*x**9/9 + c**2*x**11/11 + x**7*(2*a*c + b**2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_76():
    f = (a*x + b*x**3 + c*x**5)**2/x
    F = a**2*x**2/2 + a*b*x**4/2 + b*c*x**8/4 + c**2*x**10/10 + x**6*(2*a*c + b**2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_77():
    f = (a*x + b*x**3 + c*x**5)**2/x**2
    F = a**2*x + 2*a*b*x**3/3 + 2*b*c*x**7/7 + c**2*x**9/9 + x**5*(2*a*c + b**2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_78():
    f = x**8/(a*x + b*x**3 + c*x**5)
    F = -b*x**2/(2*c**2) + b*(-3*a*c + b**2)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**3*sqrt(-4*a*c + b**2)) + x**4/(4*c) + (-a*c + b**2)*log(a + b*x**2 + c*x**4)/(4*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_79():
    f = x**7/(a*x + b*x**3 + c*x**5)
    F = -b*x/c**2 + x**3/(3*c) + sqrt(2)*(-a*c + b**2 + b*(-3*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(5)/2)*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*(-a*c + b**2 - b*(-3*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(5)/2)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_80():
    f = x**6/(a*x + b*x**3 + c*x**5)
    F = -b*log(a + b*x**2 + c*x**4)/(4*c**2) + x**2/(2*c) - (-2*a*c + b**2)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_81():
    f = x**5/(a*x + b*x**3 + c*x**5)
    F = x/c - sqrt(2)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2))) - sqrt(2)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_82():
    f = x**4/(a*x + b*x**3 + c*x**5)
    F = b*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c*sqrt(-4*a*c + b**2)) + log(a + b*x**2 + c*x**4)/(4*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_83():
    f = x**3/(a*x + b*x**3 + c*x**5)
    F = -sqrt(2)*sqrt(b - sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(-4*a*c + b**2)) + sqrt(2)*sqrt(b + sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_84():
    f = x**2/(a*x + b*x**3 + c*x**5)
    F = -atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_85():
    f = x/(a*x + b*x**3 + c*x**5)
    F = -sqrt(2)*sqrt(c)*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)) + sqrt(2)*sqrt(c)*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_86():
    f = 1/(a*x + b*x**3 + c*x**5)
    F = b*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a*sqrt(-4*a*c + b**2)) + log(x)/a - log(a + b*x**2 + c*x**4)/(4*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_87():
    f = 1/(x*(a*x + b*x**3 + c*x**5))
    F = -sqrt(2)*sqrt(c)*(-b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*a*sqrt(b + sqrt(-4*a*c + b**2))) - sqrt(2)*sqrt(c)*(b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*a*sqrt(b - sqrt(-4*a*c + b**2))) - 1/(a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_88():
    f = 1/(x**2*(a*x + b*x**3 + c*x**5))
    F = -1/(2*a*x**2) - b*log(x)/a**2 + b*log(a + b*x**2 + c*x**4)/(4*a**2) - (-2*a*c + b**2)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_89():
    f = x**11/(a*x + b*x**3 + c*x**5)**2
    F = -b*x**4/(2*c*(-4*a*c + b**2)) - b*log(a + b*x**2 + c*x**4)/(2*c**3) + x**6*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) + x**2*(-3*a*c + b**2)/(c**2*(-4*a*c + b**2)) - (6*a**2*c**2 - 6*a*b**2*c + b**4)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(c**3*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_90():
    f = x**10/(a*x + b*x**3 + c*x**5)**2
    F = -b*x**3/(2*c*(-4*a*c + b**2)) + x**5*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) + x*(-10*a*c + 3*b**2)/(2*c**2*(-4*a*c + b**2)) - sqrt(2)*(-13*a*b*c + 3*b**3 + (20*a**2*c**2 - 19*a*b**2*c + 3*b**4)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(5)/2)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) - sqrt(2)*(-13*a*b*c + 3*b**3 - (20*a**2*c**2 - 19*a*b**2*c + 3*b**4)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(5)/2)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_91():
    f = x**9/(a*x + b*x**3 + c*x**5)**2
    F = -b*x**2/(2*c*(-4*a*c + b**2)) + b*(-6*a*c + b**2)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*c**2*(-4*a*c + b**2)**(sympy.S(3)/2)) + x**4*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) + log(a + b*x**2 + c*x**4)/(4*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_92():
    f = x**8/(a*x + b*x**3 + c*x**5)**2
    F = -b*x/(2*c*(-4*a*c + b**2)) + x**3*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*(-6*a*c + b**2 + b*(-8*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*(-6*a*c + b**2 - b*(-8*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_93():
    f = x**7/(a*x + b*x**3 + c*x**5)**2
    F = 2*a*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) + x**2*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_94():
    f = x**6/(a*x + b*x**3 + c*x**5)**2
    F = x*(2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*(b - (4*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*(4*a*c + b**2 + b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_95():
    f = x**5/(a*x + b*x**3 + c*x**5)**2
    F = -b*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) + (2*a + b*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_96():
    f = x**4/(a*x + b*x**3 + c*x**5)**2
    F = -sqrt(2)*sqrt(c)*(2*b + sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + sqrt(2)*sqrt(c)*(2*b - sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) - x*(b + 2*c*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_97():
    f = x**3/(a*x + b*x**3 + c*x**5)**2
    F = 2*c*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) - (b + 2*c*x**2)/((-8*a*c + 2*b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_98():
    f = x**2/(a*x + b*x**3 + c*x**5)**2
    F = -sqrt(2)*sqrt(c)*(-12*a*c + b**2 - b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + sqrt(2)*sqrt(c)*(-12*a*c + b**2 + b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + x*(-2*a*c + b**2 + b*c*x**2)/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_99():
    f = x/(a*x + b*x**3 + c*x**5)**2
    F = (-2*a*c + b**2 + b*c*x**2)/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + b*(-6*a*c + b**2)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**2*(-4*a*c + b**2)**(sympy.S(3)/2)) + log(x)/a**2 - log(a + b*x**2 + c*x**4)/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_100():
    f = (a*x + b*x**3 + c*x**5)**(-2)
    F = (-2*a*c + b**2 + b*c*x**2)/(2*a*x*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*sqrt(c)*(-16*a*b*c + 3*b**3 - (-10*a*c + 3*b**2)*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a**2*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) - sqrt(2)*sqrt(c)*(-16*a*b*c + 3*b**3 + (-10*a*c + 3*b**2)*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a**2*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) - (-10*a*c + 3*b**2)/(2*a**2*x*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_101():
    f = 1/(x*(a*x + b*x**3 + c*x**5)**2)
    F = (-2*a*c + b**2 + b*c*x**2)/(2*a*x**2*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) - (-3*a*c + b**2)/(a**2*x**2*(-4*a*c + b**2)) - 2*b*log(x)/a**3 + b*log(a + b*x**2 + c*x**4)/(2*a**3) - (6*a**2*c**2 - 6*a*b**2*c + b**4)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(a**3*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_102():
    f = 1/(x**2*(a*x + b*x**3 + c*x**5)**2)
    F = (-2*a*c + b**2 + b*c*x**2)/(2*a*x**3*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) - (-14*a*c + 5*b**2)/(6*a**2*x**3*(-4*a*c + b**2)) + b*(-19*a*c + 5*b**2)/(2*a**3*x*(-4*a*c + b**2)) - sqrt(2)*sqrt(c)*(28*a**2*c**2 - 29*a*b**2*c + 5*b**4 - b*(-19*a*c + 5*b**2)*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a**3*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + sqrt(2)*sqrt(c)*(28*a**2*c**2 - 29*a*b**2*c + 5*b**4 + b*(-19*a*c + 5*b**2)*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a**3*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_103():
    f = 1/(x**3*(a*x + b*x**3 + c*x**5)**2)
    F = (-2*a*c + b**2 + b*c*x**2)/(2*a*x**4*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) - (-8*a*c + 3*b**2)/(4*a**2*x**4*(-4*a*c + b**2)) + b*(-11*a*c + 3*b**2)/(2*a**3*x**2*(-4*a*c + b**2)) + b*(30*a**2*c**2 - 20*a*b**2*c + 3*b**4)*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/(2*a**4*(-4*a*c + b**2)**(sympy.S(3)/2)) + (-2*a*c + 3*b**2)*log(x)/a**4 - (-2*a*c + 3*b**2)*log(a + b*x**2 + c*x**4)/(4*a**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_104():
    f = x/sqrt(a*x + b*x**3 + c*x**5)
    F = 2*x**2*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(3)/4, sympy.S.Half, sympy.S.Half, sympy.S(7)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(3*sqrt(a*x + b*x**3 + c*x**5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_105():
    f = x**(sympy.S(3)/2)*sqrt(a*x + b*x**3 + c*x**5)
    F = 2*a**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-3*a*c + b**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(15*c**(sympy.S(7)/4)*sqrt(a*x + b*x**3 + c*x**5)) - a**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(sqrt(a)*b*sqrt(c) - 6*a*c + 2*b**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(30*c**(sympy.S(7)/4)*sqrt(a*x + b*x**3 + c*x**5)) + sqrt(x)*(b + 3*c*x**2)*sqrt(a*x + b*x**3 + c*x**5)/(15*c) - x**(sympy.S(3)/2)*(-6*a*c + 2*b**2)*(a + b*x**2 + c*x**4)/(15*c**(sympy.S(3)/2)*(sqrt(a) + sqrt(c)*x**2)*sqrt(a*x + b*x**3 + c*x**5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_106():
    f = sqrt(x)*sqrt(a*x + b*x**3 + c*x**5)
    F = (b + 2*c*x**2)*sqrt(a*x + b*x**3 + c*x**5)/(8*c*sqrt(x)) - sqrt(x)*(-4*a*c + b**2)*sqrt(a + b*x**2 + c*x**4)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(16*c**(sympy.S(3)/2)*sqrt(a*x + b*x**3 + c*x**5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_107():
    f = sqrt(a*x + b*x**3 + c*x**5)/sqrt(x)
    F = -a**(sympy.S(1)/4)*b*sqrt(x)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(3*c**(sympy.S(3)/4)*sqrt(a*x + b*x**3 + c*x**5)) + a**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(2*sqrt(a)*sqrt(c) + b)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(6*c**(sympy.S(3)/4)*sqrt(a*x + b*x**3 + c*x**5)) + b*x**(sympy.S(3)/2)*(a + b*x**2 + c*x**4)/(3*sqrt(c)*(sqrt(a) + sqrt(c)*x**2)*sqrt(a*x + b*x**3 + c*x**5)) + sqrt(x)*sqrt(a*x + b*x**3 + c*x**5)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_108():
    f = sqrt(a*x + b*x**3 + c*x**5)/x**(sympy.S(3)/2)
    F = -sqrt(a)*sqrt(x)*sqrt(a + b*x**2 + c*x**4)*atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(2*sqrt(a*x + b*x**3 + c*x**5)) + b*sqrt(x)*sqrt(a + b*x**2 + c*x**4)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(4*sqrt(c)*sqrt(a*x + b*x**3 + c*x**5)) + sqrt(a*x + b*x**3 + c*x**5)/(2*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_109():
    f = x**(sympy.S(3)/2)*(a*x + b*x**3 + c*x**5)**(sympy.S(3)/2)
    F = -3*b*sqrt(x)*(-4*a*c + b**2)**2*sqrt(a + b*x**2 + c*x**4)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(512*c**(sympy.S(7)/2)*sqrt(a*x + b*x**3 + c*x**5)) + sqrt(x)*(3*b + 8*c*x**2)*(a*x + b*x**3 + c*x**5)**(sympy.S(3)/2)/(80*c) - x**(sympy.S(3)/2)*(b*(-4*a*c + 5*b**2) + 4*c*x**2*(-16*a*c + 5*b**2))*sqrt(a*x + b*x**3 + c*x**5)/(640*c**2) + sqrt(a*x + b*x**3 + c*x**5)*(128*a**2*c**2 - 100*a*b**2*c + 15*b**4)/(1280*c**3*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_110():
    f = sqrt(x)*(a*x + b*x**3 + c*x**5)**(sympy.S(3)/2)
    F = -a**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(84*a**2*c**2 - 57*a*b**2*c + 8*b**4)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(315*c**(sympy.S(11)/4)*sqrt(a*x + b*x**3 + c*x**5)) + a**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(4*sqrt(a)*b*sqrt(c)*(-6*a*c + b**2) + 84*a**2*c**2 - 57*a*b**2*c + 8*b**4)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(630*c**(sympy.S(11)/4)*sqrt(a*x + b*x**3 + c*x**5)) + (3*b + 7*c*x**2)*(a*x + b*x**3 + c*x**5)**(sympy.S(3)/2)/(63*c*sqrt(x)) - sqrt(x)*(b*(-9*a*c + 4*b**2) + 6*c*x**2*(-7*a*c + 2*b**2))*sqrt(a*x + b*x**3 + c*x**5)/(315*c**2) + x**(sympy.S(3)/2)*(a + b*x**2 + c*x**4)*(84*a**2*c**2 - 57*a*b**2*c + 8*b**4)/(315*c**(sympy.S(5)/2)*(sqrt(a) + sqrt(c)*x**2)*sqrt(a*x + b*x**3 + c*x**5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_111():
    f = (a*x + b*x**3 + c*x**5)**(sympy.S(3)/2)/sqrt(x)
    F = (b + 2*c*x**2)*(a*x + b*x**3 + c*x**5)**(sympy.S(3)/2)/(16*c*x**(sympy.S(3)/2)) - (b + 2*c*x**2)*(-12*a*c + 3*b**2)*sqrt(a*x + b*x**3 + c*x**5)/(128*c**2*sqrt(x)) + 3*sqrt(x)*(-4*a*c + b**2)**2*sqrt(a + b*x**2 + c*x**4)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(256*c**(sympy.S(5)/2)*sqrt(a*x + b*x**3 + c*x**5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_112():
    f = (a*x + b*x**3 + c*x**5)**(sympy.S(3)/2)/x**(sympy.S(3)/2)
    F = 2*a**(sympy.S(1)/4)*b*sqrt(x)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-8*a*c + b**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(35*c**(sympy.S(7)/4)*sqrt(a*x + b*x**3 + c*x**5)) - a**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(sqrt(a)*sqrt(c)*(-20*a*c + b**2) + 2*b*(-8*a*c + b**2))*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(70*c**(sympy.S(7)/4)*sqrt(a*x + b*x**3 + c*x**5)) - 2*b*x**(sympy.S(3)/2)*(-8*a*c + b**2)*(a + b*x**2 + c*x**4)/(35*c**(sympy.S(3)/2)*(sqrt(a) + sqrt(c)*x**2)*sqrt(a*x + b*x**3 + c*x**5)) + (a*x + b*x**3 + c*x**5)**(sympy.S(3)/2)/(7*sqrt(x)) + sqrt(x)*(10*a*c + b**2 + 3*b*c*x**2)*sqrt(a*x + b*x**3 + c*x**5)/(35*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_113():
    f = x**(sympy.S(3)/2)/sqrt(a*x + b*x**3 + c*x**5)
    F = sqrt(x)*sqrt(a + b*x**2 + c*x**4)*atanh((b + 2*c*x**2)/(2*sqrt(c)*sqrt(a + b*x**2 + c*x**4)))/(2*sqrt(c)*sqrt(a*x + b*x**3 + c*x**5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_114():
    f = sqrt(x)/sqrt(a*x + b*x**3 + c*x**5)
    F = sqrt(x)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(2*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(a*x + b*x**3 + c*x**5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_115():
    f = 1/(sqrt(x)*sqrt(a*x + b*x**3 + c*x**5))
    F = -atanh(sqrt(x)*(2*a + b*x**2)/(2*sqrt(a)*sqrt(a*x + b*x**3 + c*x**5)))/(2*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_116():
    f = 1/(x**(sympy.S(3)/2)*sqrt(a*x + b*x**3 + c*x**5))
    F = sqrt(c)*x**(sympy.S(3)/2)*(a + b*x**2 + c*x**4)/(a*(sqrt(a) + sqrt(c)*x**2)*sqrt(a*x + b*x**3 + c*x**5)) - sqrt(a*x + b*x**3 + c*x**5)/(a*x**(sympy.S(3)/2)) - c**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(a**(sympy.S(3)/4)*sqrt(a*x + b*x**3 + c*x**5)) + c**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(2*a**(sympy.S(3)/4)*sqrt(a*x + b*x**3 + c*x**5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_117():
    f = x**(sympy.S(3)/2)/(a*x + b*x**3 + c*x**5)**(sympy.S(3)/2)
    F = -b*sqrt(c)*x**(sympy.S(3)/2)*(a + b*x**2 + c*x**4)/(a*(sqrt(a) + sqrt(c)*x**2)*(-4*a*c + b**2)*sqrt(a*x + b*x**3 + c*x**5)) + x**(sympy.S(3)/2)*(-2*a*c + b**2 + b*c*x**2)/(a*(-4*a*c + b**2)*sqrt(a*x + b*x**3 + c*x**5)) + b*c**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(a**(sympy.S(3)/4)*(-4*a*c + b**2)*sqrt(a*x + b*x**3 + c*x**5)) - c**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(2*a**(sympy.S(3)/4)*(-2*sqrt(a)*sqrt(c) + b)*sqrt(a*x + b*x**3 + c*x**5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_118():
    f = sqrt(x)/(a*x + b*x**3 + c*x**5)**(sympy.S(3)/2)
    F = sqrt(x)*(-2*a*c + b**2 + b*c*x**2)/(a*(-4*a*c + b**2)*sqrt(a*x + b*x**3 + c*x**5)) - atanh(sqrt(x)*(2*a + b*x**2)/(2*sqrt(a)*sqrt(a*x + b*x**3 + c*x**5)))/(2*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_119():
    f = 1/(sqrt(x)*(a*x + b*x**3 + c*x**5)**(sympy.S(3)/2))
    F = (-2*a*c + b**2 + b*c*x**2)/(a*sqrt(x)*(-4*a*c + b**2)*sqrt(a*x + b*x**3 + c*x**5)) + 2*sqrt(c)*x**(sympy.S(3)/2)*(-3*a*c + b**2)*(a + b*x**2 + c*x**4)/(a**2*(sqrt(a) + sqrt(c)*x**2)*(-4*a*c + b**2)*sqrt(a*x + b*x**3 + c*x**5)) - (-6*a*c + 2*b**2)*sqrt(a*x + b*x**3 + c*x**5)/(a**2*x**(sympy.S(3)/2)*(-4*a*c + b**2)) - 2*c**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-3*a*c + b**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(a**(sympy.S(7)/4)*(-4*a*c + b**2)*sqrt(a*x + b*x**3 + c*x**5)) + c**(sympy.S(1)/4)*sqrt(x)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(sqrt(a)*b*sqrt(c) - 6*a*c + 2*b**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(2*a**(sympy.S(7)/4)*(-4*a*c + b**2)*sqrt(a*x + b*x**3 + c*x**5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_120():
    f = 1/(x**(sympy.S(3)/2)*(a*x + b*x**3 + c*x**5)**(sympy.S(3)/2))
    F = (-2*a*c + b**2 + b*c*x**2)/(a*x**(sympy.S(3)/2)*(-4*a*c + b**2)*sqrt(a*x + b*x**3 + c*x**5)) - (-8*a*c + 3*b**2)*sqrt(a*x + b*x**3 + c*x**5)/(2*a**2*x**(sympy.S(5)/2)*(-4*a*c + b**2)) + 3*b*atanh(sqrt(x)*(2*a + b*x**2)/(2*sqrt(a)*sqrt(a*x + b*x**3 + c*x**5)))/(4*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_121():
    f = x**(3*n/2 + sympy.S(-3)/2)/(a*x**(n - 1) + b*x**n + c*x**(n + 1))**(sympy.S(3)/2)
    F = -2*x**(n/2 + sympy.S(-1)/2)*(b + 2*c*x)/((-4*a*c + b**2)*sqrt(a*x**(n - 1) + b*x**n + c*x**(n + 1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_122():
    f = x*(d + e*x**2)/sqrt(a*x + b*x**3 + c*x**5)
    F = 2*d*x**2*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(3)/4, sympy.S.Half, sympy.S.Half, sympy.S(7)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(3*sqrt(a*x + b*x**3 + c*x**5)) + 2*e*x**4*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(7)/4, sympy.S.Half, sympy.S.Half, sympy.S(11)/4, -2*c*x**2/(b - sqrt(-4*a*c + b**2)), -2*c*x**2/(b + sqrt(-4*a*c + b**2)))/(7*sqrt(a*x + b*x**3 + c*x**5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_123():
    f = 1/sqrt(x**6 - 3*x**4 + 3*x**2)
    F = -sqrt(3)*atanh(sqrt(3)*x*(6 - 3*x**2)/(6*sqrt(x**6 - 3*x**4 + 3*x**2)))/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_124():
    f = 1/sqrt(x**2*(x**4 - 3*x**2 + 3))
    F = -sqrt(3)*atanh(sqrt(3)*x*(6 - 3*x**2)/(6*sqrt(x**6 - 3*x**4 + 3*x**2)))/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_125():
    f = 1/sqrt(1 - (1 - x**2)**3)
    F = -sqrt(3)*atanh(sqrt(3)*x*(6 - 3*x**2)/(6*sqrt(x**6 - 3*x**4 + 3*x**2)))/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_126():
    f = sqrt(x**6 - 3*x**4 + 3*x**2)
    F = -(3 - 2*x**2)*sqrt(x**6 - 3*x**4 + 3*x**2)/(8*x) - 3*sqrt(x**6 - 3*x**4 + 3*x**2)*asinh(sqrt(3)*(3 - 2*x**2)/3)/(16*x*sqrt(x**4 - 3*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_127():
    f = sqrt(x**2*(x**4 - 3*x**2 + 3))
    F = -(3 - 2*x**2)*sqrt(x**6 - 3*x**4 + 3*x**2)/(8*x) - 3*sqrt(x**6 - 3*x**4 + 3*x**2)*asinh(sqrt(3)*(3 - 2*x**2)/3)/(16*x*sqrt(x**4 - 3*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_128():
    f = sqrt(1 - (1 - x**2)**3)
    F = -(3 - 2*x**2)*sqrt(x**6 - 3*x**4 + 3*x**2)/(8*x) - 3*sqrt(x**6 - 3*x**4 + 3*x**2)*asinh(sqrt(3)*(3 - 2*x**2)/3)/(16*x*sqrt(x**4 - 3*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_129():
    f = 1/(x*sqrt(a + b*x + c*x**2))
    F = -atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_130():
    f = 1/sqrt(x**2*(a + b*x + c*x**2))
    F = -atanh(x*(2*a + b*x)/(2*sqrt(a)*sqrt(a*x**2 + b*x**3 + c*x**4)))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_131():
    f = 1/(sqrt(x)*sqrt(x*(a + b*x + c*x**2)))
    F = -atanh(sqrt(x)*(2*a + b*x)/(2*sqrt(a)*sqrt(a*x + b*x**2 + c*x**3)))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_132():
    f = sqrt(x)/sqrt(x**3*(a + b*x + c*x**2))
    F = -atanh(x**(sympy.S(3)/2)*(2*a + b*x)/(2*sqrt(a)*sqrt(a*x**3 + b*x**4 + c*x**5)))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_133():
    f = 1/(x*sqrt(a + b*x**2 + c*x**4))
    F = -atanh((2*a + b*x**2)/(2*sqrt(a)*sqrt(a + b*x**2 + c*x**4)))/(2*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_134():
    f = 1/sqrt(x**2*(a + b*x**2 + c*x**4))
    F = -atanh(x*(2*a + b*x**2)/(2*sqrt(a)*sqrt(a*x**2 + b*x**4 + c*x**6)))/(2*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_135():
    f = 1/(sqrt(x)*sqrt(x*(a + b*x**2 + c*x**4)))
    F = -atanh(sqrt(x)*(2*a + b*x**2)/(2*sqrt(a)*sqrt(a*x + b*x**3 + c*x**5)))/(2*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_136():
    f = sqrt(x)/sqrt(x**3*(a + b*x**2 + c*x**4))
    F = -atanh(x**(sympy.S(3)/2)*(2*a + b*x**2)/(2*sqrt(a)*sqrt(a*x**3 + b*x**5 + c*x**7)))/(2*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_137():
    f = 1/(x*sqrt(x**4 - 3*x**2 + 3))
    F = -sqrt(3)*atanh(sqrt(3)*(2 - x**2)/(2*sqrt(x**4 - 3*x**2 + 3)))/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_138():
    f = 1/sqrt(x**2*(x**4 - 3*x**2 + 3))
    F = -sqrt(3)*atanh(sqrt(3)*x*(6 - 3*x**2)/(6*sqrt(x**6 - 3*x**4 + 3*x**2)))/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_139():
    f = 1/(sqrt(x)*sqrt(x*(x**2 - 3*x + 3)))
    F = -sqrt(3)*atanh(sqrt(3)*sqrt(x)*(2 - x)/(2*sqrt(x**3 - 3*x**2 + 3*x)))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_4_Improper_1_2_4_2_d_x_pow_m_a_x_pow_q_plus_b_x_pow_n_plus_c_x_pow_2_n_minus_q_pow_p_140():
    f = x**(q/2 - 1)/sqrt(a*x**q + b*x**n + c*x**(2*n - q))
    F = -atanh(x**(q/2)*(2*a + b*x**(n - q))/(2*sqrt(a)*sqrt(a*x**q + b*x**n + c*x**(2*n - q))))/(sqrt(a)*(n - q))
    assert integrate(f, x) == F

