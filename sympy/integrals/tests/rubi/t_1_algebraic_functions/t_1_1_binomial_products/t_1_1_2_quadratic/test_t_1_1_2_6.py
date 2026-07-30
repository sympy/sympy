"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.1 Binomial products/1.1.2 Quadratic/1.1.2.6 (g x)^m (a+b x^2)^p (c+d x^2)^q (e+f x^2)^r.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, a, b, c, d, e, m, p = symbols('A B a b c d e m p')

def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_1():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)**3*(c + d*x**2)
    F = A*a**3*c*(e*x)**(m + 1)/(e*(m + 1)) + B*b**3*d*(e*x)**(m + 11)/(e**11*(m + 11)) + a**2*(e*x)**(m + 3)*(A*a*d + 3*A*b*c + B*a*c)/(e**3*(m + 3)) + a*(e*x)**(m + 5)*(3*A*b*(a*d + b*c) + B*a*(a*d + 3*b*c))/(e**5*(m + 5)) + b**2*(e*x)**(m + 9)*(A*b*d + 3*B*a*d + B*b*c)/(e**9*(m + 9)) + b*(e*x)**(m + 7)*(A*b*(3*a*d + b*c) + 3*B*a*(a*d + b*c))/(e**7*(m + 7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_2():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)**2*(c + d*x**2)
    F = A*a**2*c*(e*x)**(m + 1)/(e*(m + 1)) + B*b**2*d*(e*x)**(m + 9)/(e**9*(m + 9)) + a*(e*x)**(m + 3)*(A*a*d + 2*A*b*c + B*a*c)/(e**3*(m + 3)) + b*(e*x)**(m + 7)*(A*b*d + 2*B*a*d + B*b*c)/(e**7*(m + 7)) + (e*x)**(m + 5)*(A*b*(2*a*d + b*c) + B*a*(a*d + 2*b*c))/(e**5*(m + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_3():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)*(c + d*x**2)
    F = A*a*c*(e*x)**(m + 1)/(e*(m + 1)) + B*b*d*(e*x)**(m + 7)/(e**7*(m + 7)) + (e*x)**(m + 3)*(A*a*d + A*b*c + B*a*c)/(e**3*(m + 3)) + (e*x)**(m + 5)*(A*b*d + B*a*d + B*b*c)/(e**5*(m + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_4():
    f = (e*x)**m*(A + B*x**2)*(c + d*x**2)
    F = A*c*(e*x)**(m + 1)/(e*(m + 1)) + B*d*(e*x)**(m + 5)/(e**5*(m + 5)) + (e*x)**(m + 3)*(A*d + B*c)/(e**3*(m + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_5():
    f = (e*x)**m*(A + B*x**2)*(c + d*x**2)/(a + b*x**2)
    F = B*d*(e*x)**(m + 3)/(b*e**3*(m + 3)) + (e*x)**(m + 1)*(A*b*d - B*a*d + B*b*c)/(b**2*e*(m + 1)) + (e*x)**(m + 1)*(A*b - B*a)*(-a*d + b*c)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a*b**2*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_6():
    f = (e*x)**m*(A + B*x**2)*(c + d*x**2)/(a + b*x**2)**2
    F = (e*x)**(m + 1)*(c + d*x**2)*(A*b - B*a)/(2*a*b*e*(a + b*x**2)) - d*(e*x)**(m + 1)*(A*b*(m + 1) - B*a*(m + 3))/(2*a*b**2*e*(m + 1)) + (e*x)**(m + 1)*(A*b*(a*d*(m + 1) + b*(-c*m + c)) + B*a*(-a*d*(m + 3) + b*c*(m + 1)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(2*a**2*b**2*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_7():
    f = (e*x)**m*(A + B*x**2)*(c + d*x**2)/(a + b*x**2)**3
    F = (e*x)**(m + 1)*(c + d*x**2)*(A*b - B*a)/(4*a*b*e*(a + b*x**2)**2) - (e*x)**(m + 1)*(A*b*(a*d*(1 - m) - b*c*(3 - m)) - B*a*(-a*d*(m + 3) + b*c*(m + 1)))/(8*a**2*b**2*e*(a + b*x**2)) + (e*x)**(m + 1)*(A*b*(1 - m)*(a*d*(m + 1) + b*c*(3 - m)) + B*a*(m + 1)*(a*d*(m + 3) + b*(-c*m + c)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(8*a**3*b**2*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_8():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)**3*(c + d*x**2)**2
    F = A*a**3*c**2*(e*x)**(m + 1)/(e*(m + 1)) + B*b**3*d**2*(e*x)**(m + 13)/(e**13*(m + 13)) + a**2*c*(e*x)**(m + 3)*(2*A*a*d + 3*A*b*c + B*a*c)/(e**3*(m + 3)) + a*(e*x)**(m + 5)*(A*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2) + B*a*c*(2*a*d + 3*b*c))/(e**5*(m + 5)) + b**2*d*(e*x)**(m + 11)*(A*b*d + 3*B*a*d + 2*B*b*c)/(e**11*(m + 11)) + b*(e*x)**(m + 9)*(3*B*a**2*d**2 + 3*a*b*d*(A*d + 2*B*c) + b**2*c*(2*A*d + B*c))/(e**9*(m + 9)) + (e*x)**(m + 7)*(A*b*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2) + B*a*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2))/(e**7*(m + 7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_9():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)**2*(c + d*x**2)**2
    F = A*a**2*c**2*(e*x)**(m + 1)/(e*(m + 1)) + B*b**2*d**2*(e*x)**(m + 11)/(e**11*(m + 11)) + a*c*(e*x)**(m + 3)*(2*A*(a*d + b*c) + B*a*c)/(e**3*(m + 3)) + b*d*(e*x)**(m + 9)*(A*b*d + 2*B*a*d + 2*B*b*c)/(e**9*(m + 9)) + (e*x)**(m + 5)*(A*(a**2*d**2 + 4*a*b*c*d + b**2*c**2) + 2*B*a*c*(a*d + b*c))/(e**5*(m + 5)) + (e*x)**(m + 7)*(B*a**2*d**2 + 2*a*b*d*(A*d + 2*B*c) + b**2*c*(2*A*d + B*c))/(e**7*(m + 7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_10():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)*(c + d*x**2)**2
    F = A*a*c**2*(e*x)**(m + 1)/(e*(m + 1)) + B*b*d**2*(e*x)**(m + 9)/(e**9*(m + 9)) + c*(e*x)**(m + 3)*(2*A*a*d + A*b*c + B*a*c)/(e**3*(m + 3)) + d*(e*x)**(m + 7)*(A*b*d + B*a*d + 2*B*b*c)/(e**7*(m + 7)) + (e*x)**(m + 5)*(a*d*(A*d + 2*B*c) + b*c*(2*A*d + B*c))/(e**5*(m + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_11():
    f = (e*x)**m*(A + B*x**2)*(c + d*x**2)**2
    F = A*c**2*(e*x)**(m + 1)/(e*(m + 1)) + B*d**2*(e*x)**(m + 7)/(e**7*(m + 7)) + c*(e*x)**(m + 3)*(2*A*d + B*c)/(e**3*(m + 3)) + d*(e*x)**(m + 5)*(A*d + 2*B*c)/(e**5*(m + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_12():
    f = (e*x)**m*(A + B*x**2)*(c + d*x**2)**2/(a + b*x**2)
    F = B*d**2*(e*x)**(m + 5)/(b*e**5*(m + 5)) + d*(e*x)**(m + 3)*(A*b*d - B*a*d + 2*B*b*c)/(b**2*e**3*(m + 3)) + (e*x)**(m + 1)*(B*a**2*d**2 - a*b*d*(A*d + 2*B*c) + b**2*c*(2*A*d + B*c))/(b**3*e*(m + 1)) + (e*x)**(m + 1)*(A*b - B*a)*(-a*d + b*c)**2*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a*b**3*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_13():
    f = (e*x)**m*(A + B*x**2)*(c + d*x**2)**2/(a + b*x**2)**2
    F = (e*x)**(m + 1)*(c + d*x**2)**2*(A*b - B*a)/(2*a*b*e*(a + b*x**2)) - d**2*(e*x)**(m + 3)*(A*b*(m + 3) - B*a*(m + 5))/(2*a*b**2*e**3*(m + 3)) - d*(e*x)**(m + 1)*(A*b*(-a*d*(m + 3) + 2*b*c*(m + 1)) - B*a*(-a*d*(m + 5) + 2*b*c*(m + 3)))/(2*a*b**3*e*(m + 1)) + (e*x)**(m + 1)*(-a*d + b*c)*(A*b*(a*d*(m + 3) + b*(-c*m + c)) + B*a*(-a*d*(m + 5) + b*c*(m + 1)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(2*a**2*b**3*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_14():
    f = (e*x)**m*(A + B*x**2)*(c + d*x**2)**2/(a + b*x**2)**3
    F = (e*x)**(m + 1)*(c + d*x**2)**2*(A*b - B*a)/(4*a*b*e*(a + b*x**2)**2) + (e*x)**(m + 1)*(-a*d + b*c)*(c*(A*b*(3 - m) + B*a*(m + 1)) - d*x**2*(A*b*(m + 1) - B*a*(m + 5)))/(8*a**2*b**2*e*(a + b*x**2)) + d*(e*x)**(m + 1)*(A*b*(m + 1) - B*a*(m + 5))*(-a*d*(m + 3) + b*c*(m + 1))/(8*a**2*b**3*e*(m + 1)) - (e*x)**(m + 1)*(a*d*(A*b*(m + 1) - B*a*(m + 5))*(-a*d*(m + 3) + b*c*(m + 1)) - b*c*(A*b*(3 - m) + B*a*(m + 1))*(a*d*(m + 1) + b*(-c*m + c)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(8*a**3*b**3*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_15():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)**3*(c + d*x**2)**3
    F = A*a**3*c**3*(e*x)**(m + 1)/(e*(m + 1)) + B*b**3*d**3*(e*x)**(m + 15)/(e**15*(m + 15)) + a**2*c**2*(e*x)**(m + 3)*(3*A*(a*d + b*c) + B*a*c)/(e**3*(m + 3)) + 3*a*c*(e*x)**(m + 5)*(A*(a**2*d**2 + 3*a*b*c*d + b**2*c**2) + B*a*c*(a*d + b*c))/(e**5*(m + 5)) + b**2*d**2*(e*x)**(m + 13)*(A*b*d + 3*B*a*d + 3*B*b*c)/(e**13*(m + 13)) + 3*b*d*(e*x)**(m + 11)*(B*a**2*d**2 + a*b*d*(A*d + 3*B*c) + b**2*c*(A*d + B*c))/(e**11*(m + 11)) + (e*x)**(m + 7)*(A*(a**3*d**3 + 9*a**2*b*c*d**2 + 9*a*b**2*c**2*d + b**3*c**3) + 3*B*a*c*(a**2*d**2 + 3*a*b*c*d + b**2*c**2))/(e**7*(m + 7)) + (e*x)**(m + 9)*(B*a**3*d**3 + 3*a**2*b*d**2*(A*d + 3*B*c) + 9*a*b**2*c*d*(A*d + B*c) + b**3*c**2*(3*A*d + B*c))/(e**9*(m + 9))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_16():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)**2*(c + d*x**2)**3
    F = A*a**2*c**3*(e*x)**(m + 1)/(e*(m + 1)) + B*b**2*d**3*(e*x)**(m + 13)/(e**13*(m + 13)) + a*c**2*(e*x)**(m + 3)*(3*A*a*d + 2*A*b*c + B*a*c)/(e**3*(m + 3)) + b*d**2*(e*x)**(m + 11)*(A*b*d + 2*B*a*d + 3*B*b*c)/(e**11*(m + 11)) + c*(e*x)**(m + 5)*(A*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2) + B*a*c*(3*a*d + 2*b*c))/(e**5*(m + 5)) + d*(e*x)**(m + 9)*(B*a**2*d**2 + 2*a*b*d*(A*d + 3*B*c) + 3*b**2*c*(A*d + B*c))/(e**9*(m + 9)) + (e*x)**(m + 7)*(a**2*d**2*(A*d + 3*B*c) + 6*a*b*c*d*(A*d + B*c) + b**2*c**2*(3*A*d + B*c))/(e**7*(m + 7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_17():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)*(c + d*x**2)**3
    F = A*a*c**3*(e*x)**(m + 1)/(e*(m + 1)) + B*b*d**3*(e*x)**(m + 11)/(e**11*(m + 11)) + c**2*(e*x)**(m + 3)*(3*A*a*d + A*b*c + B*a*c)/(e**3*(m + 3)) + c*(e*x)**(m + 5)*(3*a*d*(A*d + B*c) + b*c*(3*A*d + B*c))/(e**5*(m + 5)) + d**2*(e*x)**(m + 9)*(A*b*d + B*a*d + 3*B*b*c)/(e**9*(m + 9)) + d*(e*x)**(m + 7)*(a*d*(A*d + 3*B*c) + 3*b*c*(A*d + B*c))/(e**7*(m + 7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_18():
    f = (e*x)**m*(A + B*x**2)*(c + d*x**2)**3
    F = A*c**3*(e*x)**(m + 1)/(e*(m + 1)) + B*d**3*(e*x)**(m + 9)/(e**9*(m + 9)) + c**2*(e*x)**(m + 3)*(3*A*d + B*c)/(e**3*(m + 3)) + 3*c*d*(e*x)**(m + 5)*(A*d + B*c)/(e**5*(m + 5)) + d**2*(e*x)**(m + 7)*(A*d + 3*B*c)/(e**7*(m + 7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_19():
    f = (e*x)**m*(A + B*x**2)*(c + d*x**2)**3/(a + b*x**2)
    F = B*d**3*(e*x)**(m + 7)/(b*e**7*(m + 7)) + d**2*(e*x)**(m + 5)*(A*b*d - B*a*d + 3*B*b*c)/(b**2*e**5*(m + 5)) + d*(e*x)**(m + 3)*(B*a**2*d**2 - a*b*d*(A*d + 3*B*c) + 3*b**2*c*(A*d + B*c))/(b**3*e**3*(m + 3)) - (e*x)**(m + 1)*(B*a**3*d**3 - a**2*b*d**2*(A*d + 3*B*c) + 3*a*b**2*c*d*(A*d + B*c) - b**3*c**2*(3*A*d + B*c))/(b**4*e*(m + 1)) + (e*x)**(m + 1)*(A*b - B*a)*(-a*d + b*c)**3*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a*b**4*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_20():
    f = (e*x)**m*(A + B*x**2)*(c + d*x**2)**3/(a + b*x**2)**2
    F = (e*x)**(m + 1)*(c + d*x**2)**3*(A*b - B*a)/(2*a*b*e*(a + b*x**2)) - d**3*(e*x)**(m + 5)*(A*b*(m + 5) - B*a*(m + 7))/(2*a*b**2*e**5*(m + 5)) - d**2*(e*x)**(m + 3)*(A*b*(-a*d*(m + 5) + 3*b*c*(m + 3)) - B*a*(-a*d*(m + 7) + 3*b*c*(m + 5)))/(2*a*b**3*e**3*(m + 3)) - d*(e*x)**(m + 1)*(A*b*(a**2*d**2*(m + 5) - 3*a*b*c*d*(m + 3) + 3*b**2*c**2*(m + 1)) - B*a*(a**2*d**2*(m + 7) - 3*a*b*c*d*(m + 5) + 3*b**2*c**2*(m + 3)))/(2*a*b**4*e*(m + 1)) + (e*x)**(m + 1)*(-a*d + b*c)**2*(A*b*(a*d*(m + 5) + b*(-c*m + c)) + B*a*(-a*d*(m + 7) + b*c*(m + 1)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(2*a**2*b**4*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_21():
    f = (e*x)**m*(A + B*x**2)*(c + d*x**2)**3/(a + b*x**2)**3
    F = (e*x)**(m + 1)*(c + d*x**2)**3*(A*b - B*a)/(4*a*b*e*(a + b*x**2)**2) + (e*x)**(m + 1)*(c + d*x**2)**2*(A*b*(a*d*(m + 3) + b*c*(3 - m)) + B*a*(-a*d*(m + 7) + b*c*(m + 1)))/(8*a**2*b**2*e*(a + b*x**2)) - d**2*(e*x)**(m + 3)*(A*b*(m + 3)*(a*d*(m + 5) + b*c*(3 - m)) + B*a*(-a*d*(m**2 + 12*m + 35) + b*c*(m**2 + 4*m + 3)))/(8*a**2*b**3*e**3*(m + 3)) - d*(e*x)**(m + 1)*(A*b*(-a**2*d**2*(m**2 + 8*m + 15) + 3*a*b*c*d*(m**2 + 4*m + 3) + 2*b**2*c**2*(-m**2 + 2*m + 3)) + B*a*(a**2*d**2*(m**2 + 12*m + 35) - 3*a*b*c*d*(m**2 + 8*m + 15) + 2*b**2*c**2*(m + 1)**2))/(8*a**2*b**4*e*(m + 1)) + (e*x)**(m + 1)*(-a*d + b*c)*(A*b*(a**2*d**2*(m**2 + 8*m + 15) + 2*a*b*c*d*(-m**2 - 2*m + 3) + b**2*c**2*(m**2 - 4*m + 3)) + B*a*(-a**2*d**2*(m**2 + 12*m + 35) + 2*a*b*c*d*(m**2 + 6*m + 5) + b**2*c**2*(1 - m**2)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(8*a**3*b**4*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_22():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)**4/(c + d*x**2)
    F = B*b**4*(e*x)**(m + 9)/(d*e**9*(m + 9)) - b**3*(e*x)**(m + 7)*(-A*b*d - 4*B*a*d + B*b*c)/(d**2*e**7*(m + 7)) + b**2*(e*x)**(m + 5)*(6*B*a**2*d**2 - 4*a*b*d*(-A*d + B*c) + b**2*c*(-A*d + B*c))/(d**3*e**5*(m + 5)) + b*(e*x)**(m + 3)*(4*B*a**3*d**3 - 6*a**2*b*d**2*(-A*d + B*c) + 4*a*b**2*c*d*(-A*d + B*c) - b**3*c**2*(-A*d + B*c))/(d**4*e**3*(m + 3)) + (e*x)**(m + 1)*(B*a**4*d**4 - 4*a**3*b*d**3*(-A*d + B*c) + 6*a**2*b**2*c*d**2*(-A*d + B*c) - 4*a*b**3*c**2*d*(-A*d + B*c) + b**4*c**3*(-A*d + B*c))/(d**5*e*(m + 1)) - (e*x)**(m + 1)*(-A*d + B*c)*(-a*d + b*c)**4*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(c*d**5*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_23():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)**3/(c + d*x**2)
    F = B*b**3*(e*x)**(m + 7)/(d*e**7*(m + 7)) - b**2*(e*x)**(m + 5)*(-A*b*d - 3*B*a*d + B*b*c)/(d**2*e**5*(m + 5)) + b*(e*x)**(m + 3)*(3*B*a**2*d**2 - 3*a*b*d*(-A*d + B*c) + b**2*c*(-A*d + B*c))/(d**3*e**3*(m + 3)) + (e*x)**(m + 1)*(B*a**3*d**3 - 3*a**2*b*d**2*(-A*d + B*c) + 3*a*b**2*c*d*(-A*d + B*c) - b**3*c**2*(-A*d + B*c))/(d**4*e*(m + 1)) + (e*x)**(m + 1)*(-A*d + B*c)*(-a*d + b*c)**3*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(c*d**4*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_24():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)**2/(c + d*x**2)
    F = B*b**2*(e*x)**(m + 5)/(d*e**5*(m + 5)) - b*(e*x)**(m + 3)*(-A*b*d - 2*B*a*d + B*b*c)/(d**2*e**3*(m + 3)) + (e*x)**(m + 1)*(B*a**2*d**2 - 2*a*b*d*(-A*d + B*c) + b**2*c*(-A*d + B*c))/(d**3*e*(m + 1)) - (e*x)**(m + 1)*(-A*d + B*c)*(-a*d + b*c)**2*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(c*d**3*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_25():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)/(c + d*x**2)
    F = B*b*(e*x)**(m + 3)/(d*e**3*(m + 3)) - (e*x)**(m + 1)*(-A*b*d - B*a*d + B*b*c)/(d**2*e*(m + 1)) + (e*x)**(m + 1)*(-A*d + B*c)*(-a*d + b*c)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(c*d**2*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_26():
    f = (e*x)**m*(A + B*x**2)/(c + d*x**2)
    F = B*(e*x)**(m + 1)/(d*e*(m + 1)) - (e*x)**(m + 1)*(-A*d + B*c)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(c*d*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_27():
    f = (e*x)**m*(A + B*x**2)/((a + b*x**2)*(c + d*x**2))
    F = (e*x)**(m + 1)*(-A*d + B*c)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(c*e*(m + 1)*(-a*d + b*c)) + (e*x)**(m + 1)*(A*b - B*a)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a*e*(m + 1)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_28():
    f = (e*x)**m*(A + B*x**2)/((a + b*x**2)**2*(c + d*x**2))
    F = -d*(e*x)**(m + 1)*(-A*d + B*c)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(c*e*(m + 1)*(-a*d + b*c)**2) + (e*x)**(m + 1)*(A*b - B*a)/(2*a*e*(a + b*x**2)*(-a*d + b*c)) + (e*x)**(m + 1)*(A*b*(-a*d*(3 - m) + b*c*(1 - m)) + B*a*(a*d*(1 - m) + b*c*(m + 1)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(2*a**2*e*(m + 1)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_29():
    f = (e*x)**m*(A + B*x**2)/((a + b*x**2)**3*(c + d*x**2))
    F = d**2*(e*x)**(m + 1)*(-A*d + B*c)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(c*e*(m + 1)*(-a*d + b*c)**3) + (e*x)**(m + 1)*(A*b - B*a)/(4*a*e*(a + b*x**2)**2*(-a*d + b*c)) + (e*x)**(m + 1)*(A*b*(-a*d*(7 - m) + b*c*(3 - m)) + B*a*(a*d*(3 - m) + b*c*(m + 1)))/(8*a**2*e*(a + b*x**2)*(-a*d + b*c)**2) + (e*x)**(m + 1)*(A*b*(a**2*d**2*(m**2 - 8*m + 15) - 2*a*b*c*d*(m**2 - 6*m + 5) + b**2*c**2*(m**2 - 4*m + 3)) + B*a*(-a**2*d**2*(m**2 - 4*m + 3) - 2*a*b*c*d*(-m**2 + 2*m + 3) + b**2*c**2*(1 - m**2)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(8*a**3*e*(m + 1)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_30():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)**3/(c + d*x**2)**2
    F = -b**3*(e*x)**(m + 5)*(A*d*(m + 5) - B*c*(m + 7))/(2*c*d**2*e**5*(m + 5)) - b**2*(e*x)**(m + 3)*(3*a*d*(A*d*(m + 3) - B*c*(m + 5)) - b*c*(A*d*(m + 5) - B*c*(m + 7)))/(2*c*d**3*e**3*(m + 3)) - b*(e*x)**(m + 1)*(3*a**2*d**2*(A*d*(m + 1) - B*c*(m + 3)) - 3*a*b*c*d*(A*d*(m + 3) - B*c*(m + 5)) + b**2*c**2*(A*d*(m + 5) - B*c*(m + 7)))/(2*c*d**4*e*(m + 1)) - (e*x)**(m + 1)*(a + b*x**2)**3*(-A*d + B*c)/(2*c*d*e*(c + d*x**2)) + (e*x)**(m + 1)*(-a*d + b*c)**2*(a*d*(A*d*(1 - m) + B*c*(m + 1)) + b*c*(A*d*(m + 5) - B*c*(m + 7)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(2*c**2*d**4*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_31():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)**2/(c + d*x**2)**2
    F = -b**2*(e*x)**(m + 3)*(A*d*(m + 3) - B*c*(m + 5))/(2*c*d**2*e**3*(m + 3)) - b*(e*x)**(m + 1)*(2*a*d*(A*d*(m + 1) - B*c*(m + 3)) - b*c*(A*d*(m + 3) - B*c*(m + 5)))/(2*c*d**3*e*(m + 1)) - (e*x)**(m + 1)*(a + b*x**2)**2*(-A*d + B*c)/(2*c*d*e*(c + d*x**2)) - (e*x)**(m + 1)*(-a*d + b*c)*(a*d*(A*d*(1 - m) + B*c*(m + 1)) + b*c*(A*d*(m + 3) - B*c*(m + 5)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(2*c**2*d**3*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_32():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)/(c + d*x**2)**2
    F = -B*(e*x)**(m + 1)*(a*d*(m + 1) - b*c*(m + 3))/(2*c*d**2*e*(m + 1)) - (e*x)**(m + 1)*(A + B*x**2)*(-a*d + b*c)/(2*c*d*e*(c + d*x**2)) + (e*x)**(m + 1)*(a*d*(A*d*(1 - m) + B*c*(m + 1)) + b*c*(A*d*(m + 1) - B*c*(m + 3)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(2*c**2*d**2*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_33():
    f = (e*x)**m*(A + B*x**2)/(c + d*x**2)**2
    F = -(e*x)**(m + 1)*(-A*d + B*c)/(2*c*d*e*(c + d*x**2)) + (e*x)**(m + 1)*(A*d*(1 - m) + B*c*(m + 1))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(2*c**2*d*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_34():
    f = (e*x)**m*(A + B*x**2)/((a + b*x**2)*(c + d*x**2)**2)
    F = (e*x)**(m + 1)*(-A*d + B*c)/(2*c*e*(c + d*x**2)*(-a*d + b*c)) + (e*x)**(m + 1)*(a*d*(A*d*(1 - m) + B*c*(m + 1)) + b*c*(-A*d*(3 - m) + B*c*(1 - m)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(2*c**2*e*(m + 1)*(-a*d + b*c)**2) + b*(e*x)**(m + 1)*(A*b - B*a)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a*e*(m + 1)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_35():
    f = (e*x)**m*(A + B*x**2)/((a + b*x**2)**2*(c + d*x**2)**2)
    F = -d*(e*x)**(m + 1)*(a*d*(A*d*(1 - m) + B*c*(m + 1)) + b*c*(-A*d*(5 - m) + B*c*(3 - m)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(2*c**2*e*(m + 1)*(-a*d + b*c)**3) + (e*x)**(m + 1)*(A*b - B*a)/(2*a*e*(a + b*x**2)*(c + d*x**2)*(-a*d + b*c)) + d*(e*x)**(m + 1)*(A*a*d + A*b*c - 2*B*a*c)/(2*a*c*e*(c + d*x**2)*(-a*d + b*c)**2) + b*(e*x)**(m + 1)*(A*b*(-a*d*(5 - m) + b*c*(1 - m)) + B*a*(a*d*(3 - m) + b*c*(m + 1)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(2*a**2*e*(m + 1)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_36():
    f = (e*x)**m*(A + B*x**2)/((a + b*x**2)**3*(c + d*x**2)**2)
    F = d**2*(e*x)**(m + 1)*(a*d*(A*d*(1 - m) + B*c*(m + 1)) + b*c*(-A*d*(7 - m) + B*c*(5 - m)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(2*c**2*e*(m + 1)*(-a*d + b*c)**4) + (e*x)**(m + 1)*(A*b - B*a)/(4*a*e*(a + b*x**2)**2*(c + d*x**2)*(-a*d + b*c)) + (e*x)**(m + 1)*(A*b*(-a*d*(9 - m) + b*c*(3 - m)) + B*a*(a*d*(5 - m) + b*c*(m + 1)))/(8*a**2*e*(a + b*x**2)*(c + d*x**2)*(-a*d + b*c)**2) - d*(e*x)**(m + 1)*(A*(4*a**2*d**2 + a*b*c*d*(11 - m) - b**2*c**2*(3 - m)) - B*a*c*(a*d*(11 - m) + b*c*(m + 1)))/(8*a**2*c*e*(c + d*x**2)*(-a*d + b*c)**3) + b*(e*x)**(m + 1)*(A*b*(a**2*d**2*(m**2 - 12*m + 35) - 2*a*b*c*d*(m**2 - 8*m + 7) + b**2*c**2*(m**2 - 4*m + 3)) + B*a*(-a**2*d**2*(m**2 - 8*m + 15) - 2*a*b*c*d*(-m**2 + 4*m + 5) + b**2*c**2*(1 - m**2)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(8*a**3*e*(m + 1)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_37():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)**3/(c + d*x**2)**3
    F = -b**2*(e*x)**(m + 3)*(a*d*(m + 3)*(A*d*(3 - m) + B*c*(m + 1)) + b*c*(m + 5)*(A*d*(m + 3) - B*c*(m + 7)))/(8*c**2*d**3*e**3*(m + 3)) - b*(e*x)**(m + 1)*(2*a**2*d**2*(m + 1)*(A*d*(3 - m) + B*c*(m + 1)) + 3*a*b*c*d*(m + 3)*(A*d*(m + 1) - B*c*(m + 5)) - b**2*c**2*(m + 5)*(A*d*(m + 3) - B*c*(m + 7)))/(8*c**2*d**4*e*(m + 1)) - (e*x)**(m + 1)*(a + b*x**2)**3*(-A*d + B*c)/(4*c*d*e*(c + d*x**2)**2) + (e*x)**(m + 1)*(a + b*x**2)**2*(a*d*(A*d*(3 - m) + B*c*(m + 1)) + b*c*(A*d*(m + 3) - B*c*(m + 7)))/(8*c**2*d**2*e*(c + d*x**2)) - (e*x)**(m + 1)*(-a*d + b*c)*(a**2*d**2*(1 - m)*(A*d*(3 - m) + B*c*(m + 1)) + 2*a*b*c*d*(A*d*(-m**2 - 2*m + 3) + B*c*(m**2 + 6*m + 5)) + b**2*c**2*(m + 5)*(A*d*(m + 3) - B*c*(m + 7)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(8*c**3*d**4*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_38():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)**2/(c + d*x**2)**3
    F = b*(e*x)**(m + 1)*(A*d*(m + 1) - B*c*(m + 5))*(a*d*(m + 1) - b*c*(m + 3))/(8*c**2*d**3*e*(m + 1)) - (e*x)**(m + 1)*(a + b*x**2)**2*(-A*d + B*c)/(4*c*d*e*(c + d*x**2)**2) - (e*x)**(m + 1)*(-a*d + b*c)*(a*(A*d*(3 - m) + B*c*(m + 1)) - b*x**2*(A*d*(m + 1) - B*c*(m + 5)))/(8*c**2*d**2*e*(c + d*x**2)) + (e*x)**(m + 1)*(a*d*(A*d*(3 - m) + B*c*(m + 1))*(a*d*(1 - m) + b*c*(m + 1)) - b*c*(A*d*(m + 1) - B*c*(m + 5))*(a*d*(m + 1) - b*c*(m + 3)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(8*c**3*d**3*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_39():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)/(c + d*x**2)**3
    F = -(e*x)**(m + 1)*(A + B*x**2)*(-a*d + b*c)/(4*c*d*e*(c + d*x**2)**2) + (e*x)**(m + 1)*(a*d*(A*d*(3 - m) - B*(-c*m + c)) + b*c*(A*d*(m + 1) - B*c*(m + 3)))/(8*c**2*d**2*e*(c + d*x**2)) + (e*x)**(m + 1)*(a*d*(1 - m)*(A*d*(3 - m) + B*c*(m + 1)) + b*c*(m + 1)*(A*d*(1 - m) + B*c*(m + 3)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(8*c**3*d**2*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_40():
    f = (e*x)**m*(A + B*x**2)/(c + d*x**2)**3
    F = -(e*x)**(m + 1)*(-A*d + B*c)/(4*c*d*e*(c + d*x**2)**2) + (e*x)**(m + 1)*(A*d*(3 - m) + B*c*(m + 1))*hyper((2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(4*c**3*d*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_41():
    f = (e*x)**m*(A + B*x**2)/((a + b*x**2)*(c + d*x**2)**3)
    F = (e*x)**(m + 1)*(-A*d + B*c)/(4*c*e*(c + d*x**2)**2*(-a*d + b*c)) + (e*x)**(m + 1)*(a*d*(A*d*(3 - m) + B*c*(m + 1)) + b*c*(-A*d*(7 - m) + B*c*(3 - m)))/(8*c**2*e*(c + d*x**2)*(-a*d + b*c)**2) + (e*x)**(m + 1)*(-a**2*d**2*(1 - m)*(A*d*(3 - m) + B*c*(m + 1)) + 2*a*b*c*d*(A*d*(m**2 - 6*m + 5) + B*c*(-m**2 + 2*m + 3)) + b**2*c**2*(3 - m)*(-A*d*(5 - m) + B*c*(1 - m)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(8*c**3*e*(m + 1)*(-a*d + b*c)**3) + b**2*(e*x)**(m + 1)*(A*b - B*a)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(a*e*(m + 1)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_42():
    f = (e*x)**m*(A + B*x**2)/((a + b*x**2)**2*(c + d*x**2)**3)
    F = -d*(e*x)**(m + 1)*(-a**2*d**2*(1 - m)*(A*d*(3 - m) + B*c*(m + 1)) + 2*a*b*c*d*(A*d*(m**2 - 8*m + 7) + B*c*(-m**2 + 4*m + 5)) + b**2*c**2*(5 - m)*(-A*d*(7 - m) + B*c*(3 - m)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(8*c**3*e*(m + 1)*(-a*d + b*c)**4) + (e*x)**(m + 1)*(A*b - B*a)/(2*a*e*(a + b*x**2)*(c + d*x**2)**2*(-a*d + b*c)) + d*(e*x)**(m + 1)*(A*a*d + 2*A*b*c - 3*B*a*c)/(4*a*c*e*(c + d*x**2)**2*(-a*d + b*c)**2) + d*(e*x)**(m + 1)*(A*(-a**2*d**2*(3 - m) + a*b*c*d*(11 - m) + 4*b**2*c**2) - B*a*c*(a*d*(m + 1) + b*c*(11 - m)))/(8*a*c**2*e*(c + d*x**2)*(-a*d + b*c)**3) + b**2*(e*x)**(m + 1)*(A*b*(-a*d*(7 - m) + b*c*(1 - m)) + B*a*(a*d*(5 - m) + b*c*(m + 1)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(2*a**2*e*(m + 1)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_43():
    f = (e*x)**m*(A + B*x**2)/((a + b*x**2)**3*(c + d*x**2)**3)
    F = d**2*(e*x)**(m + 1)*(-a**2*d**2*(1 - m)*(A*d*(3 - m) + B*c*(m + 1)) + 2*a*b*c*d*(A*d*(m**2 - 10*m + 9) + B*c*(-m**2 + 6*m + 7)) + b**2*c**2*(7 - m)*(-A*d*(9 - m) + B*c*(5 - m)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -d*x**2/c)/(8*c**3*e*(m + 1)*(-a*d + b*c)**5) + (e*x)**(m + 1)*(A*b - B*a)/(4*a*e*(a + b*x**2)**2*(c + d*x**2)**2*(-a*d + b*c)) + (e*x)**(m + 1)*(A*b*(-a*d*(11 - m) + b*c*(3 - m)) + B*a*(a*d*(7 - m) + b*c*(m + 1)))/(8*a**2*e*(a + b*x**2)*(c + d*x**2)**2*(-a*d + b*c)**2) - d*(e*x)**(m + 1)*(A*(2*a**2*d**2 + a*b*c*d*(13 - m) - b**2*c**2*(3 - m)) - B*a*c*(a*d*(11 - m) + b*c*(m + 1)))/(8*a**2*c*e*(c + d*x**2)**2*(-a*d + b*c)**3) + d*(e*x)**(m + 1)*(A*(a*d + b*c)*(a**2*d**2*(3 - m) - 2*a*b*c*d*(9 - m) + b**2*c**2*(3 - m)) + B*a*c*(a**2*d**2*(m + 1) + 2*a*b*c*d*(11 - m) + b**2*c**2*(m + 1)))/(8*a**2*c**2*e*(c + d*x**2)*(-a*d + b*c)**4) + b**2*(e*x)**(m + 1)*(A*b*(a**2*d**2*(m**2 - 16*m + 63) - 2*a*b*c*d*(m**2 - 10*m + 9) + b**2*c**2*(m**2 - 4*m + 3)) + B*a*(-a**2*d**2*(m**2 - 12*m + 35) - 2*a*b*c*d*(-m**2 + 6*m + 7) + b**2*c**2*(1 - m**2)))*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(8*a**3*e*(m + 1)*(-a*d + b*c)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_44():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)**p*(c + d*x**2)**3
    F = B*(e*x)**(m + 1)*(a + b*x**2)**(p + 1)*(c + d*x**2)**3/(b*e*(m + 2*p + 9)) - (e*x)**(m + 1)*(a + b*x**2)**(p + 1)*(c + d*x**2)**2*(B*a*d*(m + 7) - b*(A*d*(m + 2*p + 9) + 6*B*c))/(b**2*e*(m + 2*p + 7)*(m + 2*p + 9)) + (e*x)**(m + 1)*(a + b*x**2)**(p + 1)*(c + d*x**2)*(B*a**2*d**2*(m**2 + 12*m + 35) - a*b*d*(A*d*(m + 5)*(m + 2*p + 9) + B*c*(m**2 + 2*m*(p + 9) + 2*p + 65)) + b**2*c*(A*d*(m**2 + 4*m*(p + 5) + 4*p**2 + 40*p + 99) + 24*B*c))/(b**3*e*(m + 2*p + 5)*(m + 2*p + 7)*(m + 2*p + 9)) - (e*x)**(m + 1)*(a + b*x**2)**(p + 1)*(B*a**3*d**3*(m**3 + 15*m**2 + 71*m + 105) - a**2*b*d**2*(m + 5)*(A*d*(m + 3)*(m + 2*p + 9) + 2*B*c*(m**2 + 2*m*p + 13*m + 2*p + 30)) + a*b**2*c*d*(2*A*d*(m**3 + 4*m**2*(p + 5) + m*(4*p**2 + 44*p + 123) + 8*p**2 + 84*p + 216) + B*c*(m**3 + m**2*(4*p + 21) + m*(4*p**2 + 44*p + 143) + 4*p**2 + 40*p + 267)) - b**3*c**2*(A*d*(m**3 + m**2*(6*p + 23) + m*(12*p**2 + 92*p + 183) + 8*p**3 + 92*p**2 + 366*p + 513) + 48*B*c))/(b**4*e*(m + 2*p + 3)*(m + 2*p + 5)*(m + 2*p + 7)*(m + 2*p + 9)) + (e*x)**(m + 1)*(a + b*x**2)**p*(a*(m + 1)*(B*a**3*d**3*(m**3 + 15*m**2 + 71*m + 105) - a**2*b*d**2*(m + 5)*(A*d*(m + 3)*(m + 2*p + 9) + 2*B*c*(m**2 + 2*m*p + 13*m + 2*p + 30)) + a*b**2*c*d*(2*A*d*(m**3 + 4*m**2*(p + 5) + m*(4*p**2 + 44*p + 123) + 8*p**2 + 84*p + 216) + B*c*(m**3 + m**2*(4*p + 21) + m*(4*p**2 + 44*p + 143) + 4*p**2 + 40*p + 267)) - b**3*c**2*(A*d*(m**3 + m**2*(6*p + 23) + m*(12*p**2 + 92*p + 183) + 8*p**3 + 92*p**2 + 366*p + 513) + 48*B*c)) - b*c*(2*b*c*(p + 2)*(2*b*c*(p + 3)*(-A*b*(m + 2*p + 9) + B*a*(m + 1)) + (m + 1)*(-a*d + b*c)*(-A*b*(m + 2*p + 9) + B*a*(m + 7))) + (m + 1)*(-a*(2*b*c*d*(p + 3)*(-A*b*(m + 2*p + 9) + B*a*(m + 1)) + d*(m + 1)*(-a*d + b*c)*(-A*b*(m + 2*p + 9) + B*a*(m + 7)) + (-4*a*d + 4*b*c)*(B*a*d*(m + 7) - b*(A*d*(m + 2*p + 9) + 6*B*c))) + b*c*(2*b*c*(p + 3)*(-A*b*(m + 2*p + 9) + B*a*(m + 1)) + (m + 1)*(-a*d + b*c)*(-A*b*(m + 2*p + 9) + B*a*(m + 7)))))*(m + 2*p + 3))*hyper((-p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(b**4*e*(1 + b*x**2/a)**p*(m + 1)*(m + 2*p + 3)*(m + 2*p + 5)*(m + 2*p + 7)*(m + 2*p + 9))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_45():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)**p*(c + d*x**2)**2
    F = B*(e*x)**(m + 1)*(a + b*x**2)**(p + 1)*(c + d*x**2)**2/(b*e*(m + 2*p + 7)) - (e*x)**(m + 1)*(a + b*x**2)**(p + 1)*(c + d*x**2)*(B*a*d*(m + 5) - b*(A*d*(m + 2*p + 7) + 4*B*c))/(b**2*e*(m + 2*p + 5)*(m + 2*p + 7)) + (e*x)**(m + 1)*(a + b*x**2)**(p + 1)*(B*a**2*d**2*(m**2 + 8*m + 15) - a*b*d*(A*d*(m + 3)*(m + 2*p + 7) + B*c*(m**2 + 2*m*(p + 6) + 2*p + 27)) + b**2*c*(A*d*(m + 2*p + 7)**2 + 8*B*c))/(b**3*e*(m + 2*p + 3)*(m + 2*p + 5)*(m + 2*p + 7)) - (e*x)**(m + 1)*(a + b*x**2)**p*(-a*(m + 1)*(2*b*c*d*(p + 2)*(-A*b*(m + 2*p + 7) + B*a*(m + 1)) + d*(m + 1)*(-a*d + b*c)*(-A*b*(m + 2*p + 7) + B*a*(m + 5)) + (-2*a*d + 2*b*c)*(B*a*d*(m + 5) - b*(A*d*(m + 2*p + 7) + 4*B*c))) + b*c*(2*b*c*(p + 2)*(-A*b*(m + 2*p + 7) + B*a*(m + 1)) + (m + 1)*(-a*d + b*c)*(-A*b*(m + 2*p + 7) + B*a*(m + 5)))*(m + 2*p + 3))*hyper((-p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(b**3*e*(1 + b*x**2/a)**p*(m + 1)*(m + 2*p + 3)*(m + 2*p + 5)*(m + 2*p + 7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_46():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)**p*(c + d*x**2)
    F = d*(e*x)**(m + 1)*(A + B*x**2)*(a + b*x**2)**(p + 1)/(b*e*(m + 2*p + 5)) - (e*x)**(m + 1)*(a + b*x**2)**(p + 1)*(B*a*d*(m + 3) - b*(2*A*d + B*c*(m + 2*p + 5)))/(b**2*e*(m + 2*p + 3)*(m + 2*p + 5)) - (e*x)**(m + 1)*(a + b*x**2)**p*(A*b*(a*d*(m + 1) - b*c*(m + 2*p + 5))*(m + 2*p + 3) - a*(m + 1)*(B*a*d*(m + 3) - b*(2*A*d + B*c*(m + 2*p + 5))))*hyper((-p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(b**2*e*(1 + b*x**2/a)**p*(m + 1)*(m + 2*p + 3)*(m + 2*p + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_47():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)**p/(c + d*x**2)
    F = B*(e*x)**(m + 1)*(a + b*x**2)**p*hyper((-p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(d*e*(1 + b*x**2/a)**p*(m + 1)) - (e*x)**(m + 1)*(a + b*x**2)**p*(-A*d + B*c)*appellf1(m/2 + sympy.S.Half, 1, -p, m/2 + sympy.S(3)/2, -d*x**2/c, -b*x**2/a)/(c*d*e*(1 + b*x**2/a)**p*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_48():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)**p/(c + d*x**2)**2
    F = -b*(e*x)**(m + 1)*(a + b*x**2)**p*(-A*d + B*c)*(m + 2*p + 1)*hyper((-p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(2*c*d*e*(1 + b*x**2/a)**p*(m + 1)*(-a*d + b*c)) + (e*x)**(m + 1)*(a + b*x**2)**(p + 1)*(-A*d + B*c)/(2*c*e*(c + d*x**2)*(-a*d + b*c)) - (e*x)**(m + 1)*(a + b*x**2)**p*(a*d*(A*d*(1 - m) + B*c*(m + 1)) - b*c*(A*d*(-m - 2*p + 1) + B*c*(m + 2*p + 1)))*appellf1(m/2 + sympy.S.Half, 1, -p, m/2 + sympy.S(3)/2, -d*x**2/c, -b*x**2/a)/(2*c**2*d*e*(1 + b*x**2/a)**p*(m + 1)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_49():
    f = (e*x)**m*(A + B*x**2)*(a + b*x**2)**p/(c + d*x**2)**3
    F = -b*(e*x)**(m + 1)*(a + b*x**2)**p*(a*d*(A*d*(3 - m) + B*c*(m + 1)) + b*c*(-A*d*(-m - 2*p + 5) + B*c*(-m - 2*p + 1)))*(m + 2*p + 1)*hyper((-p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -b*x**2/a)/(8*c**2*d*e*(1 + b*x**2/a)**p*(m + 1)*(-a*d + b*c)**2) + (e*x)**(m + 1)*(a + b*x**2)**(p + 1)*(-A*d + B*c)/(4*c*e*(c + d*x**2)**2*(-a*d + b*c)) + (e*x)**(m + 1)*(a + b*x**2)**(p + 1)*(a*d*(A*d*(3 - m) + B*c*(m + 1)) + b*c*(-A*d*(-m - 2*p + 5) + B*c*(-m - 2*p + 1)))/(8*c**2*e*(c + d*x**2)*(-a*d + b*c)**2) + (e*x)**(m + 1)*(a + b*x**2)**p*(a**2*d**2*(1 - m)*(A*d*(3 - m) + B*c*(m + 1)) - 2*a*b*c*d*(A*d*(1 - m)*(-m - 2*p + 3) + B*c*(m + 1)*(-m - 2*p + 1)) + b**2*c**2*(A*d*(-m - 2*p + 3) + B*c*(m + 2*p + 1))*(-m - 2*p + 1))*appellf1(m/2 + sympy.S.Half, 1, -p, m/2 + sympy.S(3)/2, -d*x**2/c, -b*x**2/a)/(8*c**3*d*e*(1 + b*x**2/a)**p*(m + 1)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_50():
    f = (A + B*x**2)*sqrt(a + b*x**2)*(c + d*x**2)/x
    F = -A*sqrt(a)*c*atanh(sqrt(a + b*x**2)/sqrt(a)) + A*c*sqrt(a + b*x**2) - (a + b*x**2)**(sympy.S(3)/2)*(2*B*a*d - 3*B*b*d*x**2 - 5*b*(A*d + B*c))/(15*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_6_g_x_pow_m_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_e_plus_f_x_pow_2_pow_r_51():
    f = (A + B*x**2)*(a + b*x**2)*sqrt(c + d*x**2)/x
    F = -A*a*sqrt(c)*atanh(sqrt(c + d*x**2)/sqrt(c)) + A*a*sqrt(c + d*x**2) - (c + d*x**2)**(sympy.S(3)/2)*(2*B*b*c - 3*B*b*d*x**2 - d*(5*A*b + 5*B*a))/(15*d**2)
    assert integrate(f, x) == F

