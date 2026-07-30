"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.1 Binomial products/1.1.3 General/1.1.3.6 (g x)^m (a+b x^n)^p (c+d x^n)^q (e+f x^n)^r.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, a, b, c, d, e, m, n, p, q = symbols('A B a b c d e m n p q')

def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_1():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)**3*(c + d*x**n)
    F = A*a**3*c*(e*x)**(m + 1)/(e*(m + 1)) + B*b**3*d*x**(5*n + 1)*(e*x)**m/(m + 5*n + 1) + a**2*x**(n + 1)*(e*x)**m*(A*a*d + 3*A*b*c + B*a*c)/(m + n + 1) + a*x**(2*n + 1)*(e*x)**m*(3*A*b*(a*d + b*c) + B*a*(a*d + 3*b*c))/(m + 2*n + 1) + b**2*x**(4*n + 1)*(e*x)**m*(A*b*d + 3*B*a*d + B*b*c)/(m + 4*n + 1) + b*x**(3*n + 1)*(e*x)**m*(A*b*(3*a*d + b*c) + 3*B*a*(a*d + b*c))/(m + 3*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_2():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)**2*(c + d*x**n)
    F = A*a**2*c*(e*x)**(m + 1)/(e*(m + 1)) + B*b**2*d*x**(4*n + 1)*(e*x)**m/(m + 4*n + 1) + a*x**(n + 1)*(e*x)**m*(A*a*d + 2*A*b*c + B*a*c)/(m + n + 1) + b*x**(3*n + 1)*(e*x)**m*(A*b*d + 2*B*a*d + B*b*c)/(m + 3*n + 1) + x**(2*n + 1)*(e*x)**m*(A*b*(2*a*d + b*c) + B*a*(a*d + 2*b*c))/(m + 2*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_3():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)*(c + d*x**n)
    F = A*a*c*(e*x)**(m + 1)/(e*(m + 1)) + B*b*d*x**(3*n + 1)*(e*x)**m/(m + 3*n + 1) + x**(n + 1)*(e*x)**m*(A*a*d + A*b*c + B*a*c)/(m + n + 1) + x**(2*n + 1)*(e*x)**m*(A*b*d + B*a*d + B*b*c)/(m + 2*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_4():
    f = (e*x)**m*(A + B*x**n)*(c + d*x**n)
    F = A*c*(e*x)**(m + 1)/(e*(m + 1)) + B*d*x**(2*n + 1)*(e*x)**m/(m + 2*n + 1) + x**(n + 1)*(e*x)**m*(A*d + B*c)/(m + n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_5():
    f = (e*x)**m*(A + B*x**n)*(c + d*x**n)/(a + b*x**n)
    F = B*d*x**(n + 1)*(e*x)**m/(b*(m + n + 1)) + (e*x)**(m + 1)*(A*b*d - B*a*d + B*b*c)/(b**2*e*(m + 1)) + (e*x)**(m + 1)*(A*b - B*a)*(-a*d + b*c)*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(a*b**2*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_6():
    f = (e*x)**m*(A + B*x**n)*(c + d*x**n)/(a + b*x**n)**2
    F = (e*x)**(m + 1)*(c + d*x**n)*(A*b - B*a)/(a*b*e*n*(a + b*x**n)) - d*(e*x)**(m + 1)*(A*b*(m + 1) - B*a*(m + n + 1))/(a*b**2*e*n*(m + 1)) + (e*x)**(m + 1)*(a*d*(A*b*(m + 1) - B*a*(m + n + 1)) + b*c*(-A*b*(m - n + 1) + B*a*(m + 1)))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(a**2*b**2*e*n*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_7():
    f = (e*x)**m*(A + B*x**n)*(c + d*x**n)/(a + b*x**n)**3
    F = (e*x)**(m + 1)*(c + d*x**n)*(A*b - B*a)/(2*a*b*e*n*(a + b*x**n)**2) - (e*x)**(m + 1)*(A*b*(-a*d*(m - n + 1) + b*c*(m - 2*n + 1)) - B*a*(-a*d*(m + n + 1) + b*c*(m + 1)))/(2*a**2*b**2*e*n**2*(a + b*x**n)) - (e*x)**(m + 1)*(a*d*(m + 1)*(A*b*(m - n + 1) - B*a*(m + n + 1)) + b*c*(-A*b*(m - 2*n + 1) + B*a*(m + 1))*(m - n + 1))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(2*a**3*b**2*e*n**2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_8():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)**3*(c + d*x**n)**2
    F = A*a**3*c**2*(e*x)**(m + 1)/(e*(m + 1)) + B*b**3*d**2*x**(6*n + 1)*(e*x)**m/(m + 6*n + 1) + a**2*c*x**(n + 1)*(e*x)**m*(2*A*a*d + 3*A*b*c + B*a*c)/(m + n + 1) + a*x**(2*n + 1)*(e*x)**m*(A*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2) + B*a*c*(2*a*d + 3*b*c))/(m + 2*n + 1) + b**2*d*x**(5*n + 1)*(e*x)**m*(A*b*d + 3*B*a*d + 2*B*b*c)/(m + 5*n + 1) + b*x**(4*n + 1)*(e*x)**m*(3*B*a**2*d**2 + 3*a*b*d*(A*d + 2*B*c) + b**2*c*(2*A*d + B*c))/(m + 4*n + 1) + x**(3*n + 1)*(e*x)**m*(A*b*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2) + B*a*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2))/(m + 3*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_9():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)**2*(c + d*x**n)**2
    F = A*a**2*c**2*(e*x)**(m + 1)/(e*(m + 1)) + B*b**2*d**2*x**(5*n + 1)*(e*x)**m/(m + 5*n + 1) + a*c*x**(n + 1)*(e*x)**m*(2*A*(a*d + b*c) + B*a*c)/(m + n + 1) + b*d*x**(4*n + 1)*(e*x)**m*(A*b*d + 2*B*a*d + 2*B*b*c)/(m + 4*n + 1) + x**(2*n + 1)*(e*x)**m*(A*(a**2*d**2 + 4*a*b*c*d + b**2*c**2) + 2*B*a*c*(a*d + b*c))/(m + 2*n + 1) + x**(3*n + 1)*(e*x)**m*(B*a**2*d**2 + 2*a*b*d*(A*d + 2*B*c) + b**2*c*(2*A*d + B*c))/(m + 3*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_10():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)*(c + d*x**n)**2
    F = A*a*c**2*(e*x)**(m + 1)/(e*(m + 1)) + B*b*d**2*x**(4*n + 1)*(e*x)**m/(m + 4*n + 1) + c*x**(n + 1)*(e*x)**m*(2*A*a*d + A*b*c + B*a*c)/(m + n + 1) + d*x**(3*n + 1)*(e*x)**m*(A*b*d + B*a*d + 2*B*b*c)/(m + 3*n + 1) + x**(2*n + 1)*(e*x)**m*(a*d*(A*d + 2*B*c) + b*c*(2*A*d + B*c))/(m + 2*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_11():
    f = (e*x)**m*(A + B*x**n)*(c + d*x**n)**2
    F = A*c**2*(e*x)**(m + 1)/(e*(m + 1)) + B*d**2*x**(3*n + 1)*(e*x)**m/(m + 3*n + 1) + c*x**(n + 1)*(e*x)**m*(2*A*d + B*c)/(m + n + 1) + d*x**(2*n + 1)*(e*x)**m*(A*d + 2*B*c)/(m + 2*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_12():
    f = (e*x)**m*(A + B*x**n)*(c + d*x**n)**2/(a + b*x**n)
    F = B*d**2*x**(2*n + 1)*(e*x)**m/(b*(m + 2*n + 1)) + d*x**(n + 1)*(e*x)**m*(A*b*d - B*a*d + 2*B*b*c)/(b**2*(m + n + 1)) + (e*x)**(m + 1)*(B*a**2*d**2 - a*b*d*(A*d + 2*B*c) + b**2*c*(2*A*d + B*c))/(b**3*e*(m + 1)) + (e*x)**(m + 1)*(A*b - B*a)*(-a*d + b*c)**2*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(a*b**3*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_13():
    f = (e*x)**m*(A + B*x**n)*(c + d*x**n)**2/(a + b*x**n)**2
    F = (e*x)**(m + 1)*(c + d*x**n)**2*(A*b - B*a)/(a*b*e*n*(a + b*x**n)) - d**2*x**(n + 1)*(e*x)**m*(A*b*(m + n + 1) - B*a*(m + 2*n + 1))/(a*b**2*n*(m + n + 1)) - d*(e*x)**(m + 1)*(A*b*(-a*d*(m + n + 1) + 2*b*c*(m + 1)) - B*a*(-a*d*(m + 2*n + 1) + 2*b*c*(m + n + 1)))/(a*b**3*e*n*(m + 1)) - (e*x)**(m + 1)*(-a*d + b*c)*(A*b*(-a*d*(m + n + 1) + b*c*(m - n + 1)) - B*a*(-a*d*(m + 2*n + 1) + b*c*(m + 1)))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(a**2*b**3*e*n*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_14():
    f = (e*x)**m*(A + B*x**n)*(c + d*x**n)**2/(a + b*x**n)**3
    F = (e*x)**(m + 1)*(c + d*x**n)**2*(A*b - B*a)/(2*a*b*e*n*(a + b*x**n)**2) + (e*x)**(m + 1)*(-a*d + b*c)*(c*(-A*b*(m - 2*n + 1) + B*a*(m + 1)) - d*x**n*(A*b*(m + 1) - B*a*(m + 2*n + 1)))/(2*a**2*b**2*e*n**2*(a + b*x**n)) + d*(e*x)**(m + 1)*(A*b*(m + 1) - B*a*(m + 2*n + 1))*(-a*d*(m + n + 1) + b*c*(m + 1))/(2*a**2*b**3*e*n**2*(m + 1)) + (e*x)**(m + 1)*(-a*d*(A*b*(m + 1) - B*a*(m + 2*n + 1))*(-a*d*(m + n + 1) + b*c*(m + 1)) + b*c*(-A*b*(m - 2*n + 1) + B*a*(m + 1))*(a*d*(m + 1) - b*c*(m - n + 1)))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(2*a**3*b**3*e*n**2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_15():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)**3*(c + d*x**n)**3
    F = A*a**3*c**3*(e*x)**(m + 1)/(e*(m + 1)) + B*b**3*d**3*x**(7*n + 1)*(e*x)**m/(m + 7*n + 1) + a**2*c**2*x**(n + 1)*(e*x)**m*(3*A*(a*d + b*c) + B*a*c)/(m + n + 1) + 3*a*c*x**(2*n + 1)*(e*x)**m*(A*(a**2*d**2 + 3*a*b*c*d + b**2*c**2) + B*a*c*(a*d + b*c))/(m + 2*n + 1) + b**2*d**2*x**(6*n + 1)*(e*x)**m*(A*b*d + 3*B*a*d + 3*B*b*c)/(m + 6*n + 1) + 3*b*d*x**(5*n + 1)*(e*x)**m*(B*a**2*d**2 + a*b*d*(A*d + 3*B*c) + b**2*c*(A*d + B*c))/(m + 5*n + 1) + x**(3*n + 1)*(e*x)**m*(A*(a**3*d**3 + 9*a**2*b*c*d**2 + 9*a*b**2*c**2*d + b**3*c**3) + 3*B*a*c*(a**2*d**2 + 3*a*b*c*d + b**2*c**2))/(m + 3*n + 1) + x**(4*n + 1)*(e*x)**m*(B*a**3*d**3 + 3*a**2*b*d**2*(A*d + 3*B*c) + 9*a*b**2*c*d*(A*d + B*c) + b**3*c**2*(3*A*d + B*c))/(m + 4*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_16():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)**2*(c + d*x**n)**3
    F = A*a**2*c**3*(e*x)**(m + 1)/(e*(m + 1)) + B*b**2*d**3*x**(6*n + 1)*(e*x)**m/(m + 6*n + 1) + a*c**2*x**(n + 1)*(e*x)**m*(3*A*a*d + 2*A*b*c + B*a*c)/(m + n + 1) + b*d**2*x**(5*n + 1)*(e*x)**m*(A*b*d + 2*B*a*d + 3*B*b*c)/(m + 5*n + 1) + c*x**(2*n + 1)*(e*x)**m*(A*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2) + B*a*c*(3*a*d + 2*b*c))/(m + 2*n + 1) + d*x**(4*n + 1)*(e*x)**m*(B*a**2*d**2 + 2*a*b*d*(A*d + 3*B*c) + 3*b**2*c*(A*d + B*c))/(m + 4*n + 1) + x**(3*n + 1)*(e*x)**m*(a**2*d**2*(A*d + 3*B*c) + 6*a*b*c*d*(A*d + B*c) + b**2*c**2*(3*A*d + B*c))/(m + 3*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_17():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)*(c + d*x**n)**3
    F = A*a*c**3*(e*x)**(m + 1)/(e*(m + 1)) + B*b*d**3*x**(5*n + 1)*(e*x)**m/(m + 5*n + 1) + c**2*x**(n + 1)*(e*x)**m*(3*A*a*d + A*b*c + B*a*c)/(m + n + 1) + c*x**(2*n + 1)*(e*x)**m*(3*a*d*(A*d + B*c) + b*c*(3*A*d + B*c))/(m + 2*n + 1) + d**2*x**(4*n + 1)*(e*x)**m*(A*b*d + B*a*d + 3*B*b*c)/(m + 4*n + 1) + d*x**(3*n + 1)*(e*x)**m*(a*d*(A*d + 3*B*c) + 3*b*c*(A*d + B*c))/(m + 3*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_18():
    f = (e*x)**m*(A + B*x**n)*(c + d*x**n)**3
    F = A*c**3*(e*x)**(m + 1)/(e*(m + 1)) + B*d**3*x**(4*n + 1)*(e*x)**m/(m + 4*n + 1) + c**2*x**(n + 1)*(e*x)**m*(3*A*d + B*c)/(m + n + 1) + 3*c*d*x**(2*n + 1)*(e*x)**m*(A*d + B*c)/(m + 2*n + 1) + d**2*x**(3*n + 1)*(e*x)**m*(A*d + 3*B*c)/(m + 3*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_19():
    f = (e*x)**m*(A + B*x**n)*(c + d*x**n)**3/(a + b*x**n)
    F = B*d**3*x**(3*n + 1)*(e*x)**m/(b*(m + 3*n + 1)) + d**2*x**(2*n + 1)*(e*x)**m*(A*b*d - B*a*d + 3*B*b*c)/(b**2*(m + 2*n + 1)) + d*x**(n + 1)*(e*x)**m*(B*a**2*d**2 - a*b*d*(A*d + 3*B*c) + 3*b**2*c*(A*d + B*c))/(b**3*(m + n + 1)) - (e*x)**(m + 1)*(B*a**3*d**3 - a**2*b*d**2*(A*d + 3*B*c) + 3*a*b**2*c*d*(A*d + B*c) - b**3*c**2*(3*A*d + B*c))/(b**4*e*(m + 1)) + (e*x)**(m + 1)*(A*b - B*a)*(-a*d + b*c)**3*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(a*b**4*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_20():
    f = (e*x)**m*(A + B*x**n)*(c + d*x**n)**3/(a + b*x**n)**2
    F = (e*x)**(m + 1)*(c + d*x**n)**3*(A*b - B*a)/(a*b*e*n*(a + b*x**n)) - d**3*x**(2*n + 1)*(e*x)**m*(A*b*(m + 2*n + 1) - B*a*(m + 3*n + 1))/(a*b**2*n*(m + 2*n + 1)) - d**2*x**(n + 1)*(e*x)**m*(A*b*(-a*d*(m + 2*n + 1) + 3*b*c*(m + n + 1)) - B*a*(-a*d*(m + 3*n + 1) + 3*b*c*(m + 2*n + 1)))/(a*b**3*n*(m + n + 1)) - d*(e*x)**(m + 1)*(A*b*(a**2*d**2*(m + 2*n + 1) - 3*a*b*c*d*(m + n + 1) + 3*b**2*c**2*(m + 1)) - B*a*(a**2*d**2*(m + 3*n + 1) - 3*a*b*c*d*(m + 2*n + 1) + 3*b**2*c**2*(m + n + 1)))/(a*b**4*e*n*(m + 1)) - (e*x)**(m + 1)*(-a*d + b*c)**2*(A*b*(-a*d*(m + 2*n + 1) + b*c*(m - n + 1)) - B*a*(-a*d*(m + 3*n + 1) + b*c*(m + 1)))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(a**2*b**4*e*n*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_21():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)**4/(c + d*x**n)
    F = B*b**4*x**(4*n + 1)*(e*x)**m/(d*(m + 4*n + 1)) - b**3*x**(3*n + 1)*(e*x)**m*(-A*b*d - 4*B*a*d + B*b*c)/(d**2*(m + 3*n + 1)) + b**2*x**(2*n + 1)*(e*x)**m*(6*B*a**2*d**2 - 4*a*b*d*(-A*d + B*c) + b**2*c*(-A*d + B*c))/(d**3*(m + 2*n + 1)) + b*x**(n + 1)*(e*x)**m*(4*B*a**3*d**3 - 6*a**2*b*d**2*(-A*d + B*c) + 4*a*b**2*c*d*(-A*d + B*c) - b**3*c**2*(-A*d + B*c))/(d**4*(m + n + 1)) + (e*x)**(m + 1)*(B*a**4*d**4 - 4*a**3*b*d**3*(-A*d + B*c) + 6*a**2*b**2*c*d**2*(-A*d + B*c) - 4*a*b**3*c**2*d*(-A*d + B*c) + b**4*c**3*(-A*d + B*c))/(d**5*e*(m + 1)) - (e*x)**(m + 1)*(-A*d + B*c)*(-a*d + b*c)**4*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(c*d**5*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_22():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)**3/(c + d*x**n)
    F = B*b**3*x**(3*n + 1)*(e*x)**m/(d*(m + 3*n + 1)) - b**2*x**(2*n + 1)*(e*x)**m*(-A*b*d - 3*B*a*d + B*b*c)/(d**2*(m + 2*n + 1)) + b*x**(n + 1)*(e*x)**m*(3*B*a**2*d**2 - 3*a*b*d*(-A*d + B*c) + b**2*c*(-A*d + B*c))/(d**3*(m + n + 1)) + (e*x)**(m + 1)*(B*a**3*d**3 - 3*a**2*b*d**2*(-A*d + B*c) + 3*a*b**2*c*d*(-A*d + B*c) - b**3*c**2*(-A*d + B*c))/(d**4*e*(m + 1)) + (e*x)**(m + 1)*(-A*d + B*c)*(-a*d + b*c)**3*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(c*d**4*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_23():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)**2/(c + d*x**n)
    F = B*b**2*x**(2*n + 1)*(e*x)**m/(d*(m + 2*n + 1)) - b*x**(n + 1)*(e*x)**m*(-A*b*d - 2*B*a*d + B*b*c)/(d**2*(m + n + 1)) + (e*x)**(m + 1)*(B*a**2*d**2 - 2*a*b*d*(-A*d + B*c) + b**2*c*(-A*d + B*c))/(d**3*e*(m + 1)) - (e*x)**(m + 1)*(-A*d + B*c)*(-a*d + b*c)**2*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(c*d**3*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_24():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)/(c + d*x**n)
    F = B*b*x**(n + 1)*(e*x)**m/(d*(m + n + 1)) - (e*x)**(m + 1)*(-A*b*d - B*a*d + B*b*c)/(d**2*e*(m + 1)) + (e*x)**(m + 1)*(-A*d + B*c)*(-a*d + b*c)*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(c*d**2*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_25():
    f = (e*x)**m*(A + B*x**n)/(c + d*x**n)
    F = B*(e*x)**(m + 1)/(d*e*(m + 1)) - (e*x)**(m + 1)*(-A*d + B*c)*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(c*d*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_26():
    f = (e*x)**m*(A + B*x**n)/((a + b*x**n)*(c + d*x**n))
    F = (e*x)**(m + 1)*(-A*d + B*c)*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(c*e*(m + 1)*(-a*d + b*c)) + (e*x)**(m + 1)*(A*b - B*a)*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(a*e*(m + 1)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_27():
    f = (e*x)**m*(A + B*x**n)/((a + b*x**n)**2*(c + d*x**n))
    F = -d*(e*x)**(m + 1)*(-A*d + B*c)*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(c*e*(m + 1)*(-a*d + b*c)**2) + (e*x)**(m + 1)*(A*b - B*a)/(a*e*n*(a + b*x**n)*(-a*d + b*c)) + (e*x)**(m + 1)*(A*b*(a*d*(m - 2*n + 1) - b*c*(m - n + 1)) + B*a*(-a*d*(m - n + 1) + b*c*(m + 1)))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(a**2*e*n*(m + 1)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_28():
    f = (e*x)**m*(A + B*x**n)/((a + b*x**n)**3*(c + d*x**n))
    F = d**2*(e*x)**(m + 1)*(-A*d + B*c)*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(c*e*(m + 1)*(-a*d + b*c)**3) + (e*x)**(m + 1)*(A*b - B*a)/(2*a*e*n*(a + b*x**n)**2*(-a*d + b*c)) + (e*x)**(m + 1)*(A*b*(a*d*(m - 4*n + 1) - b*c*(m - 2*n + 1)) + B*a*(-a*d*(m - 2*n + 1) + b*c*(m + 1)))/(2*a**2*e*n**2*(a + b*x**n)*(-a*d + b*c)**2) + (e*x)**(m + 1)*(A*b*(a**2*d**2*(m**2 + m*(2 - 5*n) + 6*n**2 - 5*n + 1) - 2*a*b*c*d*(m**2 + m*(2 - 4*n) + 3*n**2 - 4*n + 1) + b**2*c**2*(m**2 + m*(2 - 3*n) + 2*n**2 - 3*n + 1)) + B*a*(-a**2*d**2*(m**2 + m*(2 - 3*n) + 2*n**2 - 3*n + 1) + 2*a*b*c*d*(m + 1)*(m - 2*n + 1) - b**2*c**2*(m + 1)*(m - n + 1)))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(2*a**3*e*n**2*(m + 1)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_29():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)**3/(c + d*x**n)**2
    F = -b**3*x**(2*n + 1)*(e*x)**m*(A*d*(m + 2*n + 1) - B*c*(m + 3*n + 1))/(c*d**2*n*(m + 2*n + 1)) - b**2*x**(n + 1)*(e*x)**m*(3*a*d*(A*d*(m + n + 1) - B*c*(m + 2*n + 1)) - b*c*(A*d*(m + 2*n + 1) - B*c*(m + 3*n + 1)))/(c*d**3*n*(m + n + 1)) - b*(e*x)**(m + 1)*(3*a**2*d**2*(A*d*(m + 1) - B*c*(m + n + 1)) - 3*a*b*c*d*(A*d*(m + n + 1) - B*c*(m + 2*n + 1)) + b**2*c**2*(A*d*(m + 2*n + 1) - B*c*(m + 3*n + 1)))/(c*d**4*e*n*(m + 1)) - (e*x)**(m + 1)*(a + b*x**n)**3*(-A*d + B*c)/(c*d*e*n*(c + d*x**n)) + (e*x)**(m + 1)*(-a*d + b*c)**2*(a*d*(-A*d*(m - n + 1) + B*c*(m + 1)) + b*c*(A*d*(m + 2*n + 1) - B*c*(m + 3*n + 1)))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(c**2*d**4*e*n*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_30():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)**2/(c + d*x**n)**2
    F = -b**2*x**(n + 1)*(e*x)**m*(A*d*(m + n + 1) - B*c*(m + 2*n + 1))/(c*d**2*n*(m + n + 1)) - b*(e*x)**(m + 1)*(2*a*d*(A*d*(m + 1) - B*c*(m + n + 1)) - b*c*(A*d*(m + n + 1) - B*c*(m + 2*n + 1)))/(c*d**3*e*n*(m + 1)) - (e*x)**(m + 1)*(a + b*x**n)**2*(-A*d + B*c)/(c*d*e*n*(c + d*x**n)) - (e*x)**(m + 1)*(-a*d + b*c)*(a*d*(-A*d*(m - n + 1) + B*c*(m + 1)) + b*c*(A*d*(m + n + 1) - B*c*(m + 2*n + 1)))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(c**2*d**3*e*n*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_31():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)/(c + d*x**n)**2
    F = -B*(e*x)**(m + 1)*(a*d*(m + 1) - b*c*(m + n + 1))/(c*d**2*e*n*(m + 1)) - (e*x)**(m + 1)*(A + B*x**n)*(-a*d + b*c)/(c*d*e*n*(c + d*x**n)) + (e*x)**(m + 1)*(A*d*(-a*d*(m - n + 1) + b*c*(m + 1)) + B*c*(a*d*(m + 1) - b*c*(m + n + 1)))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(c**2*d**2*e*n*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_32():
    f = (e*x)**m*(A + B*x**n)/(c + d*x**n)**2
    F = -(e*x)**(m + 1)*(-A*d + B*c)/(c*d*e*n*(c + d*x**n)) + (e*x)**(m + 1)*(-A*d*(m - n + 1) + B*c*(m + 1))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(c**2*d*e*n*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_33():
    f = (e*x)**m*(A + B*x**n)/((a + b*x**n)*(c + d*x**n)**2)
    F = (e*x)**(m + 1)*(-A*d + B*c)/(c*e*n*(c + d*x**n)*(-a*d + b*c)) + (e*x)**(m + 1)*(a*d*(-A*d*(m - n + 1) + B*c*(m + 1)) + b*c*(A*d*(m - 2*n + 1) - B*c*(m - n + 1)))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(c**2*e*n*(m + 1)*(-a*d + b*c)**2) + b*(e*x)**(m + 1)*(A*b - B*a)*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(a*e*(m + 1)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_34():
    f = (e*x)**m*(A + B*x**n)/((a + b*x**n)**2*(c + d*x**n)**2)
    F = -d*(e*x)**(m + 1)*(a*d*(-A*d*(m - n + 1) + B*c*(m + 1)) + b*c*(A*d*(m - 3*n + 1) - B*c*(m - 2*n + 1)))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(c**2*e*n*(m + 1)*(-a*d + b*c)**3) + (e*x)**(m + 1)*(A*b - B*a)/(a*e*n*(a + b*x**n)*(c + d*x**n)*(-a*d + b*c)) + d*(e*x)**(m + 1)*(A*a*d + A*b*c - 2*B*a*c)/(a*c*e*n*(c + d*x**n)*(-a*d + b*c)**2) + b*(e*x)**(m + 1)*(A*b*(a*d*(m - 3*n + 1) - b*c*(m - n + 1)) + B*a*(-a*d*(m - 2*n + 1) + b*c*(m + 1)))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(a**2*e*n*(m + 1)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_35():
    f = (e*x)**m*(A + B*x**n)/((a + b*x**n)**3*(c + d*x**n)**2)
    F = d**2*(e*x)**(m + 1)*(a*d*(-A*d*(m - n + 1) + B*c*(m + 1)) + b*c*(A*d*(m - 4*n + 1) - B*c*(m - 3*n + 1)))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(c**2*e*n*(m + 1)*(-a*d + b*c)**4) + (e*x)**(m + 1)*(A*b - B*a)/(2*a*e*n*(a + b*x**n)**2*(c + d*x**n)*(-a*d + b*c)) + (e*x)**(m + 1)*(A*b*(a*d*(m - 5*n + 1) - b*c*(m - 2*n + 1)) + B*a*(-a*d*(m - 3*n + 1) + b*c*(m + 1)))/(2*a**2*e*n**2*(a + b*x**n)*(c + d*x**n)*(-a*d + b*c)**2) + d*(e*x)**(m + 1)*(A*(-2*a**2*d**2*n + a*b*c*d*(m - 6*n + 1) - b**2*c**2*(m - 2*n + 1)) + B*a*c*(-a*d*(m - 6*n + 1) + b*c*(m + 1)))/(2*a**2*c*e*n**2*(c + d*x**n)*(-a*d + b*c)**3) + b*(e*x)**(m + 1)*(A*b*(a**2*d**2*(m**2 + m*(2 - 7*n) + 12*n**2 - 7*n + 1) - 2*a*b*c*d*(m**2 + m*(2 - 5*n) + 4*n**2 - 5*n + 1) + b**2*c**2*(m**2 + m*(2 - 3*n) + 2*n**2 - 3*n + 1)) + B*a*(-a**2*d**2*(m**2 + m*(2 - 5*n) + 6*n**2 - 5*n + 1) + 2*a*b*c*d*(m + 1)*(m - 3*n + 1) - b**2*c**2*(m + 1)*(m - n + 1)))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(2*a**3*e*n**2*(m + 1)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_36():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)**2/(c + d*x**n)**3
    F = b*(e*x)**(m + 1)*(A*d*(m + 1) - B*c*(m + 2*n + 1))*(a*d*(m + 1) - b*c*(m + n + 1))/(2*c**2*d**3*e*n**2*(m + 1)) - (e*x)**(m + 1)*(a + b*x**n)**2*(-A*d + B*c)/(2*c*d*e*n*(c + d*x**n)**2) - (e*x)**(m + 1)*(-a*d + b*c)*(a*(-A*d*(m - 2*n + 1) + B*c*(m + 1)) - b*x**n*(A*d*(m + 1) - B*c*(m + 2*n + 1)))/(2*c**2*d**2*e*n**2*(c + d*x**n)) + (e*x)**(m + 1)*(a*d*(-A*d*(m - 2*n + 1) + B*c*(m + 1))*(-a*d*(m - n + 1) + b*c*(m + 1)) - b*c*(A*d*(m + 1) - B*c*(m + 2*n + 1))*(a*d*(m + 1) - b*c*(m + n + 1)))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(2*c**3*d**3*e*n**2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_37():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)/(c + d*x**n)**3
    F = -(e*x)**(m + 1)*(A + B*x**n)*(-a*d + b*c)/(2*c*d*e*n*(c + d*x**n)**2) - (e*x)**(m + 1)*(a*d*(A*d*(m - 2*n + 1) - B*c*(m - n + 1)) - b*c*(A*d*(m + 1) - B*c*(m + n + 1)))/(2*c**2*d**2*e*n**2*(c + d*x**n)) - (e*x)**(m + 1)*(A*d*(-a*d*(m - 2*n + 1) + b*c*(m + 1))*(m - n + 1) + B*c*(m + 1)*(a*d*(m - n + 1) - b*c*(m + n + 1)))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(2*c**3*d**2*e*n**2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_38():
    f = (e*x)**m*(A + B*x**n)/(c + d*x**n)**3
    F = -(e*x)**(m + 1)*(-A*d + B*c)/(2*c*d*e*n*(c + d*x**n)**2) + (e*x)**(m + 1)*(-A*d*(m - 2*n + 1) + B*c*(m + 1))*hyper((2, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(2*c**3*d*e*n*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_39():
    f = (e*x)**m*(A + B*x**n)/((a + b*x**n)*(c + d*x**n)**3)
    F = (e*x)**(m + 1)*(-A*d + B*c)/(2*c*e*n*(c + d*x**n)**2*(-a*d + b*c)) + (e*x)**(m + 1)*(a*d*(-A*d*(m - 2*n + 1) + B*c*(m + 1)) + b*c*(A*d*(m - 4*n + 1) - B*c*(m - 2*n + 1)))/(2*c**2*e*n**2*(c + d*x**n)*(-a*d + b*c)**2) - (e*x)**(m + 1)*(-a**2*d**2*(-A*d*(m - 2*n + 1) + B*c*(m + 1))*(m - n + 1) + 2*a*b*c*d*(-A*d*(m**2 + m*(2 - 4*n) + 3*n**2 - 4*n + 1) + B*c*(m + 1)*(m - 2*n + 1)) + b**2*c**2*(A*d*(m - 3*n + 1) - B*c*(m - n + 1))*(m - 2*n + 1))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(2*c**3*e*n**2*(m + 1)*(-a*d + b*c)**3) + b**2*(e*x)**(m + 1)*(A*b - B*a)*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(a*e*(m + 1)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_40():
    f = (e*x)**m*(A + B*x**n)/((a + b*x**n)**2*(c + d*x**n)**3)
    F = d*(e*x)**(m + 1)*(-a**2*d**2*(-A*d*(m - 2*n + 1) + B*c*(m + 1))*(m - n + 1) + 2*a*b*c*d*(-A*d*(m**2 + m*(2 - 5*n) + 4*n**2 - 5*n + 1) + B*c*(m + 1)*(m - 3*n + 1)) + b**2*c**2*(A*d*(m - 4*n + 1) - B*c*(m - 2*n + 1))*(m - 3*n + 1))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(2*c**3*e*n**2*(m + 1)*(-a*d + b*c)**4) + (e*x)**(m + 1)*(A*b - B*a)/(a*e*n*(a + b*x**n)*(c + d*x**n)**2*(-a*d + b*c)) + d*(e*x)**(m + 1)*(A*a*d + 2*A*b*c - 3*B*a*c)/(2*a*c*e*n*(c + d*x**n)**2*(-a*d + b*c)**2) - d*(e*x)**(m + 1)*(-2*A*b**2*c**2*n + a**2*d*(-A*d*(m - 2*n + 1) + B*c*(m + 1)) - a*b*c*(-A*d + B*c)*(m - 6*n + 1))/(2*a*c**2*e*n**2*(c + d*x**n)*(-a*d + b*c)**3) + b**2*(e*x)**(m + 1)*(A*b*(a*d*(m - 4*n + 1) - b*c*(m - n + 1)) + B*a*(-a*d*(m - 3*n + 1) + b*c*(m + 1)))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(a**2*e*n*(m + 1)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_41():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)**p*(c + d*x**n)**q
    F = A*(e*x)**(m + 1)*(a + b*x**n)**p*(c + d*x**n)**q*appellf1((m + 1)/n, -p, -q, (m + n + 1)/n, -b*x**n/a, -d*x**n/c)/(e*(1 + b*x**n/a)**p*(1 + d*x**n/c)**q*(m + 1)) + B*x**(n + 1)*(e*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*appellf1((m + n + 1)/n, -p, -q, (m + 2*n + 1)/n, -b*x**n/a, -d*x**n/c)/((1 + b*x**n/a)**p*(1 + d*x**n/c)**q*(m + n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_42():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)**p*(c + d*x**n)
    F = d*(e*x)**(m + 1)*(A + B*x**n)*(a + b*x**n)**(p + 1)/(b*e*(m + n*(p + 2) + 1)) - (e*x)**(m + 1)*(a + b*x**n)**(p + 1)*(B*a*d*(m + n + 1) - b*(A*d*n + B*c*(m + n*(p + 2) + 1)))/(b**2*e*(m + n*(p + 2) + 1)*(m + n*p + n + 1)) - (e*x)**(m + 1)*(a + b*x**n)**p*(A*b*(a*d*(m + 1) - b*c*(m + n*(p + 2) + 1))*(m + n*p + n + 1) - a*(m + 1)*(B*a*d*(m + n + 1) - b*(A*d*n + B*c*(m + n*(p + 2) + 1))))*hyper((-p, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(b**2*e*(1 + b*x**n/a)**p*(m + 1)*(m + n*(p + 2) + 1)*(m + n*p + n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_43():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)**p/(c + d*x**n)
    F = B*(e*x)**(m + 1)*(a + b*x**n)**p*hyper((-p, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(d*e*(1 + b*x**n/a)**p*(m + 1)) - (e*x)**(m + 1)*(a + b*x**n)**p*(-A*d + B*c)*appellf1((m + 1)/n, 1, -p, (m + n + 1)/n, -d*x**n/c, -b*x**n/a)/(c*d*e*(1 + b*x**n/a)**p*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_44():
    f = (e*x)**m*(A + B*x**n)*(a + b*x**n)**p/(c + d*x**n)**2
    F = -b*(e*x)**(m + 1)*(a + b*x**n)**p*(-A*d + B*c)*(m + n*p + 1)*hyper((-p, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(c*d*e*n*(1 + b*x**n/a)**p*(m + 1)*(-a*d + b*c)) + (e*x)**(m + 1)*(a + b*x**n)**(p + 1)*(-A*d + B*c)/(c*e*n*(c + d*x**n)*(-a*d + b*c)) - (e*x)**(m + 1)*(a + b*x**n)**p*(a*d*(-A*d*(m - n + 1) + B*c*(m + 1)) + b*c*(A*d*(m - n*(1 - p) + 1) - B*c*(m + n*p + 1)))*appellf1((m + 1)/n, 1, -p, (m + n + 1)/n, -d*x**n/c, -b*x**n/a)/(c**2*d*e*n*(1 + b*x**n/a)**p*(m + 1)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_6_g_x_pow_m_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_e_plus_f_x_pow_n_pow_r_45():
    f = (-a + b*x**(n/2))**(-1 + 1/n)*(a + b*x**(n/2))**(-1 + 1/n)*(c + d*x**n)/x**2
    F = (-a + b*x**(n/2))**(1/n)*(a + b*x**(n/2))**(1/n)*(d/b**2 + c/a**2)/x - d*(-a + b*x**(n/2))**(1/n)*(a + b*x**(n/2))**(1/n)*hyper((-1/n, -1/n), (-(1 - n)/n,), b**2*x**n/a**2)/(b**2*x*(1 - b**2*x**n/a**2)**(1/n))
    assert integrate(f, x) == F

