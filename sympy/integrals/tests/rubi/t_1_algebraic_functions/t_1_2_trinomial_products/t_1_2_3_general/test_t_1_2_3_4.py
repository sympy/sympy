"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.2 Trinomial products/1.2.3 General/1.2.3.4 (f x)^m (d+e x^n)^q (a+b x^n+c x^(2 n))^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n, p, q = symbols('a b c d e f m n p q')

def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_1():
    f = (d + e*x**3)**5*(a + b*x**3 + c*x**6)
    F = a*d**5*x + c*e**5*x**22/22 + d**4*x**4*(5*a*e + b*d)/4 + d**3*x**7*(c*d**2 + 5*e*(2*a*e + b*d))/7 + d**2*e*x**10*(c*d**2 + 2*e*(a*e + b*d))/2 + 5*d*e**2*x**13*(2*c*d**2 + e*(a*e + 2*b*d))/13 + e**4*x**19*(b*e + 5*c*d)/19 + e**3*x**16*(10*c*d**2 + e*(a*e + 5*b*d))/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_2():
    f = (d + e*x**3)**4*(a + b*x**3 + c*x**6)
    F = a*d**4*x + c*e**4*x**19/19 + d**3*x**4*(4*a*e + b*d)/4 + d**2*x**7*(6*a*e**2 + 4*b*d*e + c*d**2)/7 + d*e*x**10*(2*c*d**2 + e*(2*a*e + 3*b*d))/5 + e**3*x**16*(b*e + 4*c*d)/16 + e**2*x**13*(6*c*d**2 + e*(a*e + 4*b*d))/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_3():
    f = (d + e*x**3)**3*(a + b*x**3 + c*x**6)
    F = a*d**3*x + c*e**3*x**16/16 + d**2*x**4*(3*a*e + b*d)/4 + d*x**7*(c*d**2 + 3*e*(a*e + b*d))/7 + e**2*x**13*(b*e + 3*c*d)/13 + e*x**10*(3*c*d**2 + e*(a*e + 3*b*d))/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_4():
    f = (d + e*x**3)**2*(a + b*x**3 + c*x**6)
    F = a*d**2*x + c*e**2*x**13/13 + d*x**4*(2*a*e + b*d)/4 + e*x**10*(b*e + 2*c*d)/10 + x**7*(c*d**2 + e*(a*e + 2*b*d))/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_5():
    f = (d + e*x**3)*(a + b*x**3 + c*x**6)
    F = a*d*x + c*e*x**10/10 + x**7*(b*e + c*d)/7 + x**4*(a*e + b*d)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_6():
    f = (a + b*x**3 + c*x**6)/(d + e*x**3)
    F = c*x**4/(4*e) - x*(-b*e + c*d)/e**2 + (a*e**2 - b*d*e + c*d**2)*log(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)/(3*d**(sympy.S(2)/3)*e**(sympy.S(7)/3)) - (a*e**2 - b*d*e + c*d**2)*log(d**(sympy.S(2)/3) - d**(sympy.S(1)/3)*e**(sympy.S(1)/3)*x + e**(sympy.S(2)/3)*x**2)/(6*d**(sympy.S(2)/3)*e**(sympy.S(7)/3)) - sqrt(3)*(a*e**2 - b*d*e + c*d**2)*atan(sqrt(3)*(d**(sympy.S(1)/3) - 2*e**(sympy.S(1)/3)*x)/(3*d**(sympy.S(1)/3)))/(3*d**(sympy.S(2)/3)*e**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_7():
    f = (a + b*x**3 + c*x**6)/(d + e*x**3)**2
    F = c*x/e**2 + x*(a*e**2 - b*d*e + c*d**2)/(3*d*e**2*(d + e*x**3)) - (4*c*d**2 - e*(2*a*e + b*d))*log(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)/(9*d**(sympy.S(5)/3)*e**(sympy.S(7)/3)) + (4*c*d**2 - e*(2*a*e + b*d))*log(d**(sympy.S(2)/3) - d**(sympy.S(1)/3)*e**(sympy.S(1)/3)*x + e**(sympy.S(2)/3)*x**2)/(18*d**(sympy.S(5)/3)*e**(sympy.S(7)/3)) + sqrt(3)*(4*c*d**2 - e*(2*a*e + b*d))*atan(sqrt(3)*(d**(sympy.S(1)/3) - 2*e**(sympy.S(1)/3)*x)/(3*d**(sympy.S(1)/3)))/(9*d**(sympy.S(5)/3)*e**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_8():
    f = (a + b*x**3 + c*x**6)/(d + e*x**3)**3
    F = x*(a*e**2 - b*d*e + c*d**2)/(6*d*e**2*(d + e*x**3)**2) - x*(7*c*d**2 - e*(5*a*e + b*d))/(18*d**2*e**2*(d + e*x**3)) + (2*c*d**2 + e*(5*a*e + b*d))*log(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)/(27*d**(sympy.S(8)/3)*e**(sympy.S(7)/3)) - (2*c*d**2 + e*(5*a*e + b*d))*log(d**(sympy.S(2)/3) - d**(sympy.S(1)/3)*e**(sympy.S(1)/3)*x + e**(sympy.S(2)/3)*x**2)/(54*d**(sympy.S(8)/3)*e**(sympy.S(7)/3)) - sqrt(3)*(2*c*d**2 + e*(5*a*e + b*d))*atan(sqrt(3)*(d**(sympy.S(1)/3) - 2*e**(sympy.S(1)/3)*x)/(3*d**(sympy.S(1)/3)))/(27*d**(sympy.S(8)/3)*e**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_9():
    f = x**8*(d + e*x**3)/(a + b*x**3 + c*x**6)
    F = e*x**6/(6*c) + x**3*(-b*e + c*d)/(3*c**2) - (a*c*e - b**2*e + b*c*d)*log(a + b*x**3 + c*x**6)/(6*c**3) - (3*a*b*c*e - 2*a*c**2*d - b**3*e + b**2*c*d)*atanh((b + 2*c*x**3)/sqrt(-4*a*c + b**2))/(3*c**3*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_10():
    f = x**5*(d + e*x**3)/(a + b*x**3 + c*x**6)
    F = e*x**3/(3*c) + (-b*e + c*d)*log(a + b*x**3 + c*x**6)/(6*c**2) + (2*a*c*e - b**2*e + b*c*d)*atanh((b + 2*c*x**3)/sqrt(-4*a*c + b**2))/(3*c**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_11():
    f = x**2*(d + e*x**3)/(a + b*x**3 + c*x**6)
    F = e*log(a + b*x**3 + c*x**6)/(6*c) - (-b*e + 2*c*d)*atanh((b + 2*c*x**3)/sqrt(-4*a*c + b**2))/(3*c*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_12():
    f = (d + e*x**3)/(x*(a + b*x**3 + c*x**6))
    F = d*log(x)/a - d*log(a + b*x**3 + c*x**6)/(6*a) + (-2*a*e + b*d)*atanh((b + 2*c*x**3)/sqrt(-4*a*c + b**2))/(3*a*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_13():
    f = (d + e*x**3)/(x**4*(a + b*x**3 + c*x**6))
    F = -d/(3*a*x**3) - (-a*e + b*d)*log(x)/a**2 + (-a*e + b*d)*log(a + b*x**3 + c*x**6)/(6*a**2) - (-a*b*e - 2*a*c*d + b**2*d)*atanh((b + 2*c*x**3)/sqrt(-4*a*c + b**2))/(3*a**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_14():
    f = x**4*(d + e*x**3)/(a + b*x**3 + c*x**6)
    F = e*x**2/(2*c) - 2**(sympy.S(1)/3)*(-b*e + c*d + (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(5)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*(-b*e + c*d + (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(5)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*sqrt(3)*(-b*e + c*d + (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(5)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*(-b*e + c*d - (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(5)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*(-b*e + c*d - (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(5)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*sqrt(3)*(-b*e + c*d - (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(5)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_15():
    f = x**3*(d + e*x**3)/(a + b*x**3 + c*x**6)
    F = e*x/c + 2**(sympy.S(2)/3)*(-b*e + c*d + (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(4)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(-b*e + c*d + (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(4)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*sqrt(3)*(-b*e + c*d + (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(4)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(-b*e + c*d - (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(4)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(-b*e + c*d - (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(4)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*sqrt(3)*(-b*e + c*d - (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(4)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_16():
    f = x*(d + e*x**3)/(a + b*x**3 + c*x**6)
    F = -2**(sympy.S(1)/3)*(e - (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(2)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*(e - (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(2)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*sqrt(3)*(e - (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(2)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*(e + (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(2)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*(e + (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(2)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*sqrt(3)*(e + (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(2)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_17():
    f = (d + e*x**3)/(a + b*x**3 + c*x**6)
    F = 2**(sympy.S(2)/3)*(e - (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(1)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(e - (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(1)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*sqrt(3)*(e - (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(1)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(e + (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(1)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(e + (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(1)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*sqrt(3)*(e + (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(1)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_18():
    f = (d + e*x**3)/(x**2*(a + b*x**3 + c*x**6))
    F = 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*(d - (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*a*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*(d - (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*a*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*sqrt(3)*c**(sympy.S(1)/3)*(d - (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*a*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*(d + (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*a*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*(d + (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*a*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*sqrt(3)*c**(sympy.S(1)/3)*(d + (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*a*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) - d/(a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_19():
    f = (d + e*x**3)/(x**3*(a + b*x**3 + c*x**6))
    F = -2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*(d - (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*a*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*(d - (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*a*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*sqrt(3)*c**(sympy.S(2)/3)*(d - (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*a*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*(d + (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*a*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*(d + (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*a*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*sqrt(3)*c**(sympy.S(2)/3)*(d + (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*a*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - d/(2*a*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_20():
    f = x**8*(1 - x**3)/(x**6 - x**3 + 1)
    F = -x**6/6 + log(x**6 - x**3 + 1)/6 - sqrt(3)*atan(sqrt(3)*(1 - 2*x**3)/3)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_21():
    f = x**5*(1 - x**3)/(x**6 - x**3 + 1)
    F = -x**3/3 - 2*sqrt(3)*atan(sqrt(3)*(1 - 2*x**3)/3)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_22():
    f = x**2*(1 - x**3)/(x**6 - x**3 + 1)
    F = -log(x**6 - x**3 + 1)/6 - sqrt(3)*atan(sqrt(3)*(1 - 2*x**3)/3)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_23():
    f = (1 - x**3)/(x*(x**6 - x**3 + 1))
    F = log(x) - log(x**6 - x**3 + 1)/6 + sqrt(3)*atan(sqrt(3)*(1 - 2*x**3)/3)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_24():
    f = (1 - x**3)/(x**4*(x**6 - x**3 + 1))
    F = 2*sqrt(3)*atan(sqrt(3)*(1 - 2*x**3)/3)/9 - 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_25():
    f = x**6*(1 - x**3)/(x**6 - x**3 + 1)
    F = -x**4/4 + 2**(sympy.S(2)/3)*(3 + sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 - sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 - sqrt(3)*I)**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(3 - sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 + sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 + sqrt(3)*I)**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(3 + sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 - 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 - sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 - sqrt(3)*I)**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(3 - sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 + 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 + sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 + sqrt(3)*I)**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half - sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 - sqrt(3)*I)**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(-sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half + sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 + sqrt(3)*I)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_26():
    f = x**4*(1 - x**3)/(x**6 - x**3 + 1)
    F = -x**2/2 + sqrt(3)*I*log(-2**(sympy.S(1)/3)*x + (1 - sqrt(3)*I)**(sympy.S(1)/3))/(9*(sympy.S.Half - sqrt(3)*I/2)**(sympy.S(1)/3)) - sqrt(3)*I*log(-2**(sympy.S(1)/3)*x + (1 + sqrt(3)*I)**(sympy.S(1)/3))/(9*(sympy.S.Half + sqrt(3)*I/2)**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*sqrt(3)*I*log(2**(sympy.S(2)/3)*x**2 + x*(2 - 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 - sqrt(3)*I)**(sympy.S(2)/3))/(18*(1 - sqrt(3)*I)**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*sqrt(3)*I*log(2**(sympy.S(2)/3)*x**2 + x*(2 + 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 + sqrt(3)*I)**(sympy.S(2)/3))/(18*(1 + sqrt(3)*I)**(sympy.S(1)/3)) + I*atan(sqrt(3)*(2*x/(sympy.S.Half - sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(3*(sympy.S.Half - sqrt(3)*I/2)**(sympy.S(1)/3)) - I*atan(sqrt(3)*(2*x/(sympy.S.Half + sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(3*(sympy.S.Half + sqrt(3)*I/2)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_27():
    f = x**3*(1 - x**3)/(x**6 - x**3 + 1)
    F = -x + sqrt(3)*I*log(-2**(sympy.S(1)/3)*x + (1 - sqrt(3)*I)**(sympy.S(1)/3))/(9*(sympy.S.Half - sqrt(3)*I/2)**(sympy.S(2)/3)) - sqrt(3)*I*log(-2**(sympy.S(1)/3)*x + (1 + sqrt(3)*I)**(sympy.S(1)/3))/(9*(sympy.S.Half + sqrt(3)*I/2)**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*sqrt(3)*I*log(2**(sympy.S(2)/3)*x**2 + x*(2 - 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 - sqrt(3)*I)**(sympy.S(2)/3))/(18*(1 - sqrt(3)*I)**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*sqrt(3)*I*log(2**(sympy.S(2)/3)*x**2 + x*(2 + 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 + sqrt(3)*I)**(sympy.S(2)/3))/(18*(1 + sqrt(3)*I)**(sympy.S(2)/3)) - I*atan(sqrt(3)*(2*x/(sympy.S.Half - sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(3*(sympy.S.Half - sqrt(3)*I/2)**(sympy.S(2)/3)) + I*atan(sqrt(3)*(2*x/(sympy.S.Half + sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(3*(sympy.S.Half + sqrt(3)*I/2)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_28():
    f = x*(1 - x**3)/(x**6 - x**3 + 1)
    F = -2**(sympy.S(1)/3)*(3 - sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 - sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 - sqrt(3)*I)**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*(3 + sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 + sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 + sqrt(3)*I)**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*(3 - sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 - 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 - sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 - sqrt(3)*I)**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*(3 + sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 + 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 + sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 + sqrt(3)*I)**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*(-sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half - sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 - sqrt(3)*I)**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*(sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half + sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 + sqrt(3)*I)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_29():
    f = (1 - x**3)/(x**6 - x**3 + 1)
    F = -2**(sympy.S(2)/3)*(3 - sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 - sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 - sqrt(3)*I)**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(3 + sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 + sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 + sqrt(3)*I)**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(3 - sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 - 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 - sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 - sqrt(3)*I)**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(3 + sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 + 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 + sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 + sqrt(3)*I)**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(-sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half - sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 - sqrt(3)*I)**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half + sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 + sqrt(3)*I)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_30():
    f = (1 - x**3)/(x**2*(x**6 - x**3 + 1))
    F = -2**(sympy.S(1)/3)*(3 + sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 - sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 - sqrt(3)*I)**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*(3 - sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 + sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 + sqrt(3)*I)**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*(3 + sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 - 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 - sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 - sqrt(3)*I)**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*(3 - sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 + 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 + sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 + sqrt(3)*I)**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*(sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half - sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 - sqrt(3)*I)**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*(-sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half + sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 + sqrt(3)*I)**(sympy.S(1)/3)) - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_31():
    f = (1 - x**3)/(x**3*(x**6 - x**3 + 1))
    F = -2**(sympy.S(2)/3)*(3 + sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 - sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 - sqrt(3)*I)**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(3 - sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 + sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 + sqrt(3)*I)**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(3 + sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 - 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 - sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 - sqrt(3)*I)**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(3 - sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 + 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 + sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 + sqrt(3)*I)**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half - sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 - sqrt(3)*I)**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(-sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half + sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 + sqrt(3)*I)**(sympy.S(2)/3)) - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_32():
    f = x**2*(x**3 - 2)/(x**6 - x**3 + 1)
    F = log(x**6 - x**3 + 1)/6 + sqrt(3)*atan(sqrt(3)*(1 - 2*x**3)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_33():
    f = (x**3 + 1)/(x*(x**6 - x**3 + 1))
    F = log(x) - log(x**6 - x**3 + 1)/6 - sqrt(3)*atan(sqrt(3)*(1 - 2*x**3)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_34():
    f = (x**3 + 1)/(x**7 - x**4 + x)
    F = log(x) - log(x**6 - x**3 + 1)/6 - sqrt(3)*atan(sqrt(3)*(1 - 2*x**3)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_35():
    f = (d + e*x**3)**(sympy.S(5)/2)*(a + b*x**3 + c*x**6)
    F = 2*c*x**4*(d + e*x**3)**(sympy.S(7)/2)/(29*e) + 54*3**(sympy.S(3)/4)*d**3*sqrt((d**(sympy.S(2)/3) - d**(sympy.S(1)/3)*e**(sympy.S(1)/3)*x + e**(sympy.S(2)/3)*x**2)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)*(667*a*e**2 - 58*b*d*e + 16*c*d**2)*elliptic_f(asin((d**(sympy.S(1)/3)*(1 - sqrt(3)) + e**(sympy.S(1)/3)*x)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(124729*e**(sympy.S(7)/3)*sqrt(d**(sympy.S(1)/3)*(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)**2)*sqrt(d + e*x**3)) + 54*d**2*x*sqrt(d + e*x**3)*(667*a*e**2 - 58*b*d*e + 16*c*d**2)/(124729*e**2) + 30*d*x*(d + e*x**3)**(sympy.S(3)/2)*(667*a*e**2 - 58*b*d*e + 16*c*d**2)/(124729*e**2) - x*(d + e*x**3)**(sympy.S(7)/2)*(-58*b*e + 16*c*d)/(667*e**2) + x*(d + e*x**3)**(sympy.S(5)/2)*(1334*a*e**2 - 116*b*d*e + 32*c*d**2)/(11339*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_36():
    f = (d + e*x**3)**(sympy.S(3)/2)*(a + b*x**3 + c*x**6)
    F = 2*c*x**4*(d + e*x**3)**(sympy.S(5)/2)/(23*e) + 18*3**(sympy.S(3)/4)*d**2*sqrt((d**(sympy.S(2)/3) - d**(sympy.S(1)/3)*e**(sympy.S(1)/3)*x + e**(sympy.S(2)/3)*x**2)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)*(391*a*e**2 - 46*b*d*e + 16*c*d**2)*elliptic_f(asin((d**(sympy.S(1)/3)*(1 - sqrt(3)) + e**(sympy.S(1)/3)*x)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(21505*e**(sympy.S(7)/3)*sqrt(d**(sympy.S(1)/3)*(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)**2)*sqrt(d + e*x**3)) + 18*d*x*sqrt(d + e*x**3)*(391*a*e**2 - 46*b*d*e + 16*c*d**2)/(21505*e**2) - x*(d + e*x**3)**(sympy.S(5)/2)*(-46*b*e + 16*c*d)/(391*e**2) + x*(d + e*x**3)**(sympy.S(3)/2)*(782*a*e**2 - 92*b*d*e + 32*c*d**2)/(4301*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_37():
    f = sqrt(d + e*x**3)*(a + b*x**3 + c*x**6)
    F = 2*c*x**4*(d + e*x**3)**(sympy.S(3)/2)/(17*e) + 2*3**(sympy.S(3)/4)*d*sqrt((d**(sympy.S(2)/3) - d**(sympy.S(1)/3)*e**(sympy.S(1)/3)*x + e**(sympy.S(2)/3)*x**2)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)*(187*a*e**2 - 34*b*d*e + 16*c*d**2)*elliptic_f(asin((d**(sympy.S(1)/3)*(1 - sqrt(3)) + e**(sympy.S(1)/3)*x)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(935*e**(sympy.S(7)/3)*sqrt(d**(sympy.S(1)/3)*(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)**2)*sqrt(d + e*x**3)) - x*(d + e*x**3)**(sympy.S(3)/2)*(-34*b*e + 16*c*d)/(187*e**2) + x*sqrt(d + e*x**3)*(374*a*e**2 - 68*b*d*e + 32*c*d**2)/(935*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_38():
    f = (a + b*x**3 + c*x**6)/sqrt(d + e*x**3)
    F = 2*c*x**4*sqrt(d + e*x**3)/(11*e) - x*sqrt(d + e*x**3)*(-22*b*e + 16*c*d)/(55*e**2) + 2*3**(sympy.S(3)/4)*sqrt((d**(sympy.S(2)/3) - d**(sympy.S(1)/3)*e**(sympy.S(1)/3)*x + e**(sympy.S(2)/3)*x**2)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)*(55*a*e**2 - 22*b*d*e + 16*c*d**2)*elliptic_f(asin((d**(sympy.S(1)/3)*(1 - sqrt(3)) + e**(sympy.S(1)/3)*x)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(165*e**(sympy.S(7)/3)*sqrt(d**(sympy.S(1)/3)*(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)**2)*sqrt(d + e*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_39():
    f = (a + b*x**3 + c*x**6)/(d + e*x**3)**(sympy.S(3)/2)
    F = 2*c*x*sqrt(d + e*x**3)/(5*e**2) + x*(2*a*e**2 - 2*b*d*e + 2*c*d**2)/(3*d*e**2*sqrt(d + e*x**3)) - 2*3**(sympy.S(3)/4)*sqrt((d**(sympy.S(2)/3) - d**(sympy.S(1)/3)*e**(sympy.S(1)/3)*x + e**(sympy.S(2)/3)*x**2)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)*(16*c*d**2 - 5*e*(a*e + 2*b*d))*elliptic_f(asin((d**(sympy.S(1)/3)*(1 - sqrt(3)) + e**(sympy.S(1)/3)*x)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(45*d*e**(sympy.S(7)/3)*sqrt(d**(sympy.S(1)/3)*(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)**2)*sqrt(d + e*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_40():
    f = (a + b*x**3 + c*x**6)/(d + e*x**3)**(sympy.S(5)/2)
    F = x*(2*a*e**2 - 2*b*d*e + 2*c*d**2)/(9*d*e**2*(d + e*x**3)**(sympy.S(3)/2)) - x*(-14*a*e**2 - 4*b*d*e + 22*c*d**2)/(27*d**2*e**2*sqrt(d + e*x**3)) + 2*3**(sympy.S(3)/4)*sqrt((d**(sympy.S(2)/3) - d**(sympy.S(1)/3)*e**(sympy.S(1)/3)*x + e**(sympy.S(2)/3)*x**2)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)*(16*c*d**2 + e*(7*a*e + 2*b*d))*elliptic_f(asin((d**(sympy.S(1)/3)*(1 - sqrt(3)) + e**(sympy.S(1)/3)*x)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(81*d**2*e**(sympy.S(7)/3)*sqrt(d**(sympy.S(1)/3)*(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)**2)*sqrt(d + e*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_41():
    f = (a + b*x**3 + c*x**6)/(d + e*x**3)**(sympy.S(7)/2)
    F = x*(2*a*e**2 - 2*b*d*e + 2*c*d**2)/(15*d*e**2*(d + e*x**3)**(sympy.S(5)/2)) - x*(-26*a*e**2 - 4*b*d*e + 34*c*d**2)/(135*d**2*e**2*(d + e*x**3)**(sympy.S(3)/2)) + x*(182*a*e**2 + 28*b*d*e + 32*c*d**2)/(405*d**3*e**2*sqrt(d + e*x**3)) + 2*3**(sympy.S(3)/4)*sqrt((d**(sympy.S(2)/3) - d**(sympy.S(1)/3)*e**(sympy.S(1)/3)*x + e**(sympy.S(2)/3)*x**2)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)*(91*a*e**2 + 14*b*d*e + 16*c*d**2)*elliptic_f(asin((d**(sympy.S(1)/3)*(1 - sqrt(3)) + e**(sympy.S(1)/3)*x)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1215*d**3*e**(sympy.S(7)/3)*sqrt(d**(sympy.S(1)/3)*(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)**2)*sqrt(d + e*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_42():
    f = (a + b*x**3 + c*x**6)/(d + e*x**3)**(sympy.S(9)/2)
    F = x*(2*a*e**2 - 2*b*d*e + 2*c*d**2)/(21*d*e**2*(d + e*x**3)**(sympy.S(7)/2)) - x*(-38*a*e**2 - 4*b*d*e + 46*c*d**2)/(315*d**2*e**2*(d + e*x**3)**(sympy.S(5)/2)) + x*(494*a*e**2 + 52*b*d*e + 32*c*d**2)/(2835*d**3*e**2*(d + e*x**3)**(sympy.S(3)/2)) + x*(494*a*e**2 + 52*b*d*e + 32*c*d**2)/(1215*d**4*e**2*sqrt(d + e*x**3)) + 2*3**(sympy.S(3)/4)*sqrt((d**(sympy.S(2)/3) - d**(sympy.S(1)/3)*e**(sympy.S(1)/3)*x + e**(sympy.S(2)/3)*x**2)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)*(247*a*e**2 + 26*b*d*e + 16*c*d**2)*elliptic_f(asin((d**(sympy.S(1)/3)*(1 - sqrt(3)) + e**(sympy.S(1)/3)*x)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3645*d**4*e**(sympy.S(7)/3)*sqrt(d**(sympy.S(1)/3)*(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)/(d**(sympy.S(1)/3)*(1 + sqrt(3)) + e**(sympy.S(1)/3)*x)**2)*sqrt(d + e*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_43():
    f = x**4*(d + e*x**4)/(a + b*x**4 + c*x**8)
    F = e*x/c - 2**(sympy.S(3)/4)*(-b*e + c*d - (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(5)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) - 2**(sympy.S(3)/4)*(-b*e + c*d - (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(5)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) - 2**(sympy.S(3)/4)*(-b*e + c*d + (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(5)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) - 2**(sympy.S(3)/4)*(-b*e + c*d + (2*a*c*e - b**2*e + b*c*d)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(5)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_44():
    f = x**3*(d + e*x**4)/(a + b*x**4 + c*x**8)
    F = e*log(a + b*x**4 + c*x**8)/(8*c) - (-b*e + 2*c*d)*atanh((b + 2*c*x**4)/sqrt(-4*a*c + b**2))/(4*c*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_45():
    f = x**2*(d + e*x**4)/(a + b*x**4 + c*x**8)
    F = 2**(sympy.S(1)/4)*(e + (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) - 2**(sympy.S(1)/4)*(e + (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) + 2**(sympy.S(1)/4)*(e - (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) - 2**(sympy.S(1)/4)*(e - (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_46():
    f = x*(d + e*x**4)/(a + b*x**4 + c*x**8)
    F = sqrt(2)*(e - (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x**2/sqrt(b + sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*(e + (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x**2/sqrt(b - sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_47():
    f = (d + e*x**4)/(a + b*x**4 + c*x**8)
    F = -2**(sympy.S(3)/4)*(e + (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(1)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) - 2**(sympy.S(3)/4)*(e + (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(1)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) - 2**(sympy.S(3)/4)*(e - (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(1)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) - 2**(sympy.S(3)/4)*(e - (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(1)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_48():
    f = (d + e*x**4)/(x*(a + b*x**4 + c*x**8))
    F = d*log(x)/a - d*log(a + b*x**4 + c*x**8)/(8*a) + (-2*a*e + b*d)*atanh((b + 2*c*x**4)/sqrt(-4*a*c + b**2))/(4*a*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_49():
    f = (d + e*x**4)/(x**2*(a + b*x**4 + c*x**8))
    F = -2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(d + (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*a*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) + 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(d + (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*a*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) - 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(d - (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*a*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) + 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(d - (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*a*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) - d/(a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_50():
    f = (d + e*x**4)/(x**3*(a + b*x**4 + c*x**8))
    F = -sqrt(2)*sqrt(c)*(d - (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x**2/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(b + sqrt(-4*a*c + b**2))) - sqrt(2)*sqrt(c)*(d + (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x**2/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*sqrt(b - sqrt(-4*a*c + b**2))) - d/(2*a*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_51():
    f = (d + e*x**4)/(x**4*(a + b*x**4 + c*x**8))
    F = 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(d + (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*a*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(d + (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*a*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(d - (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*a*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(d - (-2*a*e + b*d)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*a*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) - d/(3*a*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_52():
    f = x**4*(1 - x**4)/(x**8 - x**4 + 1)
    F = -x - sqrt(6)*log(x**2 - x*sqrt(2 - sqrt(3)) + 1)/24 + sqrt(6)*log(x**2 + x*sqrt(2 - sqrt(3)) + 1)/24 - sqrt(6)*log(x**2 - x*sqrt(sqrt(3) + 2) + 1)/24 + sqrt(6)*log(x**2 + x*sqrt(sqrt(3) + 2) + 1)/24 - sqrt(6)*atan((-2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/12 + sqrt(6)*atan((2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/12 - sqrt(6)*atan((-2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/12 + sqrt(6)*atan((2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_53():
    f = x**3*(1 - x**4)/(x**8 - x**4 + 1)
    F = -log(x**8 - x**4 + 1)/8 - sqrt(3)*atan(sqrt(3)*(1 - 2*x**4)/3)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_54():
    f = x**2*(1 - x**4)/(x**8 - x**4 + 1)
    F = sqrt(sympy.S(2)/3 - sqrt(3)/3)*log(x**2 - x*sqrt(2 - sqrt(3)) + 1)/8 - sqrt(sympy.S(2)/3 - sqrt(3)/3)*log(x**2 + x*sqrt(2 - sqrt(3)) + 1)/8 - sqrt(sqrt(3)/3 + sympy.S(2)/3)*log(x**2 - x*sqrt(sqrt(3) + 2) + 1)/8 + sqrt(sqrt(3)/3 + sympy.S(2)/3)*log(x**2 + x*sqrt(sqrt(3) + 2) + 1)/8 - atan((-2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/(4*sqrt(3*sqrt(3) + 6)) + atan((2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/(4*sqrt(3*sqrt(3) + 6)) + atan((-2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/(4*sqrt(6 - 3*sqrt(3))) - atan((2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/(4*sqrt(6 - 3*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_55():
    f = x*(1 - x**4)/(x**8 - x**4 + 1)
    F = -sqrt(3)*log(x**4 - sqrt(3)*x**2 + 1)/12 + sqrt(3)*log(x**4 + sqrt(3)*x**2 + 1)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_56():
    f = (1 - x**4)/(x**8 - x**4 + 1)
    F = sqrt(sympy.S(2)/3 - sqrt(3)/3)*log(x**2 - x*sqrt(2 - sqrt(3)) + 1)/8 - sqrt(sympy.S(2)/3 - sqrt(3)/3)*log(x**2 + x*sqrt(2 - sqrt(3)) + 1)/8 - sqrt(sqrt(3)/3 + sympy.S(2)/3)*log(x**2 - x*sqrt(sqrt(3) + 2) + 1)/8 + sqrt(sqrt(3)/3 + sympy.S(2)/3)*log(x**2 + x*sqrt(sqrt(3) + 2) + 1)/8 + atan((-2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/(4*sqrt(3*sqrt(3) + 6)) - atan((2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/(4*sqrt(3*sqrt(3) + 6)) - atan((-2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/(4*sqrt(6 - 3*sqrt(3))) + atan((2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/(4*sqrt(6 - 3*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_57():
    f = (1 - x**4)/(x*(x**8 - x**4 + 1))
    F = log(x) - log(x**8 - x**4 + 1)/8 + sqrt(3)*atan(sqrt(3)*(1 - 2*x**4)/3)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_58():
    f = (1 - x**4)/(x**2*(x**8 - x**4 + 1))
    F = -sqrt(6)*log(x**2 - x*sqrt(2 - sqrt(3)) + 1)/24 + sqrt(6)*log(x**2 + x*sqrt(2 - sqrt(3)) + 1)/24 - sqrt(6)*log(x**2 - x*sqrt(sqrt(3) + 2) + 1)/24 + sqrt(6)*log(x**2 + x*sqrt(sqrt(3) + 2) + 1)/24 + sqrt(6)*atan((-2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/12 - sqrt(6)*atan((2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/12 + sqrt(6)*atan((-2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/12 - sqrt(6)*atan((2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/12 - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_59():
    f = (1 - x**4)/(x**3*(x**8 - x**4 + 1))
    F = -sqrt(3)*log(x**4 - sqrt(3)*x**2 + 1)/24 + sqrt(3)*log(x**4 + sqrt(3)*x**2 + 1)/24 - atan(2*x**2 - sqrt(3))/4 - atan(2*x**2 + sqrt(3))/4 - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_60():
    f = (1 - x**4)/(x**4*(x**8 - x**4 + 1))
    F = sqrt(sqrt(3)/3 + sympy.S(2)/3)*log(x**2 - x*sqrt(2 - sqrt(3)) + 1)/8 - sqrt(sqrt(3)/3 + sympy.S(2)/3)*log(x**2 + x*sqrt(2 - sqrt(3)) + 1)/8 - sqrt(sympy.S(2)/3 - sqrt(3)/3)*log(x**2 - x*sqrt(sqrt(3) + 2) + 1)/8 + sqrt(sympy.S(2)/3 - sqrt(3)/3)*log(x**2 + x*sqrt(sqrt(3) + 2) + 1)/8 + sqrt(sqrt(3)/3 + sympy.S(2)/3)*atan((-2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/4 - sqrt(sqrt(3)/3 + sympy.S(2)/3)*atan((2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/4 - sqrt(sympy.S(2)/3 - sqrt(3)/3)*atan((-2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/4 + sqrt(sympy.S(2)/3 - sqrt(3)/3)*atan((2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/4 - 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_61():
    f = x**3/((d + e*x)*(a + b/x + c/x**2))
    F = -d**5*log(d + e*x)/(e**4*(a*d**2 - e*(b*d - c*e))) + x**3/(3*a*e) - x**2*(a*d + b*e)/(2*a**2*e**2) + x*(a**2*d**2 + a*e*(b*d - c*e) + b**2*e**2)/(a**3*e**3) + (a**2*c**2*d - 3*a*b**2*c*d + 2*a*b*c**2*e + b**4*d - b**3*c*e)*log(a*x**2 + b*x + c)/(2*a**4*(a*d**2 - e*(b*d - c*e))) + (5*a**2*b*c**2*d - 2*a**2*c**3*e - 5*a*b**3*c*d + 4*a*b**2*c**2*e + b**5*d - b**4*c*e)*atanh((2*a*x + b)/sqrt(-4*a*c + b**2))/(a**4*sqrt(-4*a*c + b**2)*(a*d**2 - e*(b*d - c*e)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_62():
    f = x**2/((d + e*x)*(a + b/x + c/x**2))
    F = d**4*log(d + e*x)/(e**3*(a*d**2 - e*(b*d - c*e))) + x**2/(2*a*e) - x*(a*d + b*e)/(a**2*e**2) - (-2*a*b*c*d + a*c**2*e + b**3*d - b**2*c*e)*log(a*x**2 + b*x + c)/(2*a**3*(a*d**2 - e*(b*d - c*e))) - (2*a**2*c**2*d - 4*a*b**2*c*d + 3*a*b*c**2*e + b**4*d - b**3*c*e)*atanh((2*a*x + b)/sqrt(-4*a*c + b**2))/(a**3*sqrt(-4*a*c + b**2)*(a*d**2 - e*(b*d - c*e)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_63():
    f = x/((d + e*x)*(a + b/x + c/x**2))
    F = -d**3*log(d + e*x)/(e**2*(a*d**2 - e*(b*d - c*e))) + x/(a*e) + (-a*c*d + b**2*d - b*c*e)*log(a*x**2 + b*x + c)/(2*a**2*(a*d**2 - e*(b*d - c*e))) + (-3*a*b*c*d + 2*a*c**2*e + b**3*d - b**2*c*e)*atanh((2*a*x + b)/sqrt(-4*a*c + b**2))/(a**2*sqrt(-4*a*c + b**2)*(a*d**2 - e*(b*d - c*e)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_64():
    f = 1/((d + e*x)*(a + b/x + c/x**2))
    F = d**2*log(d + e*x)/(e*(a*d**2 - b*d*e + c*e**2)) - (b*d - c*e)*log(a*x**2 + b*x + c)/(2*a*(a*d**2 - e*(b*d - c*e))) - (-2*a*c*d + b**2*d - b*c*e)*atanh((2*a*x + b)/sqrt(-4*a*c + b**2))/(a*sqrt(-4*a*c + b**2)*(a*d**2 - e*(b*d - c*e)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_65():
    f = 1/(x*(d + e*x)*(a + b/x + c/x**2))
    F = d*log(a*x**2 + b*x + c)/(2*a*d**2 - 2*e*(b*d - c*e)) - d*log(d + e*x)/(a*d**2 - e*(b*d - c*e)) + (b*d - 2*c*e)*atanh((2*a*x + b)/sqrt(-4*a*c + b**2))/(sqrt(-4*a*c + b**2)*(a*d**2 - e*(b*d - c*e)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_66():
    f = 1/(x**2*(d + e*x)*(a + b/x + c/x**2))
    F = -e*log(a*x**2 + b*x + c)/(2*a*d**2 - 2*b*d*e + 2*c*e**2) + e*log(d + e*x)/(a*d**2 - b*d*e + c*e**2) - (2*a*d - b*e)*atanh((2*a*x + b)/sqrt(-4*a*c + b**2))/(sqrt(-4*a*c + b**2)*(a*d**2 - e*(b*d - c*e)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_67():
    f = 1/(x**3*(d + e*x)*(a + b/x + c/x**2))
    F = -e**2*log(d + e*x)/(d*(a*d**2 - b*d*e + c*e**2)) - (a*d - b*e)*log(a*x**2 + b*x + c)/(2*c*(a*d**2 - e*(b*d - c*e))) + (a*b*d + 2*a*c*e - b**2*e)*atanh((2*a*x + b)/sqrt(-4*a*c + b**2))/(c*sqrt(-4*a*c + b**2)*(a*d**2 - e*(b*d - c*e))) + log(x)/(c*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_68():
    f = 1/(x**4*(d + e*x)*(a + b/x + c/x**2))
    F = e**3*log(d + e*x)/(d**2*(a*d**2 - e*(b*d - c*e))) - 1/(c*d*x) + (a*b*d + a*c*e - b**2*e)*log(a*x**2 + b*x + c)/(2*c**2*(a*d**2 - e*(b*d - c*e))) + (2*a**2*c*d - a*b*(b*d + 3*c*e) + b**3*e)*atanh((2*a*x + b)/sqrt(-4*a*c + b**2))/(c**2*sqrt(-4*a*c + b**2)*(a*d**2 - e*(b*d - c*e))) - (b*d + c*e)*log(x)/(c**2*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_69():
    f = 1/(x**5*(d + e*x)*(a + b/x + c/x**2))
    F = -e**4*log(d + e*x)/(d**3*(a*d**2 - e*(b*d - c*e))) - 1/(2*c*d*x**2) + (b*d + c*e)/(c**2*d**2*x) + (a**2*c*d - a*b*(b*d + 2*c*e) + b**3*e)*log(a*x**2 + b*x + c)/(2*c**3*(a*d**2 - e*(b*d - c*e))) - (a**2*c*(3*b*d + 2*c*e) - a*b**2*(b*d + 4*c*e) + b**4*e)*atanh((2*a*x + b)/sqrt(-4*a*c + b**2))/(c**3*sqrt(-4*a*c + b**2)*(a*d**2 - e*(b*d - c*e))) + (b**2*d**2 + b*c*d*e - c*(a*d**2 - c*e**2))*log(x)/(c**3*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_70():
    f = x**3/((d + e*x)**2*(a + b/x + c/x**2))
    F = d**5/(e**4*(d + e*x)*(a*d**2 - e*(b*d - c*e))) + d**4*(3*a*d**2 - e*(4*b*d - 5*c*e))*log(d + e*x)/(e**4*(a*d**2 - e*(b*d - c*e))**2) + x**2/(2*a*e**2) - x*(2*a*d + b*e)/(a**2*e**3) + (4*a*b*c**2*d*e + a*c**2*(a*d**2 - c*e**2) + b**4*d**2 - 2*b**3*c*d*e - b**2*c*(3*a*d**2 - c*e**2))*log(a*x**2 + b*x + c)/(2*a**3*(a*d**2 - e*(b*d - c*e))**2) + (-4*a**2*c**3*d*e + 8*a*b**2*c**2*d*e + a*b*c**2*(5*a*d**2 - 3*c*e**2) + b**5*d**2 - 2*b**4*c*d*e - b**3*c*(5*a*d**2 - c*e**2))*atanh((2*a*x + b)/sqrt(-4*a*c + b**2))/(a**3*sqrt(-4*a*c + b**2)*(a*d**2 - e*(b*d - c*e))**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_71():
    f = x**2/((d + e*x)**2*(a + b/x + c/x**2))
    F = -d**4/(e**3*(d + e*x)*(a*d**2 - e*(b*d - c*e))) - d**3*(2*a*d**2 - e*(3*b*d - 4*c*e))*log(d + e*x)/(e**3*(a*d**2 - e*(b*d - c*e))**2) + x/(a*e**2) - (b*d - c*e)*(-2*a*c*d + b**2*d - b*c*e)*log(a*x**2 + b*x + c)/(2*a**2*(a*d**2 - e*(b*d - c*e))**2) - (6*a*b*c**2*d*e + 2*a*c**2*(a*d**2 - c*e**2) + b**4*d**2 - 2*b**3*c*d*e - b**2*c*(4*a*d**2 - c*e**2))*atanh((2*a*x + b)/sqrt(-4*a*c + b**2))/(a**2*sqrt(-4*a*c + b**2)*(a*d**2 - e*(b*d - c*e))**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_72():
    f = x/((d + e*x)**2*(a + b/x + c/x**2))
    F = d**3/(e**2*(d + e*x)*(a*d**2 - e*(b*d - c*e))) + d**2*(a*d**2 - e*(2*b*d - 3*c*e))*log(d + e*x)/(e**2*(a*d**2 - e*(b*d - c*e))**2) + (b**2*d**2 - 2*b*c*d*e - c*(a*d**2 - c*e**2))*log(a*x**2 + b*x + c)/(2*a*(a*d**2 - e*(b*d - c*e))**2) + (4*a*c**2*d*e + b**3*d**2 - 2*b**2*c*d*e - b*c*(3*a*d**2 - c*e**2))*atanh((2*a*x + b)/sqrt(-4*a*c + b**2))/(a*sqrt(-4*a*c + b**2)*(a*d**2 - e*(b*d - c*e))**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_73():
    f = 1/((d + e*x)**2*(a + b/x + c/x**2))
    F = -d**2/(e*(d + e*x)*(a*d**2 - b*d*e + c*e**2)) + d*(b*d - 2*c*e)*log(d + e*x)/(a*d**2 - e*(b*d - c*e))**2 - d*(b*d - 2*c*e)*log(a*x**2 + b*x + c)/(2*(a*d**2 - e*(b*d - c*e))**2) - (b**2*d**2 - 2*b*c*d*e - 2*c*(a*d**2 - c*e**2))*atanh((2*a*x + b)/sqrt(-4*a*c + b**2))/(sqrt(-4*a*c + b**2)*(a*d**2 - e*(b*d - c*e))**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_74():
    f = 1/(x*(d + e*x)**2*(a + b/x + c/x**2))
    F = d/((d + e*x)*(a*d**2 - b*d*e + c*e**2)) - (a*d**2 - c*e**2)*log(d + e*x)/(a*d**2 - e*(b*d - c*e))**2 + (a*d**2 - c*e**2)*log(a*x**2 + b*x + c)/(2*(a*d**2 - e*(b*d - c*e))**2) + (a*d*(b*d - 4*c*e) + b*c*e**2)*atanh((2*a*x + b)/sqrt(-4*a*c + b**2))/(sqrt(-4*a*c + b**2)*(a*d**2 - e*(b*d - c*e))**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_75():
    f = 1/(x**2*(d + e*x)**2*(a + b/x + c/x**2))
    F = e*(2*a*d - b*e)*log(d + e*x)/(a*d**2 - e*(b*d - c*e))**2 - e*(2*a*d - b*e)*log(a*x**2 + b*x + c)/(2*(a*d**2 - e*(b*d - c*e))**2) - e/((d + e*x)*(a*d**2 - b*d*e + c*e**2)) - (2*a**2*d**2 - 2*a*e*(b*d + c*e) + b**2*e**2)*atanh((2*a*x + b)/sqrt(-4*a*c + b**2))/(sqrt(-4*a*c + b**2)*(a*d**2 - e*(b*d - c*e))**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_76():
    f = 1/(x**3*(d + e*x)**2*(a + b/x + c/x**2))
    F = e**2/(d*(d + e*x)*(a*d**2 - b*d*e + c*e**2)) - e**2*(3*a*d**2 - e*(2*b*d - c*e))*log(d + e*x)/(d**2*(a*d**2 - e*(b*d - c*e))**2) - (a**2*d**2 - a*e*(2*b*d + c*e) + b**2*e**2)*log(a*x**2 + b*x + c)/(2*c*(a*d**2 - e*(b*d - c*e))**2) + (a**2*d*(b*d + 4*c*e) - a*b*e*(2*b*d + 3*c*e) + b**3*e**2)*atanh((2*a*x + b)/sqrt(-4*a*c + b**2))/(c*sqrt(-4*a*c + b**2)*(a*d**2 - e*(b*d - c*e))**2) + log(x)/(c*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_77():
    f = 1/(x**4*(d + e*x)**2*(a + b/x + c/x**2))
    F = -e**3/(d**2*(d + e*x)*(a*d**2 - e*(b*d - c*e))) + e**3*(4*a*d**2 - e*(3*b*d - 2*c*e))*log(d + e*x)/(d**3*(a*d**2 - e*(b*d - c*e))**2) - 1/(c*d**2*x) + (a*d - b*e)*(a*b*d + 2*a*c*e - b**2*e)*log(a*x**2 + b*x + c)/(2*c**2*(a*d**2 - e*(b*d - c*e))**2) + (2*a**3*c*d**2 - a**2*(b**2*d**2 + 6*b*c*d*e + 2*c**2*e**2) + 2*a*b**2*e*(b*d + 2*c*e) - b**4*e**2)*atanh((2*a*x + b)/sqrt(-4*a*c + b**2))/(c**2*sqrt(-4*a*c + b**2)*(a*d**2 - e*(b*d - c*e))**2) - (b*d + 2*c*e)*log(x)/(c**2*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_78():
    f = 1/(x**5*(d + e*x)**2*(a + b/x + c/x**2))
    F = e**4/(d**3*(d + e*x)*(a*d**2 - e*(b*d - c*e))) - e**4*(5*a*d**2 - e*(4*b*d - 3*c*e))*log(d + e*x)/(d**4*(a*d**2 - e*(b*d - c*e))**2) - 1/(2*c*d**2*x**2) + (b*d + 2*c*e)/(c**2*d**3*x) + (a**3*c*d**2 - a**2*(b**2*d**2 + 4*b*c*d*e + c**2*e**2) + a*b**2*e*(2*b*d + 3*c*e) - b**4*e**2)*log(a*x**2 + b*x + c)/(2*c**3*(a*d**2 - e*(b*d - c*e))**2) + (-a**3*c*d*(3*b*d + 4*c*e) + a**2*b*(b**2*d**2 + 8*b*c*d*e + 5*c**2*e**2) - a*b**3*e*(2*b*d + 5*c*e) + b**5*e**2)*atanh((2*a*x + b)/sqrt(-4*a*c + b**2))/(c**3*sqrt(-4*a*c + b**2)*(a*d**2 - e*(b*d - c*e))**2) + (b**2*d**2 + 2*b*c*d*e - c*(a*d**2 - 3*c*e**2))*log(x)/(c**3*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_79():
    f = x**4*sqrt(d + e*x)*sqrt(a + b/x + c/x**2)
    F = 2*x**5*sqrt(d + e*x)*sqrt(a + b/x + c/x**2)/11 + x*(d + e*x)**(sympy.S(7)/2)*(2*a*d + 2*b*e)*sqrt(a + b/x + c/x**2)/(99*a*e**4) - x*(d + e*x)**(sympy.S(5)/2)*sqrt(a + b/x + c/x**2)*(58*a**2*d**2 + 2*a*e*(19*b*d - 18*c*e) + 16*b**2*e**2)/(693*a**2*e**4) + x*(d + e*x)**(sympy.S(3)/2)*sqrt(a + b/x + c/x**2)*(466*a**3*d**3 + 8*a**2*d*e*(18*b*d - 37*c*e) + 2*a*b*e**2*(67*b*d - 157*c*e) + 96*b**3*e**3)/(3465*a**3*e**4) - x*sqrt(d + e*x)*sqrt(a + b/x + c/x**2)*(374*a**4*d**4 - 8*a**3*d**2*e*(2*b*d + 3*c*e) + 6*a**2*e**2*(3*b**2*d**2 - 29*b*c*d*e + 50*c**2*e**2) + 8*a*b**2*e**3*(7*b*d - 69*c*e) + 128*b**4*e**4)/(3465*a**4*e**4) - 2*sqrt(2)*x*sqrt(a*(d + e*x)/(2*a*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(-a*(a*x**2 + b*x + c)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(a*d**2 - e*(b*d - c*e))*sqrt(a + b/x + c/x**2)*(128*a**4*d**4 + 4*a**3*d**2*e*(2*b*d + 3*c*e) - 3*a**2*e**2*(3*b**2*d**2 - 29*b*c*d*e + 50*c**2*e**2) - 4*a*b**2*e**3*(7*b*d - 69*c*e) - 64*b**4*e**4)*elliptic_f(asin(sqrt(2)*sqrt((2*a*x + b + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*a*d - e*(b + sqrt(-4*a*c + b**2))))/(3465*a**5*e**5*sqrt(d + e*x)*(a*x**2 + b*x + c)) + sqrt(2)*x*sqrt(-a*(a*x**2 + b*x + c)/(-4*a*c + b**2))*sqrt(d + e*x)*sqrt(-4*a*c + b**2)*sqrt(a + b/x + c/x**2)*(128*a**5*d**5 - 4*a**4*d**3*e*(14*b*d - 27*c*e) - a**3*d*e**2*(37*b**2*d**2 - 135*b*c*d*e + 156*c**2*e**2) - a**2*b*e**3*(37*b**2*d**2 - 258*b*c*d*e - 771*c**2*e**2) - 8*a*b**3*e**4*(7*b*d + 87*c*e) + 128*b**5*e**5)*elliptic_e(asin(sqrt(2)*sqrt((2*a*x + b + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*a*d - e*(b + sqrt(-4*a*c + b**2))))/(3465*a**5*e**5*sqrt(a*(d + e*x)/(2*a*d - e*(b + sqrt(-4*a*c + b**2))))*(a*x**2 + b*x + c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_80():
    f = x**3*sqrt(d + e*x)*sqrt(a + b/x + c/x**2)
    F = 2*x**4*sqrt(d + e*x)*sqrt(a + b/x + c/x**2)/9 + x*(d + e*x)**(sympy.S(5)/2)*(2*a*d + 2*b*e)*sqrt(a + b/x + c/x**2)/(63*a*e**3) - x*(d + e*x)**(sympy.S(3)/2)*sqrt(a + b/x + c/x**2)*(32*a**2*d**2 + 4*a*e*(4*b*d - 7*c*e) + 12*b**2*e**2)/(315*a**2*e**3) + x*sqrt(d + e*x)*sqrt(a + b/x + c/x**2)*(38*a**3*d**3 - 12*a**2*c*d*e**2 + 6*a*b*e**2*(b*d - 9*c*e) + 16*b**3*e**3)/(315*a**3*e**3) + 2*sqrt(2)*x*sqrt(a*(d + e*x)/(2*a*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(-a*(a*x**2 + b*x + c)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(a*d**2 - e*(b*d - c*e))*sqrt(a + b/x + c/x**2)*(16*a**3*d**3 + 6*a**2*c*d*e**2 - 3*a*b*e**2*(b*d - 9*c*e) - 8*b**3*e**3)*elliptic_f(asin(sqrt(2)*sqrt((2*a*x + b + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*a*d - e*(b + sqrt(-4*a*c + b**2))))/(315*a**4*e**4*sqrt(d + e*x)*(a*x**2 + b*x + c)) - 2*sqrt(2)*x*sqrt(-a*(a*x**2 + b*x + c)/(-4*a*c + b**2))*sqrt(d + e*x)*sqrt(-4*a*c + b**2)*sqrt(a + b/x + c/x**2)*(8*a**4*d**4 - a**3*d**2*e*(4*b*d - 9*c*e) - 3*a**2*e**2*(b**2*d**2 - 5*b*c*d*e - 7*c**2*e**2) - 4*a*b**2*e**3*(b*d + 9*c*e) + 8*b**4*e**4)*elliptic_e(asin(sqrt(2)*sqrt((2*a*x + b + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*a*d - e*(b + sqrt(-4*a*c + b**2))))/(315*a**4*e**4*sqrt(a*(d + e*x)/(2*a*d - e*(b + sqrt(-4*a*c + b**2))))*(a*x**2 + b*x + c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_81():
    f = x**2*sqrt(d + e*x)*sqrt(a + b/x + c/x**2)
    F = 2*x*sqrt(d + e*x)*sqrt(a + b/x + c/x**2)*(a*x**2 + b*x + c)/(7*a) - 2*x*sqrt(d + e*x)*sqrt(a + b/x + c/x**2)*(4*a**2*d**2 - 3*a*e*x*(a*d - 4*b*e) - a*e*(2*b*d - 5*c*e) + 4*b**2*e**2)/(105*a**2*e**2) - 2*sqrt(2)*x*sqrt(a*(d + e*x)/(2*a*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(-a*(a*x**2 + b*x + c)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(a*d**2 - e*(b*d - c*e))*sqrt(a + b/x + c/x**2)*(8*a**2*d**2 - a*e*(b*d - 10*c*e) - 4*b**2*e**2)*elliptic_f(asin(sqrt(2)*sqrt((2*a*x + b + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*a*d - e*(b + sqrt(-4*a*c + b**2))))/(105*a**3*e**3*sqrt(d + e*x)*(a*x**2 + b*x + c)) + sqrt(2)*x*sqrt(-a*(a*x**2 + b*x + c)/(-4*a*c + b**2))*sqrt(d + e*x)*sqrt(-4*a*c + b**2)*sqrt(a + b/x + c/x**2)*(8*a**3*d**3 - a**2*d*e*(5*b*d - 16*c*e) - a*b*e**2*(5*b*d + 29*c*e) + 8*b**3*e**3)*elliptic_e(asin(sqrt(2)*sqrt((2*a*x + b + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*a*d - e*(b + sqrt(-4*a*c + b**2))))/(105*a**3*e**3*sqrt(a*(d + e*x)/(2*a*d - e*(b + sqrt(-4*a*c + b**2))))*(a*x**2 + b*x + c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_82():
    f = x*sqrt(d + e*x)*sqrt(a + b/x + c/x**2)
    F = 2*x*(d + e*x)**(sympy.S(3)/2)*sqrt(a + b/x + c/x**2)/(5*e) - x*sqrt(d + e*x)*(4*a*d - 2*b*e)*sqrt(a + b/x + c/x**2)/(15*a*e) + 2*sqrt(2)*x*sqrt(a*(d + e*x)/(2*a*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(-a*(a*x**2 + b*x + c)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(2*a*d - b*e)*(a*d**2 - e*(b*d - c*e))*sqrt(a + b/x + c/x**2)*elliptic_f(asin(sqrt(2)*sqrt((2*a*x + b + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*a*d - e*(b + sqrt(-4*a*c + b**2))))/(15*a**2*e**2*sqrt(d + e*x)*(a*x**2 + b*x + c)) - 2*sqrt(2)*x*sqrt(-a*(a*x**2 + b*x + c)/(-4*a*c + b**2))*sqrt(d + e*x)*sqrt(-4*a*c + b**2)*sqrt(a + b/x + c/x**2)*(a**2*d**2 - a*e*(b*d + 3*c*e) + b**2*e**2)*elliptic_e(asin(sqrt(2)*sqrt((2*a*x + b + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*a*d - e*(b + sqrt(-4*a*c + b**2))))/(15*a**2*e**2*sqrt(a*(d + e*x)/(2*a*d - e*(b + sqrt(-4*a*c + b**2))))*(a*x**2 + b*x + c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_83():
    f = sqrt(d + e*x)*sqrt(a + b/x + c/x**2)
    F = ((Integer(2) * (Integer(3))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('c') * ((x)**(Integer(2)))**(Integer(-1))) + (Symbol('b') * (x)**(Integer(-1))))) * x * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) + ((sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * ((Symbol('a') * Symbol('d')) + (Symbol('b') * Symbol('e'))) * sympy.sqrt((Symbol('a') + (Symbol('c') * ((x)**(Integer(2)))**(Integer(-1))) + (Symbol('b') * (x)**(Integer(-1))))) * x * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.sqrt((Integer(-1) * ((Symbol('a') * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('a') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e')) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))))) * ((Integer(3) * Symbol('a') * Symbol('e') * sympy.sqrt(((Symbol('a') * (Symbol('d') + (Symbol('e') * x))) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))) * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('d') * ((Symbol('a') * Symbol('d')) + (Symbol('b') * Symbol('e'))) * sympy.sqrt((Symbol('a') + (Symbol('c') * ((x)**(Integer(2)))**(Integer(-1))) + (Symbol('b') * (x)**(Integer(-1))))) * x * sympy.sqrt(((Symbol('a') * (Symbol('d') + (Symbol('e') * x))) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('a') * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('a') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e')) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))))) * ((Integer(3) * Symbol('a') * Symbol('e') * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)))) + ((Integer(4) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * ((Symbol('b') * Symbol('d')) + (Symbol('c') * Symbol('e'))) * sympy.sqrt((Symbol('a') + (Symbol('c') * ((x)**(Integer(2)))**(Integer(-1))) + (Symbol('b') * (x)**(Integer(-1))))) * x * sympy.sqrt(((Symbol('a') * (Symbol('d') + (Symbol('e') * x))) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('a') * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('a') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e')) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))))) * ((Integer(3) * Symbol('a') * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * (((sympy.sqrt(Symbol('a')) * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)) * (sympy.sqrt(Integer(2)) * Symbol('c') * sympy.sqrt(((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e'))))) * sympy.sqrt((Symbol('a') + (Symbol('c') * ((x)**(Integer(2)))**(Integer(-1))) + (Symbol('b') * (x)**(Integer(-1))))) * x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('d') + (Symbol('e') * x))) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e')))))**(Integer(-1)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('d') + (Symbol('e') * x))) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))))) * sympy.Function('EllipticPi')((((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e'))) + (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e'))) * ((Integer(2) * Symbol('a') * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * (sympy.sqrt(((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e'))))))**(Integer(-1)))), ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_84():
    f = sqrt(d + e*x)*sqrt(a + b/x + c/x**2)/x
    F = ((Integer(-1) * sympy.sqrt((Symbol('a') + (Symbol('c') * ((x)**(Integer(2)))**(Integer(-1))) + (Symbol('b') * (x)**(Integer(-1)))))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) + ((Integer(3) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * sympy.sqrt((Symbol('a') + (Symbol('c') * ((x)**(Integer(2)))**(Integer(-1))) + (Symbol('b') * (x)**(Integer(-1))))) * x * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.sqrt((Integer(-1) * ((Symbol('a') * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('a') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e')) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('a') * (Symbol('d') + (Symbol('e') * x))) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))) * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('c') * ((x)**(Integer(2)))**(Integer(-1))) + (Symbol('b') * (x)**(Integer(-1))))) * x * sympy.sqrt(((Symbol('a') * (Symbol('d') + (Symbol('e') * x))) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('a') * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('a') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e')) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))))) * ((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * ((Symbol('a') * Symbol('d')) + (Symbol('b') * Symbol('e'))) * sympy.sqrt((Symbol('a') + (Symbol('c') * ((x)**(Integer(2)))**(Integer(-1))) + (Symbol('b') * (x)**(Integer(-1))))) * x * sympy.sqrt(((Symbol('a') * (Symbol('d') + (Symbol('e') * x))) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('a') * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('a') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e')) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('d')) + (Symbol('c') * Symbol('e'))) * sympy.sqrt(((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e'))))) * sympy.sqrt((Symbol('a') + (Symbol('c') * ((x)**(Integer(2)))**(Integer(-1))) + (Symbol('b') * (x)**(Integer(-1))))) * x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('d') + (Symbol('e') * x))) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e')))))**(Integer(-1)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('d') + (Symbol('e') * x))) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))))) * sympy.Function('EllipticPi')((((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e'))) + (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e'))) * ((Integer(2) * Symbol('a') * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * (sympy.sqrt(((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e'))))))**(Integer(-1)))), ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('a')) * Symbol('d') * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_85():
    f = sqrt(d + e*x)*sqrt(a + b/x + c/x**2)/x**2
    F = (Integer(-1) * ((((Symbol('b') * Symbol('d')) + (Symbol('c') * Symbol('e'))) * sympy.sqrt((Symbol('a') + (Symbol('c') * ((x)**(Integer(2)))**(Integer(-1))) + (Symbol('b') * (x)**(Integer(-1))))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * ((Integer(4) * Symbol('c') * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + (Symbol('c') * ((x)**(Integer(2)))**(Integer(-1))) + (Symbol('b') * (x)**(Integer(-1))))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * ((Integer(2) * x))**(Integer(-1)))) + ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * ((Symbol('b') * Symbol('d')) + (Symbol('c') * Symbol('e'))) * sympy.sqrt((Symbol('a') + (Symbol('c') * ((x)**(Integer(2)))**(Integer(-1))) + (Symbol('b') * (x)**(Integer(-1))))) * x * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.sqrt((Integer(-1) * ((Symbol('a') * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('a') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e')) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Integer(2)) * Symbol('c') * Symbol('d') * sympy.sqrt(((Symbol('a') * (Symbol('d') + (Symbol('e') * x))) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))) * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(3) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e') * sympy.sqrt((Symbol('a') + (Symbol('c') * ((x)**(Integer(2)))**(Integer(-1))) + (Symbol('b') * (x)**(Integer(-1))))) * x * sympy.sqrt(((Symbol('a') * (Symbol('d') + (Symbol('e') * x))) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('a') * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('a') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e')) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * ((Symbol('b') * Symbol('d')) + (Symbol('c') * Symbol('e'))) * sympy.sqrt((Symbol('a') + (Symbol('c') * ((x)**(Integer(2)))**(Integer(-1))) + (Symbol('b') * (x)**(Integer(-1))))) * x * sympy.sqrt(((Symbol('a') * (Symbol('d') + (Symbol('e') * x))) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('a') * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('a') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e')) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * Symbol('c') * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('a') * Symbol('d')) + (Symbol('b') * Symbol('e'))) * sympy.sqrt(((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e'))))) * sympy.sqrt((Symbol('a') + (Symbol('c') * ((x)**(Integer(2)))**(Integer(-1))) + (Symbol('b') * (x)**(Integer(-1))))) * x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('d') + (Symbol('e') * x))) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e')))))**(Integer(-1)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('d') + (Symbol('e') * x))) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))))) * sympy.Function('EllipticPi')((((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e'))) + (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e'))) * ((Integer(2) * Symbol('a') * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * (sympy.sqrt(((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e'))))))**(Integer(-1)))), ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('a')) * Symbol('d') * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)))) + (((((Symbol('b') * Symbol('d')) + (Symbol('c') * Symbol('e'))))**(Integer(2)) * sympy.sqrt(((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e'))))) * sympy.sqrt((Symbol('a') + (Symbol('c') * ((x)**(Integer(2)))**(Integer(-1))) + (Symbol('b') * (x)**(Integer(-1))))) * x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('d') + (Symbol('e') * x))) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e')))))**(Integer(-1)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('d') + (Symbol('e') * x))) * (((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1)))))) * sympy.Function('EllipticPi')((((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e'))) + (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e'))) * ((Integer(2) * Symbol('a') * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * (sympy.sqrt(((Integer(2) * Symbol('a') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e'))))))**(Integer(-1)))), ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('d')) * (Symbol('e'))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('a')) * Symbol('c') * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('b') * x) + (Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_86():
    f = (f*x)**m*(a + c*x**(2*n))**p*(d + e*x**n)**q
    F = sympy.Function('Unintegrable')((((Symbol('f') * x))**(Symbol('m')) * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('n')))))**(Symbol('q')) * ((Symbol('a') + (Symbol('c') * (x)**((Integer(2) * Symbol('n'))))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_87():
    f = (f*x)**m*(a + c*x**(2*n))**p*(d + e*x**n)**3
    F = d**3*(f*x)**(m + 1)*(a + c*x**(2*n))**p*hyper((-p, (m + 1)/(2*n)), (1 + (m + 1)/(2*n),), -c*x**(2*n)/a)/(f*(1 + c*x**(2*n)/a)**p*(m + 1)) + 3*d**2*e*x**(n + 1)*(f*x)**m*(a + c*x**(2*n))**p*hyper((-p, (m + n + 1)/(2*n)), ((m + 3*n + 1)/(2*n),), -c*x**(2*n)/a)/((1 + c*x**(2*n)/a)**p*(m + n + 1)) + 3*d*e**2*x**(2*n + 1)*(f*x)**m*(a + c*x**(2*n))**p*hyper((-p, (m + 2*n + 1)/(2*n)), ((m + 4*n + 1)/(2*n),), -c*x**(2*n)/a)/((1 + c*x**(2*n)/a)**p*(m + 2*n + 1)) + e**3*x**(3*n + 1)*(f*x)**m*(a + c*x**(2*n))**p*hyper((-p, (m + 3*n + 1)/(2*n)), ((m + 5*n + 1)/(2*n),), -c*x**(2*n)/a)/((1 + c*x**(2*n)/a)**p*(m + 3*n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_88():
    f = (f*x)**m*(a + c*x**(2*n))**p*(d + e*x**n)**2
    F = d**2*(f*x)**(m + 1)*(a + c*x**(2*n))**p*hyper((-p, (m + 1)/(2*n)), (1 + (m + 1)/(2*n),), -c*x**(2*n)/a)/(f*(1 + c*x**(2*n)/a)**p*(m + 1)) + 2*d*e*x**(n + 1)*(f*x)**m*(a + c*x**(2*n))**p*hyper((-p, (m + n + 1)/(2*n)), ((m + 3*n + 1)/(2*n),), -c*x**(2*n)/a)/((1 + c*x**(2*n)/a)**p*(m + n + 1)) + e**2*x**(2*n + 1)*(f*x)**m*(a + c*x**(2*n))**p*hyper((-p, (m + 2*n + 1)/(2*n)), ((m + 4*n + 1)/(2*n),), -c*x**(2*n)/a)/((1 + c*x**(2*n)/a)**p*(m + 2*n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_89():
    f = (f*x)**m*(a + c*x**(2*n))**p*(d + e*x**n)
    F = d*(f*x)**(m + 1)*(a + c*x**(2*n))**p*hyper((-p, (m + 1)/(2*n)), (1 + (m + 1)/(2*n),), -c*x**(2*n)/a)/(f*(1 + c*x**(2*n)/a)**p*(m + 1)) + e*x**(n + 1)*(f*x)**m*(a + c*x**(2*n))**p*hyper((-p, (m + n + 1)/(2*n)), ((m + 3*n + 1)/(2*n),), -c*x**(2*n)/a)/((1 + c*x**(2*n)/a)**p*(m + n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_90():
    f = (f*x)**m*(a + c*x**(2*n))**p/(d + e*x**n)
    F = x*(f*x)**m*(a + c*x**(2*n))**p*appellf1((m + 1)/(2*n), 1, -p, 1 + (m + 1)/(2*n), e**2*x**(2*n)/d**2, -c*x**(2*n)/a)/(d*(1 + c*x**(2*n)/a)**p*(m + 1)) - e*x**(n + 1)*(f*x)**m*(a + c*x**(2*n))**p*appellf1((m + n + 1)/(2*n), 1, -p, (m + 3*n + 1)/(2*n), e**2*x**(2*n)/d**2, -c*x**(2*n)/a)/(d**2*(1 + c*x**(2*n)/a)**p*(m + n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_91():
    f = (f*x)**m*(a + c*x**(2*n))**p/(d + e*x**n)**2
    F = x*(f*x)**m*(a + c*x**(2*n))**p*appellf1((m + 1)/(2*n), 2, -p, 1 + (m + 1)/(2*n), e**2*x**(2*n)/d**2, -c*x**(2*n)/a)/(d**2*(1 + c*x**(2*n)/a)**p*(m + 1)) - 2*e*x**(n + 1)*(f*x)**m*(a + c*x**(2*n))**p*appellf1((m + n + 1)/(2*n), 2, -p, (m + 3*n + 1)/(2*n), e**2*x**(2*n)/d**2, -c*x**(2*n)/a)/(d**3*(1 + c*x**(2*n)/a)**p*(m + n + 1)) + e**2*x**(2*n + 1)*(f*x)**m*(a + c*x**(2*n))**p*appellf1((m + 2*n + 1)/(2*n), 2, -p, (m + 4*n + 1)/(2*n), e**2*x**(2*n)/d**2, -c*x**(2*n)/a)/(d**4*(1 + c*x**(2*n)/a)**p*(m + 2*n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_92():
    f = (f*x)**m*(a + c*x**(2*n))**p/(d + e*x**n)**3
    F = x*(f*x)**m*(a + c*x**(2*n))**p*appellf1((m + 1)/(2*n), 3, -p, 1 + (m + 1)/(2*n), e**2*x**(2*n)/d**2, -c*x**(2*n)/a)/(d**3*(1 + c*x**(2*n)/a)**p*(m + 1)) - 3*e*x**(n + 1)*(f*x)**m*(a + c*x**(2*n))**p*appellf1((m + n + 1)/(2*n), 3, -p, (m + 3*n + 1)/(2*n), e**2*x**(2*n)/d**2, -c*x**(2*n)/a)/(d**4*(1 + c*x**(2*n)/a)**p*(m + n + 1)) + 3*e**2*x**(2*n + 1)*(f*x)**m*(a + c*x**(2*n))**p*appellf1((m + 2*n + 1)/(2*n), 3, -p, (m + 4*n + 1)/(2*n), e**2*x**(2*n)/d**2, -c*x**(2*n)/a)/(d**5*(1 + c*x**(2*n)/a)**p*(m + 2*n + 1)) - e**3*x**(3*n + 1)*(f*x)**m*(a + c*x**(2*n))**p*appellf1((m + 3*n + 1)/(2*n), 3, -p, (m + 5*n + 1)/(2*n), e**2*x**(2*n)/d**2, -c*x**(2*n)/a)/(d**6*(1 + c*x**(2*n)/a)**p*(m + 3*n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_93():
    f = (b + 2*c*x)*(a + b*x + c*x**2)**13
    F = (a + b*x + c*x**2)**14/14
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_94():
    f = x*(b + 2*c*x**2)*(a + b*x**2 + c*x**4)**13
    F = (a + b*x**2 + c*x**4)**14/28
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_95():
    f = x**2*(b + 2*c*x**3)*(a + b*x**3 + c*x**6)**13
    F = (a + b*x**3 + c*x**6)**14/42
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_96():
    f = x**(n - 1)*(b + 2*c*x**n)*(a + b*x**n + c*x**(2*n))**13
    F = (a + b*x**n + c*x**(2*n))**14/(14*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_97():
    f = (b + 2*c*x)*(-a + b*x + c*x**2)**13
    F = (a - b*x - c*x**2)**14/14
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_98():
    f = x*(b + 2*c*x**2)*(-a + b*x**2 + c*x**4)**13
    F = (a - b*x**2 - c*x**4)**14/28
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_99():
    f = x**2*(b + 2*c*x**3)*(-a + b*x**3 + c*x**6)**13
    F = (a - b*x**3 - c*x**6)**14/42
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_100():
    f = x**(n - 1)*(b + 2*c*x**n)*(-a + b*x**n + c*x**(2*n))**13
    F = (a - b*x**n - c*x**(2*n))**14/(14*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_101():
    f = (b + 2*c*x)*(b*x + c*x**2)**13
    F = (b*x + c*x**2)**14/14
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_102():
    f = x*(b + 2*c*x**2)*(b*x**2 + c*x**4)**13
    F = x**28*(b + c*x**2)**14/28
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_103():
    f = x**2*(b + 2*c*x**3)*(b*x**3 + c*x**6)**13
    F = x**42*(b + c*x**3)**14/42
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_104():
    f = x**(n - 1)*(b + 2*c*x**n)*(b*x**n + c*x**(2*n))**13
    F = x**(14*n)*(b + c*x**n)**14/(14*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_105():
    f = (b + 2*c*x)/(a + b*x + c*x**2)
    F = log(a + b*x + c*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_106():
    f = x*(b + 2*c*x**2)/(a + b*x**2 + c*x**4)
    F = log(a + b*x**2 + c*x**4)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_107():
    f = x**2*(b + 2*c*x**3)/(a + b*x**3 + c*x**6)
    F = log(a + b*x**3 + c*x**6)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_108():
    f = x**(n - 1)*(b + 2*c*x**n)/(a + b*x**n + c*x**(2*n))
    F = log(a + b*x**n + c*x**(2*n))/n
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_109():
    f = (b + 2*c*x)/(a + b*x + c*x**2)**8
    F = -1/(7*(a + b*x + c*x**2)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_110():
    f = x*(b + 2*c*x**2)/(a + b*x**2 + c*x**4)**8
    F = -1/(14*(a + b*x**2 + c*x**4)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_111():
    f = x**2*(b + 2*c*x**3)/(a + b*x**3 + c*x**6)**8
    F = -1/(21*(a + b*x**3 + c*x**6)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_112():
    f = x**(n - 1)*(b + 2*c*x**n)/(a + b*x**n + c*x**(2*n))**8
    F = -1/(7*n*(a + b*x**n + c*x**(2*n))**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_113():
    f = (b + 2*c*x)/(-a + b*x + c*x**2)
    F = log(a - b*x - c*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_114():
    f = x*(b + 2*c*x**2)/(-a + b*x**2 + c*x**4)
    F = log(a - b*x**2 - c*x**4)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_115():
    f = x**2*(b + 2*c*x**3)/(-a + b*x**3 + c*x**6)
    F = log(a - b*x**3 - c*x**6)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_116():
    f = x**(n - 1)*(b + 2*c*x**n)/(-a + b*x**n + c*x**(2*n))
    F = log(a - b*x**n - c*x**(2*n))/n
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_117():
    f = (b + 2*c*x)/(-a + b*x + c*x**2)**8
    F = 1/(7*(a - b*x - c*x**2)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_118():
    f = x*(b + 2*c*x**2)/(-a + b*x**2 + c*x**4)**8
    F = 1/(14*(a - b*x**2 - c*x**4)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_119():
    f = x**2*(b + 2*c*x**3)/(-a + b*x**3 + c*x**6)**8
    F = 1/(21*(a - b*x**3 - c*x**6)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_120():
    f = x**(n - 1)*(b + 2*c*x**n)/(-a + b*x**n + c*x**(2*n))**8
    F = 1/(7*n*(a - b*x**n - c*x**(2*n))**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_121():
    f = (b + 2*c*x)/(b*x + c*x**2)
    F = log(b*x + c*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_122():
    f = x**(n - 1)*(b + 2*c*x**n)/(b*x**n + c*x**(2*n))
    F = log(x) + log(b + c*x**n)/n
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_123():
    f = (b + 2*c*x)/(b*x + c*x**2)**8
    F = -1/(7*(b*x + c*x**2)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_124():
    f = x*(b + 2*c*x**2)/(b*x**2 + c*x**4)**8
    F = -1/(14*x**14*(b + c*x**2)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_125():
    f = x**2*(b + 2*c*x**3)/(b*x**3 + c*x**6)**8
    F = -1/(21*x**21*(b + c*x**3)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_126():
    f = x**(n - 1)*(b + 2*c*x**n)/(b*x**n + c*x**(2*n))**8
    F = -1/(7*n*x**(7*n)*(b + c*x**n)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_127():
    f = (b + 2*c*x)*(a + b*x + c*x**2)**p
    F = (a + b*x + c*x**2)**(p + 1)/(p + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_128():
    f = x*(b + 2*c*x**2)*(a + b*x**2 + c*x**4)**p
    F = (a + b*x**2 + c*x**4)**(p + 1)/(2*p + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_129():
    f = x**2*(b + 2*c*x**3)*(a + b*x**3 + c*x**6)**p
    F = (a + b*x**3 + c*x**6)**(p + 1)/(3*p + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_130():
    f = x**(n - 1)*(b + 2*c*x**n)*(a + b*x**n + c*x**(2*n))**p
    F = (a + b*x**n + c*x**(2*n))**(p + 1)/(n*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_131():
    f = (b + 2*c*x)*(-a + b*x + c*x**2)**p
    F = (-a + b*x + c*x**2)**(p + 1)/(p + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_132():
    f = x*(b + 2*c*x**2)*(-a + b*x**2 + c*x**4)**p
    F = (-a + b*x**2 + c*x**4)**(p + 1)/(2*p + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_133():
    f = x**2*(b + 2*c*x**3)*(-a + b*x**3 + c*x**6)**p
    F = (-a + b*x**3 + c*x**6)**(p + 1)/(3*p + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_134():
    f = x**(n - 1)*(b + 2*c*x**n)*(-a + b*x**n + c*x**(2*n))**p
    F = (-a + b*x**n + c*x**(2*n))**(p + 1)/(n*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_135():
    f = (b + 2*c*x)*(b*x + c*x**2)**p
    F = (b*x + c*x**2)**(p + 1)/(p + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_136():
    f = x*(b + 2*c*x**2)*(b*x**2 + c*x**4)**p
    F = (b*x**2 + c*x**4)**(p + 1)/(2*p + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_137():
    f = x**2*(b + 2*c*x**3)*(b*x**3 + c*x**6)**p
    F = (b*x**3 + c*x**6)**(p + 1)/(3*p + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_138():
    f = x**(n - 1)*(b + 2*c*x**n)*(b*x**n + c*x**(2*n))**p
    F = (b*x**n + c*x**(2*n))**(p + 1)/(n*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_139():
    f = (f*x)**m*(d + e*x**n)/(a + b*x**n + c*x**(2*n))
    F = (f*x)**(m + 1)*(e - (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(f*(b + sqrt(-4*a*c + b**2))*(m + 1)) + (f*x)**(m + 1)*(e + (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(f*(b - sqrt(-4*a*c + b**2))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_140():
    f = (f*x)**m*(d + e*x**n)/(a + b*x**n + c*x**(2*n))**2
    F = -c*(f*x)**(m + 1)*((-2*a*e + b*d)*(m - n + 1) + (2*a*b*e*n + 4*a*c*d*(m - 2*n + 1) - b**2*d*(m - n + 1))/sqrt(-4*a*c + b**2))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*f*n*(b + sqrt(-4*a*c + b**2))*(m + 1)*(-4*a*c + b**2)) - c*(f*x)**(m + 1)*((-2*a*e + b*d)*(m - n + 1) - (2*a*b*e*n + 4*a*c*d*(m - 2*n + 1) - b**2*d*(m - n + 1))/sqrt(-4*a*c + b**2))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*f*n*(b - sqrt(-4*a*c + b**2))*(m + 1)*(-4*a*c + b**2)) + (f*x)**(m + 1)*(-a*b*e - 2*a*c*d + b**2*d + c*x**n*(-2*a*e + b*d))/(a*f*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_141():
    f = (f*x)**m*(d + e*x**n)/(a + b*x**n + c*x**(2*n))**3
    F = (f*x)**(m + 1)*(-a*b*e - 2*a*c*d + b**2*d + c*x**n*(-2*a*e + b*d))/(2*a*f*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))**2) - c*(f*x)**(m + 1)*((m - n + 1)*(-4*a**2*c*e*(m - 3*n + 1) + a*b**2*e*(m + 1) + 2*a*b*c*d*(2*m - 7*n + 2) - b**3*d*(m - 2*n + 1)) - (-4*a**2*b*c*e*(m**2 + m*(2 - n) - 3*n**2 - n + 1) - 8*a**2*c**2*d*(m**2 + m*(2 - 6*n) + 8*n**2 - 6*n + 1) + a*b**3*e*(m + 1)*(m - n + 1) + 6*a*b**2*c*d*(m**2 + m*(2 - 4*n) + 3*n**2 - 4*n + 1) - b**4*d*(m**2 + m*(2 - 3*n) + 2*n**2 - 3*n + 1))/sqrt(-4*a*c + b**2))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(2*a**2*f*n**2*(b + sqrt(-4*a*c + b**2))*(m + 1)*(-4*a*c + b**2)**2) - c*(f*x)**(m + 1)*((m - n + 1)*(-4*a**2*c*e*(m - 3*n + 1) + a*b**2*e*(m + 1) + 2*a*b*c*d*(2*m - 7*n + 2) - b**3*d*(m - 2*n + 1)) + (-4*a**2*b*c*e*(m**2 + m*(2 - n) - 3*n**2 - n + 1) - 8*a**2*c**2*d*(m**2 + m*(2 - 6*n) + 8*n**2 - 6*n + 1) + a*b**3*e*(m + 1)*(m - n + 1) + 6*a*b**2*c*d*(m**2 + m*(2 - 4*n) + 3*n**2 - 4*n + 1) - b**4*d*(m**2 + m*(2 - 3*n) + 2*n**2 - 3*n + 1))/sqrt(-4*a*c + b**2))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(2*a**2*f*n**2*(b - sqrt(-4*a*c + b**2))*(m + 1)*(-4*a*c + b**2)**2) + (f*x)**(m + 1)*(a*b*c*(-2*a*e + b*d)*(m - 3*n + 1) + c*x**n*(-4*a**2*c*e*(m - 3*n + 1) + a*b**2*e*(m + 1) + 2*a*b*c*d*(2*m - 7*n + 2) - b**3*d*(m - 2*n + 1)) + (-2*a*c + b**2)*(a*b*e*(m + 1) + 2*a*c*d*(m - 4*n + 1) - b**2*d*(m - 2*n + 1)))/(2*a**2*f*n**2*(-4*a*c + b**2)**2*(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_142():
    f = (c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x**(sympy.S(1)/3))/(-c**(sympy.S(2)/3)*d**(sympy.S(2)/3)*x + c**(sympy.S(1)/3)*d*x**(sympy.S(4)/3) + c*d**(sympy.S(1)/3)*x**(sympy.S(2)/3))
    F = -3*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x**(sympy.S(1)/3) + d**(sympy.S(2)/3)*x**(sympy.S(2)/3))/(c**(sympy.S(1)/3)*d**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_143():
    f = (f*x)**m*(d + e*x**n)**q/(a + b*x**n + c*x**(2*n))
    F = -2*c*(f*x)**(m + 1)*(d + e*x**n)**q*appellf1((m + 1)/n, 1, -q, (m + n + 1)/n, -2*c*x**n/(b + sqrt(-4*a*c + b**2)), -e*x**n/d)/(f*(1 + e*x**n/d)**q*(b + sqrt(-4*a*c + b**2))*(m + 1)*sqrt(-4*a*c + b**2)) + 2*c*(f*x)**(m + 1)*(d + e*x**n)**q*appellf1((m + 1)/n, 1, -q, (m + n + 1)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -e*x**n/d)/(f*(1 + e*x**n/d)**q*(b - sqrt(-4*a*c + b**2))*(m + 1)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_144():
    f = x**2*(d + e*x**n)**q/(a + b*x**n + c*x**(2*n))
    F = -2*c*x**3*(d + e*x**n)**q*appellf1(3/n, 1, -q, (n + 3)/n, -2*c*x**n/(b + sqrt(-4*a*c + b**2)), -e*x**n/d)/((1 + e*x**n/d)**q*(-12*a*c + 3*b**2 + 3*b*sqrt(-4*a*c + b**2))) - 2*c*x**3*(d + e*x**n)**q*appellf1(3/n, 1, -q, (n + 3)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -e*x**n/d)/((1 + e*x**n/d)**q*(-12*a*c + 3*b**2 - 3*b*sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_145():
    f = x*(d + e*x**n)**q/(a + b*x**n + c*x**(2*n))
    F = -c*x**2*(d + e*x**n)**q*appellf1(2/n, 1, -q, (n + 2)/n, -2*c*x**n/(b + sqrt(-4*a*c + b**2)), -e*x**n/d)/((1 + e*x**n/d)**q*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))) - c*x**2*(d + e*x**n)**q*appellf1(2/n, 1, -q, (n + 2)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -e*x**n/d)/((1 + e*x**n/d)**q*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_146():
    f = (d + e*x**n)**q/(a + b*x**n + c*x**(2*n))
    F = -2*c*x*(d + e*x**n)**q*appellf1(1/n, 1, -q, 1 + 1/n, -2*c*x**n/(b + sqrt(-4*a*c + b**2)), -e*x**n/d)/((1 + e*x**n/d)**q*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))) - 2*c*x*(d + e*x**n)**q*appellf1(1/n, 1, -q, 1 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -e*x**n/d)/((1 + e*x**n/d)**q*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_147():
    f = (d + e*x**n)**q/(x*(a + b*x**n + c*x**(2*n)))
    F = c*(d + e*x**n)**(q + 1)*(-b/sqrt(-4*a*c + b**2) + 1)*hyper((1, q + 1), (q + 2,), 2*c*(d + e*x**n)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(a*n*(q + 1)*(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + c*(d + e*x**n)**(q + 1)*(b/sqrt(-4*a*c + b**2) + 1)*hyper((1, q + 1), (q + 2,), 2*c*(d + e*x**n)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(a*n*(q + 1)*(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) - (d + e*x**n)**(q + 1)*hyper((1, q + 1), (q + 2,), 1 + e*x**n/d)/(a*d*n*(q + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_148():
    f = (d + e*x**n)**q/(x**2*(a + b*x**n + c*x**(2*n)))
    F = 2*c*(d + e*x**n)**q*appellf1(-1/n, 1, -q, -(1 - n)/n, -2*c*x**n/(b + sqrt(-4*a*c + b**2)), -e*x**n/d)/(x*(1 + e*x**n/d)**q*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))) + 2*c*(d + e*x**n)**q*appellf1(-1/n, 1, -q, -(1 - n)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -e*x**n/d)/(x*(1 + e*x**n/d)**q*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_149():
    f = (d + e*x**n)**q/(x**3*(a + b*x**n + c*x**(2*n)))
    F = c*(d + e*x**n)**q*appellf1(-2/n, 1, -q, -(2 - n)/n, -2*c*x**n/(b + sqrt(-4*a*c + b**2)), -e*x**n/d)/(x**2*(1 + e*x**n/d)**q*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))) + c*(d + e*x**n)**q*appellf1(-2/n, 1, -q, -(2 - n)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -e*x**n/d)/(x**2*(1 + e*x**n/d)**q*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_150():
    f = (f*x)**m*(d + e*x**n)**2*(a + b*x**n + c*x**(2*n))**p
    F = d**2*(f*x)**(m + 1)*(a + b*x**n + c*x**(2*n))**p*appellf1((m + 1)/n, -p, -p, (m + n + 1)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(f*(m + 1)*(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)**p) + 2*d*e*x**(n + 1)*(f*x)**m*(a + b*x**n + c*x**(2*n))**p*appellf1((m + n + 1)/n, -p, -p, (m + 2*n + 1)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)**p*(m + n + 1)) + e**2*x**(2*n + 1)*(f*x)**m*(a + b*x**n + c*x**(2*n))**p*appellf1((m + 2*n + 1)/n, -p, -p, (m + 3*n + 1)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)**p*(m + 2*n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_151():
    f = (f*x)**m*(d + e*x**n)*(a + b*x**n + c*x**(2*n))**p
    F = d*(f*x)**(m + 1)*(a + b*x**n + c*x**(2*n))**p*appellf1((m + 1)/n, -p, -p, (m + n + 1)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(f*(m + 1)*(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)**p) + e*x**(n + 1)*(f*x)**m*(a + b*x**n + c*x**(2*n))**p*appellf1((m + n + 1)/n, -p, -p, (m + 2*n + 1)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)**p*(m + n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_152():
    f = (f*x)**m*(a + b*x**n + c*x**(2*n))**p
    F = (f*x)**(m + 1)*(a + b*x**n + c*x**(2*n))**p*appellf1((m + 1)/n, -p, -p, (m + n + 1)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(f*(m + 1)*(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_153():
    f = (f*x)**m*(a + b*x**n + c*x**(2*n))**p/(d + e*x**n)
    F = sympy.Function('Unintegrable')(((((Symbol('f') * x))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * (x)**(Symbol('n'))) + (Symbol('c') * (x)**((Integer(2) * Symbol('n'))))))**(Symbol('p'))) * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('n')))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_4_f_x_pow_m_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_154():
    f = (f*x)**m*(a + b*x**n + c*x**(2*n))**p/(d + e*x**n)**2
    F = sympy.Function('Unintegrable')(((((Symbol('f') * x))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * (x)**(Symbol('n'))) + (Symbol('c') * (x)**((Integer(2) * Symbol('n'))))))**(Symbol('p'))) * (((Symbol('d') + (Symbol('e') * (x)**(Symbol('n')))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F

