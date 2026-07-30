"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.2 Trinomial products/1.2.3 General/1.2.3.5 P(x) (d x)^m (a+b x^n+c x^(2 n))^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, C, D, a, b, c, d, e, f, g, h, j, k, l, m, n, p = symbols('A B C D a b c d e f g h j k l m n p')

def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_5_P_x_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_1():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4 + j*x**5 + k*x**6 + l*x**7 + m*x**8)/(a + b*x**3 + c*x**6)
    F = k*x/c + l*x**2/(2*c) + m*x**3/(3*c) + (-b*m + c*j)*log(a + b*x**3 + c*x**6)/(6*c**2) - (-2*a*c*m + b**2*m - b*c*j + 2*c**2*f)*atanh((b + 2*c*x**3)/sqrt(-4*a*c + b**2))/(3*c**2*sqrt(-4*a*c + b**2)) + 2**(sympy.S(2)/3)*(-b*k/c + g - (-2*a*c*k + b**2*k - b*c*g + 2*c**2*d)/(c*sqrt(-4*a*c + b**2)))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(1)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(-b*k/c + g - (-2*a*c*k + b**2*k - b*c*g + 2*c**2*d)/(c*sqrt(-4*a*c + b**2)))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(1)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*sqrt(3)*(-b*k/c + g - (-2*a*c*k + b**2*k - b*c*g + 2*c**2*d)/(c*sqrt(-4*a*c + b**2)))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(1)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(-b*k/c + g + (b**2*k + 2*c**2*d - c*(2*a*k + b*g))/(c*sqrt(-4*a*c + b**2)))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(1)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(-b*k/c + g + (b**2*k + 2*c**2*d - c*(2*a*k + b*g))/(c*sqrt(-4*a*c + b**2)))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(1)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*sqrt(3)*(-b*k/c + g + (b**2*k + 2*c**2*d - c*(2*a*k + b*g))/(c*sqrt(-4*a*c + b**2)))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(1)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(1)/3)*(-b*l/c + h - (-2*a*c*l + b**2*l - b*c*h + 2*c**2*e)/(c*sqrt(-4*a*c + b**2)))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(2)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*(-b*l/c + h - (-2*a*c*l + b**2*l - b*c*h + 2*c**2*e)/(c*sqrt(-4*a*c + b**2)))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(2)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*sqrt(3)*(-b*l/c + h - (-2*a*c*l + b**2*l - b*c*h + 2*c**2*e)/(c*sqrt(-4*a*c + b**2)))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(2)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*(-b*l/c + h + (b**2*l + 2*c**2*e - c*(2*a*l + b*h))/(c*sqrt(-4*a*c + b**2)))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(2)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*(-b*l/c + h + (b**2*l + 2*c**2*e - c*(2*a*l + b*h))/(c*sqrt(-4*a*c + b**2)))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(2)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*sqrt(3)*(-b*l/c + h + (b**2*l + 2*c**2*e - c*(2*a*l + b*h))/(c*sqrt(-4*a*c + b**2)))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(2)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_5_P_x_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_2():
    f = 1/(a + b*x**n + c*x**(2*n))
    F = -2*c*x*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2)) - 2*c*x*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_5_P_x_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_3():
    f = (d + e*x)/(a + b*x**n + c*x**(2*n))
    F = -2*c*d*x*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2)) - 2*c*d*x*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2)) - c*e*x**2*hyper((1, 2/n), ((n + 2)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2)) - c*e*x**2*hyper((1, 2/n), ((n + 2)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_5_P_x_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_4():
    f = (d + e*x + f*x**2)/(a + b*x**n + c*x**(2*n))
    F = -2*c*d*x*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2)) - 2*c*d*x*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2)) - c*e*x**2*hyper((1, 2/n), ((n + 2)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2)) - c*e*x**2*hyper((1, 2/n), ((n + 2)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2)) - 2*c*f*x**3*hyper((1, 3/n), ((n + 3)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(-12*a*c + 3*b**2 + 3*b*sqrt(-4*a*c + b**2)) - 2*c*f*x**3*hyper((1, 3/n), ((n + 3)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(-12*a*c + 3*b**2 - 3*b*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_5_P_x_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_5():
    f = (d + e*x + f*x**2 + g*x**3)/(a + b*x**n + c*x**(2*n))
    F = -2*c*d*x*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2)) - 2*c*d*x*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2)) - c*e*x**2*hyper((1, 2/n), ((n + 2)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2)) - c*e*x**2*hyper((1, 2/n), ((n + 2)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2)) - 2*c*f*x**3*hyper((1, 3/n), ((n + 3)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(-12*a*c + 3*b**2 + 3*b*sqrt(-4*a*c + b**2)) - 2*c*f*x**3*hyper((1, 3/n), ((n + 3)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(-12*a*c + 3*b**2 - 3*b*sqrt(-4*a*c + b**2)) - c*g*x**4*hyper((1, 4/n), ((n + 4)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(-8*a*c + 2*b**2 + 2*b*sqrt(-4*a*c + b**2)) - c*g*x**4*hyper((1, 4/n), ((n + 4)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(-8*a*c + 2*b**2 - 2*b*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_5_P_x_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_6():
    f = (a + b*x**n + c*x**(2*n))**(-2)
    F = -c*x*(4*a*c*(1 - 2*n) - b**2*(1 - n) + b*(1 - n)*sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))) - c*x*(4*a*c*(1 - 2*n) - b**2*(1 - n) - b*(1 - n)*sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))) + x*(-2*a*c + b**2 + b*c*x**n)/(a*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_5_P_x_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_7():
    f = (d + e*x)/(a + b*x**n + c*x**(2*n))**2
    F = 2*b*c**2*e*x**(n + 2)*(2 - n)*hyper((1, (n + 2)/n), (2 + 2/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*n*(b + sqrt(-4*a*c + b**2))*(n + 2)*(-4*a*c + b**2)**(sympy.S(3)/2)) - 2*b*c**2*e*x**(n + 2)*(2 - n)*hyper((1, (n + 2)/n), (2 + 2/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*n*(b - sqrt(-4*a*c + b**2))*(n + 2)*(-4*a*c + b**2)**(sympy.S(3)/2)) - c*d*x*(4*a*c*(1 - 2*n) - b**2*(1 - n) + b*(1 - n)*sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))) - c*d*x*(4*a*c*(1 - 2*n) - b**2*(1 - n) - b*(1 - n)*sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))) - c*e*x**2*(4*a*c*(1 - n) - b**2*(2 - n))*hyper((1, 2/n), ((n + 2)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))) - c*e*x**2*(4*a*c*(1 - n) - b**2*(2 - n))*hyper((1, 2/n), ((n + 2)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))) + d*x*(-2*a*c + b**2 + b*c*x**n)/(a*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))) + e*x**2*(-2*a*c + b**2 + b*c*x**n)/(a*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_5_P_x_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_8():
    f = (d + e*x + f*x**2)/(a + b*x**n + c*x**(2*n))**2
    F = 2*b*c**2*e*x**(n + 2)*(2 - n)*hyper((1, (n + 2)/n), (2 + 2/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*n*(b + sqrt(-4*a*c + b**2))*(n + 2)*(-4*a*c + b**2)**(sympy.S(3)/2)) - 2*b*c**2*e*x**(n + 2)*(2 - n)*hyper((1, (n + 2)/n), (2 + 2/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*n*(b - sqrt(-4*a*c + b**2))*(n + 2)*(-4*a*c + b**2)**(sympy.S(3)/2)) + 2*b*c**2*f*x**(n + 3)*(3 - n)*hyper((1, (n + 3)/n), (2 + 3/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*n*(b + sqrt(-4*a*c + b**2))*(n + 3)*(-4*a*c + b**2)**(sympy.S(3)/2)) - 2*b*c**2*f*x**(n + 3)*(3 - n)*hyper((1, (n + 3)/n), (2 + 3/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*n*(b - sqrt(-4*a*c + b**2))*(n + 3)*(-4*a*c + b**2)**(sympy.S(3)/2)) - c*d*x*(4*a*c*(1 - 2*n) - b**2*(1 - n) + b*(1 - n)*sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))) - c*d*x*(4*a*c*(1 - 2*n) - b**2*(1 - n) - b*(1 - n)*sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))) - c*e*x**2*(4*a*c*(1 - n) - b**2*(2 - n))*hyper((1, 2/n), ((n + 2)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))) - c*e*x**2*(4*a*c*(1 - n) - b**2*(2 - n))*hyper((1, 2/n), ((n + 2)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))) - 2*c*f*x**3*(2*a*c*(3 - 2*n) - b**2*(3 - n))*hyper((1, 3/n), ((n + 3)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(3*a*n*(-4*a*c + b**2)*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))) - 2*c*f*x**3*(2*a*c*(3 - 2*n) - b**2*(3 - n))*hyper((1, 3/n), ((n + 3)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(3*a*n*(-4*a*c + b**2)*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))) + d*x*(-2*a*c + b**2 + b*c*x**n)/(a*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))) + e*x**2*(-2*a*c + b**2 + b*c*x**n)/(a*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))) + f*x**3*(-2*a*c + b**2 + b*c*x**n)/(a*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_5_P_x_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_9():
    f = (d + e*x + f*x**2 + g*x**3)/(a + b*x**n + c*x**(2*n))**2
    F = 2*b*c**2*e*x**(n + 2)*(2 - n)*hyper((1, (n + 2)/n), (2 + 2/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*n*(b + sqrt(-4*a*c + b**2))*(n + 2)*(-4*a*c + b**2)**(sympy.S(3)/2)) - 2*b*c**2*e*x**(n + 2)*(2 - n)*hyper((1, (n + 2)/n), (2 + 2/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*n*(b - sqrt(-4*a*c + b**2))*(n + 2)*(-4*a*c + b**2)**(sympy.S(3)/2)) + 2*b*c**2*f*x**(n + 3)*(3 - n)*hyper((1, (n + 3)/n), (2 + 3/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*n*(b + sqrt(-4*a*c + b**2))*(n + 3)*(-4*a*c + b**2)**(sympy.S(3)/2)) - 2*b*c**2*f*x**(n + 3)*(3 - n)*hyper((1, (n + 3)/n), (2 + 3/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*n*(b - sqrt(-4*a*c + b**2))*(n + 3)*(-4*a*c + b**2)**(sympy.S(3)/2)) + 2*b*c**2*g*x**(n + 4)*(4 - n)*hyper((1, (n + 4)/n), (2 + 4/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*n*(b + sqrt(-4*a*c + b**2))*(n + 4)*(-4*a*c + b**2)**(sympy.S(3)/2)) - 2*b*c**2*g*x**(n + 4)*(4 - n)*hyper((1, (n + 4)/n), (2 + 4/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*n*(b - sqrt(-4*a*c + b**2))*(n + 4)*(-4*a*c + b**2)**(sympy.S(3)/2)) - c*d*x*(4*a*c*(1 - 2*n) - b**2*(1 - n) + b*(1 - n)*sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))) - c*d*x*(4*a*c*(1 - 2*n) - b**2*(1 - n) - b*(1 - n)*sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))) - c*e*x**2*(4*a*c*(1 - n) - b**2*(2 - n))*hyper((1, 2/n), ((n + 2)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))) - c*e*x**2*(4*a*c*(1 - n) - b**2*(2 - n))*hyper((1, 2/n), ((n + 2)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))) - 2*c*f*x**3*(2*a*c*(3 - 2*n) - b**2*(3 - n))*hyper((1, 3/n), ((n + 3)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(3*a*n*(-4*a*c + b**2)*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))) - 2*c*f*x**3*(2*a*c*(3 - 2*n) - b**2*(3 - n))*hyper((1, 3/n), ((n + 3)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(3*a*n*(-4*a*c + b**2)*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))) - c*g*x**4*(4*a*c*(2 - n) - b**2*(4 - n))*hyper((1, 4/n), ((n + 4)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(2*a*n*(-4*a*c + b**2)*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))) - c*g*x**4*(4*a*c*(2 - n) - b**2*(4 - n))*hyper((1, 4/n), ((n + 4)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(2*a*n*(-4*a*c + b**2)*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))) + d*x*(-2*a*c + b**2 + b*c*x**n)/(a*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))) + e*x**2*(-2*a*c + b**2 + b*c*x**n)/(a*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))) + f*x**3*(-2*a*c + b**2 + b*c*x**n)/(a*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))) + g*x**4*(-2*a*c + b**2 + b*c*x**n)/(a*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_5_P_x_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_10():
    f = (-a*h*x**(n/2 - 1) + c*f*x**(n - 1) + c*g*x**(2*n - 1) + c*h*x**(5*n/2 - 1))/(a + b*x**n + c*x**(2*n))**(sympy.S(3)/2)
    F = -(2*c*x**n*(-b*g + 2*c*f) + 2*c*(-2*a*g + b*f) + 2*h*x**(n/2)*(-4*a*c + b**2))/(n*(-4*a*c + b**2)*sqrt(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_5_P_x_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_11():
    f = (a + b*x**n + c*x**(2*n))**p*(a + b*x**n*(n*p + n + 1) + c*x**(2*n)*(2*n*(p + 1) + 1))
    F = x*(a + b*x**n + c*x**(2*n))**(p + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_5_P_x_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_12():
    f = x**(n/4 - 1)*(-a*h + c*f*x**(n/4) + c*g*x**(3*n/4) + c*h*x**n)/(a + c*x**n)**(sympy.S(3)/2)
    F = -(2*a*g + 4*a*h*x**(n/4) - 2*c*f*x**(n/2))/(a*n*sqrt(a + c*x**n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_5_P_x_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_13():
    f = (d*x)**(n/4 - 1)*(-a*h + c*f*x**(n/4) + c*g*x**(3*n/4) + c*h*x**n)/(a + c*x**n)**(sympy.S(3)/2)
    F = -2*x**(1 - n/4)*(d*x)**(n/4 - 1)*(a*g + 2*a*h*x**(n/4) - c*f*x**(n/2))/(a*n*sqrt(a + c*x**n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_5_P_x_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_14():
    f = x**(n/2 - 1)*(-a*h + c*f*x**(n/2) + c*g*x**(3*n/2) + c*h*x**(2*n))/(a + b*x**n + c*x**(2*n))**(sympy.S(3)/2)
    F = -(2*c*x**n*(-b*g + 2*c*f) + 2*c*(-2*a*g + b*f) + 2*h*x**(n/2)*(-4*a*c + b**2))/(n*(-4*a*c + b**2)*sqrt(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_5_P_x_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_15():
    f = (d*x)**(n/2 - 1)*(-a*h + c*f*x**(n/2) + c*g*x**(3*n/2) + c*h*x**(2*n))/(a + b*x**n + c*x**(2*n))**(sympy.S(3)/2)
    F = -2*x**(1 - n/2)*(d*x)**(n/2 - 1)*(c*x**n*(-b*g + 2*c*f) + c*(-2*a*g + b*f) + h*x**(n/2)*(-4*a*c + b**2))/(n*(-4*a*c + b**2)*sqrt(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_5_P_x_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_16():
    f = (g*x)**m*(a + b*x**n + c*x**(2*n))**p*(a*(m + 1) + b*x**n*(m + n*p + n + 1) + c*x**(2*n)*(m + 2*n*(p + 1) + 1))
    F = (g*x)**(m + 1)*(a + b*x**n + c*x**(2*n))**(p + 1)/g
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_5_P_x_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_17():
    f = (A + B*x**n + C*x**(2*n) + D*x**(3*n))/(a + b*x**n + c*x**(2*n))**2
    F = x*(A*c*(-2*a*c + b**2) - a*(B*b*c - 2*C*a*c + D*a*b) + x**n*(-D*a*b**2 - 2*a*c*(B*c - D*a) + b*c*(A*c + C*a)))/(a*c*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))) + x*(D*a*b**2 + 2*a*c*(B*c*(1 - n) - D*a*(n + 1)) - b*c*(1 - n)*(A*c + C*a) - (A*c**2*(4*a*c*(1 - 2*n) - b**2*(1 - n)) - a*(4*C*a*c**2 - C*b**2*c*(1 - n) + D*b**3 - 2*b*c*(B*c*n + D*a*(n + 2))))/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*c*n*(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + x*(D*a*b**2 + 2*a*c*(B*c*(1 - n) - D*a*(n + 1)) - b*c*(1 - n)*(A*c + C*a) + (A*c**2*(4*a*c*(1 - 2*n) - b**2*(1 - n)) - a*(4*C*a*c**2 - C*b**2*c*(1 - n) + D*b**3 - 2*b*c*(B*c*n + D*a*(n + 2))))/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*c*n*(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F

