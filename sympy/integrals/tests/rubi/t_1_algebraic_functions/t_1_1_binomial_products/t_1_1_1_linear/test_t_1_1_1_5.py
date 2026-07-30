"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.1 Binomial products/1.1.1 Linear/1.1.1.5 P(x) (a+b x)^m (c+d x)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, C, D, a, b, c, d, m, n = symbols('A B C D a b c d m n')

def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_1():
    f = (a + b*x)**3*(A + B*x + C*x**2 + D*x**3)/sqrt(c + d*x)
    F = 2*D*b**3*(c + d*x)**(sympy.S(13)/2)/(13*d**7) + 2*b**2*(c + d*x)**(sympy.S(11)/2)*(C*b*d + 3*D*a*d - 6*D*b*c)/(11*d**7) + 2*b*(c + d*x)**(sympy.S(9)/2)*(3*D*a**2*d**2 + 3*a*b*d*(C*d - 5*D*c) - b**2*(-B*d**2 + 5*C*c*d - 15*D*c**2))/(9*d**7) + (c + d*x)**(sympy.S(7)/2)*(2*D*a**3*d**3 + 6*a**2*b*d**2*(C*d - 4*D*c) - 6*a*b**2*d*(-B*d**2 + 4*C*c*d - 10*D*c**2) + 2*b**3*(A*d**3 - 4*B*c*d**2 + 10*C*c**2*d - 20*D*c**3))/(7*d**7) - (c + d*x)**(sympy.S(5)/2)*(-2*a*d + 2*b*c)*(a**2*d**2*(C*d - 3*D*c) - a*b*d*(-3*B*d**2 + 8*C*c*d - 15*D*c**2) + b**2*(3*A*d**3 - 6*B*c*d**2 + 10*C*c**2*d - 15*D*c**3))/(5*d**7) - 2*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**2*(a*d*(-B*d**2 + 2*C*c*d - 3*D*c**2) - b*(3*A*d**3 - 4*B*c*d**2 + 5*C*c**2*d - 6*D*c**3))/(3*d**7) - 2*sqrt(c + d*x)*(-a*d + b*c)**3*(A*d**3 - B*c*d**2 + C*c**2*d - D*c**3)/d**7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_2():
    f = (a + b*x)**2*(A + B*x + C*x**2 + D*x**3)/sqrt(c + d*x)
    F = 2*D*b**2*(c + d*x)**(sympy.S(11)/2)/(11*d**6) + 2*b*(c + d*x)**(sympy.S(9)/2)*(C*b*d + 2*D*a*d - 5*D*b*c)/(9*d**6) + (c + d*x)**(sympy.S(7)/2)*(2*D*a**2*d**2 + 4*a*b*d*(C*d - 4*D*c) - 2*b**2*(-B*d**2 + 4*C*c*d - 10*D*c**2))/(7*d**6) + (c + d*x)**(sympy.S(5)/2)*(2*a**2*d**2*(C*d - 3*D*c) - 4*a*b*d*(-B*d**2 + 3*C*c*d - 6*D*c**2) + 2*b**2*(A*d**3 - 3*B*c*d**2 + 6*C*c**2*d - 10*D*c**3))/(5*d**6) + (c + d*x)**(sympy.S(3)/2)*(-2*a*d + 2*b*c)*(a*d*(-B*d**2 + 2*C*c*d - 3*D*c**2) - b*(2*A*d**3 - 3*B*c*d**2 + 4*C*c**2*d - 5*D*c**3))/(3*d**6) + 2*sqrt(c + d*x)*(-a*d + b*c)**2*(A*d**3 - B*c*d**2 + C*c**2*d - D*c**3)/d**6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_3():
    f = (a + b*x)*(A + B*x + C*x**2 + D*x**3)/sqrt(c + d*x)
    F = 2*D*b*(c + d*x)**(sympy.S(9)/2)/(9*d**5) + (c + d*x)**(sympy.S(7)/2)*(2*C*b*d + 2*D*a*d - 8*D*b*c)/(7*d**5) + (c + d*x)**(sympy.S(5)/2)*(2*a*d*(C*d - 3*D*c) - 2*b*(-B*d**2 + 3*C*c*d - 6*D*c**2))/(5*d**5) - (c + d*x)**(sympy.S(3)/2)*(2*a*d*(-B*d**2 + 2*C*c*d - 3*D*c**2) - 2*b*(A*d**3 - 2*B*c*d**2 + 3*C*c**2*d - 4*D*c**3))/(3*d**5) - sqrt(c + d*x)*(-2*a*d + 2*b*c)*(A*d**3 - B*c*d**2 + C*c**2*d - D*c**3)/d**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_4():
    f = (A + B*x + C*x**2 + D*x**3)/sqrt(c + d*x)
    F = 2*D*(c + d*x)**(sympy.S(7)/2)/(7*d**4) + (c + d*x)**(sympy.S(5)/2)*(2*C*d - 6*D*c)/(5*d**4) - (c + d*x)**(sympy.S(3)/2)*(-2*B*d**2 + 4*C*c*d - 6*D*c**2)/(3*d**4) + sqrt(c + d*x)*(2*A*d**3 - 2*B*c*d**2 + 2*C*c**2*d - 2*D*c**3)/d**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_5():
    f = (A + B*x + C*x**2 + D*x**3)/((a + b*x)*sqrt(c + d*x))
    F = 2*D*(c + d*x)**(sympy.S(5)/2)/(5*b*d**3) + (c + d*x)**(sympy.S(3)/2)*(2*C*b*d - 2*D*a*d - 4*D*b*c)/(3*b**2*d**3) + sqrt(c + d*x)*(2*D*a**2*d**2 - 2*a*b*d*(C*d - D*c) - 2*b**2*(-B*d**2 + C*c*d - D*c**2))/(b**3*d**3) - (2*A*b**3 - 2*a*(B*b**2 - C*a*b + D*a**2))*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(b**(sympy.S(7)/2)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_6():
    f = (A + B*x + C*x**2 + D*x**3)/((a + b*x)**2*sqrt(c + d*x))
    F = 2*D*(c + d*x)**(sympy.S(3)/2)/(3*b**2*d**2) - (A - a*(B*b**2 - C*a*b + D*a**2)/b**3)*sqrt(c + d*x)/((a + b*x)*(-a*d + b*c)) + sqrt(c + d*x)*(2*C*b*d - 4*D*a*d - 2*D*b*c)/(b**3*d**2) - (-5*D*a**3*d + 3*a**2*b*(C*d + 2*D*c) - a*b**2*(B*d + 4*C*c) + b**3*(-A*d + 2*B*c))*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(b**(sympy.S(7)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_7():
    f = (A + B*x + C*x**2 + D*x**3)/((a + b*x)**3*sqrt(c + d*x))
    F = 2*D*sqrt(c + d*x)/(b**3*d) - sqrt(c + d*x)*(-9*D*a**3*d + a**2*b*(5*C*d + 12*D*c) - a*b**2*(B*d + 8*C*c) + b**3*(-3*A*d + 4*B*c))/(4*b**3*(a + b*x)*(-a*d + b*c)**2) - sqrt(c + d*x)*(A*b**3 - a*(B*b**2 - C*a*b + D*a**2))/(2*b**3*(a + b*x)**2*(-a*d + b*c)) - (-15*D*a**3*d**2 + 3*a**2*b*d*(C*d + 12*D*c) - a*b**2*(-B*d**2 + 8*C*c*d + 24*D*c**2) + b**3*(3*A*d**2 - 4*B*c*d + 8*C*c**2))*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(4*b**(sympy.S(7)/2)*(-a*d + b*c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_8():
    f = (A + B*x + C*x**2 + D*x**3)/((a + b*x)**4*sqrt(c + d*x))
    F = -sqrt(c + d*x)*(-11*D*a**3*d**2 + a**2*b*d*(C*d + 30*D*c) - a*b**2*(-B*d**2 + 4*C*c*d + 24*D*c**2) + b**3*(5*A*d**2 - 6*B*c*d + 8*C*c**2))/(8*b**3*(a + b*x)*(-a*d + b*c)**3) - sqrt(c + d*x)*(-13*D*a**3*d + a**2*b*(7*C*d + 18*D*c) - a*b**2*(B*d + 12*C*c) + b**3*(-5*A*d + 6*B*c))/(12*b**3*(a + b*x)**2*(-a*d + b*c)**2) - sqrt(c + d*x)*(A*b**3 - a*(B*b**2 - C*a*b + D*a**2))/(3*b**3*(a + b*x)**3*(-a*d + b*c)) + (5*D*a**3*d**3 + a**2*b*d**2*(C*d - 18*D*c) - a*b**2*d*(-B*d**2 + 4*C*c*d - 24*D*c**2) + b**3*(5*A*d**3 - 6*B*c*d**2 + 8*C*c**2*d - 16*D*c**3))*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(8*b**(sympy.S(7)/2)*(-a*d + b*c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_9():
    f = (A + B*x + C*x**2 + D*x**3)/((a + b*x)**5*sqrt(c + d*x))
    F = sqrt(c + d*x)*(5*D*a**3*d**3 + 3*a**2*b*d**2*(C*d - 8*D*c) - a*b**2*d*(-5*B*d**2 + 16*C*c*d - 48*D*c**2) + b**3*(35*A*d**3 - 40*B*c*d**2 + 48*C*c**2*d - 64*D*c**3))/(64*b**3*(a + b*x)*(-a*d + b*c)**4) - sqrt(c + d*x)*(-59*D*a**3*d**2 + 3*a**2*b*d*(C*d + 56*D*c) - a*b**2*(-5*B*d**2 + 16*C*c*d + 144*D*c**2) + b**3*(35*A*d**2 - 40*B*c*d + 48*C*c**2))/(96*b**3*(a + b*x)**2*(-a*d + b*c)**3) - sqrt(c + d*x)*(-17*D*a**3*d + 3*a**2*b*(3*C*d + 8*D*c) - a*b**2*(B*d + 16*C*c) + b**3*(-7*A*d + 8*B*c))/(24*b**3*(a + b*x)**3*(-a*d + b*c)**2) - sqrt(c + d*x)*(A*b**3 - a*(B*b**2 - C*a*b + D*a**2))/(4*b**3*(a + b*x)**4*(-a*d + b*c)) - d*(5*D*a**3*d**3 + 3*a**2*b*d**2*(C*d - 8*D*c) - a*b**2*d*(-5*B*d**2 + 16*C*c*d - 48*D*c**2) + b**3*(35*A*d**3 - 40*B*c*d**2 + 48*C*c**2*d - 64*D*c**3))*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(64*b**(sympy.S(7)/2)*(-a*d + b*c)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_10():
    f = (a + b*x)**3*(A + B*x + C*x**2 + D*x**3)/(c + d*x)**(sympy.S(3)/2)
    F = 2*D*b**3*(c + d*x)**(sympy.S(11)/2)/(11*d**7) + 2*b**2*(c + d*x)**(sympy.S(9)/2)*(C*b*d + 3*D*a*d - 6*D*b*c)/(9*d**7) + 2*b*(c + d*x)**(sympy.S(7)/2)*(3*D*a**2*d**2 + 3*a*b*d*(C*d - 5*D*c) - b**2*(-B*d**2 + 5*C*c*d - 15*D*c**2))/(7*d**7) + (c + d*x)**(sympy.S(5)/2)*(2*D*a**3*d**3 + 6*a**2*b*d**2*(C*d - 4*D*c) - 6*a*b**2*d*(-B*d**2 + 4*C*c*d - 10*D*c**2) + 2*b**3*(A*d**3 - 4*B*c*d**2 + 10*C*c**2*d - 20*D*c**3))/(5*d**7) - (c + d*x)**(sympy.S(3)/2)*(-2*a*d + 2*b*c)*(a**2*d**2*(C*d - 3*D*c) - a*b*d*(-3*B*d**2 + 8*C*c*d - 15*D*c**2) + b**2*(3*A*d**3 - 6*B*c*d**2 + 10*C*c**2*d - 15*D*c**3))/(3*d**7) - 2*sqrt(c + d*x)*(-a*d + b*c)**2*(a*d*(-B*d**2 + 2*C*c*d - 3*D*c**2) - b*(3*A*d**3 - 4*B*c*d**2 + 5*C*c**2*d - 6*D*c**3))/d**7 + 2*(-a*d + b*c)**3*(A*d**3 - B*c*d**2 + C*c**2*d - D*c**3)/(d**7*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_11():
    f = (a + b*x)**2*(A + B*x + C*x**2 + D*x**3)/(c + d*x)**(sympy.S(3)/2)
    F = 2*D*b**2*(c + d*x)**(sympy.S(9)/2)/(9*d**6) + 2*b*(c + d*x)**(sympy.S(7)/2)*(C*b*d + 2*D*a*d - 5*D*b*c)/(7*d**6) + (c + d*x)**(sympy.S(5)/2)*(2*D*a**2*d**2 + 4*a*b*d*(C*d - 4*D*c) - 2*b**2*(-B*d**2 + 4*C*c*d - 10*D*c**2))/(5*d**6) + (c + d*x)**(sympy.S(3)/2)*(2*a**2*d**2*(C*d - 3*D*c) - 4*a*b*d*(-B*d**2 + 3*C*c*d - 6*D*c**2) + 2*b**2*(A*d**3 - 3*B*c*d**2 + 6*C*c**2*d - 10*D*c**3))/(3*d**6) + sqrt(c + d*x)*(-2*a*d + 2*b*c)*(a*d*(-B*d**2 + 2*C*c*d - 3*D*c**2) - b*(2*A*d**3 - 3*B*c*d**2 + 4*C*c**2*d - 5*D*c**3))/d**6 - 2*(-a*d + b*c)**2*(A*d**3 - B*c*d**2 + C*c**2*d - D*c**3)/(d**6*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_12():
    f = (a + b*x)*(A + B*x + C*x**2 + D*x**3)/(c + d*x)**(sympy.S(3)/2)
    F = 2*D*b*(c + d*x)**(sympy.S(7)/2)/(7*d**5) + (c + d*x)**(sympy.S(5)/2)*(2*C*b*d + 2*D*a*d - 8*D*b*c)/(5*d**5) + (c + d*x)**(sympy.S(3)/2)*(2*a*d*(C*d - 3*D*c) - 2*b*(-B*d**2 + 3*C*c*d - 6*D*c**2))/(3*d**5) - sqrt(c + d*x)*(2*a*d*(-B*d**2 + 2*C*c*d - 3*D*c**2) - 2*b*(A*d**3 - 2*B*c*d**2 + 3*C*c**2*d - 4*D*c**3))/d**5 + (-2*a*d + 2*b*c)*(A*d**3 - B*c*d**2 + C*c**2*d - D*c**3)/(d**5*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_13():
    f = (A + B*x + C*x**2 + D*x**3)/(c + d*x)**(sympy.S(3)/2)
    F = 2*D*(c + d*x)**(sympy.S(5)/2)/(5*d**4) + (c + d*x)**(sympy.S(3)/2)*(2*C*d - 6*D*c)/(3*d**4) - sqrt(c + d*x)*(-2*B*d**2 + 4*C*c*d - 6*D*c**2)/d**4 - (2*A*d**3 - 2*B*c*d**2 + 2*C*c**2*d - 2*D*c**3)/(d**4*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_14():
    f = (A + B*x + C*x**2 + D*x**3)/((a + b*x)*(c + d*x)**(sympy.S(3)/2))
    F = -2*D*c*sqrt(c + d*x)/(b*d**3) + 2*D*(c + d*x)**(sympy.S(3)/2)/(3*b*d**3) + (2*A*d**3 - 2*B*c*d**2 + 2*C*c**2*d - 2*D*c**3)/(d**3*sqrt(c + d*x)*(-a*d + b*c)) + sqrt(c + d*x)*(2*C*b*d - 2*D*a*d - 2*D*b*c)/(b**2*d**3) - (2*A*b**3 - 2*a*(B*b**2 - C*a*b + D*a**2))*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(b**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_15():
    f = (A + B*x + C*x**2 + D*x**3)/((a + b*x)**2*(c + d*x)**(sympy.S(3)/2))
    F = 2*D*sqrt(c + d*x)/(b**2*d**2) - (A - a*(B*b**2 - C*a*b + D*a**2)/b**3)/((a + b*x)*sqrt(c + d*x)*(-a*d + b*c)) + (B*a*b**2*d**3 - C*a**2*b*d**3 + D*a**3*d**3 - b**3*(3*A*d**3 - 2*B*c*d**2 + 2*C*c**2*d - 2*D*c**3))/(b**3*d**2*sqrt(c + d*x)*(-a*d + b*c)**2) - (-3*D*a**3*d + a**2*b*(C*d + 6*D*c) - a*b**2*(-B*d + 4*C*c) + b**3*(-3*A*d + 2*B*c))*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(b**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_16():
    f = (A + B*x + C*x**2 + D*x**3)/((a + b*x)**3*(c + d*x)**(sympy.S(3)/2))
    F = -sqrt(c + d*x)*(-7*D*a**3*d + 3*a**2*b*(C*d + 4*D*c) - a*b**2*(-B*d + 8*C*c) + b**3*(-5*A*d + 4*B*c))/(4*b**2*(a + b*x)*(-a*d + b*c)**3) - (A*b**3 - a*(B*b**2 - C*a*b + D*a**2))/(2*b**3*(a + b*x)**2*sqrt(c + d*x)*(-a*d + b*c)) - (B*a*b**2*d**3 - C*a**2*b*d**3 + D*a**3*d**3 - b**3*(5*A*d**3 - 4*B*c*d**2 + 4*C*c**2*d - 4*D*c**3))/(2*b**3*d*sqrt(c + d*x)*(-a*d + b*c)**3) - (-3*D*a**3*d**2 - a**2*b*d*(C*d - 12*D*c) + a*b**2*(-3*B*d**2 + 8*C*c*d - 24*D*c**2) + b**3*(15*A*d**2 - 12*B*c*d + 8*C*c**2))*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(4*b**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_17():
    f = (A + B*x + C*x**2 + D*x**3)/((a + b*x)**4*(c + d*x)**(sympy.S(3)/2))
    F = -sqrt(c + d*x)*(5*D*a**3*d**2 - a**2*b*d*(11*C*d - 18*D*c) + a*b**2*(-7*B*d**2 + 36*C*c*d - 72*D*c**2) + b**3*(49*A*d**2 - 42*B*c*d + 24*C*c**2))/(24*b**2*(a + b*x)*(-a*d + b*c)**4) - sqrt(c + d*x)*(-11*D*a**3*d + a**2*b*(5*C*d + 18*D*c) - a*b**2*(-B*d + 12*C*c) + b**3*(-7*A*d + 6*B*c))/(12*b**2*(a + b*x)**2*(-a*d + b*c)**3) + (B*a*b**2*d**3 - C*a**2*b*d**3 + D*a**3*d**3 - b**3*(7*A*d**3 - 6*B*c*d**2 + 6*C*c**2*d - 6*D*c**3))/(3*b**3*sqrt(c + d*x)*(-a*d + b*c)**4) - (A*b**3 - a*(B*b**2 - C*a*b + D*a**2))/(3*b**3*(a + b*x)**3*sqrt(c + d*x)*(-a*d + b*c)) - (D*a**3*d**3 + a**2*b*d**2*(C*d - 6*D*c) - a*b**2*d*(-5*B*d**2 + 12*C*c*d - 24*D*c**2) - b**3*(35*A*d**3 - 30*B*c*d**2 + 24*C*c**2*d - 16*D*c**3))*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(8*b**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_18():
    f = (a + b*x)**3*(A + B*x + C*x**2 + D*x**3)/(c + d*x)**(sympy.S(5)/2)
    F = 2*D*b**3*(c + d*x)**(sympy.S(9)/2)/(9*d**7) + 2*b**2*(c + d*x)**(sympy.S(7)/2)*(C*b*d + 3*D*a*d - 6*D*b*c)/(7*d**7) + 2*b*(c + d*x)**(sympy.S(5)/2)*(3*D*a**2*d**2 + 3*a*b*d*(C*d - 5*D*c) - b**2*(-B*d**2 + 5*C*c*d - 15*D*c**2))/(5*d**7) + (c + d*x)**(sympy.S(3)/2)*(2*D*a**3*d**3 + 6*a**2*b*d**2*(C*d - 4*D*c) - 6*a*b**2*d*(-B*d**2 + 4*C*c*d - 10*D*c**2) + 2*b**3*(A*d**3 - 4*B*c*d**2 + 10*C*c**2*d - 20*D*c**3))/(3*d**7) - sqrt(c + d*x)*(-2*a*d + 2*b*c)*(a**2*d**2*(C*d - 3*D*c) - a*b*d*(-3*B*d**2 + 8*C*c*d - 15*D*c**2) + b**2*(3*A*d**3 - 6*B*c*d**2 + 10*C*c**2*d - 15*D*c**3))/d**7 + 2*(-a*d + b*c)**2*(a*d*(-B*d**2 + 2*C*c*d - 3*D*c**2) - b*(3*A*d**3 - 4*B*c*d**2 + 5*C*c**2*d - 6*D*c**3))/(d**7*sqrt(c + d*x)) + 2*(-a*d + b*c)**3*(A*d**3 - B*c*d**2 + C*c**2*d - D*c**3)/(3*d**7*(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_19():
    f = (a + b*x)**2*(A + B*x + C*x**2 + D*x**3)/(c + d*x)**(sympy.S(5)/2)
    F = 2*D*b**2*(c + d*x)**(sympy.S(7)/2)/(7*d**6) + 2*b*(c + d*x)**(sympy.S(5)/2)*(C*b*d + 2*D*a*d - 5*D*b*c)/(5*d**6) + (c + d*x)**(sympy.S(3)/2)*(2*D*a**2*d**2 + 4*a*b*d*(C*d - 4*D*c) - 2*b**2*(-B*d**2 + 4*C*c*d - 10*D*c**2))/(3*d**6) + sqrt(c + d*x)*(2*a**2*d**2*(C*d - 3*D*c) - 4*a*b*d*(-B*d**2 + 3*C*c*d - 6*D*c**2) + 2*b**2*(A*d**3 - 3*B*c*d**2 + 6*C*c**2*d - 10*D*c**3))/d**6 - (-2*a*d + 2*b*c)*(a*d*(-B*d**2 + 2*C*c*d - 3*D*c**2) - b*(2*A*d**3 - 3*B*c*d**2 + 4*C*c**2*d - 5*D*c**3))/(d**6*sqrt(c + d*x)) - 2*(-a*d + b*c)**2*(A*d**3 - B*c*d**2 + C*c**2*d - D*c**3)/(3*d**6*(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_20():
    f = (a + b*x)*(A + B*x + C*x**2 + D*x**3)/(c + d*x)**(sympy.S(5)/2)
    F = 2*D*b*(c + d*x)**(sympy.S(5)/2)/(5*d**5) + (c + d*x)**(sympy.S(3)/2)*(2*C*b*d + 2*D*a*d - 8*D*b*c)/(3*d**5) + sqrt(c + d*x)*(2*a*d*(C*d - 3*D*c) - 2*b*(-B*d**2 + 3*C*c*d - 6*D*c**2))/d**5 + (2*a*d*(-B*d**2 + 2*C*c*d - 3*D*c**2) - 2*b*(A*d**3 - 2*B*c*d**2 + 3*C*c**2*d - 4*D*c**3))/(d**5*sqrt(c + d*x)) + (-2*a*d + 2*b*c)*(A*d**3 - B*c*d**2 + C*c**2*d - D*c**3)/(3*d**5*(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_21():
    f = (A + B*x + C*x**2 + D*x**3)/(c + d*x)**(sympy.S(5)/2)
    F = 2*D*(c + d*x)**(sympy.S(3)/2)/(3*d**4) + sqrt(c + d*x)*(2*C*d - 6*D*c)/d**4 + (-2*B*d**2 + 4*C*c*d - 6*D*c**2)/(d**4*sqrt(c + d*x)) - (2*A*d**3 - 2*B*c*d**2 + 2*C*c**2*d - 2*D*c**3)/(3*d**4*(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_22():
    f = (A + B*x + C*x**2 + D*x**3)/((a + b*x)*(c + d*x)**(sympy.S(5)/2))
    F = 2*D*sqrt(c + d*x)/(b*d**3) + (2*a*d*(-B*d**2 + 2*C*c*d - 3*D*c**2) - 2*b*(-A*d**3 + C*c**2*d - 2*D*c**3))/(d**3*sqrt(c + d*x)*(-a*d + b*c)**2) + (2*A*d**3 - 2*B*c*d**2 + 2*C*c**2*d - 2*D*c**3)/(3*d**3*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)) - (2*A*b**3 - 2*a*(B*b**2 - C*a*b + D*a**2))*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(b**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_23():
    f = (A + B*x + C*x**2 + D*x**3)/((a + b*x)**2*(c + d*x)**(sympy.S(5)/2))
    F = -(A - a*(B*b**2 - C*a*b + D*a**2)/b**3)/((a + b*x)*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)) - (C*a**2*b*d**3 - D*a**3*d**3 + a*b**2*d*(-3*B*d**2 + 4*C*c*d - 6*D*c**2) - b**3*(-5*A*d**3 + 2*B*c*d**2 - 2*D*c**3))/(b**2*d**2*sqrt(c + d*x)*(-a*d + b*c)**3) + (3*B*a*b**2*d**3 - 3*C*a**2*b*d**3 + 3*D*a**3*d**3 - b**3*(5*A*d**3 - 2*B*c*d**2 + 2*C*c**2*d - 2*D*c**3))/(3*b**3*d**2*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**2) - (-D*a**3*d - a**2*b*(C*d - 6*D*c) - a*b**2*(-3*B*d + 4*C*c) + b**3*(-5*A*d + 2*B*c))*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(b**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_24():
    f = (A + B*x + C*x**2 + D*x**3)/((a + b*x)**3*(c + d*x)**(sympy.S(5)/2))
    F = -sqrt(c + d*x)*(-5*D*a**3*d + a**2*b*(C*d + 12*D*c) - a*b**2*(-3*B*d + 8*C*c) + b**3*(-7*A*d + 4*B*c))/(4*b*(a + b*x)*(-a*d + b*c)**4) + (C*a**2*b*d**2 - D*a**3*d**2 + a*b**2*(-3*B*d**2 + 4*C*c*d - 6*D*c**2) + b**3*(7*A*d**2 - 4*B*c*d + 2*C*c**2))/(b**2*sqrt(c + d*x)*(-a*d + b*c)**4) - (A*b**3 - a*(B*b**2 - C*a*b + D*a**2))/(2*b**3*(a + b*x)**2*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)) - (3*B*a*b**2*d**3 - 3*C*a**2*b*d**3 + 3*D*a**3*d**3 - b**3*(7*A*d**3 - 4*B*c*d**2 + 4*C*c**2*d - 4*D*c**3))/(6*b**3*d*(c + d*x)**(sympy.S(3)/2)*(-a*d + b*c)**3) - (D*a**3*d**2 + 3*a**2*b*d*(C*d - 4*D*c) + 3*a*b**2*(-5*B*d**2 + 8*C*c*d - 8*D*c**2) + b**3*(35*A*d**2 - 20*B*c*d + 8*C*c**2))*atanh(sqrt(b)*sqrt(c + d*x)/sqrt(-a*d + b*c))/(4*b**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_25():
    f = (a + b*x)**3*(c + d*x)**n*(A + B*x + C*x**2 + D*x**3)
    F = D*b**3*(c + d*x)**(n + 7)/(d**7*(n + 7)) + b**2*(c + d*x)**(n + 6)*(C*b*d + 3*D*a*d - 6*D*b*c)/(d**7*(n + 6)) + b*(c + d*x)**(n + 5)*(3*D*a**2*d**2 + 3*a*b*d*(C*d - 5*D*c) - b**2*(-B*d**2 + 5*C*c*d - 15*D*c**2))/(d**7*(n + 5)) - (c + d*x)**(n + 1)*(-a*d + b*c)**3*(A*d**3 - B*c*d**2 + C*c**2*d - D*c**3)/(d**7*(n + 1)) - (c + d*x)**(n + 2)*(-a*d + b*c)**2*(a*d*(-B*d**2 + 2*C*c*d - 3*D*c**2) - b*(3*A*d**3 - 4*B*c*d**2 + 5*C*c**2*d - 6*D*c**3))/(d**7*(n + 2)) - (c + d*x)**(n + 3)*(-a*d + b*c)*(a**2*d**2*(C*d - 3*D*c) - a*b*d*(-3*B*d**2 + 8*C*c*d - 15*D*c**2) + b**2*(3*A*d**3 - 6*B*c*d**2 + 10*C*c**2*d - 15*D*c**3))/(d**7*(n + 3)) + (c + d*x)**(n + 4)*(D*a**3*d**3 + 3*a**2*b*d**2*(C*d - 4*D*c) - 3*a*b**2*d*(-B*d**2 + 4*C*c*d - 10*D*c**2) + b**3*(A*d**3 - 4*B*c*d**2 + 10*C*c**2*d - 20*D*c**3))/(d**7*(n + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_26():
    f = (a + b*x)**2*(c + d*x)**n*(A + B*x + C*x**2 + D*x**3)
    F = D*b**2*(c + d*x)**(n + 6)/(d**6*(n + 6)) + b*(c + d*x)**(n + 5)*(C*b*d + 2*D*a*d - 5*D*b*c)/(d**6*(n + 5)) + (c + d*x)**(n + 1)*(-a*d + b*c)**2*(A*d**3 - B*c*d**2 + C*c**2*d - D*c**3)/(d**6*(n + 1)) + (c + d*x)**(n + 2)*(-a*d + b*c)*(a*d*(-B*d**2 + 2*C*c*d - 3*D*c**2) - b*(2*A*d**3 - 3*B*c*d**2 + 4*C*c**2*d - 5*D*c**3))/(d**6*(n + 2)) + (c + d*x)**(n + 3)*(a**2*d**2*(C*d - 3*D*c) - 2*a*b*d*(-B*d**2 + 3*C*c*d - 6*D*c**2) + b**2*(A*d**3 - 3*B*c*d**2 + 6*C*c**2*d - 10*D*c**3))/(d**6*(n + 3)) + (c + d*x)**(n + 4)*(D*a**2*d**2 + 2*a*b*d*(C*d - 4*D*c) - b**2*(-B*d**2 + 4*C*c*d - 10*D*c**2))/(d**6*(n + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_27():
    f = (a + b*x)*(c + d*x)**n*(A + B*x + C*x**2 + D*x**3)
    F = D*b*(c + d*x)**(n + 5)/(d**5*(n + 5)) - (c + d*x)**(n + 1)*(-a*d + b*c)*(A*d**3 - B*c*d**2 + C*c**2*d - D*c**3)/(d**5*(n + 1)) - (c + d*x)**(n + 2)*(a*d*(-B*d**2 + 2*C*c*d - 3*D*c**2) - b*(A*d**3 - 2*B*c*d**2 + 3*C*c**2*d - 4*D*c**3))/(d**5*(n + 2)) + (c + d*x)**(n + 3)*(a*d*(C*d - 3*D*c) - b*(-B*d**2 + 3*C*c*d - 6*D*c**2))/(d**5*(n + 3)) + (c + d*x)**(n + 4)*(C*b*d + D*a*d - 4*D*b*c)/(d**5*(n + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_28():
    f = (c + d*x)**n*(A + B*x + C*x**2 + D*x**3)
    F = D*(c + d*x)**(n + 4)/(d**4*(n + 4)) + (c + d*x)**(n + 1)*(A*d**3 - B*c*d**2 + C*c**2*d - D*c**3)/(d**4*(n + 1)) - (c + d*x)**(n + 2)*(-B*d**2 + 2*C*c*d - 3*D*c**2)/(d**4*(n + 2)) + (c + d*x)**(n + 3)*(C*d - 3*D*c)/(d**4*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_29():
    f = (c + d*x)**n*(A + B*x + C*x**2 + D*x**3)/(a + b*x)
    F = D*(c + d*x)**(n + 3)/(b*d**3*(n + 3)) + (c + d*x)**(n + 2)*(C*b*d - D*a*d - 2*D*b*c)/(b**2*d**3*(n + 2)) - (c + d*x)**(n + 1)*(A*b**3 - a*(B*b**2 - C*a*b + D*a**2))*hyper((1, n + 1), (n + 2,), b*(c + d*x)/(-a*d + b*c))/(b**3*(n + 1)*(-a*d + b*c)) + (c + d*x)**(n + 1)*(D*a**2*d**2 - a*b*d*(C*d - D*c) - b**2*(-B*d**2 + C*c*d - D*c**2))/(b**3*d**3*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_30():
    f = (c + d*x)**n*(A + B*x + C*x**2 + D*x**3)/(a + b*x)**2
    F = D*(c + d*x)**(n + 2)/(b**2*d**2*(n + 2)) - (A - a*(B*b**2 - C*a*b + D*a**2)/b**3)*(c + d*x)**(n + 1)/((a + b*x)*(-a*d + b*c)) + (c + d*x)**(n + 1)*(D*a**3*d*(n + 3) - a**2*b*(C*d*(n + 2) + 3*D*c) + a*b**2*(B*d*(n + 1) + 2*C*c) - b**3*(A*d*n + B*c))*hyper((1, n + 1), (n + 2,), b*(c + d*x)/(-a*d + b*c))/(b**3*(n + 1)*(-a*d + b*c)**2) + (c + d*x)**(n + 1)*(C*b*d - 2*D*a*d - D*b*c)/(b**3*d**2*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_31():
    f = (c + d*x)**n*(A + B*x + C*x**2 + D*x**3)/(a + b*x)**3
    F = D*(c + d*x)**(n + 1)/(b**3*d*(n + 1)) - (c + d*x)**(n + 1)*(-D*a**3*d**2*(n**2 + 5*n + 6) + a**2*b*d*(n + 2)*(C*d*(n + 1) + 6*D*c) - a*b**2*(B*d**2*n*(n + 1) + 4*C*c*d*(n + 1) + 6*D*c**2) + b**3*(-A*d**2*n*(1 - n) + 2*B*c*d*n + 2*C*c**2))*hyper((1, n + 1), (n + 2,), b*(c + d*x)/(-a*d + b*c))/(2*b**3*(n + 1)*(-a*d + b*c)**3) - (c + d*x)**(n + 1)*(-D*a**3*d*(n + 5) + a**2*b*(C*d*(n + 3) + 6*D*c) - a*b**2*(B*d*(n + 1) + 4*C*c) + b**3*(-A*d*(1 - n) + 2*B*c))/(2*b**3*(a + b*x)*(-a*d + b*c)**2) - (c + d*x)**(n + 1)*(A*b**3 - a*(B*b**2 - C*a*b + D*a**2))/(2*b**3*(a + b*x)**2*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_32():
    f = (A + B*x)*(a + b*x)**m*(c + d*x)**n
    F = B*(a + b*x)**(m + 1)*(c + d*x)**(n + 1)/(b*d*(m + n + 2)) + (a + b*x)**(m + 1)*(c + d*x)**n*(A*b*d*(m + n + 2) - B*(a*d*(n + 1) + b*c*(m + 1)))*hyper((-n, m + 1), (m + 2,), -d*(a + b*x)/(-a*d + b*c))/(b**2*d*(b*(c + d*x)/(-a*d + b*c))**n*(m + 1)*(m + n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_33():
    f = (a + b*x)**m*(c + d*x)**n*(A + B*x + C*x**2)
    F = C*(a + b*x)**(m + 2)*(c + d*x)**(n + 1)/(b**2*d*(m + n + 3)) - (a + b*x)**(m + 1)*(c + d*x)**(n + 1)*(C*a*d*(m + 2*n + 4) + b*(-B*d*(m + n + 3) + C*c*(m + 2)))/(b**2*d**2*(m + n + 2)*(m + n + 3)) - (a + b*x)**(m + 1)*(c + d*x)**n*(d*(m + n + 2)*(-A*b**2*d*(m + n + 3) + C*a**2*d*(n + 1) + C*a*b*c*(m + 2)) - (a*d*(n + 1) + b*c*(m + 1))*(C*a*d*(m + 2*n + 4) + b*(-B*d*(m + n + 3) + C*c*(m + 2))))*hyper((-n, m + 1), (m + 2,), -d*(a + b*x)/(-a*d + b*c))/(b**3*d**2*(b*(c + d*x)/(-a*d + b*c))**n*(m + 1)*(m + n + 2)*(m + n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_1_Linear_1_1_1_5_P_x_a_plus_b_x_pow_m_c_plus_d_x_pow_n_34():
    f = (a + b*x)**m*(c + d*x)**n*(A + B*x + C*x**2 + D*x**3)
    F = D*(a + b*x)**(m + 3)*(c + d*x)**(n + 1)/(b**3*d*(m + n + 4)) - (a + b*x)**(m + 2)*(c + d*x)**(n + 1)*(D*a*d*(2*m + 3*n + 9) + b*(-C*d*(m + n + 4) + D*c*(m + 3)))/(b**3*d**2*(m + n + 3)*(m + n + 4)) + (a + b*x)**(m + 1)*(c + d*x)**(n + 1)*(D*a**2*d**2*(m**2 + m*(3*n + 8) + 3*n**2 + 15*n + 18) + a*b*d*(-C*d*(m**2 + m*(3*n + 8) + 2*n**2 + 12*n + 16) + D*c*(m + 2)*(m + 3*n + 6)) + b**2*(B*d**2*(m**2 + m*(2*n + 7) + n**2 + 7*n + 12) - C*c*d*(m + 2)*(m + n + 4) + D*c**2*(m**2 + 5*m + 6)))/(b**3*d**3*(m + n + 2)*(m + n + 3)*(m + n + 4)) + (a + b*x)**(m + 1)*(c + d*x)**n*(d*(m + n + 2)*(A*b**3*d**2*(m**2 + m*(2*n + 7) + n**2 + 7*n + 12) + D*a**3*d**2*(n + 1)*(m + 2*n + 6) - a**2*b*d*(C*d*(n + 1)*(m + n + 4) - D*c*(m + 2)*(m + 3*n + 6)) + a*b**2*c*(m + 2)*(-C*d*(m + n + 4) + D*c*(m + 3))) - (a*d*(n + 1) + b*c*(m + 1))*(D*a**2*d**2*(m**2 + m*(3*n + 8) + 3*n**2 + 15*n + 18) + a*b*d*(-C*d*(m**2 + m*(3*n + 8) + 2*n**2 + 12*n + 16) + D*c*(m + 2)*(m + 3*n + 6)) + b**2*(B*d**2*(m**2 + m*(2*n + 7) + n**2 + 7*n + 12) - C*c*d*(m + 2)*(m + n + 4) + D*c**2*(m**2 + 5*m + 6))))*hyper((-n, m + 1), (m + 2,), -d*(a + b*x)/(-a*d + b*c))/(b**4*d**3*(b*(c + d*x)/(-a*d + b*c))**n*(m + 1)*(m + n + 2)*(m + n + 3)*(m + n + 4))
    assert integrate(f, x) == F

