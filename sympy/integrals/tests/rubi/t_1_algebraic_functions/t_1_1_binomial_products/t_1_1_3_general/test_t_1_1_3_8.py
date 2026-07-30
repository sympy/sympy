"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.1 Binomial products/1.1.3 General/1.1.3.8 P(x) (c x)^m (a+b x^n)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, C, a, b, c, d, e, f, g, h, i, j, m, n, p = symbols('A B C a b c d e f g h i j m n p')

def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_1():
    f = (c + d*x + e*x**2)/sqrt(a + b*x)
    F = 2*e*(a + b*x)**(sympy.S(5)/2)/(5*b**3) + (a + b*x)**(sympy.S(3)/2)*(-4*a*e + 2*b*d)/(3*b**3) + sqrt(a + b*x)*(2*a**2*e - 2*a*b*d + 2*b**2*c)/b**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_2():
    f = (c + d*x + e*x**2)**2/sqrt(a + b*x)
    F = 2*e**2*(a + b*x)**(sympy.S(9)/2)/(9*b**5) + 4*e*(a + b*x)**(sympy.S(7)/2)*(-2*a*e + b*d)/(7*b**5) - (a + b*x)**(sympy.S(5)/2)*(-12*a**2*e**2 + 12*a*b*d*e - 2*b**2*(2*c*e + d**2))/(5*b**5) + (a + b*x)**(sympy.S(3)/2)*(-8*a*e + 4*b*d)*(a**2*e - a*b*d + b**2*c)/(3*b**5) + 2*sqrt(a + b*x)*(a**2*e - a*b*d + b**2*c)**2/b**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_3():
    f = (c + d*x + e*x**2)**3/sqrt(a + b*x)
    F = 2*e**3*(a + b*x)**(sympy.S(13)/2)/(13*b**7) + 6*e**2*(a + b*x)**(sympy.S(11)/2)*(-2*a*e + b*d)/(11*b**7) - 2*e*(a + b*x)**(sympy.S(9)/2)*(-5*a**2*e**2 + 5*a*b*d*e - b**2*(c*e + d**2))/(3*b**7) - (a + b*x)**(sympy.S(7)/2)*(-4*a*e + 2*b*d)*(-10*a**2*e**2 + 10*a*b*d*e - b**2*(6*c*e + d**2))/(7*b**7) - (a + b*x)**(sympy.S(5)/2)*(6*a**2*e - 6*a*b*d + 6*b**2*c)*(-5*a**2*e**2 + 5*a*b*d*e - b**2*(c*e + d**2))/(5*b**7) + (a + b*x)**(sympy.S(3)/2)*(-4*a*e + 2*b*d)*(a**2*e - a*b*d + b**2*c)**2/b**7 + 2*sqrt(a + b*x)*(a**2*e - a*b*d + b**2*c)**3/b**7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_4():
    f = (c + d*x + e*x**2 + f*x**3)/sqrt(a + b*x)
    F = 2*f*(a + b*x)**(sympy.S(7)/2)/(7*b**4) + (a + b*x)**(sympy.S(5)/2)*(-6*a*f + 2*b*e)/(5*b**4) + (a + b*x)**(sympy.S(3)/2)*(6*a**2*f - 4*a*b*e + 2*b**2*d)/(3*b**4) + sqrt(a + b*x)*(-2*a**3*f + 2*a**2*b*e - 2*a*b**2*d + 2*b**3*c)/b**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_5():
    f = (c + d*x + e*x**2 + f*x**3)**2/sqrt(a + b*x)
    F = 2*f**2*(a + b*x)**(sympy.S(13)/2)/(13*b**7) + 4*f*(a + b*x)**(sympy.S(11)/2)*(-3*a*f + b*e)/(11*b**7) - (a + b*x)**(sympy.S(9)/2)*(-30*a**2*f**2 + 20*a*b*e*f - 2*b**2*(2*d*f + e**2))/(9*b**7) + (a + b*x)**(sympy.S(7)/2)*(-40*a**3*f**2 + 40*a**2*b*e*f - 8*a*b**2*(2*d*f + e**2) + 4*b**3*(c*f + d*e))/(7*b**7) + (a + b*x)**(sympy.S(5)/2)*(30*a**4*f**2 - 40*a**3*b*e*f + 12*a**2*b**2*(2*d*f + e**2) - 12*a*b**3*(c*f + d*e) + 2*b**4*(2*c*e + d**2))/(5*b**7) + (a + b*x)**(sympy.S(3)/2)*(12*a**2*f - 8*a*b*e + 4*b**2*d)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*b**7) + 2*sqrt(a + b*x)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)**2/b**7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_6():
    f = (c + d*x + e*x**2 + f*x**3)**3/sqrt(a + b*x)
    F = 2*f**3*(a + b*x)**(sympy.S(19)/2)/(19*b**10) + 6*f**2*(a + b*x)**(sympy.S(17)/2)*(-3*a*f + b*e)/(17*b**10) - 2*f*(a + b*x)**(sympy.S(15)/2)*(-12*a**2*f**2 + 8*a*b*e*f - b**2*(d*f + e**2))/(5*b**10) + (a + b*x)**(sympy.S(13)/2)*(-168*a**3*f**3 + 168*a**2*b*e*f**2 - 42*a*b**2*f*(d*f + e**2) + 2*b**3*(3*c*f**2 + 6*d*e*f + e**3))/(13*b**10) - (a + b*x)**(sympy.S(11)/2)*(-252*a**4*f**3 + 336*a**3*b*e*f**2 - 126*a**2*b**2*f*(d*f + e**2) + 12*a*b**3*(3*c*f**2 + 6*d*e*f + e**3) - 6*b**4*(2*c*e*f + d**2*f + d*e**2))/(11*b**10) + (a + b*x)**(sympy.S(9)/2)*(-84*a**5*f**3 + 140*a**4*b*e*f**2 - 70*a**3*b**2*f*(d*f + e**2) + 10*a**2*b**3*(3*c*f**2 + 6*d*e*f + e**3) - 10*a*b**4*(2*c*e*f + d**2*f + d*e**2) + 2*b**5*(2*c*d*f + c*e**2 + d**2*e))/(3*b**10) - (a + b*x)**(sympy.S(7)/2)*(-168*a**6*f**3 + 336*a**5*b*e*f**2 - 210*a**4*b**2*f*(d*f + e**2) + 40*a**3*b**3*(3*c*f**2 + 6*d*e*f + e**3) - 60*a**2*b**4*(2*c*e*f + d**2*f + d*e**2) + 24*a*b**5*(2*c*d*f + c*e**2 + d**2*e) - 2*b**6*(3*c**2*f + 6*c*d*e + d**3))/(7*b**10) + (a + b*x)**(sympy.S(5)/2)*(-6*a**3*f + 6*a**2*b*e - 6*a*b**2*d + 6*b**3*c)*(12*a**4*f**2 - 16*a**3*b*e*f + a**2*b**2*(9*d*f + 5*e**2) - a*b**3*(3*c*f + 5*d*e) + b**4*(c*e + d**2))/(5*b**10) + (a + b*x)**(sympy.S(3)/2)*(6*a**2*f - 4*a*b*e + 2*b**2*d)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)**2/b**10 + 2*sqrt(a + b*x)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)**3/b**10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_7():
    f = (c + d*x)/(a + b*x**3)
    F = -(-a**(sympy.S(1)/3)*d/b**(sympy.S(1)/3) + c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)) + (-a**(sympy.S(1)/3)*d + b**(sympy.S(1)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(2)/3)*b**(sympy.S(2)/3)) - sqrt(3)*(a**(sympy.S(1)/3)*d + b**(sympy.S(1)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_8():
    f = (c + d*x)/(a + b*x**3)**2
    F = x*(c + d*x)/(3*a*(a + b*x**3)) + (-a**(sympy.S(1)/3)*d + 2*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(5)/3)*b**(sympy.S(2)/3)) - (-a**(sympy.S(1)/3)*d + 2*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(5)/3)*b**(sympy.S(2)/3)) - sqrt(3)*(a**(sympy.S(1)/3)*d + 2*b**(sympy.S(1)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_9():
    f = (c + d*x)/(a + b*x**3)**3
    F = x*(c + d*x)/(6*a*(a + b*x**3)**2) + x*(5*c + 4*d*x)/(18*a**2*(a + b*x**3)) + (-2*a**(sympy.S(1)/3)*d + 5*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(8)/3)*b**(sympy.S(2)/3)) - (-2*a**(sympy.S(1)/3)*d + 5*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(8)/3)*b**(sympy.S(2)/3)) - sqrt(3)*(2*a**(sympy.S(1)/3)*d + 5*b**(sympy.S(1)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(8)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_10():
    f = (c + d*x)/(a + b*x**3)**4
    F = x*(c + d*x)/(9*a*(a + b*x**3)**3) + x*(8*c + 7*d*x)/(54*a**2*(a + b*x**3)**2) + 2*x*(10*c + 7*d*x)/(81*a**3*(a + b*x**3)) + (-14*a**(sympy.S(1)/3)*d + 40*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(243*a**(sympy.S(11)/3)*b**(sympy.S(2)/3)) - (-7*a**(sympy.S(1)/3)*d + 20*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(243*a**(sympy.S(11)/3)*b**(sympy.S(2)/3)) - sqrt(3)*(14*a**(sympy.S(1)/3)*d + 40*b**(sympy.S(1)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(243*a**(sympy.S(11)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_11():
    f = (a + b*x)/(d + e*x**3)
    F = -(a - b*d**(sympy.S(1)/3)/e**(sympy.S(1)/3))*log(d**(sympy.S(2)/3) - d**(sympy.S(1)/3)*e**(sympy.S(1)/3)*x + e**(sympy.S(2)/3)*x**2)/(6*d**(sympy.S(2)/3)*e**(sympy.S(1)/3)) - (-a*e**(sympy.S(1)/3) + b*d**(sympy.S(1)/3))*log(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)/(3*d**(sympy.S(2)/3)*e**(sympy.S(2)/3)) - sqrt(3)*(a*e**(sympy.S(1)/3) + b*d**(sympy.S(1)/3))*atan(sqrt(3)*(d**(sympy.S(1)/3) - 2*e**(sympy.S(1)/3)*x)/(3*d**(sympy.S(1)/3)))/(3*d**(sympy.S(2)/3)*e**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_12():
    f = (a + b*x)/(d - e*x**3)
    F = -sqrt(3)*(-a*e**(sympy.S(1)/3) + b*d**(sympy.S(1)/3))*atan(sqrt(3)*(d**(sympy.S(1)/3) + 2*e**(sympy.S(1)/3)*x)/(3*d**(sympy.S(1)/3)))/(3*d**(sympy.S(2)/3)*e**(sympy.S(2)/3)) - (a*e**(sympy.S(1)/3) + b*d**(sympy.S(1)/3))*log(d**(sympy.S(1)/3) - e**(sympy.S(1)/3)*x)/(3*d**(sympy.S(2)/3)*e**(sympy.S(2)/3)) + (a*e**(sympy.S(1)/3) + b*d**(sympy.S(1)/3))*log(d**(sympy.S(2)/3) + d**(sympy.S(1)/3)*e**(sympy.S(1)/3)*x + e**(sympy.S(2)/3)*x**2)/(6*d**(sympy.S(2)/3)*e**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_13():
    f = (x + 1)/(x**3 + 1)
    F = -2*sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_14():
    f = (1 - x)/(1 - x**3)
    F = 2*sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_15():
    f = (x + 1)/(1 - x**3)
    F = -2*log(1 - x)/3 + log(x**2 + x + 1)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_16():
    f = (1 - x)/(x**3 + 1)
    F = 2*log(x + 1)/3 - log(x**2 - x + 1)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_17():
    f = (3 - x)/(1 - x**3)
    F = -2*log(1 - x)/3 + log(x**2 + x + 1)/3 + 4*sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_18():
    f = (c + d*x)/(c**3 + d**3*x**3)
    F = -2*sqrt(3)*atan(sqrt(3)*(c - 2*d*x)/(3*c))/(3*c*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_19():
    f = (c - d*x)/(c**3 - d**3*x**3)
    F = 2*sqrt(3)*atan(sqrt(3)*(c + 2*d*x)/(3*c))/(3*c*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_20():
    f = (B*a**(sympy.S(1)/3)*b**(sympy.S(1)/3) + B*b**(sympy.S(2)/3)*x)/(a + b*x**3)
    F = -2*sqrt(3)*B*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_21():
    f = (B*a**(sympy.S(1)/3)*(-b)**(sympy.S(1)/3) - B*x*(-b)**(sympy.S(2)/3))/(a + b*x**3)
    F = 2*sqrt(3)*B*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*x*(-b)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_22():
    f = -C*x**2/(a + b*x**3) + (B*x + C*x**2)/(a + b*x**3)
    F = -B*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)) + B*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)) - sqrt(3)*B*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_23():
    f = -C*x**2/(a + b*x**3) + (A + C*x**2)/(a + b*x**3)
    F = A*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)) - A*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)) - sqrt(3)*A*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_24():
    f = -C*x**2/(a + b*x**3) + (A + B*x + C*x**2)/(a + b*x**3)
    F = -(A - B*a**(sympy.S(1)/3)/b**(sympy.S(1)/3))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)) + (A*b**(sympy.S(1)/3) - B*a**(sympy.S(1)/3))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(2)/3)*b**(sympy.S(2)/3)) - sqrt(3)*(A*b**(sympy.S(1)/3) + B*a**(sympy.S(1)/3))*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_25():
    f = (b*x + c*x**2)/(d + e*x**3)
    F = -b*log(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)/(3*d**(sympy.S(1)/3)*e**(sympy.S(2)/3)) + b*log(d**(sympy.S(2)/3) - d**(sympy.S(1)/3)*e**(sympy.S(1)/3)*x + e**(sympy.S(2)/3)*x**2)/(6*d**(sympy.S(1)/3)*e**(sympy.S(2)/3)) - sqrt(3)*b*atan(sqrt(3)*(d**(sympy.S(1)/3) - 2*e**(sympy.S(1)/3)*x)/(3*d**(sympy.S(1)/3)))/(3*d**(sympy.S(1)/3)*e**(sympy.S(2)/3)) + c*log(d + e*x**3)/(3*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_26():
    f = (a + c*x**2)/(d - e*x**3)
    F = -a*log(d**(sympy.S(1)/3) - e**(sympy.S(1)/3)*x)/(3*d**(sympy.S(2)/3)*e**(sympy.S(1)/3)) + a*log(d**(sympy.S(2)/3) + d**(sympy.S(1)/3)*e**(sympy.S(1)/3)*x + e**(sympy.S(2)/3)*x**2)/(6*d**(sympy.S(2)/3)*e**(sympy.S(1)/3)) + sqrt(3)*a*atan(sqrt(3)*(d**(sympy.S(1)/3) + 2*e**(sympy.S(1)/3)*x)/(3*d**(sympy.S(1)/3)))/(3*d**(sympy.S(2)/3)*e**(sympy.S(1)/3)) - c*log(d - e*x**3)/(3*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_27():
    f = (2*a**2 + b**2*x**2)/(a**3 + b**3*x**3)
    F = log(a + b*x)/b - 2*sqrt(3)*atan(sqrt(3)*(a - 2*b*x)/(3*a))/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_28():
    f = (2*a**2 + b**2*x**2)/(a**3 - b**3*x**3)
    F = -log(a - b*x)/b + 2*sqrt(3)*atan(sqrt(3)*(a + 2*b*x)/(3*a))/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_29():
    f = (C*b**(sympy.S(2)/3)*x**2 + 8*C)/(b*x**3 + 8)
    F = C*log(b**(sympy.S(1)/3)*x + 2)/b**(sympy.S(1)/3) - 2*sqrt(3)*C*atan(sqrt(3)*(-b**(sympy.S(1)/3)*x + 1)/3)/(3*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_30():
    f = (C*a**(sympy.S(2)/3) + 2*C*x**2)/(a + 8*x**3)
    F = C*log(a**(sympy.S(1)/3) + 2*x)/4 - sqrt(3)*C*atan(sqrt(3)*(a**(sympy.S(1)/3) - 4*x)/(3*a**(sympy.S(1)/3)))/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_31():
    f = (C*x**2*(-b)**(sympy.S(2)/3) + 8*C)/(b*x**3 - 8)
    F = -C*log(x*(-b)**(sympy.S(1)/3) + 2)/(-b)**(sympy.S(1)/3) + 2*sqrt(3)*C*atan(sqrt(3)*(-x*(-b)**(sympy.S(1)/3) + 1)/3)/(3*(-b)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_32():
    f = (2*C*x**2 + C*(-a)**(sympy.S(2)/3))/(a - 8*x**3)
    F = -C*log(2*x + (-a)**(sympy.S(1)/3))/4 + sqrt(3)*C*atan(sqrt(3)*(-4*x/(-a)**(sympy.S(1)/3) + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_33():
    f = (C*x**2 + 2*C*(a/b)**(sympy.S(2)/3))/(a + b*x**3)
    F = C*log(x + (a/b)**(sympy.S(1)/3))/b - 2*sqrt(3)*C*atan(sqrt(3)*(-2*x/(a/b)**(sympy.S(1)/3) + 1)/3)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_34():
    f = (C*x**2 + 2*C*(-a/b)**(sympy.S(2)/3))/(a - b*x**3)
    F = -C*log(x + (-a/b)**(sympy.S(1)/3))/b + 2*sqrt(3)*C*atan(sqrt(3)*(-2*x/(-a/b)**(sympy.S(1)/3) + 1)/3)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_35():
    f = (C*x**2 + 2*C*(-a/b)**(sympy.S(2)/3))/(a + b*x**3)
    F = C*log(-x + (-a/b)**(sympy.S(1)/3))/b - 2*sqrt(3)*C*atan(sqrt(3)*(2*x/(-a/b)**(sympy.S(1)/3) + 1)/3)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_36():
    f = (C*x**2 + 2*C*(a/b)**(sympy.S(2)/3))/(a - b*x**3)
    F = -C*log(-x + (a/b)**(sympy.S(1)/3))/b + 2*sqrt(3)*C*atan(sqrt(3)*(2*x/(a/b)**(sympy.S(1)/3) + 1)/3)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_37():
    f = (2*C*a**(sympy.S(2)/3) + C*b**(sympy.S(2)/3)*x**2)/(a + b*x**3)
    F = C*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/b**(sympy.S(1)/3) - 2*sqrt(3)*C*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_38():
    f = (-2*C*a**(sympy.S(2)/3) - C*x**2*(-b)**(sympy.S(2)/3))/(a + b*x**3)
    F = C*log(a**(sympy.S(1)/3) - x*(-b)**(sympy.S(1)/3))/(-b)**(sympy.S(1)/3) - 2*sqrt(3)*C*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*x*(-b)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*(-b)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_39():
    f = (x**2 - 3)/(x**3 - 1)
    F = -2*log(1 - x)/3 + 5*log(x**2 + x + 1)/6 + sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_40():
    f = (B*a**(sympy.S(1)/3)*b**(sympy.S(1)/3) + B*b**(sympy.S(2)/3)*x + 2*C*a**(sympy.S(2)/3) + C*b**(sympy.S(2)/3)*x**2)/(a + b*x**3)
    F = C*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/b**(sympy.S(1)/3) - sqrt(3)*(2*B/a**(sympy.S(1)/3) + 2*C/b**(sympy.S(1)/3))*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_41():
    f = (B*a**(sympy.S(1)/3)*(-b)**(sympy.S(1)/3) - B*x*(-b)**(sympy.S(2)/3) - 2*C*a**(sympy.S(2)/3) - C*x**2*(-b)**(sympy.S(2)/3))/(a + b*x**3)
    F = C*log(a**(sympy.S(1)/3) - x*(-b)**(sympy.S(1)/3))/(-b)**(sympy.S(1)/3) + sqrt(3)*(2*B*b + 2*C*a**(sympy.S(1)/3)*(-b)**(sympy.S(2)/3))*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2*x*(-b)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3)*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_42():
    f = (B**2 + B*C*x + C**2*x**2)/(-B**3 + C**3*x**3)
    F = log(B - C*x)/C
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_43():
    f = (C*a**(sympy.S(2)/3) - C*a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + C*b**(sympy.S(2)/3)*x**2)/(a + b*x**3)
    F = C*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/b**(sympy.S(1)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_44():
    f = (B*x + B*(a/b)**(sympy.S(1)/3) + C*x**2 + 2*C*(a/b)**(sympy.S(2)/3))/(a + b*x**3)
    F = C*log(x + (a/b)**(sympy.S(1)/3))/b - 2*sqrt(3)*(a/b)**(sympy.S(2)/3)*(B + C*(a/b)**(sympy.S(1)/3))*atan(sqrt(3)*(-2*x/(a/b)**(sympy.S(1)/3) + 1)/3)/(3*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_45():
    f = (B*x + B*(-a/b)**(sympy.S(1)/3) + C*x**2 + 2*C*(-a/b)**(sympy.S(2)/3))/(a - b*x**3)
    F = -C*log(x + (-a/b)**(sympy.S(1)/3))/b + sqrt(3)*(2*B + 2*C*(-a/b)**(sympy.S(1)/3))*atan(sqrt(3)*(-2*x/(-a/b)**(sympy.S(1)/3) + 1)/3)/(3*b*(-a/b)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_46():
    f = (B*x - B*(-a/b)**(sympy.S(1)/3) + C*x**2 + 2*C*(-a/b)**(sympy.S(2)/3))/(a + b*x**3)
    F = C*log(-x + (-a/b)**(sympy.S(1)/3))/b + sqrt(3)*(2*B - 2*C*(-a/b)**(sympy.S(1)/3))*atan(sqrt(3)*(2*x/(-a/b)**(sympy.S(1)/3) + 1)/3)/(3*b*(-a/b)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_47():
    f = (B*x - B*(a/b)**(sympy.S(1)/3) + C*x**2 + 2*C*(a/b)**(sympy.S(2)/3))/(a - b*x**3)
    F = -C*log(-x + (a/b)**(sympy.S(1)/3))/b - 2*sqrt(3)*(a/b)**(sympy.S(2)/3)*(B - C*(a/b)**(sympy.S(1)/3))*atan(sqrt(3)*(2*x/(a/b)**(sympy.S(1)/3) + 1)/3)/(3*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_48():
    f = (a*x + a + c*x**2)/(1 - x**3)
    F = (-2*a/3 - c/3)*log(1 - x) + (a/3 - c/3)*log(x**2 + x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_49():
    f = (a + b*x + c*x**2)/(1 - x**3)
    F = sqrt(3)*(a - b)*atan(sqrt(3)*(2*x + 1)/3)/3 + (a/6 + b/6 - c/3)*log(x**2 + x + 1) - (a/3 + b/3 + c/3)*log(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_50():
    f = (x**2 + x + 1)/(1 - x**3)
    F = -log(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_51():
    f = (3*x**2 - x + 1)/(1 - x**3)
    F = -log(1 - x**3) + 2*sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_52():
    f = (4*x**2 + x + 1)/(1 - x**3)
    F = -2*log(1 - x) - log(x**2 + x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_53():
    f = (a + b*x**3)**3*(a*c + a*d*x + b*c*x**3 + b*d*x**4)
    F = a**4*c*x + a**4*d*x**2/2 + a**3*b*c*x**4 + 4*a**3*b*d*x**5/5 + 6*a**2*b**2*c*x**7/7 + 3*a**2*b**2*d*x**8/4 + 2*a*b**3*c*x**10/5 + 4*a*b**3*d*x**11/11 + b**4*c*x**13/13 + b**4*d*x**14/14
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_54():
    f = (a + b*x**3)**2*(a*c + a*d*x + b*c*x**3 + b*d*x**4)
    F = a**3*c*x + a**3*d*x**2/2 + 3*a**2*b*c*x**4/4 + 3*a**2*b*d*x**5/5 + 3*a*b**2*c*x**7/7 + 3*a*b**2*d*x**8/8 + b**3*c*x**10/10 + b**3*d*x**11/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_55():
    f = (a + b*x**3)*(a*c + a*d*x + b*c*x**3 + b*d*x**4)
    F = a**2*c*x + a**2*d*x**2/2 + a*b*c*x**4/2 + 2*a*b*d*x**5/5 + b**2*c*x**7/7 + b**2*d*x**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_56():
    f = (a*c + a*d*x + b*c*x**3 + b*d*x**4)/(a + b*x**3)
    F = c*x + d*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_57():
    f = (a*c + a*d*x + b*c*x**3 + b*d*x**4)/(a + b*x**3)**2
    F = -(-a**(sympy.S(1)/3)*d/b**(sympy.S(1)/3) + c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)) + (-a**(sympy.S(1)/3)*d + b**(sympy.S(1)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(2)/3)*b**(sympy.S(2)/3)) - sqrt(3)*(a**(sympy.S(1)/3)*d + b**(sympy.S(1)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_58():
    f = (a*c + a*d*x + b*c*x**3 + b*d*x**4)/(a + b*x**3)**3
    F = x*(c + d*x)/(3*a*(a + b*x**3)) + (-a**(sympy.S(1)/3)*d + 2*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(5)/3)*b**(sympy.S(2)/3)) - (-a**(sympy.S(1)/3)*d + 2*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(5)/3)*b**(sympy.S(2)/3)) - sqrt(3)*(a**(sympy.S(1)/3)*d + 2*b**(sympy.S(1)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_59():
    f = (a + b*x**3)**(sympy.S(3)/2)*(a*c + a*d*x + b*c*x**3 + b*d*x**4)
    F = -405*3**(sympy.S(1)/4)*a**(sympy.S(10)/3)*d*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1729*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 810*a**3*d*sqrt(a + b*x**3)/(1729*b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 54*3**(sympy.S(3)/4)*a**3*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*d*(935 - 935*sqrt(3)) + 1729*b**(sympy.S(1)/3)*c)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(323323*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 54*a**2*sqrt(a + b*x**3)*(1729*c*x + 935*d*x**2)/323323 + 30*a*(a + b*x**3)**(sympy.S(3)/2)*(247*c*x + 187*d*x**2)/46189 + (a + b*x**3)**(sympy.S(5)/2)*(2*c*x/17 + 2*d*x**2/19)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_60():
    f = sqrt(a + b*x**3)*(a*c + a*d*x + b*c*x**3 + b*d*x**4)
    F = -27*3**(sympy.S(1)/4)*a**(sympy.S(7)/3)*d*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(91*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 54*a**2*d*sqrt(a + b*x**3)/(91*b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 18*3**(sympy.S(3)/4)*a**2*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*d*(55 - 55*sqrt(3)) + 91*b**(sympy.S(1)/3)*c)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(5005*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 18*a*sqrt(a + b*x**3)*(91*c*x + 55*d*x**2)/5005 + (a + b*x**3)**(sympy.S(3)/2)*(2*c*x/11 + 2*d*x**2/13)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_61():
    f = (a*c + a*d*x + b*c*x**3 + b*d*x**4)/sqrt(a + b*x**3)
    F = -3*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*d*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(7*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 6*a*d*sqrt(a + b*x**3)/(7*b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*3**(sympy.S(3)/4)*a*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*d*(5 - 5*sqrt(3)) + 7*b**(sympy.S(1)/3)*c)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(35*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + sqrt(a + b*x**3)*(2*c*x/5 + 2*d*x**2/7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_62():
    f = (a*c + a*d*x + b*c*x**3 + b*d*x**4)/(a + b*x**3)**(sympy.S(3)/2)
    F = -3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*d*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*d*sqrt(a + b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*d*(1 - sqrt(3)) + b**(sympy.S(1)/3)*c)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_63():
    f = (a*c + a*d*x + b*c*x**3 + b*d*x**4)/(a + b*x**3)**(sympy.S(5)/2)
    F = 2*x*(c + d*x)/(3*a*sqrt(a + b*x**3)) - 2*d*sqrt(a + b*x**3)/(3*a*b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*d*(1 - sqrt(3)) + b**(sympy.S(1)/3)*c)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*a*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 3**(sympy.S(1)/4)*d*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*a**(sympy.S(2)/3)*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_64():
    f = (a*c + a*d*x + b*c*x**3 + b*d*x**4)/(a + b*x**3)**(sympy.S(7)/2)
    F = 2*x*(c + d*x)/(9*a*(a + b*x**3)**(sympy.S(3)/2)) + 2*x*(7*c + 5*d*x)/(27*a**2*sqrt(a + b*x**3)) - 10*d*sqrt(a + b*x**3)/(27*a**2*b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*d*(5 - 5*sqrt(3)) + 7*b**(sympy.S(1)/3)*c)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(81*a**2*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 5*3**(sympy.S(1)/4)*d*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(27*a**(sympy.S(5)/3)*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_65():
    f = (a*c + a*d*x + b*c*x**3 + b*d*x**4)/(a + b*x**3)**(sympy.S(9)/2)
    F = 2*x*(c + d*x)/(15*a*(a + b*x**3)**(sympy.S(5)/2)) + 2*x*(13*c + 11*d*x)/(135*a**2*(a + b*x**3)**(sympy.S(3)/2)) + 2*x*(91*c + 55*d*x)/(405*a**3*sqrt(a + b*x**3)) - 22*d*sqrt(a + b*x**3)/(81*a**3*b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*d*(55 - 55*sqrt(3)) + 91*b**(sympy.S(1)/3)*c)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1215*a**3*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 11*3**(sympy.S(1)/4)*d*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(81*a**(sympy.S(8)/3)*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_66():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4)/sqrt(a + b*x**3)
    F = -3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-4*a*g + 7*b*d)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(7*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*e*sqrt(a + b*x**3)/(3*b) + 2*f*x*sqrt(a + b*x**3)/(5*b) + 2*g*x**2*sqrt(a + b*x**3)/(7*b) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*(5 - 5*sqrt(3))*(-4*a*g + 7*b*d) + 7*b**(sympy.S(1)/3)*(-2*a*f + 5*b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(105*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + sqrt(a + b*x**3)*(-8*a*g + 14*b*d)/(7*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_67():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4)/(a + b*x**3)**(sympy.S(3)/2)
    F = -2*e*sqrt(a + b*x**3)/(3*a*b) + 2*x*(-a*f + b*c + b*e*x**2 + x*(-a*g + b*d))/(3*a*b*sqrt(a + b*x**3)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*(1 - sqrt(3))*(-4*a*g + b*d) + b**(sympy.S(1)/3)*(2*a*f + b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*a*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - sqrt(a + b*x**3)*(-8*a*g + 2*b*d)/(3*a*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-4*a*g + b*d)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*a**(sympy.S(2)/3)*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_68():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4)/(a + b*x**3)**(sympy.S(5)/2)
    F = 2*x*(-a*f + b*c + b*e*x**2 + x*(-a*g + b*d))/(9*a*b*(a + b*x**3)**(sympy.S(3)/2)) - (6*a*e - 2*x*(2*a*f + 7*b*c + x*(4*a*g + 5*b*d)))/(27*a**2*b*sqrt(a + b*x**3)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*(1 - sqrt(3))*(4*a*g + 5*b*d) + b**(sympy.S(1)/3)*(2*a*f + 7*b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(81*a**2*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - sqrt(a + b*x**3)*(8*a*g + 10*b*d)/(27*a**2*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(4*a*g + 5*b*d)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(27*a**(sympy.S(5)/3)*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_69():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4)/(a + b*x**3)**(sympy.S(7)/2)
    F = 2*x*(-a*f + b*c + b*e*x**2 + x*(-a*g + b*d))/(15*a*b*(a + b*x**3)**(sympy.S(5)/2)) - (18*a*e - 2*x*(2*a*f + 13*b*c + x*(4*a*g + 11*b*d)))/(135*a**2*b*(a + b*x**3)**(sympy.S(3)/2)) + 2*x*(14*a*f + 91*b*c + x*(20*a*g + 55*b*d))/(405*a**3*b*sqrt(a + b*x**3)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*(5 - 5*sqrt(3))*(4*a*g + 11*b*d) + 7*b**(sympy.S(1)/3)*(2*a*f + 13*b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1215*a**3*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - sqrt(a + b*x**3)*(8*a*g + 22*b*d)/(81*a**3*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(4*a*g + 11*b*d)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(81*a**(sympy.S(8)/3)*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_70():
    f = (a + b*x)**2/(c + d*x**3)
    F = -a*(-a*d**(sympy.S(1)/3) + 2*b*c**(sympy.S(1)/3))*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(2)/3)*d**(sympy.S(2)/3)) + a*(-a*d**(sympy.S(1)/3) + 2*b*c**(sympy.S(1)/3))*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(6*c**(sympy.S(2)/3)*d**(sympy.S(2)/3)) - sqrt(3)*a*(a*d**(sympy.S(1)/3) + 2*b*c**(sympy.S(1)/3))*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(3*c**(sympy.S(2)/3)*d**(sympy.S(2)/3)) + b**2*log(c + d*x**3)/(3*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_71():
    f = (a + b*x)**3/(c + d*x**3)
    F = a*b**2*log(c + d*x**3)/d + b**3*x/d + sqrt(3)*(-a**3*d - 3*a**2*b*c**(sympy.S(1)/3)*d**(sympy.S(2)/3) + b**3*c)*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(3*c**(sympy.S(2)/3)*d**(sympy.S(4)/3)) - (-a**3*d + 3*a**2*b*c**(sympy.S(1)/3)*d**(sympy.S(2)/3) + b**3*c)*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(2)/3)*d**(sympy.S(4)/3)) + (-a**3*d + 3*a**2*b*c**(sympy.S(1)/3)*d**(sympy.S(2)/3) + b**3*c)*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(6*c**(sympy.S(2)/3)*d**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_72():
    f = (a + b*x)**4/(c + d*x**3)
    F = 2*a**2*b**2*log(c + d*x**3)/d + 4*a*b**3*x/d + b**4*x**2/(2*d) + (b*c**(sympy.S(1)/3)*(-4*a**3*d + b**3*c) - d**(sympy.S(1)/3)*(-a**4*d + 4*a*b**3*c))*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(2)/3)*d**(sympy.S(5)/3)) - (b*c**(sympy.S(1)/3)*(-4*a**3*d + b**3*c) - d**(sympy.S(1)/3)*(-a**4*d + 4*a*b**3*c))*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(6*c**(sympy.S(2)/3)*d**(sympy.S(5)/3)) + sqrt(3)*(-a**4*d**(sympy.S(4)/3) - 4*a**3*b*c**(sympy.S(1)/3)*d + 4*a*b**3*c*d**(sympy.S(1)/3) + b**4*c**(sympy.S(4)/3))*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(3*c**(sympy.S(2)/3)*d**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_73():
    f = (a + b*x + c*x**2)**2/(d + e*x**3)
    F = 2*b*c*x/e + c**2*x**2/(2*e) + (2*a*c + b**2)*log(d + e*x**3)/(3*e) - (-d**(sympy.S(1)/3)*(-2*a*b*e + c**2*d) + e**(sympy.S(1)/3)*(-a**2*e + 2*b*c*d))*log(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)/(3*d**(sympy.S(2)/3)*e**(sympy.S(5)/3)) + (-d**(sympy.S(1)/3)*(-2*a*b*e + c**2*d) + e**(sympy.S(1)/3)*(-a**2*e + 2*b*c*d))*log(d**(sympy.S(2)/3) - d**(sympy.S(1)/3)*e**(sympy.S(1)/3)*x + e**(sympy.S(2)/3)*x**2)/(6*d**(sympy.S(2)/3)*e**(sympy.S(5)/3)) + sqrt(3)*(-a*e*(a*e**(sympy.S(1)/3) + 2*b*d**(sympy.S(1)/3)) + 2*b*c*d*e**(sympy.S(1)/3) + c**2*d**(sympy.S(4)/3))*atan(sqrt(3)*(d**(sympy.S(1)/3) - 2*e**(sympy.S(1)/3)*x)/(3*d**(sympy.S(1)/3)))/(3*d**(sympy.S(2)/3)*e**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_74():
    f = (a + b*x + c*x**2)**3/(d + e*x**3)
    F = b*c**2*x**3/e + c**3*x**4/(4*e) + 3*c*x**2*(a*c + b**2)/(2*e) - x*(-6*a*b*c*e - b**3*e + c**3*d)/e**2 - (-a**2*c*e - a*b**2*e + b*c**2*d)*log(d + e*x**3)/e**2 + (-6*a*b*c*d*e + c**3*d**2 + 3*d**(sympy.S(1)/3)*e**(sympy.S(2)/3)*(-a**2*b*e + a*c**2*d + b**2*c*d) - e*(-a**3*e + b**3*d))*log(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)/(3*d**(sympy.S(2)/3)*e**(sympy.S(7)/3)) - (-6*a*b*c*d*e + c**3*d**2 + 3*d**(sympy.S(1)/3)*e**(sympy.S(2)/3)*(-a**2*b*e + a*c**2*d + b**2*c*d) - e*(-a**3*e + b**3*d))*log(d**(sympy.S(2)/3) - d**(sympy.S(1)/3)*e**(sympy.S(1)/3)*x + e**(sympy.S(2)/3)*x**2)/(6*d**(sympy.S(2)/3)*e**(sympy.S(7)/3)) - sqrt(3)*(a**3*e**2 + 3*a**2*b*d**(sympy.S(1)/3)*e**(sympy.S(5)/3) - 6*a*b*c*d*e - 3*a*c**2*d**(sympy.S(4)/3)*e**(sympy.S(2)/3) - b**3*d*e - 3*b**2*c*d**(sympy.S(4)/3)*e**(sympy.S(2)/3) + c**3*d**2)*atan(sqrt(3)*(d**(sympy.S(1)/3) - 2*e**(sympy.S(1)/3)*x)/(3*d**(sympy.S(1)/3)))/(3*d**(sympy.S(2)/3)*e**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_75():
    f = (a + b*x + c*x**2)**4/(d + e*x**3)
    F = 4*b*c**3*x**5/(5*e) + c**4*x**6/(6*e) + c**2*x**4*(2*a*c + 3*b**2)/(2*e) - c*x**3*(-12*a*b*c*e - 4*b**3*e + c**3*d)/(3*e**2) - x**2*(-6*a**2*c**2*e - 12*a*b**2*c*e - b**4*e + 4*b*c**3*d)/(2*e**2) - x*(-12*a**2*b*c*e - 4*a*b**3*e + 4*a*c**3*d + 6*b**2*c**2*d)/e**2 + (6*a**2*b**2*e**2 - 12*a*b*c**2*d*e + c**4*d**2 - 4*c*e*(-a**3*e + b**3*d))*log(d + e*x**3)/(3*e**3) - sqrt(3)*(a*e**(sympy.S(1)/3) + b*d**(sympy.S(1)/3))*(-12*a*b*c*d*e + 4*c**3*d**2 + 6*c**2*(-a*d**(sympy.S(4)/3)*e**(sympy.S(2)/3) + b*d**(sympy.S(5)/3)*e**(sympy.S(1)/3)) - e*(-a**3*e - 3*a**2*b*d**(sympy.S(1)/3)*e**(sympy.S(2)/3) + 3*a*b**2*d**(sympy.S(2)/3)*e**(sympy.S(1)/3) + b**3*d))*atan(sqrt(3)*(d**(sympy.S(1)/3) - 2*e**(sympy.S(1)/3)*x)/(3*d**(sympy.S(1)/3)))/(3*d**(sympy.S(2)/3)*e**(sympy.S(8)/3)) + (d**(sympy.S(1)/3)*(6*a**2*c**2*d*e + 12*a*b**2*c*d*e + b**4*d*e - 4*b*(a**3*e**2 + c**3*d**2)) + e**(sympy.S(1)/3)*(a**4*e**2 - 12*a**2*b*c*d*e - 4*a*b**3*d*e + 4*a*c**3*d**2 + 6*b**2*c**2*d**2))*log(d**(sympy.S(1)/3) + e**(sympy.S(1)/3)*x)/(3*d**(sympy.S(2)/3)*e**(sympy.S(8)/3)) - (d**(sympy.S(1)/3)*(6*a**2*c**2*d*e + 12*a*b**2*c*d*e + b**4*d*e - 4*b*(a**3*e**2 + c**3*d**2)) + e**(sympy.S(1)/3)*(a**4*e**2 - 12*a**2*b*c*d*e - 4*a*b**3*d*e + 4*a*c**3*d**2 + 6*b**2*c**2*d**2))*log(d**(sympy.S(2)/3) - d**(sympy.S(1)/3)*e**(sympy.S(1)/3)*x + e**(sympy.S(2)/3)*x**2)/(6*d**(sympy.S(2)/3)*e**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_76():
    f = (x**4 + 2*x**2)/(x**3 + 1)
    F = x**2/2 + log(x + 1) + log(x**2 - x + 1)/2 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_77():
    f = (x**4 + 2*x**2)/(1 - x**3)
    F = -x**2/2 - log(1 - x) - log(x**2 + x + 1)/2 - sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_78():
    f = (4*x**3 - x + 1)/(x**3 + 1)
    F = 4*x - 2*log(x + 1)/3 + log(x**2 - x + 1)/3 + 4*sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_79():
    f = (x + 1 + sqrt(3))/sqrt(x**3 + 1)
    F = 2*sqrt(x**3 + 1)/(x + 1 + sqrt(3)) - 3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(2 - sqrt(3))*(x + 1)*elliptic_e(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1)) + 4*3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_80():
    f = (-x + 1 + sqrt(3))/sqrt(1 - x**3)
    F = -2*sqrt(1 - x**3)/(-x + 1 + sqrt(3)) + 3**(sympy.S(1)/4)*sqrt((x**2 + x + 1)/(-x + 1 + sqrt(3))**2)*(1 - x)*sqrt(2 - sqrt(3))*elliptic_e(asin((-x - sqrt(3) + 1)/(-x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(sqrt((1 - x)/(-x + 1 + sqrt(3))**2)*sqrt(1 - x**3)) - 4*3**(sympy.S(1)/4)*sqrt((x**2 + x + 1)/(-x + 1 + sqrt(3))**2)*(1 - x)*sqrt(sqrt(3) + 2)*elliptic_f(asin((-x - sqrt(3) + 1)/(-x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(sqrt((1 - x)/(-x + 1 + sqrt(3))**2)*sqrt(1 - x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_81():
    f = (-x + 1 + sqrt(3))/sqrt(x**3 - 1)
    F = -3**(sympy.S(1)/4)*sqrt((x**2 + x + 1)/(-x - sqrt(3) + 1)**2)*(1 - x)*sqrt(sqrt(3) + 2)*elliptic_e(asin((-x + 1 + sqrt(3))/(-x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(sqrt(-(1 - x)/(-x - sqrt(3) + 1)**2)*sqrt(x**3 - 1)) + 2*sqrt(x**3 - 1)/(-x - sqrt(3) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_82():
    f = (x + 1 + sqrt(3))/sqrt(-x**3 - 1)
    F = 3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x - sqrt(3) + 1)**2)*sqrt(sqrt(3) + 2)*(x + 1)*elliptic_e(asin((x + 1 + sqrt(3))/(x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(sqrt(-(x + 1)/(x - sqrt(3) + 1)**2)*sqrt(-x**3 - 1)) - 2*sqrt(-x**3 - 1)/(x - sqrt(3) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_83():
    f = (a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/sqrt(a + b*x**3)
    F = -3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(b**(sympy.S(1)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 4*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(b**(sympy.S(1)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*sqrt(a + b*x**3)/(b**(sympy.S(1)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_84():
    f = (a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/sqrt(a - b*x**3)
    F = 3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(b**(sympy.S(1)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(a - b*x**3)) - 4*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(b**(sympy.S(1)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(a - b*x**3)) - 2*sqrt(a - b*x**3)/(b**(sympy.S(1)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_85():
    f = (a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/sqrt(-a + b*x**3)
    F = -3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(b**(sympy.S(1)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(-a + b*x**3)) + 2*sqrt(-a + b*x**3)/(b**(sympy.S(1)/3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_86():
    f = (a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/sqrt(-a - b*x**3)
    F = 3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(b**(sympy.S(1)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(-a - b*x**3)) - 2*sqrt(-a - b*x**3)/(b**(sympy.S(1)/3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_87():
    f = (x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3))/sqrt(a + b*x**3)
    F = -3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*(b/a)**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*(b/a)**(sympy.S(1)/3)*sqrt(a + b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*(b/a)**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*(1 + sqrt(3)))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_88():
    f = (-x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3))/sqrt(a - b*x**3)
    F = 3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*(b/a)**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(a - b*x**3)) - 2*(b/a)**(sympy.S(1)/3)*sqrt(a - b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)) - 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*(b/a)**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*(1 + sqrt(3)))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(a - b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_89():
    f = (-x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3))/sqrt(-a + b*x**3)
    F = -3**(sympy.S(1)/4)*sqrt((x**2*(b/a)**(sympy.S(2)/3) + x*(b/a)**(sympy.S(1)/3) + 1)/(-x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*sqrt(sqrt(3) + 2)*(-x*(b/a)**(sympy.S(1)/3) + 1)*elliptic_e(asin((-x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3))/(-x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/((b/a)**(sympy.S(1)/3)*sqrt(-(-x*(b/a)**(sympy.S(1)/3) + 1)/(-x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*sqrt(-a + b*x**3)) + 2*(b/a)**(sympy.S(2)/3)*sqrt(-a + b*x**3)/(b*(-x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_90():
    f = (x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3))/sqrt(-a - b*x**3)
    F = 3**(sympy.S(1)/4)*sqrt((x**2*(b/a)**(sympy.S(2)/3) - x*(b/a)**(sympy.S(1)/3) + 1)/(x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*sqrt(sqrt(3) + 2)*(x*(b/a)**(sympy.S(1)/3) + 1)*elliptic_e(asin((x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3))/(x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1)), -7 + 4*sqrt(3))/((b/a)**(sympy.S(1)/3)*sqrt(-(x*(b/a)**(sympy.S(1)/3) + 1)/(x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1)**2)*sqrt(-a - b*x**3)) - 2*(b/a)**(sympy.S(2)/3)*sqrt(-a - b*x**3)/(b*(x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_91():
    f = (x - sqrt(3) + 1)/sqrt(x**3 + 1)
    F = 2*sqrt(x**3 + 1)/(x + 1 + sqrt(3)) - 3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(2 - sqrt(3))*(x + 1)*elliptic_e(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_92():
    f = (-x - sqrt(3) + 1)/sqrt(1 - x**3)
    F = -2*sqrt(1 - x**3)/(-x + 1 + sqrt(3)) + 3**(sympy.S(1)/4)*sqrt((x**2 + x + 1)/(-x + 1 + sqrt(3))**2)*(1 - x)*sqrt(2 - sqrt(3))*elliptic_e(asin((-x - sqrt(3) + 1)/(-x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(sqrt((1 - x)/(-x + 1 + sqrt(3))**2)*sqrt(1 - x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_93():
    f = (-x - sqrt(3) + 1)/sqrt(x**3 - 1)
    F = -3**(sympy.S(1)/4)*sqrt((x**2 + x + 1)/(-x - sqrt(3) + 1)**2)*(1 - x)*sqrt(sqrt(3) + 2)*elliptic_e(asin((-x + 1 + sqrt(3))/(-x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(sqrt(-(1 - x)/(-x - sqrt(3) + 1)**2)*sqrt(x**3 - 1)) + 4*3**(sympy.S(1)/4)*sqrt((x**2 + x + 1)/(-x - sqrt(3) + 1)**2)*(1 - x)*sqrt(2 - sqrt(3))*elliptic_f(asin((-x + 1 + sqrt(3))/(-x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(sqrt(-(1 - x)/(-x - sqrt(3) + 1)**2)*sqrt(x**3 - 1)) + 2*sqrt(x**3 - 1)/(-x - sqrt(3) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_94():
    f = (x - sqrt(3) + 1)/sqrt(-x**3 - 1)
    F = 3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x - sqrt(3) + 1)**2)*sqrt(sqrt(3) + 2)*(x + 1)*elliptic_e(asin((x + 1 + sqrt(3))/(x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(sqrt(-(x + 1)/(x - sqrt(3) + 1)**2)*sqrt(-x**3 - 1)) - 4*3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x - sqrt(3) + 1)**2)*sqrt(2 - sqrt(3))*(x + 1)*elliptic_f(asin((x + 1 + sqrt(3))/(x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(sqrt(-(x + 1)/(x - sqrt(3) + 1)**2)*sqrt(-x**3 - 1)) - 2*sqrt(-x**3 - 1)/(x - sqrt(3) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_95():
    f = (-x - 1 + sqrt(3))/sqrt(x**3 + 1)
    F = -2*sqrt(x**3 + 1)/(x + 1 + sqrt(3)) + 3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(2 - sqrt(3))*(x + 1)*elliptic_e(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_96():
    f = (x - 1 + sqrt(3))/sqrt(1 - x**3)
    F = 2*sqrt(1 - x**3)/(-x + 1 + sqrt(3)) - 3**(sympy.S(1)/4)*sqrt((x**2 + x + 1)/(-x + 1 + sqrt(3))**2)*(1 - x)*sqrt(2 - sqrt(3))*elliptic_e(asin((-x - sqrt(3) + 1)/(-x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(sqrt((1 - x)/(-x + 1 + sqrt(3))**2)*sqrt(1 - x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_97():
    f = (x - 1 + sqrt(3))/sqrt(x**3 - 1)
    F = 3**(sympy.S(1)/4)*sqrt((x**2 + x + 1)/(-x - sqrt(3) + 1)**2)*(1 - x)*sqrt(sqrt(3) + 2)*elliptic_e(asin((-x + 1 + sqrt(3))/(-x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(sqrt(-(1 - x)/(-x - sqrt(3) + 1)**2)*sqrt(x**3 - 1)) - 4*3**(sympy.S(1)/4)*sqrt((x**2 + x + 1)/(-x - sqrt(3) + 1)**2)*(1 - x)*sqrt(2 - sqrt(3))*elliptic_f(asin((-x + 1 + sqrt(3))/(-x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(sqrt(-(1 - x)/(-x - sqrt(3) + 1)**2)*sqrt(x**3 - 1)) - 2*sqrt(x**3 - 1)/(-x - sqrt(3) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_98():
    f = (-x - 1 + sqrt(3))/sqrt(-x**3 - 1)
    F = -3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x - sqrt(3) + 1)**2)*sqrt(sqrt(3) + 2)*(x + 1)*elliptic_e(asin((x + 1 + sqrt(3))/(x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(sqrt(-(x + 1)/(x - sqrt(3) + 1)**2)*sqrt(-x**3 - 1)) + 4*3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x - sqrt(3) + 1)**2)*sqrt(2 - sqrt(3))*(x + 1)*elliptic_f(asin((x + 1 + sqrt(3))/(x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(sqrt(-(x + 1)/(x - sqrt(3) + 1)**2)*sqrt(-x**3 - 1)) + 2*sqrt(-x**3 - 1)/(x - sqrt(3) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_99():
    f = (a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/sqrt(a + b*x**3)
    F = -3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(b**(sympy.S(1)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*sqrt(a + b*x**3)/(b**(sympy.S(1)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_100():
    f = (a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/sqrt(a - b*x**3)
    F = 3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(b**(sympy.S(1)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(a - b*x**3)) - 2*sqrt(a - b*x**3)/(b**(sympy.S(1)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_101():
    f = (a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/sqrt(-a + b*x**3)
    F = -3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(b**(sympy.S(1)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(-a + b*x**3)) + 4*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(b**(sympy.S(1)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(-a + b*x**3)) + 2*sqrt(-a + b*x**3)/(b**(sympy.S(1)/3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_102():
    f = (a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/sqrt(-a - b*x**3)
    F = 3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(b**(sympy.S(1)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(-a - b*x**3)) - 4*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(b**(sympy.S(1)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(-a - b*x**3)) - 2*sqrt(-a - b*x**3)/(b**(sympy.S(1)/3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_103():
    f = (x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1)/sqrt(a + b*x**3)
    F = -3**(sympy.S(1)/4)*sqrt((x**2*(b/a)**(sympy.S(2)/3) - x*(b/a)**(sympy.S(1)/3) + 1)/(x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*sqrt(2 - sqrt(3))*(x*(b/a)**(sympy.S(1)/3) + 1)*elliptic_e(asin((x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1)/(x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3))), -7 - 4*sqrt(3))/((b/a)**(sympy.S(1)/3)*sqrt((x*(b/a)**(sympy.S(1)/3) + 1)/(x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*sqrt(a + b*x**3)) + 2*(b/a)**(sympy.S(2)/3)*sqrt(a + b*x**3)/(b*(x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_104():
    f = (-x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1)/sqrt(a - b*x**3)
    F = 3**(sympy.S(1)/4)*sqrt((x**2*(b/a)**(sympy.S(2)/3) + x*(b/a)**(sympy.S(1)/3) + 1)/(-x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*sqrt(2 - sqrt(3))*(-x*(b/a)**(sympy.S(1)/3) + 1)*elliptic_e(asin((-x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1)/(-x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3))), -7 - 4*sqrt(3))/((b/a)**(sympy.S(1)/3)*sqrt((-x*(b/a)**(sympy.S(1)/3) + 1)/(-x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*sqrt(a - b*x**3)) - 2*(b/a)**(sympy.S(2)/3)*sqrt(a - b*x**3)/(b*(-x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_105():
    f = (-x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1)/sqrt(-a + b*x**3)
    F = -3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*(b/a)**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(-a + b*x**3)) + 2*(b/a)**(sympy.S(1)/3)*sqrt(-a + b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)) - 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*(b/a)**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*(1 - sqrt(3)))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(-a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_106():
    f = (x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1)/sqrt(-a - b*x**3)
    F = 3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*(b/a)**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(-a - b*x**3)) - 2*(b/a)**(sympy.S(1)/3)*sqrt(-a - b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*(b/a)**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*(1 - sqrt(3)))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(-a - b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_107():
    f = (c + d*x)/sqrt(a + b*x**3)
    F = -3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*d*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*d*sqrt(a + b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*d*(1 - sqrt(3)) + b**(sympy.S(1)/3)*c)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_108():
    f = (c + d*x)/sqrt(a - b*x**3)
    F = -3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*d*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(a - b*x**3)) + 2*d*sqrt(a - b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)) - 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*d*(1 - sqrt(3)) + b**(sympy.S(1)/3)*c)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(a - b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_109():
    f = (c + d*x)/sqrt(-a + b*x**3)
    F = 3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*d*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(-a + b*x**3)) - 2*d*sqrt(-a + b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)) - 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*d*(1 + sqrt(3)) + b**(sympy.S(1)/3)*c)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(-a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_110():
    f = (c + d*x)/sqrt(-a - b*x**3)
    F = 3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*d*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(-a - b*x**3)) - 2*d*sqrt(-a - b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*d*(1 + sqrt(3)) + b**(sympy.S(1)/3)*c)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(-a - b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_111():
    f = (c + d*x)/sqrt(x**3 + 1)
    F = 2*d*sqrt(x**3 + 1)/(x + 1 + sqrt(3)) - 3**(sympy.S(1)/4)*d*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(2 - sqrt(3))*(x + 1)*elliptic_e(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1)) + 2*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*(c - d*(1 - sqrt(3)))*(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_112():
    f = (c + d*x)/sqrt(1 - x**3)
    F = 2*d*sqrt(1 - x**3)/(-x + 1 + sqrt(3)) - 3**(sympy.S(1)/4)*d*sqrt((x**2 + x + 1)/(-x + 1 + sqrt(3))**2)*(1 - x)*sqrt(2 - sqrt(3))*elliptic_e(asin((-x - sqrt(3) + 1)/(-x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(sqrt((1 - x)/(-x + 1 + sqrt(3))**2)*sqrt(1 - x**3)) - 2*3**(sympy.S(3)/4)*sqrt((x**2 + x + 1)/(-x + 1 + sqrt(3))**2)*(1 - x)*sqrt(sqrt(3) + 2)*(c - sqrt(3)*d + d)*elliptic_f(asin((-x - sqrt(3) + 1)/(-x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((1 - x)/(-x + 1 + sqrt(3))**2)*sqrt(1 - x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_113():
    f = (c + d*x)/sqrt(x**3 - 1)
    F = 3**(sympy.S(1)/4)*d*sqrt((x**2 + x + 1)/(-x - sqrt(3) + 1)**2)*(1 - x)*sqrt(sqrt(3) + 2)*elliptic_e(asin((-x + 1 + sqrt(3))/(-x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(sqrt(-(1 - x)/(-x - sqrt(3) + 1)**2)*sqrt(x**3 - 1)) - 2*d*sqrt(x**3 - 1)/(-x - sqrt(3) + 1) - 2*3**(sympy.S(3)/4)*sqrt((x**2 + x + 1)/(-x - sqrt(3) + 1)**2)*(1 - x)*sqrt(2 - sqrt(3))*(c + d + sqrt(3)*d)*elliptic_f(asin((-x + 1 + sqrt(3))/(-x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(3*sqrt(-(1 - x)/(-x - sqrt(3) + 1)**2)*sqrt(x**3 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_114():
    f = (c + d*x)/sqrt(-x**3 - 1)
    F = 3**(sympy.S(1)/4)*d*sqrt((x**2 - x + 1)/(x - sqrt(3) + 1)**2)*sqrt(sqrt(3) + 2)*(x + 1)*elliptic_e(asin((x + 1 + sqrt(3))/(x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(sqrt(-(x + 1)/(x - sqrt(3) + 1)**2)*sqrt(-x**3 - 1)) - 2*d*sqrt(-x**3 - 1)/(x - sqrt(3) + 1) + 2*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x - sqrt(3) + 1)**2)*sqrt(2 - sqrt(3))*(c - d*(1 + sqrt(3)))*(x + 1)*elliptic_f(asin((x + 1 + sqrt(3))/(x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(3*sqrt(-(x + 1)/(x - sqrt(3) + 1)**2)*sqrt(-x**3 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_115():
    f = (c + d*x)/(a - b*x**4)
    F = d*atanh(sqrt(b)*x**2/sqrt(a))/(2*sqrt(a)*sqrt(b)) + c*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)) + c*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_116():
    f = (c + d*x)/(a + b*x**4)
    F = d*atan(sqrt(b)*x**2/sqrt(a))/(2*sqrt(a)*sqrt(b)) - sqrt(2)*c*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)) + sqrt(2)*c*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)) - sqrt(2)*c*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)) + sqrt(2)*c*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_117():
    f = (c + d*x)/(a - b*x**4)**2
    F = x*(c + d*x)/(4*a*(a - b*x**4)) + d*atanh(sqrt(b)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*sqrt(b)) + 3*c*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*b**(sympy.S(1)/4)) + 3*c*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*b**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_118():
    f = (c + d*x)/(a + b*x**4)**2
    F = x*(c + d*x)/(4*a*(a + b*x**4)) + d*atan(sqrt(b)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*sqrt(b)) - 3*sqrt(2)*c*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(1)/4)) + 3*sqrt(2)*c*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(1)/4)) - 3*sqrt(2)*c*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(1)/4)) + 3*sqrt(2)*c*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_119():
    f = (c + d*x)/(a - b*x**4)**3
    F = x*(c + d*x)/(8*a*(a - b*x**4)**2) + x*(7*c + 6*d*x)/(32*a**2*(a - b*x**4)) + 3*d*atanh(sqrt(b)*x**2/sqrt(a))/(16*a**(sympy.S(5)/2)*sqrt(b)) + 21*c*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(64*a**(sympy.S(11)/4)*b**(sympy.S(1)/4)) + 21*c*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(64*a**(sympy.S(11)/4)*b**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_120():
    f = (c + d*x)/(a + b*x**4)**3
    F = x*(c + d*x)/(8*a*(a + b*x**4)**2) + x*(7*c + 6*d*x)/(32*a**2*(a + b*x**4)) + 3*d*atan(sqrt(b)*x**2/sqrt(a))/(16*a**(sympy.S(5)/2)*sqrt(b)) - 21*sqrt(2)*c*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(256*a**(sympy.S(11)/4)*b**(sympy.S(1)/4)) + 21*sqrt(2)*c*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(256*a**(sympy.S(11)/4)*b**(sympy.S(1)/4)) - 21*sqrt(2)*c*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*b**(sympy.S(1)/4)) + 21*sqrt(2)*c*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*b**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_121():
    f = (c + d*x)/(a - b*x**4)**4
    F = x*(c + d*x)/(12*a*(a - b*x**4)**3) + x*(11*c + 10*d*x)/(96*a**2*(a - b*x**4)**2) + x*(77*c + 60*d*x)/(384*a**3*(a - b*x**4)) + 5*d*atanh(sqrt(b)*x**2/sqrt(a))/(32*a**(sympy.S(7)/2)*sqrt(b)) + 77*c*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(256*a**(sympy.S(15)/4)*b**(sympy.S(1)/4)) + 77*c*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(256*a**(sympy.S(15)/4)*b**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_122():
    f = (c + d*x)/(a + b*x**4)**4
    F = x*(c + d*x)/(12*a*(a + b*x**4)**3) + x*(11*c + 10*d*x)/(96*a**2*(a + b*x**4)**2) + x*(77*c + 60*d*x)/(384*a**3*(a + b*x**4)) + 5*d*atan(sqrt(b)*x**2/sqrt(a))/(32*a**(sympy.S(7)/2)*sqrt(b)) - 77*sqrt(2)*c*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(1024*a**(sympy.S(15)/4)*b**(sympy.S(1)/4)) + 77*sqrt(2)*c*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(1024*a**(sympy.S(15)/4)*b**(sympy.S(1)/4)) - 77*sqrt(2)*c*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(512*a**(sympy.S(15)/4)*b**(sympy.S(1)/4)) + 77*sqrt(2)*c*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(512*a**(sympy.S(15)/4)*b**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_123():
    f = (c + d*x)/(1 - x**4)
    F = c*atan(x)/2 + c*atanh(x)/2 + d*atanh(x**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_124():
    f = (c + d*x)/(x**4 + 1)
    F = -sqrt(2)*c*log(x**2 - sqrt(2)*x + 1)/8 + sqrt(2)*c*log(x**2 + sqrt(2)*x + 1)/8 + sqrt(2)*c*atan(sqrt(2)*x - 1)/4 + sqrt(2)*c*atan(sqrt(2)*x + 1)/4 + d*atan(x**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_125():
    f = (c + d*x + e*x**2)/(a - b*x**4)
    F = d*atanh(sqrt(b)*x**2/sqrt(a))/(2*sqrt(a)*sqrt(b)) + (-sqrt(a)*e + sqrt(b)*c)*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)) + (sqrt(a)*e + sqrt(b)*c)*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_126():
    f = (c + d*x + e*x**2)/(a + b*x**4)
    F = d*atan(sqrt(b)*x**2/sqrt(a))/(2*sqrt(a)*sqrt(b)) - sqrt(2)*(-sqrt(a)*e + sqrt(b)*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(-sqrt(a)*e + sqrt(b)*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)) - sqrt(2)*(sqrt(a)*e + sqrt(b)*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(sqrt(a)*e + sqrt(b)*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_127():
    f = (c + d*x + e*x**2)/(a - b*x**4)**2
    F = x*(c + d*x + e*x**2)/(4*a*(a - b*x**4)) + d*atanh(sqrt(b)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*sqrt(b)) + (-sqrt(a)*e + 3*sqrt(b)*c)*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*b**(sympy.S(3)/4)) + (sqrt(a)*e + 3*sqrt(b)*c)*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_128():
    f = (c + d*x + e*x**2)/(a + b*x**4)**2
    F = x*(c + d*x + e*x**2)/(4*a*(a + b*x**4)) + d*atan(sqrt(b)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*sqrt(b)) - sqrt(2)*(-sqrt(a)*e + 3*sqrt(b)*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(-sqrt(a)*e + 3*sqrt(b)*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(3)/4)) - sqrt(2)*(sqrt(a)*e + 3*sqrt(b)*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(sqrt(a)*e + 3*sqrt(b)*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_129():
    f = (c + d*x + e*x**2)/(a - b*x**4)**3
    F = x*(c + d*x + e*x**2)/(8*a*(a - b*x**4)**2) + x*(7*c + 6*d*x + 5*e*x**2)/(32*a**2*(a - b*x**4)) + 3*d*atanh(sqrt(b)*x**2/sqrt(a))/(16*a**(sympy.S(5)/2)*sqrt(b)) + (-5*sqrt(a)*e + 21*sqrt(b)*c)*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(64*a**(sympy.S(11)/4)*b**(sympy.S(3)/4)) + (5*sqrt(a)*e + 21*sqrt(b)*c)*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(64*a**(sympy.S(11)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_130():
    f = (c + d*x + e*x**2)/(a + b*x**4)**3
    F = x*(c + d*x + e*x**2)/(8*a*(a + b*x**4)**2) + x*(7*c + 6*d*x + 5*e*x**2)/(32*a**2*(a + b*x**4)) + 3*d*atan(sqrt(b)*x**2/sqrt(a))/(16*a**(sympy.S(5)/2)*sqrt(b)) - sqrt(2)*(-5*sqrt(a)*e + 21*sqrt(b)*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(256*a**(sympy.S(11)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(-5*sqrt(a)*e + 21*sqrt(b)*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(256*a**(sympy.S(11)/4)*b**(sympy.S(3)/4)) - sqrt(2)*(5*sqrt(a)*e + 21*sqrt(b)*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(5*sqrt(a)*e + 21*sqrt(b)*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_131():
    f = (c + d*x + e*x**2)/(a - b*x**4)**4
    F = x*(c + d*x + e*x**2)/(12*a*(a - b*x**4)**3) + x*(11*c + 10*d*x + 9*e*x**2)/(96*a**2*(a - b*x**4)**2) + x*(77*c + 60*d*x + 45*e*x**2)/(384*a**3*(a - b*x**4)) + 5*d*atanh(sqrt(b)*x**2/sqrt(a))/(32*a**(sympy.S(7)/2)*sqrt(b)) + (-15*sqrt(a)*e + 77*sqrt(b)*c)*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(256*a**(sympy.S(15)/4)*b**(sympy.S(3)/4)) + (15*sqrt(a)*e + 77*sqrt(b)*c)*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(256*a**(sympy.S(15)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_132():
    f = (c + d*x + e*x**2)/(a + b*x**4)**4
    F = x*(c + d*x + e*x**2)/(12*a*(a + b*x**4)**3) + x*(11*c + 10*d*x + 9*e*x**2)/(96*a**2*(a + b*x**4)**2) + x*(77*c + 60*d*x + 45*e*x**2)/(384*a**3*(a + b*x**4)) + 5*d*atan(sqrt(b)*x**2/sqrt(a))/(32*a**(sympy.S(7)/2)*sqrt(b)) - sqrt(2)*(-15*sqrt(a)*e + 77*sqrt(b)*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(1024*a**(sympy.S(15)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(-15*sqrt(a)*e + 77*sqrt(b)*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(1024*a**(sympy.S(15)/4)*b**(sympy.S(3)/4)) - sqrt(2)*(15*sqrt(a)*e + 77*sqrt(b)*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(512*a**(sympy.S(15)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(15*sqrt(a)*e + 77*sqrt(b)*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(512*a**(sympy.S(15)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_133():
    f = a*(e + f*x**4)**2
    F = a*e**2*x + 2*a*e*f*x**5/5 + a*f**2*x**9/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_134():
    f = b*x*(e + f*x**4)**2
    F = b*e**2*x**2/2 + b*e*f*x**6/3 + b*f**2*x**10/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_135():
    f = (a + b*x)*(e + f*x**4)**2
    F = a*e**2*x + 2*a*e*f*x**5/5 + a*f**2*x**9/9 + b*e**2*x**2/2 + b*e*f*x**6/3 + b*f**2*x**10/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_136():
    f = c*x**2*(e + f*x**4)**2
    F = c*e**2*x**3/3 + 2*c*e*f*x**7/7 + c*f**2*x**11/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_137():
    f = (a + c*x**2)*(e + f*x**4)**2
    F = a*e**2*x + 2*a*e*f*x**5/5 + a*f**2*x**9/9 + c*e**2*x**3/3 + 2*c*e*f*x**7/7 + c*f**2*x**11/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_138():
    f = (e + f*x**4)**2*(b*x + c*x**2)
    F = b*e**2*x**2/2 + b*e*f*x**6/3 + b*f**2*x**10/10 + c*e**2*x**3/3 + 2*c*e*f*x**7/7 + c*f**2*x**11/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_139():
    f = (e + f*x**4)**2*(a + b*x + c*x**2)
    F = a*e**2*x + 2*a*e*f*x**5/5 + a*f**2*x**9/9 + b*e**2*x**2/2 + b*e*f*x**6/3 + b*f**2*x**10/10 + c*e**2*x**3/3 + 2*c*e*f*x**7/7 + c*f**2*x**11/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_140():
    f = d*x**3*(e + f*x**4)**2
    F = d*(e + f*x**4)**3/(12*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_141():
    f = (a + d*x**3)*(e + f*x**4)**2
    F = a*e**2*x + 2*a*e*f*x**5/5 + a*f**2*x**9/9 + d*(e + f*x**4)**3/(12*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_142():
    f = (e + f*x**4)**2*(b*x + d*x**3)
    F = b*e**2*x**2/2 + b*e*f*x**6/3 + b*f**2*x**10/10 + d*(e + f*x**4)**3/(12*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_143():
    f = (e + f*x**4)**2*(a + b*x + d*x**3)
    F = a*e**2*x + 2*a*e*f*x**5/5 + a*f**2*x**9/9 + b*e**2*x**2/2 + b*e*f*x**6/3 + b*f**2*x**10/10 + d*(e + f*x**4)**3/(12*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_144():
    f = (e + f*x**4)**2*(c*x**2 + d*x**3)
    F = c*e**2*x**3/3 + 2*c*e*f*x**7/7 + c*f**2*x**11/11 + d*(e + f*x**4)**3/(12*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_145():
    f = (e + f*x**4)**2*(a + c*x**2 + d*x**3)
    F = a*e**2*x + 2*a*e*f*x**5/5 + a*f**2*x**9/9 + c*e**2*x**3/3 + 2*c*e*f*x**7/7 + c*f**2*x**11/11 + d*(e + f*x**4)**3/(12*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_146():
    f = (e + f*x**4)**2*(b*x + c*x**2 + d*x**3)
    F = b*e**2*x**2/2 + b*e*f*x**6/3 + b*f**2*x**10/10 + c*e**2*x**3/3 + 2*c*e*f*x**7/7 + c*f**2*x**11/11 + d*(e + f*x**4)**3/(12*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_147():
    f = (a + b*x**4)**2*(c + d*x + e*x**2 + f*x**3)
    F = a**2*c*x + a**2*d*x**2/2 + a**2*e*x**3/3 + 2*a*b*c*x**5/5 + a*b*d*x**6/3 + 2*a*b*e*x**7/7 + b**2*c*x**9/9 + b**2*d*x**10/10 + b**2*e*x**11/11 + f*(a + b*x**4)**3/(12*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_148():
    f = (a + b*x**4)**3*(c + d*x + e*x**2 + f*x**3)
    F = a**3*c*x + a**3*d*x**2/2 + a**3*e*x**3/3 + 3*a**2*b*c*x**5/5 + a**2*b*d*x**6/2 + 3*a**2*b*e*x**7/7 + a*b**2*c*x**9/3 + 3*a*b**2*d*x**10/10 + 3*a*b**2*e*x**11/11 + b**3*c*x**13/13 + b**3*d*x**14/14 + b**3*e*x**15/15 + f*(a + b*x**4)**4/(16*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_149():
    f = (c + d*x + e*x**2 + f*x**3)/(a - b*x**4)**2
    F = (a*f + b*x*(c + d*x + e*x**2))/(4*a*b*(a - b*x**4)) + d*atanh(sqrt(b)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*sqrt(b)) + (-sqrt(a)*e + 3*sqrt(b)*c)*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*b**(sympy.S(3)/4)) + (sqrt(a)*e + 3*sqrt(b)*c)*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_150():
    f = (c + d*x + e*x**2 + f*x**3)/(a - b*x**4)**3
    F = (a*f + b*x*(c + d*x + e*x**2))/(8*a*b*(a - b*x**4)**2) + x*(7*c + 6*d*x + 5*e*x**2)/(32*a**2*(a - b*x**4)) + 3*d*atanh(sqrt(b)*x**2/sqrt(a))/(16*a**(sympy.S(5)/2)*sqrt(b)) + (-5*sqrt(a)*e + 21*sqrt(b)*c)*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(64*a**(sympy.S(11)/4)*b**(sympy.S(3)/4)) + (5*sqrt(a)*e + 21*sqrt(b)*c)*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(64*a**(sympy.S(11)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_151():
    f = (c + d*x + e*x**2 + f*x**3)/(a - b*x**4)**4
    F = (a*f + b*x*(c + d*x + e*x**2))/(12*a*b*(a - b*x**4)**3) + x*(11*c + 10*d*x + 9*e*x**2)/(96*a**2*(a - b*x**4)**2) + x*(77*c + 60*d*x + 45*e*x**2)/(384*a**3*(a - b*x**4)) + 5*d*atanh(sqrt(b)*x**2/sqrt(a))/(32*a**(sympy.S(7)/2)*sqrt(b)) + (-15*sqrt(a)*e + 77*sqrt(b)*c)*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(256*a**(sympy.S(15)/4)*b**(sympy.S(3)/4)) + (15*sqrt(a)*e + 77*sqrt(b)*c)*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(256*a**(sympy.S(15)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_152():
    f = a/(3*x**4 + 2)
    F = -6**(sympy.S(3)/4)*a*log(3*x**2 - 6**(sympy.S(3)/4)*x + sqrt(6))/48 + 6**(sympy.S(3)/4)*a*log(3*x**2 + 6**(sympy.S(3)/4)*x + sqrt(6))/48 + 6**(sympy.S(3)/4)*a*atan(6**(sympy.S(1)/4)*x - 1)/24 + 6**(sympy.S(3)/4)*a*atan(6**(sympy.S(1)/4)*x + 1)/24
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_153():
    f = b*x/(3*x**4 + 2)
    F = sqrt(6)*b*atan(sqrt(6)*x**2/2)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_154():
    f = (a + b*x)/(3*x**4 + 2)
    F = -6**(sympy.S(3)/4)*a*log(3*x**2 - 6**(sympy.S(3)/4)*x + sqrt(6))/48 + 6**(sympy.S(3)/4)*a*log(3*x**2 + 6**(sympy.S(3)/4)*x + sqrt(6))/48 + 6**(sympy.S(3)/4)*a*atan(6**(sympy.S(1)/4)*x - 1)/24 + 6**(sympy.S(3)/4)*a*atan(6**(sympy.S(1)/4)*x + 1)/24 + sqrt(6)*b*atan(sqrt(6)*x**2/2)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_155():
    f = c*x**2/(3*x**4 + 2)
    F = 6**(sympy.S(1)/4)*c*log(3*x**2 - 6**(sympy.S(3)/4)*x + sqrt(6))/24 - 6**(sympy.S(1)/4)*c*log(3*x**2 + 6**(sympy.S(3)/4)*x + sqrt(6))/24 + 6**(sympy.S(1)/4)*c*atan(6**(sympy.S(1)/4)*x - 1)/12 + 6**(sympy.S(1)/4)*c*atan(6**(sympy.S(1)/4)*x + 1)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_156():
    f = (a + c*x**2)/(3*x**4 + 2)
    F = -6**(sympy.S(1)/4)*(sqrt(6)*a - 2*c)*log(3*x**2 - 6**(sympy.S(3)/4)*x + sqrt(6))/48 + 6**(sympy.S(1)/4)*(sqrt(6)*a - 2*c)*log(3*x**2 + 6**(sympy.S(3)/4)*x + sqrt(6))/48 + 6**(sympy.S(1)/4)*(sqrt(6)*a + 2*c)*atan(6**(sympy.S(1)/4)*x - 1)/24 + 6**(sympy.S(1)/4)*(sqrt(6)*a + 2*c)*atan(6**(sympy.S(1)/4)*x + 1)/24
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_157():
    f = (b*x + c*x**2)/(3*x**4 + 2)
    F = sqrt(6)*b*atan(sqrt(6)*x**2/2)/12 + 6**(sympy.S(1)/4)*c*log(3*x**2 - 6**(sympy.S(3)/4)*x + sqrt(6))/24 - 6**(sympy.S(1)/4)*c*log(3*x**2 + 6**(sympy.S(3)/4)*x + sqrt(6))/24 + 6**(sympy.S(1)/4)*c*atan(6**(sympy.S(1)/4)*x - 1)/12 + 6**(sympy.S(1)/4)*c*atan(6**(sympy.S(1)/4)*x + 1)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_158():
    f = (a + b*x + c*x**2)/(3*x**4 + 2)
    F = sqrt(6)*b*atan(sqrt(6)*x**2/2)/12 - 6**(sympy.S(1)/4)*(sqrt(6)*a - 2*c)*log(3*x**2 - 6**(sympy.S(3)/4)*x + sqrt(6))/48 + 6**(sympy.S(1)/4)*(sqrt(6)*a - 2*c)*log(3*x**2 + 6**(sympy.S(3)/4)*x + sqrt(6))/48 + 6**(sympy.S(1)/4)*(sqrt(6)*a + 2*c)*atan(6**(sympy.S(1)/4)*x - 1)/24 + 6**(sympy.S(1)/4)*(sqrt(6)*a + 2*c)*atan(6**(sympy.S(1)/4)*x + 1)/24
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_159():
    f = d*x**3/(3*x**4 + 2)
    F = d*log(3*x**4 + 2)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_160():
    f = (a + d*x**3)/(3*x**4 + 2)
    F = -6**(sympy.S(3)/4)*a*log(3*x**2 - 6**(sympy.S(3)/4)*x + sqrt(6))/48 + 6**(sympy.S(3)/4)*a*log(3*x**2 + 6**(sympy.S(3)/4)*x + sqrt(6))/48 + 6**(sympy.S(3)/4)*a*atan(6**(sympy.S(1)/4)*x - 1)/24 + 6**(sympy.S(3)/4)*a*atan(6**(sympy.S(1)/4)*x + 1)/24 + d*log(3*x**4 + 2)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_161():
    f = (b*x + d*x**3)/(3*x**4 + 2)
    F = sqrt(6)*b*atan(sqrt(6)*x**2/2)/12 + d*log(3*x**4 + 2)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_162():
    f = (a + b*x + d*x**3)/(3*x**4 + 2)
    F = -6**(sympy.S(3)/4)*a*log(3*x**2 - 6**(sympy.S(3)/4)*x + sqrt(6))/48 + 6**(sympy.S(3)/4)*a*log(3*x**2 + 6**(sympy.S(3)/4)*x + sqrt(6))/48 + 6**(sympy.S(3)/4)*a*atan(6**(sympy.S(1)/4)*x - 1)/24 + 6**(sympy.S(3)/4)*a*atan(6**(sympy.S(1)/4)*x + 1)/24 + sqrt(6)*b*atan(sqrt(6)*x**2/2)/12 + d*log(3*x**4 + 2)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_163():
    f = (c*x**2 + d*x**3)/(3*x**4 + 2)
    F = 6**(sympy.S(1)/4)*c*log(3*x**2 - 6**(sympy.S(3)/4)*x + sqrt(6))/24 - 6**(sympy.S(1)/4)*c*log(3*x**2 + 6**(sympy.S(3)/4)*x + sqrt(6))/24 + 6**(sympy.S(1)/4)*c*atan(6**(sympy.S(1)/4)*x - 1)/12 + 6**(sympy.S(1)/4)*c*atan(6**(sympy.S(1)/4)*x + 1)/12 + d*log(3*x**4 + 2)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_164():
    f = (a + c*x**2 + d*x**3)/(3*x**4 + 2)
    F = d*log(3*x**4 + 2)/12 - 6**(sympy.S(1)/4)*(sqrt(6)*a - 2*c)*log(3*x**2 - 6**(sympy.S(3)/4)*x + sqrt(6))/48 + 6**(sympy.S(1)/4)*(sqrt(6)*a - 2*c)*log(3*x**2 + 6**(sympy.S(3)/4)*x + sqrt(6))/48 + 6**(sympy.S(1)/4)*(sqrt(6)*a + 2*c)*atan(6**(sympy.S(1)/4)*x - 1)/24 + 6**(sympy.S(1)/4)*(sqrt(6)*a + 2*c)*atan(6**(sympy.S(1)/4)*x + 1)/24
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_165():
    f = (b*x + c*x**2 + d*x**3)/(3*x**4 + 2)
    F = sqrt(6)*b*atan(sqrt(6)*x**2/2)/12 + 6**(sympy.S(1)/4)*c*log(3*x**2 - 6**(sympy.S(3)/4)*x + sqrt(6))/24 - 6**(sympy.S(1)/4)*c*log(3*x**2 + 6**(sympy.S(3)/4)*x + sqrt(6))/24 + 6**(sympy.S(1)/4)*c*atan(6**(sympy.S(1)/4)*x - 1)/12 + 6**(sympy.S(1)/4)*c*atan(6**(sympy.S(1)/4)*x + 1)/12 + d*log(3*x**4 + 2)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_166():
    f = (a + b*x + c*x**2 + d*x**3)/(3*x**4 + 2)
    F = sqrt(6)*b*atan(sqrt(6)*x**2/2)/12 + d*log(3*x**4 + 2)/12 - 6**(sympy.S(1)/4)*(sqrt(6)*a - 2*c)*log(3*x**2 - 6**(sympy.S(3)/4)*x + sqrt(6))/48 + 6**(sympy.S(1)/4)*(sqrt(6)*a - 2*c)*log(3*x**2 + 6**(sympy.S(3)/4)*x + sqrt(6))/48 + 6**(sympy.S(1)/4)*(sqrt(6)*a + 2*c)*atan(6**(sympy.S(1)/4)*x - 1)/24 + 6**(sympy.S(1)/4)*(sqrt(6)*a + 2*c)*atan(6**(sympy.S(1)/4)*x + 1)/24
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_167():
    f = (x**3 + x**2 + x + 1)/(1 - x**4)
    F = -log(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_168():
    f = (x**3 + x**2 + x + 1)/(x**4 + 1)
    F = log(x**4 + 1)/4 + atan(x**2)/2 + sqrt(2)*atan(sqrt(2)*x - 1)/2 + sqrt(2)*atan(sqrt(2)*x + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_169():
    f = (x**3 + x**2 + x + 1)/(a - b*x**4)
    F = -log(a - b*x**4)/(4*b) + atanh(sqrt(b)*x**2/sqrt(a))/(2*sqrt(a)*sqrt(b)) - (sqrt(a) - sqrt(b))*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)) + (sqrt(a) + sqrt(b))*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_170():
    f = (x**3 + x**2 + x + 1)/(a + b*x**4)
    F = log(a + b*x**4)/(4*b) + atan(sqrt(b)*x**2/sqrt(a))/(2*sqrt(a)*sqrt(b)) + sqrt(2)*(sqrt(a) - sqrt(b))*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)) - sqrt(2)*(sqrt(a) - sqrt(b))*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)) - sqrt(2)*(sqrt(a) + sqrt(b))*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(sqrt(a) + sqrt(b))*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_171():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4)/(a - b*x**4)
    F = -f*log(a - b*x**4)/(4*b) - g*x/b + d*atanh(sqrt(b)*x**2/sqrt(a))/(2*sqrt(a)*sqrt(b)) + (-sqrt(a)*sqrt(b)*e + a*g + b*c)*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)) + (sqrt(a)*sqrt(b)*e + a*g + b*c)*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_172():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4)/(a - b*x**4)**2
    F = x*(a*g + b*c + b*d*x + b*e*x**2 + b*f*x**3)/(4*a*b*(a - b*x**4)) + d*atanh(sqrt(b)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*sqrt(b)) + (-sqrt(a)*sqrt(b)*e - a*g + 3*b*c)*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*b**(sympy.S(5)/4)) + (sqrt(a)*sqrt(b)*e - a*g + 3*b*c)*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_173():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4)/(a - b*x**4)**3
    F = x*(a*g + b*c + b*d*x + b*e*x**2 + b*f*x**3)/(8*a*b*(a - b*x**4)**2) + (4*a*f + x*(-a*g + 7*b*c + 6*b*d*x + 5*b*e*x**2))/(32*a**2*b*(a - b*x**4)) + 3*d*atanh(sqrt(b)*x**2/sqrt(a))/(16*a**(sympy.S(5)/2)*sqrt(b)) + (-5*sqrt(a)*sqrt(b)*e - 3*a*g + 21*b*c)*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(64*a**(sympy.S(11)/4)*b**(sympy.S(5)/4)) + (5*sqrt(a)*sqrt(b)*e - 3*a*g + 21*b*c)*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(64*a**(sympy.S(11)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_174():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4)/(a - b*x**4)**4
    F = x*(a*g + b*c + b*d*x + b*e*x**2 + b*f*x**3)/(12*a*b*(a - b*x**4)**3) + (8*a*f + x*(-a*g + 11*b*c + 10*b*d*x + 9*b*e*x**2))/(96*a**2*b*(a - b*x**4)**2) + x*(-7*a*g + 77*b*c + 60*b*d*x + 45*b*e*x**2)/(384*a**3*b*(a - b*x**4)) + 5*d*atanh(sqrt(b)*x**2/sqrt(a))/(32*a**(sympy.S(7)/2)*sqrt(b)) + (-15*sqrt(a)*sqrt(b)*e - 7*a*g + 77*b*c)*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(256*a**(sympy.S(15)/4)*b**(sympy.S(5)/4)) + (15*sqrt(a)*sqrt(b)*e - 7*a*g + 77*b*c)*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(256*a**(sympy.S(15)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_175():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4)/(a + b*x**4)
    F = f*log(a + b*x**4)/(4*b) + g*x/b + d*atan(sqrt(b)*x**2/sqrt(a))/(2*sqrt(a)*sqrt(b)) - sqrt(2)*(-sqrt(a)*sqrt(b)*e - a*g + b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(-sqrt(a)*sqrt(b)*e - a*g + b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)) - sqrt(2)*(sqrt(a)*sqrt(b)*e - a*g + b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(sqrt(a)*sqrt(b)*e - a*g + b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_176():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4)/(a + b*x**4)**2
    F = x*(-a*g + b*c + b*d*x + b*e*x**2 + b*f*x**3)/(4*a*b*(a + b*x**4)) + d*atan(sqrt(b)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*sqrt(b)) - sqrt(2)*(-sqrt(a)*sqrt(b)*e + a*g + 3*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(-sqrt(a)*sqrt(b)*e + a*g + 3*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(5)/4)) - sqrt(2)*(sqrt(a)*sqrt(b)*e + a*g + 3*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(sqrt(a)*sqrt(b)*e + a*g + 3*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_177():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4)/(a + b*x**4)**3
    F = x*(-a*g + b*c + b*d*x + b*e*x**2 + b*f*x**3)/(8*a*b*(a + b*x**4)**2) - (4*a*f - x*(a*g + 7*b*c + 6*b*d*x + 5*b*e*x**2))/(32*a**2*b*(a + b*x**4)) + 3*d*atan(sqrt(b)*x**2/sqrt(a))/(16*a**(sympy.S(5)/2)*sqrt(b)) - sqrt(2)*(-5*sqrt(a)*sqrt(b)*e + 3*a*g + 21*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(256*a**(sympy.S(11)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(-5*sqrt(a)*sqrt(b)*e + 3*a*g + 21*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(256*a**(sympy.S(11)/4)*b**(sympy.S(5)/4)) - sqrt(2)*(5*sqrt(a)*sqrt(b)*e + 3*a*g + 21*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(5*sqrt(a)*sqrt(b)*e + 3*a*g + 21*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_178():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4)/(a + b*x**4)**4
    F = x*(-a*g + b*c + b*d*x + b*e*x**2 + b*f*x**3)/(12*a*b*(a + b*x**4)**3) - (8*a*f - x*(a*g + 11*b*c + 10*b*d*x + 9*b*e*x**2))/(96*a**2*b*(a + b*x**4)**2) + x*(7*a*g + 77*b*c + 60*b*d*x + 45*b*e*x**2)/(384*a**3*b*(a + b*x**4)) + 5*d*atan(sqrt(b)*x**2/sqrt(a))/(32*a**(sympy.S(7)/2)*sqrt(b)) - sqrt(2)*(-15*sqrt(a)*sqrt(b)*e + 7*a*g + 77*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(1024*a**(sympy.S(15)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(-15*sqrt(a)*sqrt(b)*e + 7*a*g + 77*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(1024*a**(sympy.S(15)/4)*b**(sympy.S(5)/4)) - sqrt(2)*(15*sqrt(a)*sqrt(b)*e + 7*a*g + 77*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(512*a**(sympy.S(15)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(15*sqrt(a)*sqrt(b)*e + 7*a*g + 77*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(512*a**(sympy.S(15)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_179():
    f = (1 - x**4)**3/(x**3 + x**2 + x + 1)**3
    F = -(1 - x)**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_180():
    f = (1 - x**4)**2/(x**3 + x**2 + x + 1)**2
    F = -(1 - x)**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_181():
    f = (1 - x**4)/(x**3 + x**2 + x + 1)
    F = -x**2/2 + x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_182():
    f = (x**3 + x**2 + x + 1)/(1 - x**4)
    F = -log(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_183():
    f = (x**3 + x**2 + x + 1)**2/(1 - x**4)**2
    F = 1/(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_184():
    f = (x**3 + x**2 + x + 1)**3/(1 - x**4)**3
    F = 1/(2*(1 - x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_185():
    f = (x**3 + x**2 + x + 1)**4/(1 - x**4)**4
    F = 1/(3*(1 - x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_186():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a - b*x**4)
    F = -f*log(a - b*x**4)/(4*b) - g*x/b - h*x**2/(2*b) + (a*h + b*d)*atanh(sqrt(b)*x**2/sqrt(a))/(2*sqrt(a)*b**(sympy.S(3)/2)) + (-sqrt(a)*sqrt(b)*e + a*g + b*c)*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)) + (sqrt(a)*sqrt(b)*e + a*g + b*c)*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_187():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5 + i*x**6)/(a - b*x**4)
    F = -f*log(a - b*x**4)/(4*b) - g*x/b - h*x**2/(2*b) - i*x**3/(3*b) + (a*h + b*d)*atanh(sqrt(b)*x**2/sqrt(a))/(2*sqrt(a)*b**(sympy.S(3)/2)) - (a*i + b*e - sqrt(b)*(a*g + b*c)/sqrt(a))*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)*b**(sympy.S(7)/4)) + (a*i + b*e + sqrt(b)*(a*g + b*c)/sqrt(a))*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_188():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5 + i*x**6 + j*x**7)/(a - b*x**4)
    F = -g*x/b - h*x**2/(2*b) - i*x**3/(3*b) - j*x**4/(4*b) - (a*j + b*f)*log(a - b*x**4)/(4*b**2) + (a*h + b*d)*atanh(sqrt(b)*x**2/sqrt(a))/(2*sqrt(a)*b**(sympy.S(3)/2)) - (a*i + b*e - sqrt(b)*(a*g + b*c)/sqrt(a))*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)*b**(sympy.S(7)/4)) + (a*i + b*e + sqrt(b)*(a*g + b*c)/sqrt(a))*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_189():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a + b*x**4)
    F = f*log(a + b*x**4)/(4*b) + g*x/b + h*x**2/(2*b) + (-a*h + b*d)*atan(sqrt(b)*x**2/sqrt(a))/(2*sqrt(a)*b**(sympy.S(3)/2)) - sqrt(2)*(-sqrt(a)*sqrt(b)*e - a*g + b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(-sqrt(a)*sqrt(b)*e - a*g + b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)) - sqrt(2)*(sqrt(a)*sqrt(b)*e - a*g + b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(sqrt(a)*sqrt(b)*e - a*g + b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_190():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5 + i*x**6)/(a + b*x**4)
    F = f*log(a + b*x**4)/(4*b) + g*x/b + h*x**2/(2*b) + i*x**3/(3*b) + (-a*h + b*d)*atan(sqrt(b)*x**2/sqrt(a))/(2*sqrt(a)*b**(sympy.S(3)/2)) - sqrt(2)*(-sqrt(a)*(-a*i + b*e) + sqrt(b)*(-a*g + b*c))*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(-sqrt(a)*(-a*i + b*e) + sqrt(b)*(-a*g + b*c))*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(7)/4)) - sqrt(2)*(sqrt(a)*(-a*i + b*e) + sqrt(b)*(-a*g + b*c))*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(sqrt(a)*(-a*i + b*e) + sqrt(b)*(-a*g + b*c))*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_191():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5 + i*x**6 + j*x**7)/(a + b*x**4)
    F = g*x/b + h*x**2/(2*b) + i*x**3/(3*b) + j*x**4/(4*b) + (-a*j + b*f)*log(a + b*x**4)/(4*b**2) + (-a*h + b*d)*atan(sqrt(b)*x**2/sqrt(a))/(2*sqrt(a)*b**(sympy.S(3)/2)) - sqrt(2)*(-sqrt(a)*(-a*i + b*e) + sqrt(b)*(-a*g + b*c))*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(-sqrt(a)*(-a*i + b*e) + sqrt(b)*(-a*g + b*c))*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(7)/4)) - sqrt(2)*(sqrt(a)*(-a*i + b*e) + sqrt(b)*(-a*g + b*c))*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(sqrt(a)*(-a*i + b*e) + sqrt(b)*(-a*g + b*c))*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_192():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a - b*x**4)**2
    F = x*(a*g + b*c + b*e*x**2 + b*f*x**3 + x*(a*h + b*d))/(4*a*b*(a - b*x**4)) + (-a*h + b*d)*atanh(sqrt(b)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*b**(sympy.S(3)/2)) + (-sqrt(a)*sqrt(b)*e - a*g + 3*b*c)*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*b**(sympy.S(5)/4)) + (sqrt(a)*sqrt(b)*e - a*g + 3*b*c)*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(8*a**(sympy.S(7)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_193():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5 + i*x**6)/(a - b*x**4)**2
    F = x*(a*g + b*c + b*f*x**3 + x**2*(a*i + b*e) + x*(a*h + b*d))/(4*a*b*(a - b*x**4)) + (-a*h + b*d)*atanh(sqrt(b)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*b**(sympy.S(3)/2)) - (-3*a*i + b*e - sqrt(b)*(-a*g + 3*b*c)/sqrt(a))*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(8*a**(sympy.S(5)/4)*b**(sympy.S(7)/4)) + (-3*a*i + b*e + sqrt(b)*(-a*g + 3*b*c)/sqrt(a))*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(8*a**(sympy.S(5)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_194():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5 + i*x**6 + j*x**7)/(a - b*x**4)**2
    F = j*log(a - b*x**4)/(4*b**2) + x*(a*g + b*c + x**3*(a*j + b*f) + x**2*(a*i + b*e) + x*(a*h + b*d))/(4*a*b*(a - b*x**4)) + (-a*h + b*d)*atanh(sqrt(b)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*b**(sympy.S(3)/2)) - (-3*a*i + b*e - sqrt(b)*(-a*g + 3*b*c)/sqrt(a))*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(8*a**(sympy.S(5)/4)*b**(sympy.S(7)/4)) + (-3*a*i + b*e + sqrt(b)*(-a*g + 3*b*c)/sqrt(a))*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(8*a**(sympy.S(5)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_195():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a + b*x**4)**2
    F = x*(-a*g + b*c + b*e*x**2 + b*f*x**3 + x*(-a*h + b*d))/(4*a*b*(a + b*x**4)) + (a*h + b*d)*atan(sqrt(b)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*b**(sympy.S(3)/2)) - sqrt(2)*(-sqrt(a)*sqrt(b)*e + a*g + 3*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(-sqrt(a)*sqrt(b)*e + a*g + 3*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(5)/4)) - sqrt(2)*(sqrt(a)*sqrt(b)*e + a*g + 3*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(sqrt(a)*sqrt(b)*e + a*g + 3*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_196():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5 + i*x**6)/(a + b*x**4)**2
    F = x*(-a*g + b*c + b*f*x**3 + x**2*(-a*i + b*e) + x*(-a*h + b*d))/(4*a*b*(a + b*x**4)) + (a*h + b*d)*atan(sqrt(b)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*b**(sympy.S(3)/2)) - sqrt(2)*(-sqrt(a)*(3*a*i + b*e) + sqrt(b)*(a*g + 3*b*c))*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(-sqrt(a)*(3*a*i + b*e) + sqrt(b)*(a*g + 3*b*c))*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(7)/4)) - sqrt(2)*(sqrt(a)*(3*a*i + b*e) + sqrt(b)*(a*g + 3*b*c))*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(sqrt(a)*(3*a*i + b*e) + sqrt(b)*(a*g + 3*b*c))*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_197():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5 + i*x**6 + j*x**7)/(a + b*x**4)**2
    F = j*log(a + b*x**4)/(4*b**2) + x*(-a*g + b*c + x**3*(-a*j + b*f) + x**2*(-a*i + b*e) + x*(-a*h + b*d))/(4*a*b*(a + b*x**4)) + (a*h + b*d)*atan(sqrt(b)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*b**(sympy.S(3)/2)) - sqrt(2)*(-sqrt(a)*(3*a*i + b*e) + sqrt(b)*(a*g + 3*b*c))*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(-sqrt(a)*(3*a*i + b*e) + sqrt(b)*(a*g + 3*b*c))*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(7)/4)) - sqrt(2)*(sqrt(a)*(3*a*i + b*e) + sqrt(b)*(a*g + 3*b*c))*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(sqrt(a)*(3*a*i + b*e) + sqrt(b)*(a*g + 3*b*c))*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_198():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a - b*x**4)**3
    F = x*(a*g + b*c + b*e*x**2 + b*f*x**3 + x*(a*h + b*d))/(8*a*b*(a - b*x**4)**2) + (4*a*f + x*(-a*g + 7*b*c + 5*b*e*x**2 + x*(-2*a*h + 6*b*d)))/(32*a**2*b*(a - b*x**4)) + (-a*h + 3*b*d)*atanh(sqrt(b)*x**2/sqrt(a))/(16*a**(sympy.S(5)/2)*b**(sympy.S(3)/2)) + (-5*sqrt(a)*sqrt(b)*e - 3*a*g + 21*b*c)*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(64*a**(sympy.S(11)/4)*b**(sympy.S(5)/4)) + (5*sqrt(a)*sqrt(b)*e - 3*a*g + 21*b*c)*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(64*a**(sympy.S(11)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_199():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5 + i*x**6)/(a - b*x**4)**3
    F = x*(a*g + b*c + b*f*x**3 + x**2*(a*i + b*e) + x*(a*h + b*d))/(8*a*b*(a - b*x**4)**2) + (4*a*f + x*(-a*g + 7*b*c + x**2*(-3*a*i + 5*b*e) + x*(-2*a*h + 6*b*d)))/(32*a**2*b*(a - b*x**4)) + (-a*h + 3*b*d)*atanh(sqrt(b)*x**2/sqrt(a))/(16*a**(sympy.S(5)/2)*b**(sympy.S(3)/2)) - (-3*a*i + 5*b*e - 3*sqrt(b)*(-a*g + 7*b*c)/sqrt(a))*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(64*a**(sympy.S(9)/4)*b**(sympy.S(7)/4)) + (-3*a*i + 5*b*e + 3*sqrt(b)*(-a*g + 7*b*c)/sqrt(a))*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(64*a**(sympy.S(9)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_200():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5 + i*x**6 + j*x**7)/(a - b*x**4)**3
    F = x*(a*g + b*c + x**3*(a*j + b*f) + x**2*(a*i + b*e) + x*(a*h + b*d))/(8*a*b*(a - b*x**4)**2) + (4*a*(-a*j + b*f) + x*(b*x**2*(-3*a*i + 5*b*e) + 2*b*x*(-a*h + 3*b*d) + b*(-a*g + 7*b*c)))/(32*a**2*b**2*(a - b*x**4)) + (-a*h + 3*b*d)*atanh(sqrt(b)*x**2/sqrt(a))/(16*a**(sympy.S(5)/2)*b**(sympy.S(3)/2)) - (-3*a*i + 5*b*e - 3*sqrt(b)*(-a*g + 7*b*c)/sqrt(a))*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(64*a**(sympy.S(9)/4)*b**(sympy.S(7)/4)) + (-3*a*i + 5*b*e + 3*sqrt(b)*(-a*g + 7*b*c)/sqrt(a))*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(64*a**(sympy.S(9)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_201():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a + b*x**4)**3
    F = x*(-a*g + b*c + b*e*x**2 + b*f*x**3 + x*(-a*h + b*d))/(8*a*b*(a + b*x**4)**2) - (4*a*f - x*(a*g + 7*b*c + 5*b*e*x**2 + x*(2*a*h + 6*b*d)))/(32*a**2*b*(a + b*x**4)) + (a*h + 3*b*d)*atan(sqrt(b)*x**2/sqrt(a))/(16*a**(sympy.S(5)/2)*b**(sympy.S(3)/2)) - sqrt(2)*(-5*sqrt(a)*sqrt(b)*e + 3*a*g + 21*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(256*a**(sympy.S(11)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(-5*sqrt(a)*sqrt(b)*e + 3*a*g + 21*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(256*a**(sympy.S(11)/4)*b**(sympy.S(5)/4)) - sqrt(2)*(5*sqrt(a)*sqrt(b)*e + 3*a*g + 21*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(5*sqrt(a)*sqrt(b)*e + 3*a*g + 21*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_202():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5 + i*x**6)/(a + b*x**4)**3
    F = x*(-a*g + b*c + b*f*x**3 + x**2*(-a*i + b*e) + x*(-a*h + b*d))/(8*a*b*(a + b*x**4)**2) - (4*a*f - x*(a*g + 7*b*c + x**2*(3*a*i + 5*b*e) + x*(2*a*h + 6*b*d)))/(32*a**2*b*(a + b*x**4)) + (a*h + 3*b*d)*atan(sqrt(b)*x**2/sqrt(a))/(16*a**(sympy.S(5)/2)*b**(sympy.S(3)/2)) - sqrt(2)*(-sqrt(a)*(3*a*i + 5*b*e) + 3*sqrt(b)*(a*g + 7*b*c))*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(256*a**(sympy.S(11)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(-sqrt(a)*(3*a*i + 5*b*e) + 3*sqrt(b)*(a*g + 7*b*c))*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(256*a**(sympy.S(11)/4)*b**(sympy.S(7)/4)) - sqrt(2)*(sqrt(a)*(3*a*i + 5*b*e) + 3*sqrt(b)*(a*g + 7*b*c))*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(sqrt(a)*(3*a*i + 5*b*e) + 3*sqrt(b)*(a*g + 7*b*c))*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_203():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5 + i*x**6 + j*x**7)/(a + b*x**4)**3
    F = x*(-a*g + b*c + x**3*(-a*j + b*f) + x**2*(-a*i + b*e) + x*(-a*h + b*d))/(8*a*b*(a + b*x**4)**2) - (4*a*(a*j + b*f) - x*(b*x**2*(3*a*i + 5*b*e) + 2*b*x*(a*h + 3*b*d) + b*(a*g + 7*b*c)))/(32*a**2*b**2*(a + b*x**4)) + (a*h + 3*b*d)*atan(sqrt(b)*x**2/sqrt(a))/(16*a**(sympy.S(5)/2)*b**(sympy.S(3)/2)) - sqrt(2)*(-sqrt(a)*(3*a*i + 5*b*e) + 3*sqrt(b)*(a*g + 7*b*c))*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(256*a**(sympy.S(11)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(-sqrt(a)*(3*a*i + 5*b*e) + 3*sqrt(b)*(a*g + 7*b*c))*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(256*a**(sympy.S(11)/4)*b**(sympy.S(7)/4)) - sqrt(2)*(sqrt(a)*(3*a*i + 5*b*e) + 3*sqrt(b)*(a*g + 7*b*c))*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(sqrt(a)*(3*a*i + 5*b*e) + 3*sqrt(b)*(a*g + 7*b*c))*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_204():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a - b*x**4)**4
    F = x*(a*g + b*c + b*e*x**2 + b*f*x**3 + x*(a*h + b*d))/(12*a*b*(a - b*x**4)**3) + (8*a*f + x*(-a*g + 11*b*c + 9*b*e*x**2 + x*(-2*a*h + 10*b*d)))/(96*a**2*b*(a - b*x**4)**2) + x*(-7*a*g + 77*b*c + 45*b*e*x**2 + x*(-12*a*h + 60*b*d))/(384*a**3*b*(a - b*x**4)) + (-a*h + 5*b*d)*atanh(sqrt(b)*x**2/sqrt(a))/(32*a**(sympy.S(7)/2)*b**(sympy.S(3)/2)) + (-15*sqrt(a)*sqrt(b)*e - 7*a*g + 77*b*c)*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(256*a**(sympy.S(15)/4)*b**(sympy.S(5)/4)) + (15*sqrt(a)*sqrt(b)*e - 7*a*g + 77*b*c)*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(256*a**(sympy.S(15)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_205():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5 + i*x**6)/(a - b*x**4)**4
    F = x*(a*g + b*c + b*f*x**3 + x**2*(a*i + b*e) + x*(a*h + b*d))/(12*a*b*(a - b*x**4)**3) + (8*a*f + x*(-a*g + 11*b*c + x**2*(-3*a*i + 9*b*e) + x*(-2*a*h + 10*b*d)))/(96*a**2*b*(a - b*x**4)**2) + x*(-7*a*g + 77*b*c + x**2*(-15*a*i + 45*b*e) + x*(-12*a*h + 60*b*d))/(384*a**3*b*(a - b*x**4)) + (-a*h + 5*b*d)*atanh(sqrt(b)*x**2/sqrt(a))/(32*a**(sympy.S(7)/2)*b**(sympy.S(3)/2)) + (-5*a*i + 15*b*e + 7*sqrt(b)*(-a*g + 11*b*c)/sqrt(a))*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(256*a**(sympy.S(13)/4)*b**(sympy.S(7)/4)) + (5*a*i - 15*b*e + 7*sqrt(b)*(-a*g + 11*b*c)/sqrt(a))*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(256*a**(sympy.S(13)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_206():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5 + i*x**6 + j*x**7)/(a - b*x**4)**4
    F = x*(a*g + b*c + x**3*(a*j + b*f) + x**2*(a*i + b*e) + x*(a*h + b*d))/(12*a*b*(a - b*x**4)**3) + (4*a*(-a*j + 2*b*f) + x*(3*b*x**2*(-a*i + 3*b*e) + 2*b*x*(-a*h + 5*b*d) + b*(-a*g + 11*b*c)))/(96*a**2*b**2*(a - b*x**4)**2) + x*(-7*a*g + 77*b*c + x**2*(-15*a*i + 45*b*e) + x*(-12*a*h + 60*b*d))/(384*a**3*b*(a - b*x**4)) + (-a*h + 5*b*d)*atanh(sqrt(b)*x**2/sqrt(a))/(32*a**(sympy.S(7)/2)*b**(sympy.S(3)/2)) + (-5*a*i + 15*b*e + 7*sqrt(b)*(-a*g + 11*b*c)/sqrt(a))*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(256*a**(sympy.S(13)/4)*b**(sympy.S(7)/4)) + (5*a*i - 15*b*e + 7*sqrt(b)*(-a*g + 11*b*c)/sqrt(a))*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(256*a**(sympy.S(13)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_207():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a + b*x**4)**4
    F = x*(-a*g + b*c + b*e*x**2 + b*f*x**3 + x*(-a*h + b*d))/(12*a*b*(a + b*x**4)**3) - (8*a*f - x*(a*g + 11*b*c + 9*b*e*x**2 + x*(2*a*h + 10*b*d)))/(96*a**2*b*(a + b*x**4)**2) + x*(7*a*g + 77*b*c + 45*b*e*x**2 + x*(12*a*h + 60*b*d))/(384*a**3*b*(a + b*x**4)) + (a*h + 5*b*d)*atan(sqrt(b)*x**2/sqrt(a))/(32*a**(sympy.S(7)/2)*b**(sympy.S(3)/2)) - sqrt(2)*(-15*sqrt(a)*sqrt(b)*e + 7*a*g + 77*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(1024*a**(sympy.S(15)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(-15*sqrt(a)*sqrt(b)*e + 7*a*g + 77*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(1024*a**(sympy.S(15)/4)*b**(sympy.S(5)/4)) - sqrt(2)*(15*sqrt(a)*sqrt(b)*e + 7*a*g + 77*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(512*a**(sympy.S(15)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(15*sqrt(a)*sqrt(b)*e + 7*a*g + 77*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(512*a**(sympy.S(15)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_208():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5 + i*x**6)/(a + b*x**4)**4
    F = x*(-a*g + b*c + b*f*x**3 + x**2*(-a*i + b*e) + x*(-a*h + b*d))/(12*a*b*(a + b*x**4)**3) - (8*a*f - x*(a*g + 11*b*c + x**2*(3*a*i + 9*b*e) + x*(2*a*h + 10*b*d)))/(96*a**2*b*(a + b*x**4)**2) + x*(7*a*g + 77*b*c + x**2*(15*a*i + 45*b*e) + x*(12*a*h + 60*b*d))/(384*a**3*b*(a + b*x**4)) + (a*h + 5*b*d)*atan(sqrt(b)*x**2/sqrt(a))/(32*a**(sympy.S(7)/2)*b**(sympy.S(3)/2)) - sqrt(2)*(-5*sqrt(a)*(a*i + 3*b*e) + 7*sqrt(b)*(a*g + 11*b*c))*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(1024*a**(sympy.S(15)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(-5*sqrt(a)*(a*i + 3*b*e) + 7*sqrt(b)*(a*g + 11*b*c))*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(1024*a**(sympy.S(15)/4)*b**(sympy.S(7)/4)) - sqrt(2)*(5*sqrt(a)*(a*i + 3*b*e) + 7*sqrt(b)*(a*g + 11*b*c))*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(512*a**(sympy.S(15)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(5*sqrt(a)*(a*i + 3*b*e) + 7*sqrt(b)*(a*g + 11*b*c))*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(512*a**(sympy.S(15)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_209():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5 + i*x**6 + j*x**7)/(a + b*x**4)**4
    F = x*(-a*g + b*c + x**3*(-a*j + b*f) + x**2*(-a*i + b*e) + x*(-a*h + b*d))/(12*a*b*(a + b*x**4)**3) - (4*a*(a*j + 2*b*f) - x*(3*b*x**2*(a*i + 3*b*e) + 2*b*x*(a*h + 5*b*d) + b*(a*g + 11*b*c)))/(96*a**2*b**2*(a + b*x**4)**2) + x*(7*a*g + 77*b*c + x**2*(15*a*i + 45*b*e) + x*(12*a*h + 60*b*d))/(384*a**3*b*(a + b*x**4)) + (a*h + 5*b*d)*atan(sqrt(b)*x**2/sqrt(a))/(32*a**(sympy.S(7)/2)*b**(sympy.S(3)/2)) - sqrt(2)*(-5*sqrt(a)*(a*i + 3*b*e) + 7*sqrt(b)*(a*g + 11*b*c))*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(1024*a**(sympy.S(15)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(-5*sqrt(a)*(a*i + 3*b*e) + 7*sqrt(b)*(a*g + 11*b*c))*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(1024*a**(sympy.S(15)/4)*b**(sympy.S(7)/4)) - sqrt(2)*(5*sqrt(a)*(a*i + 3*b*e) + 7*sqrt(b)*(a*g + 11*b*c))*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(512*a**(sympy.S(15)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(5*sqrt(a)*(a*i + 3*b*e) + 7*sqrt(b)*(a*g + 11*b*c))*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(512*a**(sympy.S(15)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_210():
    f = (c + d*x)/sqrt(a + b*x**4)
    F = d*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(2*sqrt(b)) + c*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_211():
    f = (c + d*x)/sqrt(a - b*x**4)
    F = a**(sympy.S(1)/4)*c*sqrt(1 - b*x**4/a)*elliptic_f(asin(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), -1)/(b**(sympy.S(1)/4)*sqrt(a - b*x**4)) + d*atan(sqrt(b)*x**2/sqrt(a - b*x**4))/(2*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_212():
    f = (c + d*x)/sqrt(-a + b*x**4)
    F = a**(sympy.S(1)/4)*c*sqrt(1 - b*x**4/a)*elliptic_f(asin(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), -1)/(b**(sympy.S(1)/4)*sqrt(-a + b*x**4)) + d*atanh(sqrt(b)*x**2/sqrt(-a + b*x**4))/(2*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_213():
    f = (c + d*x)/sqrt(-a - b*x**4)
    F = d*atan(sqrt(b)*x**2/sqrt(-a - b*x**4))/(2*sqrt(b)) + c*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*sqrt(-a - b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_214():
    f = (c + d*x + e*x**2)/sqrt(a + b*x**4)
    F = -a**(sympy.S(1)/4)*e*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(b**(sympy.S(3)/4)*sqrt(a + b*x**4)) + a**(sympy.S(1)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(e + sqrt(b)*c/sqrt(a))*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*b**(sympy.S(3)/4)*sqrt(a + b*x**4)) + d*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(2*sqrt(b)) + e*x*sqrt(a + b*x**4)/(sqrt(b)*(sqrt(a) + sqrt(b)*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_215():
    f = (a*g - b*g*x**4)/(a + b*x**4)**(sympy.S(3)/2)
    F = g*x/sqrt(a + b*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_216():
    f = (a*g - b*g*x**4 + e*x)/(a + b*x**4)**(sympy.S(3)/2)
    F = (2*a*g*x + e*x**2)/(2*a*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_217():
    f = (a*g - b*g*x**4 + f*x**3)/(a + b*x**4)**(sympy.S(3)/2)
    F = -(-2*b*g*x + f)/(2*b*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_218():
    f = (a*g - b*g*x**4 + e*x + f*x**3)/(a + b*x**4)**(sympy.S(3)/2)
    F = -(-2*a*b*g*x + a*f - b*e*x**2)/(2*a*b*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_219():
    f = (x**4 - 1)/(x**4 + 1)**(sympy.S(3)/2)
    F = -x/sqrt(x**4 + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_220():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5 + i*x**6)/sqrt(a + b*x**4)
    F = -a**(sympy.S(1)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(-3*a*i + 5*b*e)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(7)/4)*sqrt(a + b*x**4)) + a**(sympy.S(1)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(-9*a*i + 15*b*e + 5*sqrt(b)*(-a*g + 3*b*c)/sqrt(a))*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(30*b**(sympy.S(7)/4)*sqrt(a + b*x**4)) + f*sqrt(a + b*x**4)/(2*b) + g*x*sqrt(a + b*x**4)/(3*b) + h*x**2*sqrt(a + b*x**4)/(4*b) + i*x**3*sqrt(a + b*x**4)/(5*b) + x*sqrt(a + b*x**4)*(-3*a*i + 5*b*e)/(5*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x**2)) + (-a*h + 2*b*d)*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(4*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_221():
    f = (x + 1)/(x**5 + 1)
    F = -(-1)**(sympy.S(1)/5)*(1 + (-1)**(sympy.S(1)/5))*log(-x + (-1)**(sympy.S(1)/5))/5 + (-1)**(sympy.S(4)/5)*(1 - (-1)**(sympy.S(4)/5))*log(-x - (-1)**(sympy.S(4)/5))/5 + (-1)**(sympy.S(2)/5)*(1 - (-1)**(sympy.S(2)/5))*log(x + (-1)**(sympy.S(2)/5))/5 - (-1)**(sympy.S(3)/5)*(1 + (-1)**(sympy.S(3)/5))*log(x - (-1)**(sympy.S(3)/5))/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_222():
    f = (1 - x)/(1 - x**5)
    F = -(-1)**(sympy.S(2)/5)*(1 - (-1)**(sympy.S(2)/5))*log(-x + (-1)**(sympy.S(2)/5))/5 + (-1)**(sympy.S(3)/5)*(1 + (-1)**(sympy.S(3)/5))*log(-x - (-1)**(sympy.S(3)/5))/5 + (-1)**(sympy.S(1)/5)*(1 + (-1)**(sympy.S(1)/5))*log(x + (-1)**(sympy.S(1)/5))/5 - (-1)**(sympy.S(4)/5)*(1 - (-1)**(sympy.S(4)/5))*log(x - (-1)**(sympy.S(4)/5))/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_223():
    f = x**11*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)
    F = -a**3*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a + b*x**3)/(3*b**7) + a**2*x**3*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*b**6) - a*x**6*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*b**5) + f*x**18/(18*b) + x**15*(-a*f + b*e)/(15*b**2) + x**12*(a**2*f - a*b*e + b**2*d)/(12*b**3) + x**9*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(9*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_224():
    f = x**8*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)
    F = a**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a + b*x**3)/(3*b**6) - a*x**3*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*b**5) + f*x**15/(15*b) + x**12*(-a*f + b*e)/(12*b**2) + x**9*(a**2*f - a*b*e + b**2*d)/(9*b**3) + x**6*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_225():
    f = x**5*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)
    F = -a*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a + b*x**3)/(3*b**5) + f*x**12/(12*b) + x**9*(-a*f + b*e)/(9*b**2) + x**6*(a**2*f - a*b*e + b**2*d)/(6*b**3) + x**3*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_226():
    f = x**2*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)
    F = f*x**9/(9*b) + x**6*(-a*f + b*e)/(6*b**2) + x**3*(a**2*f - a*b*e + b**2*d)/(3*b**3) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a + b*x**3)/(3*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_227():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x*(a + b*x**3))
    F = f*x**6/(6*b) + x**3*(-a*f + b*e)/(3*b**2) + c*log(x)/a - (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a + b*x**3)/(3*a*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_228():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**4*(a + b*x**3))
    F = f*x**3/(3*b) - c/(3*a*x**3) - (-a*d + b*c)*log(x)/a**2 + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a + b*x**3)/(3*a**2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_229():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**7*(a + b*x**3))
    F = -c/(6*a*x**6) + (-a*d + b*c)/(3*a**2*x**3) + (a**2*e - a*b*d + b**2*c)*log(x)/a**3 - (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a + b*x**3)/(3*a**3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_230():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**10*(a + b*x**3))
    F = -c/(9*a*x**9) + (-a*d + b*c)/(6*a**2*x**6) - (a**2*e - a*b*d + b**2*c)/(3*a**3*x**3) - (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(x)/a**4 + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a + b*x**3)/(3*a**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_231():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**13*(a + b*x**3))
    F = -c/(12*a*x**12) + (-a*d + b*c)/(9*a**2*x**9) - (a**2*e - a*b*d + b**2*c)/(6*a**3*x**6) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*a**4*x**3) + b*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(x)/a**5 - b*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a + b*x**3)/(3*a**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_232():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**16*(a + b*x**3))
    F = -c/(15*a*x**15) + (-a*d + b*c)/(12*a**2*x**12) - (a**2*e - a*b*d + b**2*c)/(9*a**3*x**9) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*a**4*x**6) - b*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*a**5*x**3) - b**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(x)/a**6 + b**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a + b*x**3)/(3*a**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_233():
    f = x**9*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)
    F = -a**(sympy.S(7)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*b**(sympy.S(19)/3)) + a**(sympy.S(7)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*b**(sympy.S(19)/3)) + sqrt(3)*a**(sympy.S(7)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(19)/3)) + a**2*x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/b**6 - a*x**4*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(4*b**5) + f*x**16/(16*b) + x**13*(-a*f + b*e)/(13*b**2) + x**10*(a**2*f - a*b*e + b**2*d)/(10*b**3) + x**7*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(7*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_234():
    f = x**7*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)
    F = -a**(sympy.S(5)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*b**(sympy.S(17)/3)) + a**(sympy.S(5)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*b**(sympy.S(17)/3)) - sqrt(3)*a**(sympy.S(5)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(17)/3)) - a*x**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(2*b**5) + f*x**14/(14*b) + x**11*(-a*f + b*e)/(11*b**2) + x**8*(a**2*f - a*b*e + b**2*d)/(8*b**3) + x**5*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(5*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_235():
    f = x**6*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)
    F = a**(sympy.S(4)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*b**(sympy.S(16)/3)) - a**(sympy.S(4)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*b**(sympy.S(16)/3)) - sqrt(3)*a**(sympy.S(4)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(16)/3)) - a*x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/b**5 + f*x**13/(13*b) + x**10*(-a*f + b*e)/(10*b**2) + x**7*(a**2*f - a*b*e + b**2*d)/(7*b**3) + x**4*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(4*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_236():
    f = x**4*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)
    F = a**(sympy.S(2)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*b**(sympy.S(14)/3)) - a**(sympy.S(2)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*b**(sympy.S(14)/3)) + sqrt(3)*a**(sympy.S(2)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(14)/3)) + f*x**11/(11*b) + x**8*(-a*f + b*e)/(8*b**2) + x**5*(a**2*f - a*b*e + b**2*d)/(5*b**3) + x**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_237():
    f = x**3*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)
    F = -a**(sympy.S(1)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*b**(sympy.S(13)/3)) + a**(sympy.S(1)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*b**(sympy.S(13)/3)) + sqrt(3)*a**(sympy.S(1)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(13)/3)) + f*x**10/(10*b) + x**7*(-a*f + b*e)/(7*b**2) + x**4*(a**2*f - a*b*e + b**2*d)/(4*b**3) + x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/b**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_238():
    f = x*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)
    F = f*x**8/(8*b) + x**5*(-a*f + b*e)/(5*b**2) + x**2*(a**2*f - a*b*e + b**2*d)/(2*b**3) - (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)*b**(sympy.S(11)/3)) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(1)/3)*b**(sympy.S(11)/3)) - sqrt(3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3)*b**(sympy.S(11)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_239():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)
    F = f*x**7/(7*b) + x**4*(-a*f + b*e)/(4*b**2) + x*(a**2*f - a*b*e + b**2*d)/b**3 + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(2)/3)*b**(sympy.S(10)/3)) - (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(10)/3)) - sqrt(3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(10)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_240():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**2*(a + b*x**3))
    F = f*x**5/(5*b) + x**2*(-a*f + b*e)/(2*b**2) - c/(a*x) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(4)/3)*b**(sympy.S(8)/3)) - (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(4)/3)*b**(sympy.S(8)/3)) + sqrt(3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(4)/3)*b**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_241():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**3*(a + b*x**3))
    F = f*x**4/(4*b) + x*(-a*f + b*e)/b**2 - c/(2*a*x**2) - (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(5)/3)*b**(sympy.S(7)/3)) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(5)/3)*b**(sympy.S(7)/3)) + sqrt(3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(5)/3)*b**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_242():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**5*(a + b*x**3))
    F = f*x**2/(2*b) - c/(4*a*x**4) + (-a*d + b*c)/(a**2*x) - (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(7)/3)*b**(sympy.S(5)/3)) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(7)/3)*b**(sympy.S(5)/3)) - sqrt(3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(7)/3)*b**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_243():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**6*(a + b*x**3))
    F = f*x/b - c/(5*a*x**5) + (-a*d + b*c)/(2*a**2*x**2) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(8)/3)*b**(sympy.S(4)/3)) - (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(8)/3)*b**(sympy.S(4)/3)) - sqrt(3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(8)/3)*b**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_244():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**8*(a + b*x**3))
    F = -c/(7*a*x**7) + (-a*d + b*c)/(4*a**2*x**4) - (a**2*e - a*b*d + b**2*c)/(a**3*x) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(10)/3)*b**(sympy.S(2)/3)) - (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(10)/3)*b**(sympy.S(2)/3)) + sqrt(3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(10)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_245():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**9*(a + b*x**3))
    F = -c/(8*a*x**8) + (-a*d + b*c)/(5*a**2*x**5) - (a**2*e - a*b*d + b**2*c)/(2*a**3*x**2) - (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(11)/3)*b**(sympy.S(1)/3)) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(11)/3)*b**(sympy.S(1)/3)) + sqrt(3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(11)/3)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_246():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**11*(a + b*x**3))
    F = -c/(10*a*x**10) + (-a*d + b*c)/(7*a**2*x**7) - (a**2*e - a*b*d + b**2*c)/(4*a**3*x**4) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(a**4*x) - b**(sympy.S(1)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(13)/3)) + b**(sympy.S(1)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(13)/3)) - sqrt(3)*b**(sympy.S(1)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(13)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_247():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**12*(a + b*x**3))
    F = -c/(11*a*x**11) + (-a*d + b*c)/(8*a**2*x**8) - (a**2*e - a*b*d + b**2*c)/(5*a**3*x**5) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(2*a**4*x**2) + b**(sympy.S(2)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(14)/3)) - b**(sympy.S(2)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(14)/3)) - sqrt(3)*b**(sympy.S(2)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(14)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_248():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**14*(a + b*x**3))
    F = -c/(13*a*x**13) + (-a*d + b*c)/(10*a**2*x**10) - (a**2*e - a*b*d + b**2*c)/(7*a**3*x**7) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(4*a**4*x**4) - b*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(a**5*x) + b**(sympy.S(4)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(16)/3)) - b**(sympy.S(4)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(16)/3)) + sqrt(3)*b**(sympy.S(4)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(16)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_249():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**15*(a + b*x**3))
    F = -c/(14*a*x**14) + (-a*d + b*c)/(11*a**2*x**11) - (a**2*e - a*b*d + b**2*c)/(8*a**3*x**8) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(5*a**4*x**5) - b*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(2*a**5*x**2) - b**(sympy.S(5)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(17)/3)) + b**(sympy.S(5)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(17)/3)) + sqrt(3)*b**(sympy.S(5)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(17)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_250():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**17*(a + b*x**3))
    F = -c/(16*a*x**16) + (-a*d + b*c)/(13*a**2*x**13) - (a**2*e - a*b*d + b**2*c)/(10*a**3*x**10) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(7*a**4*x**7) - b*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(4*a**5*x**4) + b**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(a**6*x) - b**(sympy.S(7)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(19)/3)) + b**(sympy.S(7)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(19)/3)) - sqrt(3)*b**(sympy.S(7)/3)*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(19)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_251():
    f = x**11*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**2
    F = a**3*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*b**7*(a + b*x**3)) + a**2*(-6*a**3*f + 5*a**2*b*e - 4*a*b**2*d + 3*b**3*c)*log(a + b*x**3)/(3*b**7) - a*x**3*(-5*a**3*f + 4*a**2*b*e - 3*a*b**2*d + 2*b**3*c)/(3*b**6) + f*x**15/(15*b**2) + x**12*(-2*a*f + b*e)/(12*b**3) + x**9*(3*a**2*f - 2*a*b*e + b**2*d)/(9*b**4) + x**6*(-4*a**3*f + 3*a**2*b*e - 2*a*b**2*d + b**3*c)/(6*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_252():
    f = x**8*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**2
    F = -a**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*b**6*(a + b*x**3)) - a*(-5*a**3*f + 4*a**2*b*e - 3*a*b**2*d + 2*b**3*c)*log(a + b*x**3)/(3*b**6) + f*x**12/(12*b**2) + x**9*(-2*a*f + b*e)/(9*b**3) + x**6*(3*a**2*f - 2*a*b*e + b**2*d)/(6*b**4) + x**3*(-4*a**3*f + 3*a**2*b*e - 2*a*b**2*d + b**3*c)/(3*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_253():
    f = x**5*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**2
    F = a*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*b**5*(a + b*x**3)) + f*x**9/(9*b**2) + x**6*(-2*a*f + b*e)/(6*b**3) + x**3*(3*a**2*f - 2*a*b*e + b**2*d)/(3*b**4) + (-4*a**3*f + 3*a**2*b*e - 2*a*b**2*d + b**3*c)*log(a + b*x**3)/(3*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_254():
    f = x**2*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**2
    F = f*x**6/(6*b**2) + x**3*(-2*a*f + b*e)/(3*b**3) + (3*a**2*f - 2*a*b*e + b**2*d)*log(a + b*x**3)/(3*b**4) - (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*b**4*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_255():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x*(a + b*x**3)**2)
    F = f*x**3/(3*b**2) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*a*b**3*(a + b*x**3)) + c*log(x)/a**2 - (2*a**3*f - a**2*b*e + b**3*c)*log(a + b*x**3)/(3*a**2*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_256():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**4*(a + b*x**3)**2)
    F = -c/(3*a**2*x**3) - (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*a**2*b**2*(a + b*x**3)) - (-a*d + 2*b*c)*log(x)/a**3 + (a**3*f - a*b**2*d + 2*b**3*c)*log(a + b*x**3)/(3*a**3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_257():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**7*(a + b*x**3)**2)
    F = -c/(6*a**2*x**6) + (-a*d + 2*b*c)/(3*a**3*x**3) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*a**3*b*(a + b*x**3)) + (a**2*e - 2*a*b*d + 3*b**2*c)*log(x)/a**4 - (a**2*e - 2*a*b*d + 3*b**2*c)*log(a + b*x**3)/(3*a**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_258():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**10*(a + b*x**3)**2)
    F = -c/(9*a**2*x**9) + (-a*d + 2*b*c)/(6*a**3*x**6) - (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*a**4*(a + b*x**3)) - (a**2*e - 2*a*b*d + 3*b**2*c)/(3*a**4*x**3) - (-a**3*f + 2*a**2*b*e - 3*a*b**2*d + 4*b**3*c)*log(x)/a**5 + (-a**3*f + 2*a**2*b*e - 3*a*b**2*d + 4*b**3*c)*log(a + b*x**3)/(3*a**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_259():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**13*(a + b*x**3)**2)
    F = -c/(12*a**2*x**12) + (-a*d + 2*b*c)/(9*a**3*x**9) - (a**2*e - 2*a*b*d + 3*b**2*c)/(6*a**4*x**6) + b*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*a**5*(a + b*x**3)) + (-a**3*f + 2*a**2*b*e - 3*a*b**2*d + 4*b**3*c)/(3*a**5*x**3) + b*(-2*a**3*f + 3*a**2*b*e - 4*a*b**2*d + 5*b**3*c)*log(x)/a**6 - b*(-2*a**3*f + 3*a**2*b*e - 4*a*b**2*d + 5*b**3*c)*log(a + b*x**3)/(3*a**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_260():
    f = x**9*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**2
    F = a**(sympy.S(4)/3)*(-16*a**3*f + 13*a**2*b*e - 10*a*b**2*d + 7*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*b**(sympy.S(19)/3)) - a**(sympy.S(4)/3)*(-16*a**3*f + 13*a**2*b*e - 10*a*b**2*d + 7*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*b**(sympy.S(19)/3)) - sqrt(3)*a**(sympy.S(4)/3)*(-16*a**3*f + 13*a**2*b*e - 10*a*b**2*d + 7*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*b**(sympy.S(19)/3)) - a**2*x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*b**6*(a + b*x**3)) - a*x*(-5*a**3*f + 4*a**2*b*e - 3*a*b**2*d + 2*b**3*c)/b**6 + f*x**13/(13*b**2) + x**10*(-2*a*f + b*e)/(10*b**3) + x**7*(3*a**2*f - 2*a*b*e + b**2*d)/(7*b**4) + x**4*(-4*a**3*f + 3*a**2*b*e - 2*a*b**2*d + b**3*c)/(4*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_261():
    f = x**7*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**2
    F = a**(sympy.S(2)/3)*(-14*a**3*f + 11*a**2*b*e - 8*a*b**2*d + 5*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*b**(sympy.S(17)/3)) - a**(sympy.S(2)/3)*(-14*a**3*f + 11*a**2*b*e - 8*a*b**2*d + 5*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*b**(sympy.S(17)/3)) + sqrt(3)*a**(sympy.S(2)/3)*(-14*a**3*f + 11*a**2*b*e - 8*a*b**2*d + 5*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*b**(sympy.S(17)/3)) + a*x**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*b**5*(a + b*x**3)) + f*x**11/(11*b**2) + x**8*(-2*a*f + b*e)/(8*b**3) + x**5*(3*a**2*f - 2*a*b*e + b**2*d)/(5*b**4) + x**2*(-4*a**3*f + 3*a**2*b*e - 2*a*b**2*d + b**3*c)/(2*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_262():
    f = x**6*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**2
    F = -a**(sympy.S(1)/3)*(-13*a**3*f + 10*a**2*b*e - 7*a*b**2*d + 4*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*b**(sympy.S(16)/3)) + a**(sympy.S(1)/3)*(-13*a**3*f + 10*a**2*b*e - 7*a*b**2*d + 4*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*b**(sympy.S(16)/3)) + sqrt(3)*a**(sympy.S(1)/3)*(-13*a**3*f + 10*a**2*b*e - 7*a*b**2*d + 4*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*b**(sympy.S(16)/3)) + a*x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*b**5*(a + b*x**3)) + f*x**10/(10*b**2) + x**7*(-2*a*f + b*e)/(7*b**3) + x**4*(3*a**2*f - 2*a*b*e + b**2*d)/(4*b**4) + x*(-4*a**3*f + 3*a**2*b*e - 2*a*b**2*d + b**3*c)/b**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_263():
    f = x**4*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**2
    F = f*x**8/(8*b**2) + x**5*(-2*a*f + b*e)/(5*b**3) + x**2*(3*a**2*f - 2*a*b*e + b**2*d)/(2*b**4) - x**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*b**4*(a + b*x**3)) - (-11*a**3*f + 8*a**2*b*e - 5*a*b**2*d + 2*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(1)/3)*b**(sympy.S(14)/3)) + (-11*a**3*f + 8*a**2*b*e - 5*a*b**2*d + 2*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(1)/3)*b**(sympy.S(14)/3)) - sqrt(3)*(-11*a**3*f + 8*a**2*b*e - 5*a*b**2*d + 2*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(1)/3)*b**(sympy.S(14)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_264():
    f = x**3*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**2
    F = f*x**7/(7*b**2) + x**4*(-2*a*f + b*e)/(4*b**3) + x*(3*a**2*f - 2*a*b*e + b**2*d)/b**4 - x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*b**4*(a + b*x**3)) + (-10*a**3*f + 7*a**2*b*e - 4*a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(2)/3)*b**(sympy.S(13)/3)) - (-10*a**3*f + 7*a**2*b*e - 4*a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(2)/3)*b**(sympy.S(13)/3)) - sqrt(3)*(-10*a**3*f + 7*a**2*b*e - 4*a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(2)/3)*b**(sympy.S(13)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_265():
    f = x*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**2
    F = f*x**5/(5*b**2) + x**2*(-2*a*f + b*e)/(2*b**3) + x**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*a*b**3*(a + b*x**3)) - (8*a**3*f - 5*a**2*b*e + 2*a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(4)/3)*b**(sympy.S(11)/3)) + (8*a**3*f - 5*a**2*b*e + 2*a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(4)/3)*b**(sympy.S(11)/3)) - sqrt(3)*(8*a**3*f - 5*a**2*b*e + 2*a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(4)/3)*b**(sympy.S(11)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_266():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**2
    F = f*x**4/(4*b**2) + x*(-2*a*f + b*e)/b**3 + x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*a*b**3*(a + b*x**3)) + (7*a**3*f - 4*a**2*b*e + a*b**2*d + 2*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(5)/3)*b**(sympy.S(10)/3)) - (7*a**3*f - 4*a**2*b*e + a*b**2*d + 2*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(5)/3)*b**(sympy.S(10)/3)) - sqrt(3)*(7*a**3*f - 4*a**2*b*e + a*b**2*d + 2*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*b**(sympy.S(10)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_267():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**2*(a + b*x**3)**2)
    F = f*x**2/(2*b**2) - c/(a**2*x) - x**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*a**2*b**2*(a + b*x**3)) + (5*a**3*f - 2*a**2*b*e - a*b**2*d + 4*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(7)/3)*b**(sympy.S(8)/3)) - (5*a**3*f - 2*a**2*b*e - a*b**2*d + 4*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(7)/3)*b**(sympy.S(8)/3)) + sqrt(3)*(5*a**3*f - 2*a**2*b*e - a*b**2*d + 4*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(7)/3)*b**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_268():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**3*(a + b*x**3)**2)
    F = f*x/b**2 - c/(2*a**2*x**2) - x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*a**2*b**2*(a + b*x**3)) - (4*a**3*f - a**2*b*e - 2*a*b**2*d + 5*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(8)/3)*b**(sympy.S(7)/3)) + (4*a**3*f - a**2*b*e - 2*a*b**2*d + 5*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(8)/3)*b**(sympy.S(7)/3)) + sqrt(3)*(4*a**3*f - a**2*b*e - 2*a*b**2*d + 5*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(8)/3)*b**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_269():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**5*(a + b*x**3)**2)
    F = -c/(4*a**2*x**4) + (-a*d + 2*b*c)/(a**3*x) + x**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*a**3*b*(a + b*x**3)) - (2*a**3*f + a**2*b*e - 4*a*b**2*d + 7*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(10)/3)*b**(sympy.S(5)/3)) + (2*a**3*f + a**2*b*e - 4*a*b**2*d + 7*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(10)/3)*b**(sympy.S(5)/3)) - sqrt(3)*(2*a**3*f + a**2*b*e - 4*a*b**2*d + 7*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(10)/3)*b**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_270():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**6*(a + b*x**3)**2)
    F = -c/(5*a**2*x**5) + (-a*d + 2*b*c)/(2*a**3*x**2) + x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*a**3*b*(a + b*x**3)) + (a**3*f + 2*a**2*b*e - 5*a*b**2*d + 8*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(11)/3)*b**(sympy.S(4)/3)) - (a**3*f + 2*a**2*b*e - 5*a*b**2*d + 8*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(11)/3)*b**(sympy.S(4)/3)) - sqrt(3)*(a**3*f + 2*a**2*b*e - 5*a*b**2*d + 8*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(11)/3)*b**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_271():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**8*(a + b*x**3)**2)
    F = -c/(7*a**2*x**7) + (-a*d + 2*b*c)/(4*a**3*x**4) - x**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*a**4*(a + b*x**3)) - (a**2*e - 2*a*b*d + 3*b**2*c)/(a**4*x) + (-a**3*f + 4*a**2*b*e - 7*a*b**2*d + 10*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(13)/3)*b**(sympy.S(2)/3)) - (-a**3*f + 4*a**2*b*e - 7*a*b**2*d + 10*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(13)/3)*b**(sympy.S(2)/3)) + sqrt(3)*(-a**3*f + 4*a**2*b*e - 7*a*b**2*d + 10*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(13)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_272():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**9*(a + b*x**3)**2)
    F = -c/(8*a**2*x**8) + (-a*d + 2*b*c)/(5*a**3*x**5) - x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*a**4*(a + b*x**3)) - (a**2*e - 2*a*b*d + 3*b**2*c)/(2*a**4*x**2) - (-2*a**3*f + 5*a**2*b*e - 8*a*b**2*d + 11*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(14)/3)*b**(sympy.S(1)/3)) + (-2*a**3*f + 5*a**2*b*e - 8*a*b**2*d + 11*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(14)/3)*b**(sympy.S(1)/3)) + sqrt(3)*(-2*a**3*f + 5*a**2*b*e - 8*a*b**2*d + 11*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(14)/3)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_273():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**11*(a + b*x**3)**2)
    F = -c/(10*a**2*x**10) + (-a*d + 2*b*c)/(7*a**3*x**7) - (a**2*e - 2*a*b*d + 3*b**2*c)/(4*a**4*x**4) + b*x**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*a**5*(a + b*x**3)) + (-a**3*f + 2*a**2*b*e - 3*a*b**2*d + 4*b**3*c)/(a**5*x) - b**(sympy.S(1)/3)*(-4*a**3*f + 7*a**2*b*e - 10*a*b**2*d + 13*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(16)/3)) + b**(sympy.S(1)/3)*(-4*a**3*f + 7*a**2*b*e - 10*a*b**2*d + 13*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(16)/3)) - sqrt(3)*b**(sympy.S(1)/3)*(-4*a**3*f + 7*a**2*b*e - 10*a*b**2*d + 13*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(16)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_274():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**12*(a + b*x**3)**2)
    F = -c/(11*a**2*x**11) + (-a*d + 2*b*c)/(8*a**3*x**8) - (a**2*e - 2*a*b*d + 3*b**2*c)/(5*a**4*x**5) + b*x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*a**5*(a + b*x**3)) + (-a**3*f + 2*a**2*b*e - 3*a*b**2*d + 4*b**3*c)/(2*a**5*x**2) + b**(sympy.S(2)/3)*(-5*a**3*f + 8*a**2*b*e - 11*a*b**2*d + 14*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(17)/3)) - b**(sympy.S(2)/3)*(-5*a**3*f + 8*a**2*b*e - 11*a*b**2*d + 14*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(17)/3)) - sqrt(3)*b**(sympy.S(2)/3)*(-5*a**3*f + 8*a**2*b*e - 11*a*b**2*d + 14*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(17)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_275():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**14*(a + b*x**3)**2)
    F = -c/(13*a**2*x**13) + (-a*d + 2*b*c)/(10*a**3*x**10) - (a**2*e - 2*a*b*d + 3*b**2*c)/(7*a**4*x**7) + (-a**3*f + 2*a**2*b*e - 3*a*b**2*d + 4*b**3*c)/(4*a**5*x**4) - b**2*x**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(3*a**6*(a + b*x**3)) - b*(-2*a**3*f + 3*a**2*b*e - 4*a*b**2*d + 5*b**3*c)/(a**6*x) + b**(sympy.S(4)/3)*(-7*a**3*f + 10*a**2*b*e - 13*a*b**2*d + 16*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(19)/3)) - b**(sympy.S(4)/3)*(-7*a**3*f + 10*a**2*b*e - 13*a*b**2*d + 16*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(19)/3)) + sqrt(3)*b**(sympy.S(4)/3)*(-7*a**3*f + 10*a**2*b*e - 13*a*b**2*d + 16*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(19)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_276():
    f = x**14*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**3
    F = -a**4*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*b**8*(a + b*x**3)**2) + a**3*(-7*a**3*f + 6*a**2*b*e - 5*a*b**2*d + 4*b**3*c)/(3*b**8*(a + b*x**3)) + a**2*(-21*a**3*f + 15*a**2*b*e - 10*a*b**2*d + 6*b**3*c)*log(a + b*x**3)/(3*b**8) - a*x**3*(-15*a**3*f + 10*a**2*b*e - 6*a*b**2*d + 3*b**3*c)/(3*b**7) + f*x**15/(15*b**3) + x**12*(-3*a*f + b*e)/(12*b**4) + x**9*(6*a**2*f - 3*a*b*e + b**2*d)/(9*b**5) + x**6*(-10*a**3*f + 6*a**2*b*e - 3*a*b**2*d + b**3*c)/(6*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_277():
    f = x**11*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**3
    F = a**3*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*b**7*(a + b*x**3)**2) - a**2*(-6*a**3*f + 5*a**2*b*e - 4*a*b**2*d + 3*b**3*c)/(3*b**7*(a + b*x**3)) - a*(-15*a**3*f + 10*a**2*b*e - 6*a*b**2*d + 3*b**3*c)*log(a + b*x**3)/(3*b**7) + f*x**12/(12*b**3) + x**9*(-3*a*f + b*e)/(9*b**4) + x**6*(6*a**2*f - 3*a*b*e + b**2*d)/(6*b**5) + x**3*(-10*a**3*f + 6*a**2*b*e - 3*a*b**2*d + b**3*c)/(3*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_278():
    f = x**8*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**3
    F = -a**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*b**6*(a + b*x**3)**2) + a*(-5*a**3*f + 4*a**2*b*e - 3*a*b**2*d + 2*b**3*c)/(3*b**6*(a + b*x**3)) + f*x**9/(9*b**3) + x**6*(-3*a*f + b*e)/(6*b**4) + x**3*(6*a**2*f - 3*a*b*e + b**2*d)/(3*b**5) + (-10*a**3*f + 6*a**2*b*e - 3*a*b**2*d + b**3*c)*log(a + b*x**3)/(3*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_279():
    f = x**5*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**3
    F = a*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*b**5*(a + b*x**3)**2) + f*x**6/(6*b**3) + x**3*(-3*a*f + b*e)/(3*b**4) + (6*a**2*f - 3*a*b*e + b**2*d)*log(a + b*x**3)/(3*b**5) - (-4*a**3*f + 3*a**2*b*e - 2*a*b**2*d + b**3*c)/(3*b**5*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_280():
    f = x**2*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**3
    F = f*x**3/(3*b**3) + (-3*a*f + b*e)*log(a + b*x**3)/(3*b**4) - (3*a**2*f - 2*a*b*e + b**2*d)/(3*b**4*(a + b*x**3)) - (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*b**4*(a + b*x**3)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_281():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x*(a + b*x**3)**3)
    F = -(-f/(3*b**3) + c/(3*a**3))*log(a + b*x**3) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*a*b**3*(a + b*x**3)**2) + (2*a**3*f - a**2*b*e + b**3*c)/(3*a**2*b**3*(a + b*x**3)) + c*log(x)/a**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_282():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**4*(a + b*x**3)**3)
    F = -(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*a**2*b**2*(a + b*x**3)**2) - c/(3*a**3*x**3) - (a**3*f - a*b**2*d + 2*b**3*c)/(3*a**3*b**2*(a + b*x**3)) - (-a*d + 3*b*c)*log(x)/a**4 + (-a*d + 3*b*c)*log(a + b*x**3)/(3*a**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_283():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**7*(a + b*x**3)**3)
    F = -c/(6*a**3*x**6) + (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*a**3*b*(a + b*x**3)**2) + (a**2*e - 2*a*b*d + 3*b**2*c)/(3*a**4*(a + b*x**3)) + (-a*d + 3*b*c)/(3*a**4*x**3) + (a**2*e - 3*a*b*d + 6*b**2*c)*log(x)/a**5 - (a**2*e - 3*a*b*d + 6*b**2*c)*log(a + b*x**3)/(3*a**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_284():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**10*(a + b*x**3)**3)
    F = -c/(9*a**3*x**9) - (-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*a**4*(a + b*x**3)**2) + (-a*d + 3*b*c)/(6*a**4*x**6) - (-a**3*f + 2*a**2*b*e - 3*a*b**2*d + 4*b**3*c)/(3*a**5*(a + b*x**3)) - (a**2*e - 3*a*b*d + 6*b**2*c)/(3*a**5*x**3) - (-a**3*f + 3*a**2*b*e - 6*a*b**2*d + 10*b**3*c)*log(x)/a**6 + (-a**3*f + 3*a**2*b*e - 6*a*b**2*d + 10*b**3*c)*log(a + b*x**3)/(3*a**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_285():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**13*(a + b*x**3)**3)
    F = -c/(12*a**3*x**12) + (-a*d + 3*b*c)/(9*a**4*x**9) + b*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*a**5*(a + b*x**3)**2) - (a**2*e - 3*a*b*d + 6*b**2*c)/(6*a**5*x**6) + b*(-2*a**3*f + 3*a**2*b*e - 4*a*b**2*d + 5*b**3*c)/(3*a**6*(a + b*x**3)) + (-a**3*f + 3*a**2*b*e - 6*a*b**2*d + 10*b**3*c)/(3*a**6*x**3) + b*(-3*a**3*f + 6*a**2*b*e - 10*a*b**2*d + 15*b**3*c)*log(x)/a**7 - b*(-3*a**3*f + 6*a**2*b*e - 10*a*b**2*d + 15*b**3*c)*log(a + b*x**3)/(3*a**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_286():
    f = x**12*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**3
    F = a**(sympy.S(4)/3)*(-152*a**3*f + 104*a**2*b*e - 65*a*b**2*d + 35*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*b**(sympy.S(22)/3)) - a**(sympy.S(4)/3)*(-152*a**3*f + 104*a**2*b*e - 65*a*b**2*d + 35*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*b**(sympy.S(22)/3)) - sqrt(3)*a**(sympy.S(4)/3)*(-152*a**3*f + 104*a**2*b*e - 65*a*b**2*d + 35*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*b**(sympy.S(22)/3)) + a**3*x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*b**7*(a + b*x**3)**2) - a**2*x*(-37*a**3*f + 31*a**2*b*e - 25*a*b**2*d + 19*b**3*c)/(18*b**7*(a + b*x**3)) - a*x*(-15*a**3*f + 10*a**2*b*e - 6*a*b**2*d + 3*b**3*c)/b**7 + f*x**13/(13*b**3) + x**10*(-3*a*f + b*e)/(10*b**4) + x**7*(6*a**2*f - 3*a*b*e + b**2*d)/(7*b**5) + x**4*(-10*a**3*f + 6*a**2*b*e - 3*a*b**2*d + b**3*c)/(4*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_287():
    f = x**10*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**3
    F = a**(sympy.S(2)/3)*(-119*a**3*f + 77*a**2*b*e - 44*a*b**2*d + 20*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*b**(sympy.S(20)/3)) - a**(sympy.S(2)/3)*(-119*a**3*f + 77*a**2*b*e - 44*a*b**2*d + 20*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*b**(sympy.S(20)/3)) + sqrt(3)*a**(sympy.S(2)/3)*(-119*a**3*f + 77*a**2*b*e - 44*a*b**2*d + 20*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*b**(sympy.S(20)/3)) - a**2*x**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*b**6*(a + b*x**3)**2) + a*x**2*(-16*a**3*f + 13*a**2*b*e - 10*a*b**2*d + 7*b**3*c)/(9*b**6*(a + b*x**3)) + f*x**11/(11*b**3) + x**8*(-3*a*f + b*e)/(8*b**4) + x**5*(6*a**2*f - 3*a*b*e + b**2*d)/(5*b**5) + x**2*(-10*a**3*f + 6*a**2*b*e - 3*a*b**2*d + b**3*c)/(2*b**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_288():
    f = x**9*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**3
    F = -a**(sympy.S(1)/3)*(-104*a**3*f + 65*a**2*b*e - 35*a*b**2*d + 14*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*b**(sympy.S(19)/3)) + a**(sympy.S(1)/3)*(-104*a**3*f + 65*a**2*b*e - 35*a*b**2*d + 14*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*b**(sympy.S(19)/3)) + sqrt(3)*a**(sympy.S(1)/3)*(-104*a**3*f + 65*a**2*b*e - 35*a*b**2*d + 14*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*b**(sympy.S(19)/3)) - a**2*x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*b**6*(a + b*x**3)**2) + a*x*(-31*a**3*f + 25*a**2*b*e - 19*a*b**2*d + 13*b**3*c)/(18*b**6*(a + b*x**3)) + f*x**10/(10*b**3) + x**7*(-3*a*f + b*e)/(7*b**4) + x**4*(6*a**2*f - 3*a*b*e + b**2*d)/(4*b**5) + x*(-10*a**3*f + 6*a**2*b*e - 3*a*b**2*d + b**3*c)/b**6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_289():
    f = x**7*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**3
    F = a*x**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*b**5*(a + b*x**3)**2) + f*x**8/(8*b**3) + x**5*(-3*a*f + b*e)/(5*b**4) + x**2*(6*a**2*f - 3*a*b*e + b**2*d)/(2*b**5) - x**2*(-13*a**3*f + 10*a**2*b*e - 7*a*b**2*d + 4*b**3*c)/(9*b**5*(a + b*x**3)) - (-77*a**3*f + 44*a**2*b*e - 20*a*b**2*d + 5*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(1)/3)*b**(sympy.S(17)/3)) + (-77*a**3*f + 44*a**2*b*e - 20*a*b**2*d + 5*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(1)/3)*b**(sympy.S(17)/3)) - sqrt(3)*(-77*a**3*f + 44*a**2*b*e - 20*a*b**2*d + 5*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(1)/3)*b**(sympy.S(17)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_290():
    f = x**6*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**3
    F = a*x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*b**5*(a + b*x**3)**2) + f*x**7/(7*b**3) + x**4*(-3*a*f + b*e)/(4*b**4) + x*(6*a**2*f - 3*a*b*e + b**2*d)/b**5 - x*(-25*a**3*f + 19*a**2*b*e - 13*a*b**2*d + 7*b**3*c)/(18*b**5*(a + b*x**3)) + (-65*a**3*f + 35*a**2*b*e - 14*a*b**2*d + 2*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(2)/3)*b**(sympy.S(16)/3)) - (-65*a**3*f + 35*a**2*b*e - 14*a*b**2*d + 2*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(2)/3)*b**(sympy.S(16)/3)) - sqrt(3)*(-65*a**3*f + 35*a**2*b*e - 14*a*b**2*d + 2*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(2)/3)*b**(sympy.S(16)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_291():
    f = x**4*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**3
    F = f*x**5/(5*b**3) + x**2*(-3*a*f + b*e)/(2*b**4) - x**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*b**4*(a + b*x**3)**2) + x**2*(-10*a**3*f + 7*a**2*b*e - 4*a*b**2*d + b**3*c)/(9*a*b**4*(a + b*x**3)) - (44*a**3*f - 20*a**2*b*e + 5*a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(4)/3)*b**(sympy.S(14)/3)) + (44*a**3*f - 20*a**2*b*e + 5*a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(4)/3)*b**(sympy.S(14)/3)) - sqrt(3)*(44*a**3*f - 20*a**2*b*e + 5*a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(4)/3)*b**(sympy.S(14)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_292():
    f = x**3*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**3
    F = f*x**4/(4*b**3) + x*(-3*a*f + b*e)/b**4 - x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*b**4*(a + b*x**3)**2) + x*(-19*a**3*f + 13*a**2*b*e - 7*a*b**2*d + b**3*c)/(18*a*b**4*(a + b*x**3)) + (35*a**3*f - 14*a**2*b*e + 2*a*b**2*d + b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(5)/3)*b**(sympy.S(13)/3)) - (35*a**3*f - 14*a**2*b*e + 2*a*b**2*d + b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(5)/3)*b**(sympy.S(13)/3)) - sqrt(3)*(35*a**3*f - 14*a**2*b*e + 2*a*b**2*d + b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(5)/3)*b**(sympy.S(13)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_293():
    f = x*(c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**3
    F = f*x**2/(2*b**3) + x**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*a*b**3*(a + b*x**3)**2) + x**2*(7*a**3*f - 4*a**2*b*e + a*b**2*d + 2*b**3*c)/(9*a**2*b**3*(a + b*x**3)) - (-20*a**3*f + 5*a**2*b*e + a*b**2*d + 2*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(7)/3)*b**(sympy.S(11)/3)) + (-20*a**3*f + 5*a**2*b*e + a*b**2*d + 2*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(7)/3)*b**(sympy.S(11)/3)) - sqrt(3)*(-20*a**3*f + 5*a**2*b*e + a*b**2*d + 2*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(7)/3)*b**(sympy.S(11)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_294():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(a + b*x**3)**3
    F = f*x/b**3 + x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*a*b**3*(a + b*x**3)**2) + x*(13*a**3*f - 7*a**2*b*e + a*b**2*d + 5*b**3*c)/(18*a**2*b**3*(a + b*x**3)) + (-14*a**3*f + 2*a**2*b*e + a*b**2*d + 5*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(8)/3)*b**(sympy.S(10)/3)) - (-14*a**3*f + 2*a**2*b*e + a*b**2*d + 5*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(8)/3)*b**(sympy.S(10)/3)) - sqrt(3)*(-14*a**3*f + 2*a**2*b*e + a*b**2*d + 5*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(8)/3)*b**(sympy.S(10)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_295():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**2*(a + b*x**3)**3)
    F = -x**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*a**2*b**2*(a + b*x**3)**2) - c/(a**3*x) - x**2*(4*a**3*f - a**2*b*e - 2*a*b**2*d + 5*b**3*c)/(9*a**3*b**2*(a + b*x**3)) + (-5*a**3*f - a**2*b*e - 2*a*b**2*d + 14*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(10)/3)*b**(sympy.S(8)/3)) - (-5*a**3*f - a**2*b*e - 2*a*b**2*d + 14*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(10)/3)*b**(sympy.S(8)/3)) + sqrt(3)*(-5*a**3*f - a**2*b*e - 2*a*b**2*d + 14*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(10)/3)*b**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_296():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**3*(a + b*x**3)**3)
    F = -x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*a**2*b**2*(a + b*x**3)**2) - c/(2*a**3*x**2) - x*(7*a**3*f - a**2*b*e - 5*a*b**2*d + 11*b**3*c)/(18*a**3*b**2*(a + b*x**3)) - (-2*a**3*f - a**2*b*e - 5*a*b**2*d + 20*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(11)/3)*b**(sympy.S(7)/3)) + (-2*a**3*f - a**2*b*e - 5*a*b**2*d + 20*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(11)/3)*b**(sympy.S(7)/3)) + sqrt(3)*(-2*a**3*f - a**2*b*e - 5*a*b**2*d + 20*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(11)/3)*b**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_297():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**5*(a + b*x**3)**3)
    F = -c/(4*a**3*x**4) + x**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*a**3*b*(a + b*x**3)**2) + (-a*d + 3*b*c)/(a**4*x) + x**2*(a**3*f + 2*a**2*b*e - 5*a*b**2*d + 8*b**3*c)/(9*a**4*b*(a + b*x**3)) - (a**3*f + 2*a**2*b*e - 14*a*b**2*d + 35*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(13)/3)*b**(sympy.S(5)/3)) + (a**3*f + 2*a**2*b*e - 14*a*b**2*d + 35*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(13)/3)*b**(sympy.S(5)/3)) - sqrt(3)*(a**3*f + 2*a**2*b*e - 14*a*b**2*d + 35*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(13)/3)*b**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_298():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**6*(a + b*x**3)**3)
    F = -c/(5*a**3*x**5) + x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*a**3*b*(a + b*x**3)**2) + (-a*d + 3*b*c)/(2*a**4*x**2) + x*(a**3*f + 5*a**2*b*e - 11*a*b**2*d + 17*b**3*c)/(18*a**4*b*(a + b*x**3)) + (a**3*f + 5*a**2*b*e - 20*a*b**2*d + 44*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(14)/3)*b**(sympy.S(4)/3)) - (a**3*f + 5*a**2*b*e - 20*a*b**2*d + 44*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(14)/3)*b**(sympy.S(4)/3)) - sqrt(3)*(a**3*f + 5*a**2*b*e - 20*a*b**2*d + 44*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(14)/3)*b**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_299():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**8*(a + b*x**3)**3)
    F = -c/(7*a**3*x**7) - x**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*a**4*(a + b*x**3)**2) + (-a*d + 3*b*c)/(4*a**4*x**4) - x**2*(-2*a**3*f + 5*a**2*b*e - 8*a*b**2*d + 11*b**3*c)/(9*a**5*(a + b*x**3)) - (a**2*e - 3*a*b*d + 6*b**2*c)/(a**5*x) + (-2*a**3*f + 14*a**2*b*e - 35*a*b**2*d + 65*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(16)/3)*b**(sympy.S(2)/3)) - (-2*a**3*f + 14*a**2*b*e - 35*a*b**2*d + 65*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(16)/3)*b**(sympy.S(2)/3)) + sqrt(3)*(-2*a**3*f + 14*a**2*b*e - 35*a*b**2*d + 65*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(16)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_300():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**9*(a + b*x**3)**3)
    F = -c/(8*a**3*x**8) - x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*a**4*(a + b*x**3)**2) + (-a*d + 3*b*c)/(5*a**4*x**5) - x*(-5*a**3*f + 11*a**2*b*e - 17*a*b**2*d + 23*b**3*c)/(18*a**5*(a + b*x**3)) - (a**2*e - 3*a*b*d + 6*b**2*c)/(2*a**5*x**2) - (-5*a**3*f + 20*a**2*b*e - 44*a*b**2*d + 77*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(17)/3)*b**(sympy.S(1)/3)) + (-5*a**3*f + 20*a**2*b*e - 44*a*b**2*d + 77*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(17)/3)*b**(sympy.S(1)/3)) + sqrt(3)*(-5*a**3*f + 20*a**2*b*e - 44*a*b**2*d + 77*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(17)/3)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_301():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**11*(a + b*x**3)**3)
    F = -c/(10*a**3*x**10) + (-a*d + 3*b*c)/(7*a**4*x**7) + b*x**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*a**5*(a + b*x**3)**2) - (a**2*e - 3*a*b*d + 6*b**2*c)/(4*a**5*x**4) + b*x**2*(-5*a**3*f + 8*a**2*b*e - 11*a*b**2*d + 14*b**3*c)/(9*a**6*(a + b*x**3)) + (-a**3*f + 3*a**2*b*e - 6*a*b**2*d + 10*b**3*c)/(a**6*x) - b**(sympy.S(1)/3)*(-14*a**3*f + 35*a**2*b*e - 65*a*b**2*d + 104*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(19)/3)) + b**(sympy.S(1)/3)*(-14*a**3*f + 35*a**2*b*e - 65*a*b**2*d + 104*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(19)/3)) - sqrt(3)*b**(sympy.S(1)/3)*(-14*a**3*f + 35*a**2*b*e - 65*a*b**2*d + 104*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(19)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_302():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**12*(a + b*x**3)**3)
    F = -c/(11*a**3*x**11) + (-a*d + 3*b*c)/(8*a**4*x**8) + b*x*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*a**5*(a + b*x**3)**2) - (a**2*e - 3*a*b*d + 6*b**2*c)/(5*a**5*x**5) + b*x*(-11*a**3*f + 17*a**2*b*e - 23*a*b**2*d + 29*b**3*c)/(18*a**6*(a + b*x**3)) + (-a**3*f + 3*a**2*b*e - 6*a*b**2*d + 10*b**3*c)/(2*a**6*x**2) + b**(sympy.S(2)/3)*(-20*a**3*f + 44*a**2*b*e - 77*a*b**2*d + 119*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(20)/3)) - b**(sympy.S(2)/3)*(-20*a**3*f + 44*a**2*b*e - 77*a*b**2*d + 119*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(20)/3)) - sqrt(3)*b**(sympy.S(2)/3)*(-20*a**3*f + 44*a**2*b*e - 77*a*b**2*d + 119*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(20)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_303():
    f = (c + d*x**3 + e*x**6 + f*x**9)/(x**14*(a + b*x**3)**3)
    F = -c/(13*a**3*x**13) + (-a*d + 3*b*c)/(10*a**4*x**10) - (a**2*e - 3*a*b*d + 6*b**2*c)/(7*a**5*x**7) - b**2*x**2*(-a**3*f + a**2*b*e - a*b**2*d + b**3*c)/(6*a**6*(a + b*x**3)**2) + (-a**3*f + 3*a**2*b*e - 6*a*b**2*d + 10*b**3*c)/(4*a**6*x**4) - b**2*x**2*(-8*a**3*f + 11*a**2*b*e - 14*a*b**2*d + 17*b**3*c)/(9*a**7*(a + b*x**3)) - b*(-3*a**3*f + 6*a**2*b*e - 10*a*b**2*d + 15*b**3*c)/(a**7*x) + b**(sympy.S(4)/3)*(-35*a**3*f + 65*a**2*b*e - 104*a*b**2*d + 152*b**3*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(22)/3)) - b**(sympy.S(4)/3)*(-35*a**3*f + 65*a**2*b*e - 104*a*b**2*d + 152*b**3*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(22)/3)) + sqrt(3)*b**(sympy.S(4)/3)*(-35*a**3*f + 65*a**2*b*e - 104*a*b**2*d + 152*b**3*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(22)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_304():
    f = x**4*(1 - x)/(x**3 + 1)
    F = -x**3/3 + x**2/2 + 2*log(x + 1)/3 + log(x**2 - x + 1)/6 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_305():
    f = x**3*(1 - x)/(x**3 + 1)
    F = -x**2/2 + x - 2*log(x + 1)/3 + log(x**2 - x + 1)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_306():
    f = x**2*(1 - x)/(x**3 + 1)
    F = -x + 2*log(x + 1)/3 + log(x**2 - x + 1)/6 - sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_307():
    f = x*(1 - x)/(x**3 + 1)
    F = -2*log(x + 1)/3 - log(x**2 - x + 1)/6 - sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_308():
    f = (1 - x)/(x*(x**3 + 1))
    F = log(x) - 2*log(x + 1)/3 - log(x**2 - x + 1)/6 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_309():
    f = (1 - x)/(x**2*(x**3 + 1))
    F = -log(x) + 2*log(x + 1)/3 + log(x**2 - x + 1)/6 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/3 - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_310():
    f = (1 - x)/(x**3*(x**3 + 1))
    F = -2*log(x + 1)/3 + log(x**2 - x + 1)/3 + 1/x - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_311():
    f = x*(2*x + 1)/(x**3 + 1)
    F = log(x + 1)/3 + 5*log(x**2 - x + 1)/6 - sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_312():
    f = x*(2*x + 1)/(1 - x**3)
    F = -log(1 - x) - log(x**2 + x + 1)/2 - sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_313():
    f = x**2*(a + b*x**3)*(c + d*x + e*x**2)
    F = a*c*x**3/3 + a*d*x**4/4 + a*e*x**5/5 + b*c*x**6/6 + b*d*x**7/7 + b*e*x**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_314():
    f = x*(a + b*x**3)*(c + d*x + e*x**2)
    F = a*c*x**2/2 + a*d*x**3/3 + a*e*x**4/4 + b*c*x**5/5 + b*d*x**6/6 + b*e*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_315():
    f = (a + b*x**3)*(c + d*x + e*x**2)
    F = a*c*x + a*d*x**2/2 + a*e*x**3/3 + b*c*x**4/4 + b*d*x**5/5 + b*e*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_316():
    f = (a + b*x**3)*(c + d*x + e*x**2)/x
    F = a*c*log(x) + a*d*x + a*e*x**2/2 + b*c*x**3/3 + b*d*x**4/4 + b*e*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_317():
    f = (a + b*x**3)*(c + d*x + e*x**2)/x**2
    F = -a*c/x + a*d*log(x) + a*e*x + b*c*x**2/2 + b*d*x**3/3 + b*e*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_318():
    f = (a + b*x**3)*(c + d*x + e*x**2)/x**3
    F = -a*c/(2*x**2) - a*d/x + a*e*log(x) + b*c*x + b*d*x**2/2 + b*e*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_319():
    f = x**2*(a + b*x**3)**2*(c + d*x + e*x**2)
    F = a**2*d*x**4/4 + a**2*e*x**5/5 + 2*a*b*d*x**7/7 + a*b*e*x**8/4 + b**2*d*x**10/10 + b**2*e*x**11/11 + c*(a + b*x**3)**3/(9*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_320():
    f = x*(a + b*x**3)**2*(c + d*x + e*x**2)
    F = a**2*c*x**2/2 + a**2*e*x**4/4 + 2*a*b*c*x**5/5 + 2*a*b*e*x**7/7 + b**2*c*x**8/8 + b**2*e*x**10/10 + d*(a + b*x**3)**3/(9*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_321():
    f = (a + b*x**3)**2*(c + d*x + e*x**2)
    F = a**2*c*x + a**2*d*x**2/2 + a*b*c*x**4/2 + 2*a*b*d*x**5/5 + b**2*c*x**7/7 + b**2*d*x**8/8 + e*(a + b*x**3)**3/(9*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_322():
    f = (a + b*x**3)**2*(c + d*x + e*x**2)/x
    F = a**2*c*log(x) + a**2*d*x + a**2*e*x**2/2 + 2*a*b*c*x**3/3 + a*b*d*x**4/2 + 2*a*b*e*x**5/5 + b**2*c*x**6/6 + b**2*d*x**7/7 + b**2*e*x**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_323():
    f = (a + b*x**3)**2*(c + d*x + e*x**2)/x**2
    F = -a**2*c/x + a**2*d*log(x) + a**2*e*x + a*b*c*x**2 + 2*a*b*d*x**3/3 + a*b*e*x**4/2 + b**2*c*x**5/5 + b**2*d*x**6/6 + b**2*e*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_324():
    f = (a + b*x**3)**2*(c + d*x + e*x**2)/x**3
    F = -a**2*c/(2*x**2) - a**2*d/x + a**2*e*log(x) + 2*a*b*c*x + a*b*d*x**2 + 2*a*b*e*x**3/3 + b**2*c*x**4/4 + b**2*d*x**5/5 + b**2*e*x**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_325():
    f = x**2*(a + b*x**3)**3*(c + d*x + e*x**2)
    F = a**3*d*x**4/4 + a**3*e*x**5/5 + 3*a**2*b*d*x**7/7 + 3*a**2*b*e*x**8/8 + 3*a*b**2*d*x**10/10 + 3*a*b**2*e*x**11/11 + b**3*d*x**13/13 + b**3*e*x**14/14 + c*(a + b*x**3)**4/(12*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_326():
    f = x*(a + b*x**3)**3*(c + d*x + e*x**2)
    F = a**3*c*x**2/2 + a**3*e*x**4/4 + 3*a**2*b*c*x**5/5 + 3*a**2*b*e*x**7/7 + 3*a*b**2*c*x**8/8 + 3*a*b**2*e*x**10/10 + b**3*c*x**11/11 + b**3*e*x**13/13 + d*(a + b*x**3)**4/(12*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_327():
    f = (a + b*x**3)**3*(c + d*x + e*x**2)
    F = a**3*c*x + a**3*d*x**2/2 + 3*a**2*b*c*x**4/4 + 3*a**2*b*d*x**5/5 + 3*a*b**2*c*x**7/7 + 3*a*b**2*d*x**8/8 + b**3*c*x**10/10 + b**3*d*x**11/11 + e*(a + b*x**3)**4/(12*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_328():
    f = (a + b*x**3)**3*(c + d*x + e*x**2)/x
    F = a**3*c*log(x) + a**3*d*x + a**3*e*x**2/2 + a**2*b*c*x**3 + 3*a**2*b*d*x**4/4 + 3*a**2*b*e*x**5/5 + a*b**2*c*x**6/2 + 3*a*b**2*d*x**7/7 + 3*a*b**2*e*x**8/8 + b**3*c*x**9/9 + b**3*d*x**10/10 + b**3*e*x**11/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_329():
    f = (a + b*x**3)**3*(c + d*x + e*x**2)/x**2
    F = -a**3*c/x + a**3*d*log(x) + a**3*e*x + 3*a**2*b*c*x**2/2 + a**2*b*d*x**3 + 3*a**2*b*e*x**4/4 + 3*a*b**2*c*x**5/5 + a*b**2*d*x**6/2 + 3*a*b**2*e*x**7/7 + b**3*c*x**8/8 + b**3*d*x**9/9 + b**3*e*x**10/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_330():
    f = (a + b*x**3)**3*(c + d*x + e*x**2)/x**3
    F = -a**3*c/(2*x**2) - a**3*d/x + a**3*e*log(x) + 3*a**2*b*c*x + 3*a**2*b*d*x**2/2 + a**2*b*e*x**3 + 3*a*b**2*c*x**4/4 + 3*a*b**2*d*x**5/5 + a*b**2*e*x**6/2 + b**3*c*x**7/7 + b**3*d*x**8/8 + b**3*e*x**9/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_331():
    f = x**2*(a + b*x**3)**4*(c + d*x + e*x**2)
    F = a**4*d*x**4/4 + a**4*e*x**5/5 + 4*a**3*b*d*x**7/7 + a**3*b*e*x**8/2 + 3*a**2*b**2*d*x**10/5 + 6*a**2*b**2*e*x**11/11 + 4*a*b**3*d*x**13/13 + 2*a*b**3*e*x**14/7 + b**4*d*x**16/16 + b**4*e*x**17/17 + c*(a + b*x**3)**5/(15*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_332():
    f = x*(a + b*x**3)**4*(c + d*x + e*x**2)
    F = a**4*c*x**2/2 + a**4*e*x**4/4 + 4*a**3*b*c*x**5/5 + 4*a**3*b*e*x**7/7 + 3*a**2*b**2*c*x**8/4 + 3*a**2*b**2*e*x**10/5 + 4*a*b**3*c*x**11/11 + 4*a*b**3*e*x**13/13 + b**4*c*x**14/14 + b**4*e*x**16/16 + d*(a + b*x**3)**5/(15*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_333():
    f = (a + b*x**3)**4*(c + d*x + e*x**2)
    F = a**4*c*x + a**4*d*x**2/2 + a**3*b*c*x**4 + 4*a**3*b*d*x**5/5 + 6*a**2*b**2*c*x**7/7 + 3*a**2*b**2*d*x**8/4 + 2*a*b**3*c*x**10/5 + 4*a*b**3*d*x**11/11 + b**4*c*x**13/13 + b**4*d*x**14/14 + e*(a + b*x**3)**5/(15*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_334():
    f = (a + b*x**3)**4*(c + d*x + e*x**2)/x
    F = a**4*c*log(x) + a**4*d*x + a**4*e*x**2/2 + 4*a**3*b*c*x**3/3 + a**3*b*d*x**4 + 4*a**3*b*e*x**5/5 + a**2*b**2*c*x**6 + 6*a**2*b**2*d*x**7/7 + 3*a**2*b**2*e*x**8/4 + 4*a*b**3*c*x**9/9 + 2*a*b**3*d*x**10/5 + 4*a*b**3*e*x**11/11 + b**4*c*x**12/12 + b**4*d*x**13/13 + b**4*e*x**14/14
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_335():
    f = (a + b*x**3)**4*(c + d*x + e*x**2)/x**2
    F = -a**4*c/x + a**4*d*log(x) + a**4*e*x + 2*a**3*b*c*x**2 + 4*a**3*b*d*x**3/3 + a**3*b*e*x**4 + 6*a**2*b**2*c*x**5/5 + a**2*b**2*d*x**6 + 6*a**2*b**2*e*x**7/7 + a*b**3*c*x**8/2 + 4*a*b**3*d*x**9/9 + 2*a*b**3*e*x**10/5 + b**4*c*x**11/11 + b**4*d*x**12/12 + b**4*e*x**13/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_336():
    f = (a + b*x**3)**4*(c + d*x + e*x**2)/x**3
    F = -a**4*c/(2*x**2) - a**4*d/x + a**4*e*log(x) + 4*a**3*b*c*x + 2*a**3*b*d*x**2 + 4*a**3*b*e*x**3/3 + 3*a**2*b**2*c*x**4/2 + 6*a**2*b**2*d*x**5/5 + a**2*b**2*e*x**6 + 4*a*b**3*c*x**7/7 + a*b**3*d*x**8/2 + 4*a*b**3*e*x**9/9 + b**4*c*x**10/10 + b**4*d*x**11/11 + b**4*e*x**12/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_337():
    f = x**3*(c + d*x + e*x**2)/(a + b*x**3)
    F = a**(sympy.S(1)/3)*(-a**(sympy.S(1)/3)*d/b**(sympy.S(1)/3) + c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*b**(sympy.S(4)/3)) - a**(sympy.S(1)/3)*(-a**(sympy.S(1)/3)*d + b**(sympy.S(1)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*b**(sympy.S(5)/3)) + sqrt(3)*a**(sympy.S(1)/3)*(a**(sympy.S(1)/3)*d + b**(sympy.S(1)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(5)/3)) - a*e*log(a + b*x**3)/(3*b**2) + c*x/b + d*x**2/(2*b) + e*x**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_338():
    f = x**2*(c + d*x + e*x**2)/(a + b*x**3)
    F = a**(sympy.S(1)/3)*(-a**(sympy.S(1)/3)*e/b**(sympy.S(1)/3) + d)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*b**(sympy.S(4)/3)) - a**(sympy.S(1)/3)*(-a**(sympy.S(1)/3)*e + b**(sympy.S(1)/3)*d)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*b**(sympy.S(5)/3)) + sqrt(3)*a**(sympy.S(1)/3)*(a**(sympy.S(1)/3)*e + b**(sympy.S(1)/3)*d)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(5)/3)) + c*log(a + b*x**3)/(3*b) + d*x/b + e*x**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_339():
    f = x*(c + d*x + e*x**2)/(a + b*x**3)
    F = d*log(a + b*x**3)/(3*b) + e*x/b - sqrt(3)*(-a**(sympy.S(2)/3)*e + b**(sympy.S(2)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3)*b**(sympy.S(4)/3)) - (a**(sympy.S(2)/3)*e + b**(sympy.S(2)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)*b**(sympy.S(4)/3)) + (a**(sympy.S(2)/3)*e + b**(sympy.S(2)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(1)/3)*b**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_340():
    f = (c + d*x + e*x**2)/(a + b*x**3)
    F = e*log(a + b*x**3)/(3*b) - (-a**(sympy.S(1)/3)*d/b**(sympy.S(1)/3) + c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)) + (-a**(sympy.S(1)/3)*d + b**(sympy.S(1)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(2)/3)*b**(sympy.S(2)/3)) - sqrt(3)*(a**(sympy.S(1)/3)*d + b**(sympy.S(1)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_341():
    f = (c + d*x + e*x**2)/(x*(a + b*x**3))
    F = c*log(x)/a - c*log(a + b*x**3)/(3*a) - (-a**(sympy.S(1)/3)*e/b**(sympy.S(1)/3) + d)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)) + (-a**(sympy.S(1)/3)*e + b**(sympy.S(1)/3)*d)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(2)/3)*b**(sympy.S(2)/3)) - sqrt(3)*(a**(sympy.S(1)/3)*e + b**(sympy.S(1)/3)*d)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_342():
    f = (c + d*x + e*x**2)/(x**2*(a + b*x**3))
    F = -c/(a*x) + d*log(x)/a - d*log(a + b*x**3)/(3*a) + sqrt(3)*(-a**(sympy.S(2)/3)*e + b**(sympy.S(2)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(4)/3)*b**(sympy.S(1)/3)) + (a**(sympy.S(2)/3)*e + b**(sympy.S(2)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(4)/3)*b**(sympy.S(1)/3)) - (a**(sympy.S(2)/3)*e + b**(sympy.S(2)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(4)/3)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_343():
    f = (c + d*x + e*x**2)/(x**3*(a + b*x**3))
    F = -c/(2*a*x**2) - d/(a*x) + e*log(x)/a - e*log(a + b*x**3)/(3*a) + b**(sympy.S(2)/3)*(-a**(sympy.S(1)/3)*d/b**(sympy.S(1)/3) + c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(5)/3)) - b**(sympy.S(1)/3)*(-a**(sympy.S(1)/3)*d + b**(sympy.S(1)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(5)/3)) + sqrt(3)*b**(sympy.S(1)/3)*(a**(sympy.S(1)/3)*d + b**(sympy.S(1)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_344():
    f = x**2*(c + d*x + e*x**2)/(a + b*x**3)**2
    F = -(c + d*x + e*x**2)/(3*b*(a + b*x**3)) - (-2*a**(sympy.S(1)/3)*e/b**(sympy.S(1)/3) + d)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(2)/3)*b**(sympy.S(4)/3)) + (-2*a**(sympy.S(1)/3)*e + b**(sympy.S(1)/3)*d)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(2)/3)*b**(sympy.S(5)/3)) - sqrt(3)*(2*a**(sympy.S(1)/3)*e + b**(sympy.S(1)/3)*d)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(2)/3)*b**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_345():
    f = x*(c + d*x + e*x**2)/(a + b*x**3)**2
    F = -x*(a*e - b*c*x - b*d*x**2)/(3*a*b*(a + b*x**3)) - (-a**(sympy.S(2)/3)*e + b**(sympy.S(2)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(4)/3)*b**(sympy.S(4)/3)) + (-a**(sympy.S(2)/3)*e + b**(sympy.S(2)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(4)/3)*b**(sympy.S(4)/3)) - sqrt(3)*(a**(sympy.S(2)/3)*e + b**(sympy.S(2)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(4)/3)*b**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_346():
    f = (c + d*x + e*x**2)/(a + b*x**3)**2
    F = -(a*e - b*x*(c + d*x))/(3*a*b*(a + b*x**3)) + (-a**(sympy.S(1)/3)*d + 2*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(5)/3)*b**(sympy.S(2)/3)) - (-a**(sympy.S(1)/3)*d + 2*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(5)/3)*b**(sympy.S(2)/3)) - sqrt(3)*(a**(sympy.S(1)/3)*d + 2*b**(sympy.S(1)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_347():
    f = (c + d*x + e*x**2)/(x*(a + b*x**3)**2)
    F = c*log(x)/a**2 - c*log(a + b*x**3)/(3*a**2) + x*(a*d + a*e*x - b*c*x**2)/(3*a**2*(a + b*x**3)) + (-a**(sympy.S(1)/3)*e + 2*b**(sympy.S(1)/3)*d)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(5)/3)*b**(sympy.S(2)/3)) - (-a**(sympy.S(1)/3)*e + 2*b**(sympy.S(1)/3)*d)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(5)/3)*b**(sympy.S(2)/3)) - sqrt(3)*(a**(sympy.S(1)/3)*e + 2*b**(sympy.S(1)/3)*d)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_348():
    f = (c + d*x + e*x**2)/(x**2*(a + b*x**3)**2)
    F = -c/(a**2*x) + d*log(x)/a**2 - d*log(a + b*x**3)/(3*a**2) + x*(a*e - b*c*x - b*d*x**2)/(3*a**2*(a + b*x**3)) + sqrt(3)*(-2*a**(sympy.S(2)/3)*e + 4*b**(sympy.S(2)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(7)/3)*b**(sympy.S(1)/3)) - (a**(sympy.S(2)/3)*e + 2*b**(sympy.S(2)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(9*a**(sympy.S(7)/3)*b**(sympy.S(1)/3)) + (2*a**(sympy.S(2)/3)*e + 4*b**(sympy.S(2)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(7)/3)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_349():
    f = (c + d*x + e*x**2)/(x**3*(a + b*x**3)**2)
    F = -c/(2*a**2*x**2) - d/(a**2*x) + e*log(x)/a**2 - e*log(a + b*x**3)/(3*a**2) - x*(b*c + b*d*x + b*e*x**2)/(3*a**2*(a + b*x**3)) - b**(sympy.S(1)/3)*(-4*a**(sympy.S(1)/3)*d + 5*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(8)/3)) + b**(sympy.S(1)/3)*(-4*a**(sympy.S(1)/3)*d + 5*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(8)/3)) + sqrt(3)*b**(sympy.S(1)/3)*(4*a**(sympy.S(1)/3)*d + 5*b**(sympy.S(1)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_350():
    f = (c + d*x + e*x**2)/(x**4*(a + b*x**3)**2)
    F = -c/(3*a**2*x**3) - d/(2*a**2*x**2) - e/(a**2*x) - x*(b*d + b*e*x - b**2*c*x**2/a)/(3*a**2*(a + b*x**3)) - 2*b*c*log(x)/a**3 + 2*b*c*log(a + b*x**3)/(3*a**3) - b**(sympy.S(1)/3)*(-4*a**(sympy.S(1)/3)*e + 5*b**(sympy.S(1)/3)*d)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(8)/3)) + b**(sympy.S(1)/3)*(-4*a**(sympy.S(1)/3)*e + 5*b**(sympy.S(1)/3)*d)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(8)/3)) + sqrt(3)*b**(sympy.S(1)/3)*(4*a**(sympy.S(1)/3)*e + 5*b**(sympy.S(1)/3)*d)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_351():
    f = x**2*(c + d*x + e*x**2)/(a + b*x**3)**3
    F = -(c + d*x + e*x**2)/(6*b*(a + b*x**3)**2) + x*(d + 2*e*x)/(18*a*b*(a + b*x**3)) - (-a**(sympy.S(1)/3)*e/b**(sympy.S(1)/3) + d)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(5)/3)*b**(sympy.S(4)/3)) + (-a**(sympy.S(1)/3)*e + b**(sympy.S(1)/3)*d)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(5)/3)*b**(sympy.S(5)/3)) - sqrt(3)*(a**(sympy.S(1)/3)*e + b**(sympy.S(1)/3)*d)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(5)/3)*b**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_352():
    f = x*(c + d*x + e*x**2)/(a + b*x**3)**3
    F = -x*(a*e - b*c*x - b*d*x**2)/(6*a*b*(a + b*x**3)**2) - (3*a*d - x*(a*e + 4*b*c*x))/(18*a**2*b*(a + b*x**3)) - (-a**(sympy.S(2)/3)*e + 2*b**(sympy.S(2)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(7)/3)*b**(sympy.S(4)/3)) + (-a**(sympy.S(2)/3)*e + 2*b**(sympy.S(2)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(7)/3)*b**(sympy.S(4)/3)) - sqrt(3)*(a**(sympy.S(2)/3)*e + 2*b**(sympy.S(2)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(7)/3)*b**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_353():
    f = (c + d*x + e*x**2)/(a + b*x**3)**3
    F = -(a*e - b*x*(c + d*x))/(6*a*b*(a + b*x**3)**2) + x*(5*c + 4*d*x)/(18*a**2*(a + b*x**3)) + (-2*a**(sympy.S(1)/3)*d + 5*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(8)/3)*b**(sympy.S(2)/3)) - (-2*a**(sympy.S(1)/3)*d + 5*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(8)/3)*b**(sympy.S(2)/3)) - sqrt(3)*(2*a**(sympy.S(1)/3)*d + 5*b**(sympy.S(1)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(8)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_354():
    f = (c + d*x + e*x**2)/(x*(a + b*x**3)**3)
    F = x*(a*d + a*e*x - b*c*x**2)/(6*a**2*(a + b*x**3)**2) + c*log(x)/a**3 - c*log(a + b*x**3)/(3*a**3) + x*(5*a*d + 4*a*e*x - 9*b*c*x**2)/(18*a**3*(a + b*x**3)) + (-2*a**(sympy.S(1)/3)*e + 5*b**(sympy.S(1)/3)*d)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(8)/3)*b**(sympy.S(2)/3)) - (-2*a**(sympy.S(1)/3)*e + 5*b**(sympy.S(1)/3)*d)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(8)/3)*b**(sympy.S(2)/3)) - sqrt(3)*(2*a**(sympy.S(1)/3)*e + 5*b**(sympy.S(1)/3)*d)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(8)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_355():
    f = (c + d*x + e*x**2)/(x**2*(a + b*x**3)**3)
    F = x*(a*e - b*c*x - b*d*x**2)/(6*a**2*(a + b*x**3)**2) - c/(a**3*x) + d*log(x)/a**3 - d*log(a + b*x**3)/(3*a**3) + x*(5*a*e - 10*b*c*x - 9*b*d*x**2)/(18*a**3*(a + b*x**3)) + sqrt(3)*(-5*a**(sympy.S(2)/3)*e + 14*b**(sympy.S(2)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(10)/3)*b**(sympy.S(1)/3)) + (5*a**(sympy.S(2)/3)*e + 14*b**(sympy.S(2)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(10)/3)*b**(sympy.S(1)/3)) - (5*a**(sympy.S(2)/3)*e + 14*b**(sympy.S(2)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(10)/3)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_356():
    f = (c + d*x + e*x**2)/(x**3*(a + b*x**3)**3)
    F = -x*(b*c + b*d*x + b*e*x**2)/(6*a**2*(a + b*x**3)**2) - c/(2*a**3*x**2) - d/(a**3*x) + e*log(x)/a**3 - e*log(a + b*x**3)/(3*a**3) - x*(11*b*c + 10*b*d*x + 9*b*e*x**2)/(18*a**3*(a + b*x**3)) - 2*b**(sympy.S(1)/3)*(-7*a**(sympy.S(1)/3)*d + 10*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(11)/3)) + b**(sympy.S(1)/3)*(-7*a**(sympy.S(1)/3)*d + 10*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(27*a**(sympy.S(11)/3)) + 2*sqrt(3)*b**(sympy.S(1)/3)*(7*a**(sympy.S(1)/3)*d + 10*b**(sympy.S(1)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(11)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_357():
    f = (c + d*x + e*x**2)/(x**4*(a + b*x**3)**3)
    F = -x*(b*d + b*e*x - b**2*c*x**2/a)/(6*a**2*(a + b*x**3)**2) - c/(3*a**3*x**3) - d/(2*a**3*x**2) - e/(a**3*x) - x*(11*b*d + 10*b*e*x - 15*b**2*c*x**2/a)/(18*a**3*(a + b*x**3)) - 3*b*c*log(x)/a**4 + b*c*log(a + b*x**3)/a**4 - 2*b**(sympy.S(1)/3)*(-7*a**(sympy.S(1)/3)*e + 10*b**(sympy.S(1)/3)*d)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(11)/3)) + b**(sympy.S(1)/3)*(-7*a**(sympy.S(1)/3)*e + 10*b**(sympy.S(1)/3)*d)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(27*a**(sympy.S(11)/3)) + 2*sqrt(3)*b**(sympy.S(1)/3)*(7*a**(sympy.S(1)/3)*e + 10*b**(sympy.S(1)/3)*d)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(11)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_358():
    f = x**2*(c + d*x + e*x**2)/(a + b*x**3)**4
    F = -(c + d*x + e*x**2)/(9*b*(a + b*x**3)**3) + x*(d + 2*e*x)/(54*a*b*(a + b*x**3)**2) + x*(5*d + 8*e*x)/(162*a**2*b*(a + b*x**3)) + (-4*a**(sympy.S(1)/3)*e + 5*b**(sympy.S(1)/3)*d)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(243*a**(sympy.S(8)/3)*b**(sympy.S(5)/3)) - (-4*a**(sympy.S(1)/3)*e + 5*b**(sympy.S(1)/3)*d)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(486*a**(sympy.S(8)/3)*b**(sympy.S(5)/3)) - sqrt(3)*(4*a**(sympy.S(1)/3)*e + 5*b**(sympy.S(1)/3)*d)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(243*a**(sympy.S(8)/3)*b**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_359():
    f = x*(c + d*x + e*x**2)/(a + b*x**3)**4
    F = -x*(a*e - b*c*x - b*d*x**2)/(9*a*b*(a + b*x**3)**3) - (6*a*d - x*(a*e + 7*b*c*x))/(54*a**2*b*(a + b*x**3)**2) + x*(5*a*e + 28*b*c*x)/(162*a**3*b*(a + b*x**3)) - (-5*a**(sympy.S(2)/3)*e + 14*b**(sympy.S(2)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(243*a**(sympy.S(10)/3)*b**(sympy.S(4)/3)) + (-5*a**(sympy.S(2)/3)*e + 14*b**(sympy.S(2)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(486*a**(sympy.S(10)/3)*b**(sympy.S(4)/3)) - sqrt(3)*(5*a**(sympy.S(2)/3)*e + 14*b**(sympy.S(2)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(243*a**(sympy.S(10)/3)*b**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_360():
    f = (c + d*x + e*x**2)/(a + b*x**3)**4
    F = -(a*e - b*x*(c + d*x))/(9*a*b*(a + b*x**3)**3) + x*(8*c + 7*d*x)/(54*a**2*(a + b*x**3)**2) + 2*x*(10*c + 7*d*x)/(81*a**3*(a + b*x**3)) + (-14*a**(sympy.S(1)/3)*d + 40*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(243*a**(sympy.S(11)/3)*b**(sympy.S(2)/3)) - (-7*a**(sympy.S(1)/3)*d + 20*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(243*a**(sympy.S(11)/3)*b**(sympy.S(2)/3)) - sqrt(3)*(14*a**(sympy.S(1)/3)*d + 40*b**(sympy.S(1)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(243*a**(sympy.S(11)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_361():
    f = (c + d*x + e*x**2)/(x*(a + b*x**3)**4)
    F = x*(a*d + a*e*x - b*c*x**2)/(9*a**2*(a + b*x**3)**3) + x*(8*a*d + 7*a*e*x - 15*b*c*x**2)/(54*a**3*(a + b*x**3)**2) + c*log(x)/a**4 - c*log(a + b*x**3)/(3*a**4) + x*(40*a*d + 28*a*e*x - 99*b*c*x**2)/(162*a**4*(a + b*x**3)) + (-14*a**(sympy.S(1)/3)*e + 40*b**(sympy.S(1)/3)*d)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(243*a**(sympy.S(11)/3)*b**(sympy.S(2)/3)) - (-7*a**(sympy.S(1)/3)*e + 20*b**(sympy.S(1)/3)*d)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(243*a**(sympy.S(11)/3)*b**(sympy.S(2)/3)) - sqrt(3)*(14*a**(sympy.S(1)/3)*e + 40*b**(sympy.S(1)/3)*d)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(243*a**(sympy.S(11)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_362():
    f = (c + d*x + e*x**2)/(x**2*(a + b*x**3)**4)
    F = x*(a*e - b*c*x - b*d*x**2)/(9*a**2*(a + b*x**3)**3) + x*(8*a*e - 16*b*c*x - 15*b*d*x**2)/(54*a**3*(a + b*x**3)**2) - c/(a**4*x) + d*log(x)/a**4 - d*log(a + b*x**3)/(3*a**4) + x*(40*a*e - 118*b*c*x - 99*b*d*x**2)/(162*a**4*(a + b*x**3)) + sqrt(3)*(-40*a**(sympy.S(2)/3)*e + 140*b**(sympy.S(2)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(243*a**(sympy.S(13)/3)*b**(sympy.S(1)/3)) - (20*a**(sympy.S(2)/3)*e + 70*b**(sympy.S(2)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(243*a**(sympy.S(13)/3)*b**(sympy.S(1)/3)) + (40*a**(sympy.S(2)/3)*e + 140*b**(sympy.S(2)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(243*a**(sympy.S(13)/3)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_363():
    f = (c + d*x + e*x**2)/(x**3*(a + b*x**3)**4)
    F = -x*(b*c + b*d*x + b*e*x**2)/(9*a**2*(a + b*x**3)**3) - x*(17*b*c + 16*b*d*x + 15*b*e*x**2)/(54*a**3*(a + b*x**3)**2) - c/(2*a**4*x**2) - d/(a**4*x) + e*log(x)/a**4 - e*log(a + b*x**3)/(3*a**4) - x*(139*b*c + 118*b*d*x + 99*b*e*x**2)/(162*a**4*(a + b*x**3)) - 20*b**(sympy.S(1)/3)*(-7*a**(sympy.S(1)/3)*d + 11*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(243*a**(sympy.S(14)/3)) + 10*b**(sympy.S(1)/3)*(-7*a**(sympy.S(1)/3)*d + 11*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(243*a**(sympy.S(14)/3)) + 20*sqrt(3)*b**(sympy.S(1)/3)*(7*a**(sympy.S(1)/3)*d + 11*b**(sympy.S(1)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(243*a**(sympy.S(14)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_364():
    f = (c + d*x + e*x**2)/(x**4*(a + b*x**3)**4)
    F = -x*(b*d + b*e*x - b**2*c*x**2/a)/(9*a**2*(a + b*x**3)**3) - x*(17*b*d + 16*b*e*x - 24*b**2*c*x**2/a)/(54*a**3*(a + b*x**3)**2) - c/(3*a**4*x**3) - d/(2*a**4*x**2) - e/(a**4*x) - x*(139*b*d + 118*b*e*x - 234*b**2*c*x**2/a)/(162*a**4*(a + b*x**3)) - 4*b*c*log(x)/a**5 + 4*b*c*log(a + b*x**3)/(3*a**5) - 20*b**(sympy.S(1)/3)*(-7*a**(sympy.S(1)/3)*e + 11*b**(sympy.S(1)/3)*d)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(243*a**(sympy.S(14)/3)) + 10*b**(sympy.S(1)/3)*(-7*a**(sympy.S(1)/3)*e + 11*b**(sympy.S(1)/3)*d)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(243*a**(sympy.S(14)/3)) + 20*sqrt(3)*b**(sympy.S(1)/3)*(7*a**(sympy.S(1)/3)*e + 11*b**(sympy.S(1)/3)*d)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(243*a**(sympy.S(14)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_365():
    f = (2*a*x - x**2)/(a**3 + x**3)
    F = -log(a + x) - 2*sqrt(3)*atan(sqrt(3)*(a - 2*x)/(3*a))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_366():
    f = x*(2*a - x)/(a**3 + x**3)
    F = -log(a + x) - 2*sqrt(3)*atan(sqrt(3)*(a - 2*x)/(3*a))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_367():
    f = (2*a*x + x**2)/(a**3 - x**3)
    F = -log(a - x) - 2*sqrt(3)*atan(sqrt(3)*(a + 2*x)/(3*a))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_368():
    f = x*(2*a + x)/(a**3 - x**3)
    F = -log(a - x) - 2*sqrt(3)*atan(sqrt(3)*(a + 2*x)/(3*a))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_369():
    f = x*(C*x - 2*C*(a/b)**(sympy.S(1)/3))/(a + b*x**3)
    F = C*log(x + (a/b)**(sympy.S(1)/3))/b + 2*sqrt(3)*C*atan(sqrt(3)*(-2*x/(a/b)**(sympy.S(1)/3) + 1)/3)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_370():
    f = x*(C*x - 2*C*(-a/b)**(sympy.S(1)/3))/(a - b*x**3)
    F = -C*log(x + (-a/b)**(sympy.S(1)/3))/b - 2*sqrt(3)*C*atan(sqrt(3)*(-2*x/(-a/b)**(sympy.S(1)/3) + 1)/3)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_371():
    f = x*(C*x + 2*C*(-a/b)**(sympy.S(1)/3))/(a + b*x**3)
    F = C*log(-x + (-a/b)**(sympy.S(1)/3))/b + 2*sqrt(3)*C*atan(sqrt(3)*(2*x/(-a/b)**(sympy.S(1)/3) + 1)/3)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_372():
    f = x*(C*x + 2*C*(a/b)**(sympy.S(1)/3))/(a - b*x**3)
    F = -C*log(-x + (a/b)**(sympy.S(1)/3))/b - 2*sqrt(3)*C*atan(sqrt(3)*(2*x/(a/b)**(sympy.S(1)/3) + 1)/3)/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_373():
    f = x**4*(a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)
    F = a*c*x**5/5 + a*d*x**6/6 + a*e*x**7/7 + b*f*x**11/11 + b*g*x**12/12 + b*h*x**13/13 + x**10*(a*h + b*e)/10 + x**9*(a*g + b*d)/9 + x**8*(a*f + b*c)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_374():
    f = x**3*(a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)
    F = a*c*x**4/4 + a*d*x**5/5 + a*e*x**6/6 + b*f*x**10/10 + b*g*x**11/11 + b*h*x**12/12 + x**9*(a*h + b*e)/9 + x**8*(a*g + b*d)/8 + x**7*(a*f + b*c)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_375():
    f = x**2*(a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)
    F = a*c*x**3/3 + a*d*x**4/4 + a*e*x**5/5 + b*f*x**9/9 + b*g*x**10/10 + b*h*x**11/11 + x**8*(a*h + b*e)/8 + x**7*(a*g + b*d)/7 + x**6*(a*f + b*c)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_376():
    f = x*(a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)
    F = a*c*x**2/2 + a*d*x**3/3 + a*e*x**4/4 + b*f*x**8/8 + b*g*x**9/9 + b*h*x**10/10 + x**7*(a*h + b*e)/7 + x**6*(a*g + b*d)/6 + x**5*(a*f + b*c)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_377():
    f = (a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)
    F = a*c*x + a*d*x**2/2 + a*e*x**3/3 + b*f*x**7/7 + b*g*x**8/8 + b*h*x**9/9 + x**6*(a*h + b*e)/6 + x**5*(a*g + b*d)/5 + x**4*(a*f + b*c)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_378():
    f = (a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/x
    F = a*c*log(x) + a*d*x + a*e*x**2/2 + b*f*x**6/6 + b*g*x**7/7 + b*h*x**8/8 + x**5*(a*h + b*e)/5 + x**4*(a*g + b*d)/4 + x**3*(a*f + b*c)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_379():
    f = (a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/x**2
    F = -a*c/x + a*d*log(x) + a*e*x + b*f*x**5/5 + b*g*x**6/6 + b*h*x**7/7 + x**4*(a*h + b*e)/4 + x**3*(a*g + b*d)/3 + x**2*(a*f + b*c)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_380():
    f = (a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/x**3
    F = -a*c/(2*x**2) - a*d/x + a*e*log(x) + b*f*x**4/4 + b*g*x**5/5 + b*h*x**6/6 + x**3*(a*h + b*e)/3 + x**2*(a*g + b*d)/2 + x*(a*f + b*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_381():
    f = (a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/x**4
    F = -a*c/(3*x**3) - a*d/(2*x**2) - a*e/x + b*f*x**3/3 + b*g*x**4/4 + b*h*x**5/5 + x**2*(a*h + b*e)/2 + x*(a*g + b*d) + (a*f + b*c)*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_382():
    f = (a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/x**5
    F = -a*c/(4*x**4) - a*d/(3*x**3) - a*e/(2*x**2) + b*f*x**2/2 + b*g*x**3/3 + b*h*x**4/4 + x*(a*h + b*e) + (a*g + b*d)*log(x) - (a*f + b*c)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_383():
    f = x**4*(a + b*x**3)**2*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)
    F = a**2*c*x**5/5 + a**2*d*x**6/6 + a**2*e*x**7/7 + a*x**10*(a*h + 2*b*e)/10 + a*x**9*(a*g + 2*b*d)/9 + a*x**8*(a*f + 2*b*c)/8 + b**2*f*x**14/14 + b**2*g*x**15/15 + b**2*h*x**16/16 + b*x**13*(2*a*h + b*e)/13 + b*x**12*(2*a*g + b*d)/12 + b*x**11*(2*a*f + b*c)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_384():
    f = x**3*(a + b*x**3)**2*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)
    F = a**2*c*x**4/4 + a**2*d*x**5/5 + a**2*e*x**6/6 + a*x**9*(a*h + 2*b*e)/9 + a*x**8*(a*g + 2*b*d)/8 + a*x**7*(a*f + 2*b*c)/7 + b**2*f*x**13/13 + b**2*g*x**14/14 + b**2*h*x**15/15 + b*x**12*(2*a*h + b*e)/12 + b*x**11*(2*a*g + b*d)/11 + b*x**10*(2*a*f + b*c)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_385():
    f = x**2*(a + b*x**3)**2*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)
    F = a**2*d*x**4/4 + a**2*e*x**5/5 + a**2*f*x**6/6 + 2*a*b*f*x**9/9 + a*x**8*(a*h + 2*b*e)/8 + a*x**7*(a*g + 2*b*d)/7 + b**2*f*x**12/12 + b**2*g*x**13/13 + b**2*h*x**14/14 + b*x**11*(2*a*h + b*e)/11 + b*x**10*(2*a*g + b*d)/10 + c*(a + b*x**3)**3/(9*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_386():
    f = x*(a + b*x**3)**2*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)
    F = a**2*c*x**2/2 + a**2*e*x**4/4 + a**2*g*x**6/6 + 2*a*b*g*x**9/9 + a*x**7*(a*h + 2*b*e)/7 + a*x**5*(a*f + 2*b*c)/5 + b**2*f*x**11/11 + b**2*g*x**12/12 + b**2*h*x**13/13 + b*x**10*(2*a*h + b*e)/10 + b*x**8*(2*a*f + b*c)/8 + d*(a + b*x**3)**3/(9*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_387():
    f = (a + b*x**3)**2*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)
    F = a**2*c*x + a**2*d*x**2/2 + a**2*h*x**6/6 + 2*a*b*h*x**9/9 + a*x**5*(a*g + 2*b*d)/5 + a*x**4*(a*f + 2*b*c)/4 + b**2*f*x**10/10 + b**2*g*x**11/11 + b**2*h*x**12/12 + b*x**8*(2*a*g + b*d)/8 + b*x**7*(2*a*f + b*c)/7 + e*(a + b*x**3)**3/(9*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_388():
    f = (a + b*x**3)**2*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/x
    F = a**2*c*log(x) + a**2*d*x + a**2*e*x**2/2 + 2*a*b*c*x**3/3 + a*x**5*(a*h + 2*b*e)/5 + a*x**4*(a*g + 2*b*d)/4 + b**2*c*x**6/6 + b**2*g*x**10/10 + b**2*h*x**11/11 + b*x**8*(2*a*h + b*e)/8 + b*x**7*(2*a*g + b*d)/7 + f*(a + b*x**3)**3/(9*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_389():
    f = (a + b*x**3)**2*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/x**2
    F = -a**2*c/x + a**2*d*log(x) + a**2*e*x + 2*a*b*d*x**3/3 + a*x**4*(a*h + 2*b*e)/4 + a*x**2*(a*f + 2*b*c)/2 + b**2*d*x**6/6 + b**2*f*x**8/8 + b**2*h*x**10/10 + b*x**7*(2*a*h + b*e)/7 + b*x**5*(2*a*f + b*c)/5 + g*(a + b*x**3)**3/(9*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_390():
    f = (a + b*x**3)**2*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/x**3
    F = -a**2*c/(2*x**2) - a**2*d/x + a**2*e*log(x) + 2*a*b*e*x**3/3 + a*x**2*(a*g + 2*b*d)/2 + a*x*(a*f + 2*b*c) + b**2*e*x**6/6 + b**2*f*x**7/7 + b**2*g*x**8/8 + b*x**5*(2*a*g + b*d)/5 + b*x**4*(2*a*f + b*c)/4 + h*(a + b*x**3)**3/(9*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_391():
    f = (a + b*x**3)**2*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/x**4
    F = -a**2*c/(3*x**3) - a**2*d/(2*x**2) - a**2*e/x + a*x**2*(a*h + 2*b*e)/2 + a*x*(a*g + 2*b*d) + a*(a*f + 2*b*c)*log(x) + b**2*f*x**6/6 + b**2*g*x**7/7 + b**2*h*x**8/8 + b*x**5*(2*a*h + b*e)/5 + b*x**4*(2*a*g + b*d)/4 + b*x**3*(2*a*f + b*c)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_392():
    f = (a + b*x**3)**2*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/x**5
    F = -a**2*c/(4*x**4) - a**2*d/(3*x**3) - a**2*e/(2*x**2) + a*x*(a*h + 2*b*e) + a*(a*g + 2*b*d)*log(x) - a*(a*f + 2*b*c)/x + b**2*f*x**5/5 + b**2*g*x**6/6 + b**2*h*x**7/7 + b*x**4*(2*a*h + b*e)/4 + b*x**3*(2*a*g + b*d)/3 + b*x**2*(2*a*f + b*c)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_393():
    f = x**4*(a + b*x**3)**3*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)
    F = a**3*c*x**5/5 + a**3*d*x**6/6 + a**3*e*x**7/7 + a**2*x**10*(a*h + 3*b*e)/10 + a**2*x**9*(a*g + 3*b*d)/9 + a**2*x**8*(a*f + 3*b*c)/8 + 3*a*b*x**13*(a*h + b*e)/13 + a*b*x**12*(a*g + b*d)/4 + 3*a*b*x**11*(a*f + b*c)/11 + b**3*f*x**17/17 + b**3*g*x**18/18 + b**3*h*x**19/19 + b**2*x**16*(3*a*h + b*e)/16 + b**2*x**15*(3*a*g + b*d)/15 + b**2*x**14*(3*a*f + b*c)/14
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_394():
    f = x**3*(a + b*x**3)**3*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)
    F = a**3*c*x**4/4 + a**3*d*x**5/5 + a**3*e*x**6/6 + a**2*x**9*(a*h + 3*b*e)/9 + a**2*x**8*(a*g + 3*b*d)/8 + a**2*x**7*(a*f + 3*b*c)/7 + a*b*x**12*(a*h + b*e)/4 + 3*a*b*x**11*(a*g + b*d)/11 + 3*a*b*x**10*(a*f + b*c)/10 + b**3*f*x**16/16 + b**3*g*x**17/17 + b**3*h*x**18/18 + b**2*x**15*(3*a*h + b*e)/15 + b**2*x**14*(3*a*g + b*d)/14 + b**2*x**13*(3*a*f + b*c)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_395():
    f = x**2*(a + b*x**3)**3*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)
    F = a**3*d*x**4/4 + a**3*e*x**5/5 + a**3*f*x**6/6 + a**2*b*f*x**9/3 + a**2*x**8*(a*h + 3*b*e)/8 + a**2*x**7*(a*g + 3*b*d)/7 + a*b**2*f*x**12/4 + 3*a*b*x**11*(a*h + b*e)/11 + 3*a*b*x**10*(a*g + b*d)/10 + b**3*f*x**15/15 + b**3*g*x**16/16 + b**3*h*x**17/17 + b**2*x**14*(3*a*h + b*e)/14 + b**2*x**13*(3*a*g + b*d)/13 + c*(a + b*x**3)**4/(12*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_396():
    f = x*(a + b*x**3)**3*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)
    F = a**3*c*x**2/2 + a**3*e*x**4/4 + a**3*g*x**6/6 + a**2*b*g*x**9/3 + a**2*x**7*(a*h + 3*b*e)/7 + a**2*x**5*(a*f + 3*b*c)/5 + a*b**2*g*x**12/4 + 3*a*b*x**10*(a*h + b*e)/10 + 3*a*b*x**8*(a*f + b*c)/8 + b**3*f*x**14/14 + b**3*g*x**15/15 + b**3*h*x**16/16 + b**2*x**13*(3*a*h + b*e)/13 + b**2*x**11*(3*a*f + b*c)/11 + d*(a + b*x**3)**4/(12*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_397():
    f = (a + b*x**3)**3*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)
    F = a**3*c*x + a**3*d*x**2/2 + a**3*h*x**6/6 + a**2*b*h*x**9/3 + a**2*x**5*(a*g + 3*b*d)/5 + a**2*x**4*(a*f + 3*b*c)/4 + a*b**2*h*x**12/4 + 3*a*b*x**8*(a*g + b*d)/8 + 3*a*b*x**7*(a*f + b*c)/7 + b**3*f*x**13/13 + b**3*g*x**14/14 + b**3*h*x**15/15 + b**2*x**11*(3*a*g + b*d)/11 + b**2*x**10*(3*a*f + b*c)/10 + e*(a + b*x**3)**4/(12*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_398():
    f = (a + b*x**3)**3*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/x
    F = a**3*c*log(x) + a**3*d*x + a**3*e*x**2/2 + a**2*b*c*x**3 + a**2*x**5*(a*h + 3*b*e)/5 + a**2*x**4*(a*g + 3*b*d)/4 + a*b**2*c*x**6/2 + 3*a*b*x**8*(a*h + b*e)/8 + 3*a*b*x**7*(a*g + b*d)/7 + b**3*c*x**9/9 + b**3*g*x**13/13 + b**3*h*x**14/14 + b**2*x**11*(3*a*h + b*e)/11 + b**2*x**10*(3*a*g + b*d)/10 + f*(a + b*x**3)**4/(12*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_399():
    f = (a + b*x**3)**3*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/x**2
    F = -a**3*c/x + a**3*d*log(x) + a**3*e*x + a**2*b*d*x**3 + a**2*x**4*(a*h + 3*b*e)/4 + a**2*x**2*(a*f + 3*b*c)/2 + a*b**2*d*x**6/2 + 3*a*b*x**7*(a*h + b*e)/7 + 3*a*b*x**5*(a*f + b*c)/5 + b**3*d*x**9/9 + b**3*f*x**11/11 + b**3*h*x**13/13 + b**2*x**10*(3*a*h + b*e)/10 + b**2*x**8*(3*a*f + b*c)/8 + g*(a + b*x**3)**4/(12*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_400():
    f = (a + b*x**3)**3*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/x**3
    F = -a**3*c/(2*x**2) - a**3*d/x + a**3*e*log(x) + a**2*b*e*x**3 + a**2*x**2*(a*g + 3*b*d)/2 + a**2*x*(a*f + 3*b*c) + a*b**2*e*x**6/2 + 3*a*b*x**5*(a*g + b*d)/5 + 3*a*b*x**4*(a*f + b*c)/4 + b**3*e*x**9/9 + b**3*f*x**10/10 + b**3*g*x**11/11 + b**2*x**8*(3*a*g + b*d)/8 + b**2*x**7*(3*a*f + b*c)/7 + h*(a + b*x**3)**4/(12*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_401():
    f = (a + b*x**3)**3*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/x**4
    F = -a**3*c/(3*x**3) - a**3*d/(2*x**2) - a**3*e/x + a**2*x**2*(a*h + 3*b*e)/2 + a**2*x*(a*g + 3*b*d) + a**2*(a*f + 3*b*c)*log(x) + 3*a*b*x**5*(a*h + b*e)/5 + 3*a*b*x**4*(a*g + b*d)/4 + a*b*x**3*(a*f + b*c) + b**3*f*x**9/9 + b**3*g*x**10/10 + b**3*h*x**11/11 + b**2*x**8*(3*a*h + b*e)/8 + b**2*x**7*(3*a*g + b*d)/7 + b**2*x**6*(3*a*f + b*c)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_402():
    f = (a + b*x**3)**3*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/x**5
    F = -a**3*c/(4*x**4) - a**3*d/(3*x**3) - a**3*e/(2*x**2) + a**2*x*(a*h + 3*b*e) + a**2*(a*g + 3*b*d)*log(x) - a**2*(a*f + 3*b*c)/x + 3*a*b*x**4*(a*h + b*e)/4 + a*b*x**3*(a*g + b*d) + 3*a*b*x**2*(a*f + b*c)/2 + b**3*f*x**8/8 + b**3*g*x**9/9 + b**3*h*x**10/10 + b**2*x**7*(3*a*h + b*e)/7 + b**2*x**6*(3*a*g + b*d)/6 + b**2*x**5*(3*a*f + b*c)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_403():
    f = x**4*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a + b*x**3)
    F = a**(sympy.S(2)/3)*(a**(sympy.S(2)/3)*(-a*h + b*e) + b**(sympy.S(2)/3)*(-a*f + b*c))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*b**(sympy.S(10)/3)) - a**(sympy.S(2)/3)*(a**(sympy.S(2)/3)*(-a*h + b*e) + b**(sympy.S(2)/3)*(-a*f + b*c))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*b**(sympy.S(10)/3)) + sqrt(3)*a**(sympy.S(2)/3)*(a**(sympy.S(5)/3)*h - a**(sympy.S(2)/3)*b*e - a*b**(sympy.S(2)/3)*f + b**(sympy.S(5)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(10)/3)) - a*x*(-a*h + b*e)/b**3 - a*(-a*g + b*d)*log(a + b*x**3)/(3*b**3) + f*x**5/(5*b) + g*x**6/(6*b) + h*x**7/(7*b) + x**4*(-a*h + b*e)/(4*b**2) + x**3*(-a*g + b*d)/(3*b**2) + x**2*(-a*f + b*c)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_404():
    f = x**3*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a + b*x**3)
    F = -a**(sympy.S(1)/3)*(-a**(sympy.S(1)/3)*(-a*g + b*d) + b**(sympy.S(1)/3)*(-a*f + b*c))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*b**(sympy.S(8)/3)) + a**(sympy.S(1)/3)*(-a**(sympy.S(1)/3)*(-a*g + b*d) + b**(sympy.S(1)/3)*(-a*f + b*c))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*b**(sympy.S(8)/3)) + sqrt(3)*a**(sympy.S(1)/3)*(-a**(sympy.S(4)/3)*g + a**(sympy.S(1)/3)*b*d - a*b**(sympy.S(1)/3)*f + b**(sympy.S(4)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(8)/3)) - a*(-a*h + b*e)*log(a + b*x**3)/(3*b**3) + f*x**4/(4*b) + g*x**5/(5*b) + h*x**6/(6*b) + x**3*(-a*h + b*e)/(3*b**2) + x**2*(-a*g + b*d)/(2*b**2) + x*(-a*f + b*c)/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_405():
    f = x**2*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a + b*x**3)
    F = -a**(sympy.S(1)/3)*(-a**(sympy.S(1)/3)*(-a*h + b*e) + b**(sympy.S(1)/3)*(-a*g + b*d))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*b**(sympy.S(8)/3)) + a**(sympy.S(1)/3)*(-a**(sympy.S(1)/3)*(-a*h + b*e) + b**(sympy.S(1)/3)*(-a*g + b*d))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*b**(sympy.S(8)/3)) + sqrt(3)*a**(sympy.S(1)/3)*(-a**(sympy.S(4)/3)*h + a**(sympy.S(1)/3)*b*e - a*b**(sympy.S(1)/3)*g + b**(sympy.S(4)/3)*d)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(8)/3)) + f*x**3/(3*b) + g*x**4/(4*b) + h*x**5/(5*b) + x**2*(-a*h + b*e)/(2*b**2) + x*(-a*g + b*d)/b**2 + (-a*f + b*c)*log(a + b*x**3)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_406():
    f = x*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a + b*x**3)
    F = f*x**2/(2*b) + g*x**3/(3*b) + h*x**4/(4*b) + x*(-a*h + b*e)/b**2 + (-a*g + b*d)*log(a + b*x**3)/(3*b**2) - (a**(sympy.S(2)/3)*(-a*h + b*e) + b**(sympy.S(2)/3)*(-a*f + b*c))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)*b**(sympy.S(7)/3)) + (a**(sympy.S(2)/3)*(-a*h + b*e) + b**(sympy.S(2)/3)*(-a*f + b*c))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(1)/3)*b**(sympy.S(7)/3)) - sqrt(3)*(a**(sympy.S(5)/3)*h - a**(sympy.S(2)/3)*b*e - a*b**(sympy.S(2)/3)*f + b**(sympy.S(5)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3)*b**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_407():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a + b*x**3)
    F = f*x/b + g*x**2/(2*b) + h*x**3/(3*b) + (-a*h + b*e)*log(a + b*x**3)/(3*b**2) + (-a**(sympy.S(1)/3)*(-a*g + b*d) + b**(sympy.S(1)/3)*(-a*f + b*c))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(2)/3)*b**(sympy.S(5)/3)) - (-a**(sympy.S(1)/3)*(-a*g + b*d) + b**(sympy.S(1)/3)*(-a*f + b*c))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(5)/3)) - sqrt(3)*(-a**(sympy.S(4)/3)*g + a**(sympy.S(1)/3)*b*d - a*b**(sympy.S(1)/3)*f + b**(sympy.S(4)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_408():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(x*(a + b*x**3))
    F = g*x/b + h*x**2/(2*b) + c*log(x)/a - (-a*f + b*c)*log(a + b*x**3)/(3*a*b) + (-a**(sympy.S(1)/3)*(-a*h + b*e) + b**(sympy.S(1)/3)*(-a*g + b*d))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(2)/3)*b**(sympy.S(5)/3)) - (-a**(sympy.S(1)/3)*(-a*h + b*e) + b**(sympy.S(1)/3)*(-a*g + b*d))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(5)/3)) - sqrt(3)*(-a**(sympy.S(4)/3)*h + a**(sympy.S(1)/3)*b*e - a*b**(sympy.S(1)/3)*g + b**(sympy.S(4)/3)*d)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_409():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(x**2*(a + b*x**3))
    F = h*x/b - c/(a*x) + d*log(x)/a - (-a*g + b*d)*log(a + b*x**3)/(3*a*b) + (a**(sympy.S(2)/3)*(-a*h + b*e) + b**(sympy.S(2)/3)*(-a*f + b*c))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(4)/3)*b**(sympy.S(4)/3)) - (a**(sympy.S(2)/3)*(-a*h + b*e) + b**(sympy.S(2)/3)*(-a*f + b*c))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(4)/3)*b**(sympy.S(4)/3)) + sqrt(3)*(a**(sympy.S(5)/3)*h - a**(sympy.S(2)/3)*b*e - a*b**(sympy.S(2)/3)*f + b**(sympy.S(5)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(4)/3)*b**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_410():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(x**3*(a + b*x**3))
    F = -c/(2*a*x**2) - d/(a*x) + e*log(x)/a - (-a*h + b*e)*log(a + b*x**3)/(3*a*b) - (-a**(sympy.S(1)/3)*(-a*g + b*d) + b**(sympy.S(1)/3)*(-a*f + b*c))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(5)/3)*b**(sympy.S(2)/3)) + (-a**(sympy.S(1)/3)*(-a*g + b*d) + b**(sympy.S(1)/3)*(-a*f + b*c))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(5)/3)*b**(sympy.S(2)/3)) + sqrt(3)*(-a**(sympy.S(4)/3)*g + a**(sympy.S(1)/3)*b*d - a*b**(sympy.S(1)/3)*f + b**(sympy.S(4)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(5)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_411():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(x**4*(a + b*x**3))
    F = -c/(3*a*x**3) - d/(2*a*x**2) - e/(a*x) - (-a*f + b*c)*log(x)/a**2 + (-a*f + b*c)*log(a + b*x**3)/(3*a**2) - (-a**(sympy.S(1)/3)*(-a*h + b*e) + b**(sympy.S(1)/3)*(-a*g + b*d))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(5)/3)*b**(sympy.S(2)/3)) + (-a**(sympy.S(1)/3)*(-a*h + b*e) + b**(sympy.S(1)/3)*(-a*g + b*d))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(5)/3)*b**(sympy.S(2)/3)) + sqrt(3)*(-a**(sympy.S(4)/3)*h + a**(sympy.S(1)/3)*b*e - a*b**(sympy.S(1)/3)*g + b**(sympy.S(4)/3)*d)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(5)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_412():
    f = x**4*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a + b*x**3)**2
    F = f*x**2/(2*b**2) + g*x**3/(3*b**2) + h*x**4/(4*b**2) + x*(-2*a*h + b*e)/b**3 + x*(a*(-a*h + b*e) - b*x**2*(-a*g + b*d) - b*x*(-a*f + b*c))/(3*b**3*(a + b*x**3)) + (-2*a*g + b*d)*log(a + b*x**3)/(3*b**3) - (a**(sympy.S(2)/3)*(-7*a*h + 4*b*e) + b**(sympy.S(2)/3)*(-5*a*f + 2*b*c))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(1)/3)*b**(sympy.S(10)/3)) + (a**(sympy.S(2)/3)*(-7*a*h + 4*b*e) + b**(sympy.S(2)/3)*(-5*a*f + 2*b*c))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(1)/3)*b**(sympy.S(10)/3)) - sqrt(3)*(7*a**(sympy.S(5)/3)*h - 4*a**(sympy.S(2)/3)*b*e - 5*a*b**(sympy.S(2)/3)*f + 2*b**(sympy.S(5)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(1)/3)*b**(sympy.S(10)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_413():
    f = x**3*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a + b*x**3)**2
    F = f*x/b**2 + g*x**2/(2*b**2) + h*x**3/(3*b**2) - x*(-a*f + b*c + x**2*(-a*h + b*e) + x*(-a*g + b*d))/(3*b**2*(a + b*x**3)) + (-2*a*h + b*e)*log(a + b*x**3)/(3*b**3) + (-a**(sympy.S(1)/3)*(-5*a*g + 2*b*d) + b**(sympy.S(1)/3)*(-4*a*f + b*c))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(2)/3)*b**(sympy.S(8)/3)) - (-a**(sympy.S(1)/3)*(-5*a*g + 2*b*d) + b**(sympy.S(1)/3)*(-4*a*f + b*c))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(2)/3)*b**(sympy.S(8)/3)) - sqrt(3)*(-5*a**(sympy.S(4)/3)*g + 2*a**(sympy.S(1)/3)*b*d - 4*a*b**(sympy.S(1)/3)*f + b**(sympy.S(4)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(2)/3)*b**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_414():
    f = x**2*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a + b*x**3)**2
    F = -(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(3*b*(a + b*x**3)) + f*log(a + b*x**3)/(3*b**2) + 4*g*x/(3*b**2) + 5*h*x**2/(6*b**2) + (-a**(sympy.S(1)/3)*(-5*a*h + 2*b*e) + b**(sympy.S(1)/3)*(-4*a*g + b*d))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(2)/3)*b**(sympy.S(8)/3)) - (-a**(sympy.S(1)/3)*(-5*a*h + 2*b*e) + b**(sympy.S(1)/3)*(-4*a*g + b*d))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(2)/3)*b**(sympy.S(8)/3)) - sqrt(3)*(-5*a**(sympy.S(4)/3)*h + 2*a**(sympy.S(1)/3)*b*e - 4*a*b**(sympy.S(1)/3)*g + b**(sympy.S(4)/3)*d)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(2)/3)*b**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_415():
    f = x*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a + b*x**3)**2
    F = g*log(a + b*x**3)/(3*b**2) + h*x/b**2 - x*(a*(-a*h + b*e) - b*x**2*(-a*g + b*d) - b*x*(-a*f + b*c))/(3*a*b**2*(a + b*x**3)) - (-a**(sympy.S(2)/3)*(-4*a*h + b*e) + b**(sympy.S(2)/3)*(2*a*f + b*c))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(4)/3)*b**(sympy.S(7)/3)) + (-a**(sympy.S(2)/3)*(-4*a*h + b*e) + b**(sympy.S(2)/3)*(2*a*f + b*c))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(4)/3)*b**(sympy.S(7)/3)) - sqrt(3)*(-4*a**(sympy.S(5)/3)*h + a**(sympy.S(2)/3)*b*e + 2*a*b**(sympy.S(2)/3)*f + b**(sympy.S(5)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(4)/3)*b**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_416():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a + b*x**3)**2
    F = h*log(a + b*x**3)/(3*b**2) + x*(-a*f + b*c + x**2*(-a*h + b*e) + x*(-a*g + b*d))/(3*a*b*(a + b*x**3)) + (-a**(sympy.S(1)/3)*(2*a*g + b*d) + b**(sympy.S(1)/3)*(a*f + 2*b*c))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(5)/3)*b**(sympy.S(5)/3)) - (-a**(sympy.S(1)/3)*(2*a*g + b*d) + b**(sympy.S(1)/3)*(a*f + 2*b*c))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(5)/3)*b**(sympy.S(5)/3)) - sqrt(3)*(2*a**(sympy.S(4)/3)*g + a**(sympy.S(1)/3)*b*d + a*b**(sympy.S(1)/3)*f + 2*b**(sympy.S(4)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*b**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_417():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(x*(a + b*x**3)**2)
    F = c*log(x)/a**2 - c*log(a + b*x**3)/(3*a**2) + x*(a*x*(-a*h + b*e) + a*(-a*g + b*d) - b*x**2*(-a*f + b*c))/(3*a**2*b*(a + b*x**3)) + (-a**(sympy.S(1)/3)*(2*a*h + b*e) + b**(sympy.S(1)/3)*(a*g + 2*b*d))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(5)/3)*b**(sympy.S(5)/3)) - (-a**(sympy.S(1)/3)*(2*a*h + b*e) + b**(sympy.S(1)/3)*(a*g + 2*b*d))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(5)/3)*b**(sympy.S(5)/3)) - sqrt(3)*(2*a**(sympy.S(4)/3)*h + a**(sympy.S(1)/3)*b*e + a*b**(sympy.S(1)/3)*g + 2*b**(sympy.S(4)/3)*d)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*b**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_418():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(x**2*(a + b*x**3)**2)
    F = -c/(a**2*x) + d*log(x)/a**2 - d*log(a + b*x**3)/(3*a**2) + x*(a*(-a*h + b*e) - b*x**2*(-a*g + b*d) - b*x*(-a*f + b*c))/(3*a**2*b*(a + b*x**3)) + (a**(sympy.S(2)/3)*(a*h + 2*b*e) + b**(sympy.S(2)/3)*(-a*f + 4*b*c))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(7)/3)*b**(sympy.S(4)/3)) - (a**(sympy.S(2)/3)*(a*h + 2*b*e) + b**(sympy.S(2)/3)*(-a*f + 4*b*c))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(7)/3)*b**(sympy.S(4)/3)) + sqrt(3)*(-a**(sympy.S(5)/3)*h - 2*a**(sympy.S(2)/3)*b*e - a*b**(sympy.S(2)/3)*f + 4*b**(sympy.S(5)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(7)/3)*b**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_419():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(x**3*(a + b*x**3)**2)
    F = -c/(2*a**2*x**2) - d/(a**2*x) + e*log(x)/a**2 - e*log(a + b*x**3)/(3*a**2) - x*(-a*f + b*c + x**2*(-a*h + b*e) + x*(-a*g + b*d))/(3*a**2*(a + b*x**3)) - (-a**(sympy.S(1)/3)*(-a*g + 4*b*d) + b**(sympy.S(1)/3)*(-2*a*f + 5*b*c))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(8)/3)*b**(sympy.S(2)/3)) + (-a**(sympy.S(1)/3)*(-a*g + 4*b*d) + b**(sympy.S(1)/3)*(-2*a*f + 5*b*c))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(8)/3)*b**(sympy.S(2)/3)) + sqrt(3)*(-a**(sympy.S(4)/3)*g + 4*a**(sympy.S(1)/3)*b*d - 2*a*b**(sympy.S(1)/3)*f + 5*b**(sympy.S(4)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(8)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_420():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(x**4*(a + b*x**3)**2)
    F = -c/(3*a**2*x**3) - d/(2*a**2*x**2) - e/(a**2*x) - x*(-a*g + b*d - b*x**2*(-f + b*c/a) + x*(-a*h + b*e))/(3*a**2*(a + b*x**3)) - (-a*f + 2*b*c)*log(x)/a**3 + (-a*f + 2*b*c)*log(a + b*x**3)/(3*a**3) - (-a**(sympy.S(1)/3)*(-a*h + 4*b*e) + b**(sympy.S(1)/3)*(-2*a*g + 5*b*d))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(8)/3)*b**(sympy.S(2)/3)) + (-a**(sympy.S(1)/3)*(-a*h + 4*b*e) + b**(sympy.S(1)/3)*(-2*a*g + 5*b*d))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(8)/3)*b**(sympy.S(2)/3)) + sqrt(3)*(-a**(sympy.S(4)/3)*h + 4*a**(sympy.S(1)/3)*b*e - 2*a*b**(sympy.S(1)/3)*g + 5*b**(sympy.S(4)/3)*d)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(8)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_421():
    f = x**4*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a + b*x**3)**3
    F = g*log(a + b*x**3)/(3*b**3) + h*x/b**3 + x*(a*(-a*h + b*e) - b*x**2*(-a*g + b*d) - b*x*(-a*f + b*c))/(6*b**3*(a + b*x**3)**2) - x*(a*(-13*a*h + 7*b*e) - 3*b*x**2*(-3*a*g + b*d) - 2*b*x*(-4*a*f + b*c))/(18*a*b**3*(a + b*x**3)) - (-2*a**(sympy.S(2)/3)*(-7*a*h + b*e) + b**(sympy.S(2)/3)*(5*a*f + b*c))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(4)/3)*b**(sympy.S(10)/3)) + (-2*a**(sympy.S(2)/3)*(-7*a*h + b*e) + b**(sympy.S(2)/3)*(5*a*f + b*c))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(4)/3)*b**(sympy.S(10)/3)) - sqrt(3)*(-14*a**(sympy.S(5)/3)*h + 2*a**(sympy.S(2)/3)*b*e + 5*a*b**(sympy.S(2)/3)*f + b**(sympy.S(5)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(4)/3)*b**(sympy.S(10)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_422():
    f = x**3*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a + b*x**3)**3
    F = -x*(-a*f + b*c + x**2*(-a*h + b*e) + x*(-a*g + b*d))/(6*b**2*(a + b*x**3)**2) + h*log(a + b*x**3)/(3*b**3) + x*(-7*a*f + b*c + x**2*(-9*a*h + 3*b*e) + x*(-8*a*g + 2*b*d))/(18*a*b**2*(a + b*x**3)) + (-a**(sympy.S(1)/3)*(5*a*g + b*d) + b**(sympy.S(1)/3)*(2*a*f + b*c))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(5)/3)*b**(sympy.S(8)/3)) - (-a**(sympy.S(1)/3)*(5*a*g + b*d) + b**(sympy.S(1)/3)*(2*a*f + b*c))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(5)/3)*b**(sympy.S(8)/3)) - sqrt(3)*(5*a**(sympy.S(4)/3)*g + a**(sympy.S(1)/3)*b*d + 2*a*b**(sympy.S(1)/3)*f + b**(sympy.S(4)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(5)/3)*b**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_423():
    f = x**2*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a + b*x**3)**3
    F = -(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(6*b*(a + b*x**3)**2) + x*(-4*a*g + b*d + 3*b*f*x**2 + x*(-5*a*h + 2*b*e))/(18*a*b**2*(a + b*x**3)) + (-a**(sympy.S(1)/3)*(5*a*h + b*e) + b**(sympy.S(1)/3)*(2*a*g + b*d))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(5)/3)*b**(sympy.S(8)/3)) - (-a**(sympy.S(1)/3)*(5*a*h + b*e) + b**(sympy.S(1)/3)*(2*a*g + b*d))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(5)/3)*b**(sympy.S(8)/3)) - sqrt(3)*(5*a**(sympy.S(4)/3)*h + a**(sympy.S(1)/3)*b*e + 2*a*b**(sympy.S(1)/3)*g + b**(sympy.S(4)/3)*d)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(5)/3)*b**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_424():
    f = x*(c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a + b*x**3)**3
    F = -x*(a*(-a*h + b*e) - b*x**2*(-a*g + b*d) - b*x*(-a*f + b*c))/(6*a*b**2*(a + b*x**3)**2) + x*(a*(-7*a*h + b*e) + 3*b*x**2*(a*g + b*d) + 2*b*x*(a*f + 2*b*c))/(18*a**2*b**2*(a + b*x**3)) - (-a**(sympy.S(2)/3)*(2*a*h + b*e) + b**(sympy.S(2)/3)*(a*f + 2*b*c))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(7)/3)*b**(sympy.S(7)/3)) + (-a**(sympy.S(2)/3)*(2*a*h + b*e) + b**(sympy.S(2)/3)*(a*f + 2*b*c))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(7)/3)*b**(sympy.S(7)/3)) - sqrt(3)*(2*a**(sympy.S(5)/3)*h + a**(sympy.S(2)/3)*b*e + a*b**(sympy.S(2)/3)*f + 2*b**(sympy.S(5)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(7)/3)*b**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_425():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(a + b*x**3)**3
    F = x*(-a*f + b*c + x**2*(-a*h + b*e) + x*(-a*g + b*d))/(6*a*b*(a + b*x**3)**2) - (3*a*(a*h + b*e) - b*x*(a*f + 5*b*c + x*(2*a*g + 4*b*d)))/(18*a**2*b**2*(a + b*x**3)) + (-a**(sympy.S(1)/3)*(a*g + 2*b*d) + b**(sympy.S(1)/3)*(a*f + 5*b*c))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(8)/3)*b**(sympy.S(5)/3)) - (-a**(sympy.S(1)/3)*(a*g + 2*b*d) + b**(sympy.S(1)/3)*(a*f + 5*b*c))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(8)/3)*b**(sympy.S(5)/3)) - sqrt(3)*(a**(sympy.S(4)/3)*g + 2*a**(sympy.S(1)/3)*b*d + a*b**(sympy.S(1)/3)*f + 5*b**(sympy.S(4)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(8)/3)*b**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_426():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(x*(a + b*x**3)**3)
    F = x*(a*x*(-a*h + b*e) + a*(-a*g + b*d) - b*x**2*(-a*f + b*c))/(6*a**2*b*(a + b*x**3)**2) + c*log(x)/a**3 - c*log(a + b*x**3)/(3*a**3) + x*(2*a*x*(a*h + 2*b*e) + a*(a*g + 5*b*d) - 3*b*x**2*(-a*f + 3*b*c))/(18*a**3*b*(a + b*x**3)) + (-a**(sympy.S(1)/3)*(a*h + 2*b*e) + b**(sympy.S(1)/3)*(a*g + 5*b*d))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(8)/3)*b**(sympy.S(5)/3)) - (-a**(sympy.S(1)/3)*(a*h + 2*b*e) + b**(sympy.S(1)/3)*(a*g + 5*b*d))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(8)/3)*b**(sympy.S(5)/3)) - sqrt(3)*(a**(sympy.S(4)/3)*h + 2*a**(sympy.S(1)/3)*b*e + a*b**(sympy.S(1)/3)*g + 5*b**(sympy.S(4)/3)*d)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(8)/3)*b**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_427():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(x**2*(a + b*x**3)**3)
    F = x*(a*(-a*h + b*e) - b*x**2*(-a*g + b*d) - b*x*(-a*f + b*c))/(6*a**2*b*(a + b*x**3)**2) - c/(a**3*x) + d*log(x)/a**3 - d*log(a + b*x**3)/(3*a**3) + x*(a*(a*h + 5*b*e) - 3*b*x**2*(-a*g + 3*b*d) - 2*b*x*(-2*a*f + 5*b*c))/(18*a**3*b*(a + b*x**3)) + (a**(sympy.S(2)/3)*(a*h + 5*b*e) + 2*b**(sympy.S(2)/3)*(-a*f + 7*b*c))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(10)/3)*b**(sympy.S(4)/3)) - (a**(sympy.S(2)/3)*(a*h + 5*b*e) + 2*b**(sympy.S(2)/3)*(-a*f + 7*b*c))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(10)/3)*b**(sympy.S(4)/3)) + sqrt(3)*(-a**(sympy.S(5)/3)*h - 5*a**(sympy.S(2)/3)*b*e - 2*a*b**(sympy.S(2)/3)*f + 14*b**(sympy.S(5)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(10)/3)*b**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_428():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(x**3*(a + b*x**3)**3)
    F = -x*(-a*f + b*c + x**2*(-a*h + b*e) + x*(-a*g + b*d))/(6*a**2*(a + b*x**3)**2) - c/(2*a**3*x**2) - d/(a**3*x) + e*log(x)/a**3 - e*log(a + b*x**3)/(3*a**3) - x*(-5*a*f + 11*b*c + x**2*(-3*a*h + 9*b*e) + x*(-4*a*g + 10*b*d))/(18*a**3*(a + b*x**3)) - (-2*a**(sympy.S(1)/3)*(-a*g + 7*b*d) + 5*b**(sympy.S(1)/3)*(-a*f + 4*b*c))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(11)/3)*b**(sympy.S(2)/3)) + (-2*a**(sympy.S(1)/3)*(-a*g + 7*b*d) + 5*b**(sympy.S(1)/3)*(-a*f + 4*b*c))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(11)/3)*b**(sympy.S(2)/3)) + sqrt(3)*(-2*a**(sympy.S(4)/3)*g + 14*a**(sympy.S(1)/3)*b*d - 5*a*b**(sympy.S(1)/3)*f + 20*b**(sympy.S(4)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(11)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_429():
    f = (c + d*x + e*x**2 + f*x**3 + g*x**4 + h*x**5)/(x**4*(a + b*x**3)**3)
    F = -x*(-a*g + b*d - b*x**2*(-f + b*c/a) + x*(-a*h + b*e))/(6*a**2*(a + b*x**3)**2) - c/(3*a**3*x**3) - d/(2*a**3*x**2) - e/(a**3*x) - x*(-5*a*g + 11*b*d - 3*b*x**2*(-3*f + 5*b*c/a) + x*(-4*a*h + 10*b*e))/(18*a**3*(a + b*x**3)) - (-a*f + 3*b*c)*log(x)/a**4 + (-a*f + 3*b*c)*log(a + b*x**3)/(3*a**4) - (-2*a**(sympy.S(1)/3)*(-a*h + 7*b*e) + 5*b**(sympy.S(1)/3)*(-a*g + 4*b*d))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(11)/3)*b**(sympy.S(2)/3)) + (-2*a**(sympy.S(1)/3)*(-a*h + 7*b*e) + 5*b**(sympy.S(1)/3)*(-a*g + 4*b*d))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(11)/3)*b**(sympy.S(2)/3)) + sqrt(3)*(-2*a**(sympy.S(4)/3)*h + 14*a**(sympy.S(1)/3)*b*e - 5*a*b**(sympy.S(1)/3)*g + 20*b**(sympy.S(4)/3)*d)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(11)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_430():
    f = x**3*(c + d*x + e*x**2)/sqrt(a + b*x**3)
    F = 4*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*d*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(7*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 4*a*e*sqrt(a + b*x**3)/(9*b**2) - 8*a*d*sqrt(a + b*x**3)/(7*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) - 4*3**(sympy.S(3)/4)*a*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*d*(10 - 10*sqrt(3)) + 7*b**(sympy.S(1)/3)*c)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(105*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*c*x*sqrt(a + b*x**3)/(5*b) + 2*d*x**2*sqrt(a + b*x**3)/(7*b) + 2*e*x**3*sqrt(a + b*x**3)/(9*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_431():
    f = x**2*(c + d*x + e*x**2)/sqrt(a + b*x**3)
    F = 4*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*e*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(7*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 8*a*e*sqrt(a + b*x**3)/(7*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) - 4*3**(sympy.S(3)/4)*a*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*e*(10 - 10*sqrt(3)) + 7*b**(sympy.S(1)/3)*d)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(105*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*c*sqrt(a + b*x**3)/(3*b) + 2*d*x*sqrt(a + b*x**3)/(5*b) + 2*e*x**2*sqrt(a + b*x**3)/(7*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_432():
    f = x*(c + d*x + e*x**2)/sqrt(a + b*x**3)
    F = -3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*c*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 2*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(2*a**(sympy.S(2)/3)*e + b**(sympy.S(2)/3)*c*(5 - 5*sqrt(3)))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(15*b**(sympy.S(4)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*d*sqrt(a + b*x**3)/(3*b) + 2*e*x*sqrt(a + b*x**3)/(5*b) + 2*c*sqrt(a + b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_433():
    f = (c + d*x + e*x**2)/sqrt(a + b*x**3)
    F = -3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*d*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*e*sqrt(a + b*x**3)/(3*b) + 2*d*sqrt(a + b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*d*(1 - sqrt(3)) + b**(sympy.S(1)/3)*c)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_434():
    f = (c + d*x + e*x**2)/(x*sqrt(a + b*x**3))
    F = -3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*e*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*e*sqrt(a + b*x**3)/(b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*e*(1 - sqrt(3)) + b**(sympy.S(1)/3)*d)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 2*c*atanh(sqrt(a + b*x**3)/sqrt(a))/(3*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_435():
    f = (c + d*x + e*x**2)/(x**2*sqrt(a + b*x**3))
    F = b**(sympy.S(1)/3)*c*sqrt(a + b*x**3)/(a*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) - c*sqrt(a + b*x**3)/(a*x) - 2*d*atanh(sqrt(a + b*x**3)/sqrt(a))/(3*sqrt(a)) - 3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*c*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(2*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-2*a**(sympy.S(2)/3)*e + b**(sympy.S(2)/3)*c*(1 - sqrt(3)))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_436():
    f = (c + d*x + e*x**2)/(x**3*sqrt(a + b*x**3))
    F = b**(sympy.S(1)/3)*d*sqrt(a + b*x**3)/(a*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) - 3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*d*(2 - 2*sqrt(3)) + b**(sympy.S(1)/3)*c)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(6*a*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - c*sqrt(a + b*x**3)/(2*a*x**2) - d*sqrt(a + b*x**3)/(a*x) - 2*e*atanh(sqrt(a + b*x**3)/sqrt(a))/(3*sqrt(a)) - 3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*d*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(2*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_437():
    f = x**5*(c + d*x + e*x**2)/(a + b*x**3)**(sympy.S(3)/2)
    F = 40*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*e*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(21*b**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 80*a*e*sqrt(a + b*x**3)/(21*b**(sympy.S(8)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) - 16*3**(sympy.S(3)/4)*a*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*e*(25 - 25*sqrt(3)) + 14*b**(sympy.S(1)/3)*d)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(315*b**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 4*c*sqrt(a + b*x**3)/(3*b**2) + 2*d*x*sqrt(a + b*x**3)/(5*b**2) + 2*e*x**2*sqrt(a + b*x**3)/(7*b**2) + 2*x*(a*d + a*e*x - b*c*x**2)/(3*b**2*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_438():
    f = x**4*(c + d*x + e*x**2)/(a + b*x**3)**(sympy.S(3)/2)
    F = -4*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*c*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 8*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(4*a**(sympy.S(2)/3)*e + b**(sympy.S(2)/3)*c*(5 - 5*sqrt(3)))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(45*b**(sympy.S(7)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 4*d*sqrt(a + b*x**3)/(3*b**2) + 2*e*x*sqrt(a + b*x**3)/(5*b**2) + 2*x*(a*e - b*c*x - b*d*x**2)/(3*b**2*sqrt(a + b*x**3)) + 8*c*sqrt(a + b*x**3)/(3*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_439():
    f = x**3*(c + d*x + e*x**2)/(a + b*x**3)**(sympy.S(3)/2)
    F = -4*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*d*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 2*x*(c + d*x + e*x**2)/(3*b*sqrt(a + b*x**3)) + 4*e*sqrt(a + b*x**3)/(3*b**2) + 8*d*sqrt(a + b*x**3)/(3*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 4*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*d*(2 - 2*sqrt(3)) + b**(sympy.S(1)/3)*c)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_440():
    f = x**2*(c + d*x + e*x**2)/(a + b*x**3)**(sympy.S(3)/2)
    F = -4*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*e*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - (2*c + 2*d*x + 2*e*x**2)/(3*b*sqrt(a + b*x**3)) + 8*e*sqrt(a + b*x**3)/(3*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 4*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*e*(2 - 2*sqrt(3)) + b**(sympy.S(1)/3)*d)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_441():
    f = x*(c + d*x + e*x**2)/(a + b*x**3)**(sympy.S(3)/2)
    F = -2*d*sqrt(a + b*x**3)/(3*a*b) - 2*x*(a*e - b*c*x - b*d*x**2)/(3*a*b*sqrt(a + b*x**3)) - 2*c*sqrt(a + b*x**3)/(3*a*b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 3**(sympy.S(1)/4)*c*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*a**(sympy.S(2)/3)*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(2*a**(sympy.S(2)/3)*e + b**(sympy.S(2)/3)*(-sqrt(3)*c + c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*a**(sympy.S(2)/3)*b**(sympy.S(4)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_442():
    f = (c + d*x + e*x**2)/(a + b*x**3)**(sympy.S(3)/2)
    F = -(2*a*e - 2*b*x*(c + d*x))/(3*a*b*sqrt(a + b*x**3)) - 2*d*sqrt(a + b*x**3)/(3*a*b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*d*(1 - sqrt(3)) + b**(sympy.S(1)/3)*c)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*a*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 3**(sympy.S(1)/4)*d*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*a**(sympy.S(2)/3)*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_443():
    f = (c + d*x + e*x**2)/(x*(a + b*x**3)**(sympy.S(3)/2))
    F = -2*e*sqrt(a + b*x**3)/(3*a*b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*e*(1 - sqrt(3)) + b**(sympy.S(1)/3)*d)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*a*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*c*sqrt(a + b*x**3)/(3*a**2) + 2*x*(a*d + a*e*x - b*c*x**2)/(3*a**2*sqrt(a + b*x**3)) - 2*c*atanh(sqrt(a + b*x**3)/sqrt(a))/(3*a**(sympy.S(3)/2)) + 3**(sympy.S(1)/4)*e*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*a**(sympy.S(2)/3)*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_444():
    f = (c + d*x + e*x**2)/(x**2*(a + b*x**3)**(sympy.S(3)/2))
    F = 5*b**(sympy.S(1)/3)*c*sqrt(a + b*x**3)/(3*a**2*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) - c*sqrt(a + b*x**3)/(a**2*x) + 2*d*sqrt(a + b*x**3)/(3*a**2) + 2*x*(a*e - b*c*x - b*d*x**2)/(3*a**2*sqrt(a + b*x**3)) - 2*d*atanh(sqrt(a + b*x**3)/sqrt(a))/(3*a**(sympy.S(3)/2)) - 5*3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*c*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(6*a**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-2*a**(sympy.S(2)/3)*e + b**(sympy.S(2)/3)*c*(5 - 5*sqrt(3)))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*a**(sympy.S(5)/3)*b**(sympy.S(1)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_445():
    f = x**3*sqrt(a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4)
    F = 12*3**(sympy.S(1)/4)*a**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-10*a*g + 19*b*d)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1729*b**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 4*a**2*e*sqrt(a + b*x**3)/(45*b**2) - 4*3**(sympy.S(3)/4)*a**2*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*(1870 - 1870*sqrt(3))*(-10*a*g + 19*b*d) + 1729*b**(sympy.S(1)/3)*(-8*a*f + 17*b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1616615*b**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 24*a**2*sqrt(a + b*x**3)*(-10*a*g + 19*b*d)/(1729*b**(sympy.S(8)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*a*e*x**3*sqrt(a + b*x**3)/(45*b) + 6*a*f*x**4*sqrt(a + b*x**3)/(187*b) + 6*a*g*x**5*sqrt(a + b*x**3)/(247*b) + 6*a*x**2*sqrt(a + b*x**3)*(-10*a*g + 19*b*d)/(1729*b**2) + 6*a*x*sqrt(a + b*x**3)*(-8*a*f + 17*b*c)/(935*b**2) + 2*x**3*sqrt(a + b*x**3)*(62985*c*x + 53295*d*x**2 + 46189*e*x**3 + 40755*f*x**4 + 36465*g*x**5)/692835
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_446():
    f = x**2*sqrt(a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4)
    F = 12*3**(sympy.S(1)/4)*a**(sympy.S(7)/3)*e*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(91*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 24*a**2*e*sqrt(a + b*x**3)/(91*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) - 4*3**(sympy.S(3)/4)*a**2*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*e*(1870 - 1870*sqrt(3)) - 728*a*g + 1547*b*d)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(85085*b**(sympy.S(7)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 6*a*e*x**2*sqrt(a + b*x**3)/(91*b) + 2*a*f*x**3*sqrt(a + b*x**3)/(45*b) + 6*a*g*x**4*sqrt(a + b*x**3)/(187*b) + 6*a*x*sqrt(a + b*x**3)*(-8*a*g + 17*b*d)/(935*b**2) + 2*a*sqrt(a + b*x**3)*(-2*a*f + 5*b*c)/(45*b**2) + 2*x**2*sqrt(a + b*x**3)*(12155*c*x + 9945*d*x**2 + 8415*e*x**3 + 7293*f*x**4 + 6435*g*x**5)/109395
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_447():
    f = x*sqrt(a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4)
    F = -3*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-4*a*f + 13*b*c)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(91*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 2*3**(sympy.S(3)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(182*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*e + (55 - 55*sqrt(3))*(-4*a*f + 13*b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(5005*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 6*a*e*x*sqrt(a + b*x**3)/(55*b) + 6*a*f*x**2*sqrt(a + b*x**3)/(91*b) + 2*a*g*x**3*sqrt(a + b*x**3)/(45*b) + 2*a*sqrt(a + b*x**3)*(-2*a*g + 5*b*d)/(45*b**2) + 6*a*sqrt(a + b*x**3)*(-4*a*f + 13*b*c)/(91*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*x*sqrt(a + b*x**3)*(6435*c*x + 5005*d*x**2 + 4095*e*x**3 + 3465*f*x**4 + 3003*g*x**5)/45045
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_448():
    f = sqrt(a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4)
    F = -3*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-4*a*g + 13*b*d)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(91*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*a*e*sqrt(a + b*x**3)/(9*b) + 6*a*f*x*sqrt(a + b*x**3)/(55*b) + 6*a*g*x**2*sqrt(a + b*x**3)/(91*b) + 2*3**(sympy.S(3)/4)*a*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*(55 - 55*sqrt(3))*(-4*a*g + 13*b*d) + 91*b**(sympy.S(1)/3)*(-2*a*f + 11*b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(5005*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 6*a*sqrt(a + b*x**3)*(-4*a*g + 13*b*d)/(91*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*sqrt(a + b*x**3)*(9009*c*x + 6435*d*x**2 + 5005*e*x**3 + 4095*f*x**4 + 3465*g*x**5)/45045
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_449():
    f = sqrt(a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4)/x
    F = -3*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*e*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(7*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 2*sqrt(a)*c*atanh(sqrt(a + b*x**3)/sqrt(a))/3 + 2*a*f*sqrt(a + b*x**3)/(9*b) + 6*a*g*x*sqrt(a + b*x**3)/(55*b) + 6*a*e*sqrt(a + b*x**3)/(7*b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*3**(sympy.S(3)/4)*a*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*e*(55 - 55*sqrt(3)) - 14*a*g + 77*b*d)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(385*b**(sympy.S(4)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*sqrt(a + b*x**3)*(1155*c*x + 693*d*x**2 + 495*e*x**3 + 385*f*x**4 + 315*g*x**5)/(3465*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_450():
    f = sqrt(a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4)/x**2
    F = -3*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(2*a*f + 7*b*c)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(14*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(14*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*e - (5 - 5*sqrt(3))*(2*a*f + 7*b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(35*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 2*sqrt(a)*d*atanh(sqrt(a + b*x**3)/sqrt(a))/3 + 2*a*g*sqrt(a + b*x**3)/(9*b) - 3*c*sqrt(a + b*x**3)/x + 2*sqrt(a + b*x**3)*(315*c*x + 105*d*x**2 + 63*e*x**3 + 45*f*x**4 + 35*g*x**5)/(315*x**2) + sqrt(a + b*x**3)*(6*a*f + 21*b*c)/(7*b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_451():
    f = sqrt(a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4)/x**3
    F = -3*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(2*a*g + 7*b*d)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(14*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 2*sqrt(a)*e*atanh(sqrt(a + b*x**3)/sqrt(a))/3 + 3*c*sqrt(a + b*x**3)/(2*x**2) - 3*d*sqrt(a + b*x**3)/x - 2*sqrt(a + b*x**3)*(105*c*x - 105*d*x**2 - 35*e*x**3 - 21*f*x**4 - 15*g*x**5)/(105*x**3) + 3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*(10 - 10*sqrt(3))*(2*a*g + 7*b*d) + 7*b**(sympy.S(1)/3)*(4*a*f + 5*b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(70*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + sqrt(a + b*x**3)*(6*a*g + 21*b*d)/(7*b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_452():
    f = sqrt(a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4)/x**4
    F = -3*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*e*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(2*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 3*b**(sympy.S(1)/3)*e*sqrt(a + b*x**3)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x) + c*sqrt(a + b*x**3)/(3*x**3) + 3*d*sqrt(a + b*x**3)/(2*x**2) - 3*e*sqrt(a + b*x**3)/x - 2*sqrt(a + b*x**3)*(5*c*x + 15*d*x**2 - 15*e*x**3 - 5*f*x**4 - 3*g*x**5)/(15*x**4) + 3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*e*(10 - 10*sqrt(3)) + 4*a*g + 5*b*d)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(10*b**(sympy.S(1)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - (2*a*f + b*c)*atanh(sqrt(a + b*x**3)/sqrt(a))/(3*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_453():
    f = sqrt(a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4)/x**5
    F = 3*c*sqrt(a + b*x**3)/(20*x**4) + d*sqrt(a + b*x**3)/(3*x**3) + 3*e*sqrt(a + b*x**3)/(2*x**2) - 2*sqrt(a + b*x**3)*(3*c*x + 5*d*x**2 + 15*e*x**3 - 15*f*x**4 - 5*g*x**5)/(15*x**5) + 3*b**(sympy.S(1)/3)*sqrt(a + b*x**3)*(8*a*f + b*c)/(8*a*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) - sqrt(a + b*x**3)*(24*a*f + 3*b*c)/(8*a*x) - (2*a*g + b*d)*atanh(sqrt(a + b*x**3)/sqrt(a))/(3*sqrt(a)) - 3*3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(8*a*f + b*c)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(16*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(4*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*e - (1 - sqrt(3))*(8*a*f + b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(8*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_454():
    f = sqrt(a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4)/x**6
    F = sqrt(a + b*x**3)*(-c/(5*x**5) - d/(4*x**4) - e/(3*x**3) - f/(2*x**2) - g/x) - 3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*(5 - 5*sqrt(3))*(8*a*g + b*d) + 2*b**(sympy.S(1)/3)*(-10*a*f + b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(40*a*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 3*b**(sympy.S(1)/3)*sqrt(a + b*x**3)*(8*a*g + b*d)/(8*a*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) - 3*b*c*sqrt(a + b*x**3)/(20*a*x**2) - 3*b*d*sqrt(a + b*x**3)/(8*a*x) - b*e*atanh(sqrt(a + b*x**3)/sqrt(a))/(3*sqrt(a)) - 3*3**(sympy.S(1)/4)*b**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(8*a*g + b*d)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(16*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_455():
    f = sqrt(a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4)/x**7
    F = sqrt(a + b*x**3)*(-c/(6*x**6) - d/(5*x**5) - e/(4*x**4) - f/(3*x**3) - g/(2*x**2)) + 3*b**(sympy.S(4)/3)*e*sqrt(a + b*x**3)/(8*a*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) - 3**(sympy.S(3)/4)*b**(sympy.S(2)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*e*(5 - 5*sqrt(3)) - 20*a*g + 2*b*d)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(40*a*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - b*c*sqrt(a + b*x**3)/(12*a*x**3) - 3*b*d*sqrt(a + b*x**3)/(20*a*x**2) - 3*b*e*sqrt(a + b*x**3)/(8*a*x) + b*(-4*a*f + b*c)*atanh(sqrt(a + b*x**3)/sqrt(a))/(12*a**(sympy.S(3)/2)) - 3*3**(sympy.S(1)/4)*b**(sympy.S(4)/3)*e*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(16*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_456():
    f = sqrt(a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4)/x**8
    F = sqrt(a + b*x**3)*(-c/(7*x**7) - d/(6*x**6) - e/(5*x**5) - f/(4*x**4) - g/(3*x**3)) - 3*b*c*sqrt(a + b*x**3)/(56*a*x**4) - b*d*sqrt(a + b*x**3)/(12*a*x**3) - 3*b*e*sqrt(a + b*x**3)/(20*a*x**2) - 3*b**(sympy.S(4)/3)*sqrt(a + b*x**3)*(-14*a*f + 5*b*c)/(112*a**2*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 3*b*sqrt(a + b*x**3)*(-14*a*f + 5*b*c)/(112*a**2*x) + b*(-4*a*g + b*d)*atanh(sqrt(a + b*x**3)/sqrt(a))/(12*a**(sympy.S(3)/2)) + 3*3**(sympy.S(1)/4)*b**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-14*a*f + 5*b*c)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(224*a**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 3**(sympy.S(3)/4)*b**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(28*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*e - (5 - 5*sqrt(3))*(-14*a*f + 5*b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(560*a**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_457():
    f = sqrt(a + b*x**3)*(c + d*x + e*x**2 + f*x**3 + g*x**4)/x**9
    F = sqrt(a + b*x**3)*(-c/(8*x**8) - d/(7*x**7) - e/(6*x**6) - f/(5*x**5) - g/(4*x**4)) - 3*b*c*sqrt(a + b*x**3)/(80*a*x**5) - 3*b*d*sqrt(a + b*x**3)/(56*a*x**4) - b*e*sqrt(a + b*x**3)/(12*a*x**3) + 3**(sympy.S(3)/4)*b**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*(20 - 20*sqrt(3))*(-14*a*g + 5*b*d) + 7*b**(sympy.S(1)/3)*(-16*a*f + 7*b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(2240*a**2*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 3*b**(sympy.S(4)/3)*sqrt(a + b*x**3)*(-14*a*g + 5*b*d)/(112*a**2*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 3*b*sqrt(a + b*x**3)*(-14*a*g + 5*b*d)/(112*a**2*x) + 3*b*sqrt(a + b*x**3)*(-16*a*f + 7*b*c)/(320*a**2*x**2) + b**2*e*atanh(sqrt(a + b*x**3)/sqrt(a))/(12*a**(sympy.S(3)/2)) + 3*3**(sympy.S(1)/4)*b**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-14*a*g + 5*b*d)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(224*a**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_458():
    f = x**3*(a + b*x**3)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3 + g*x**4)
    F = 108*3**(sympy.S(1)/4)*a**(sympy.S(10)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-2*a*g + 5*b*d)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(8645*b**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 4*a**3*e*sqrt(a + b*x**3)/(105*b**2) - 36*3**(sympy.S(3)/4)*a**3*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*(8602 - 8602*sqrt(3))*(-2*a*g + 5*b*d) + 1729*b**(sympy.S(1)/3)*(-8*a*f + 23*b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(37182145*b**(sympy.S(8)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 216*a**3*sqrt(a + b*x**3)*(-2*a*g + 5*b*d)/(8645*b**(sympy.S(8)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*a**2*e*x**3*sqrt(a + b*x**3)/(105*b) + 54*a**2*f*x**4*sqrt(a + b*x**3)/(4301*b) + 54*a**2*g*x**5*sqrt(a + b*x**3)/(6175*b) + 54*a**2*x**2*sqrt(a + b*x**3)*(-2*a*g + 5*b*d)/(8645*b**2) + 54*a**2*x*sqrt(a + b*x**3)*(-8*a*f + 23*b*c)/(21505*b**2) + 2*a*x**3*sqrt(a + b*x**3)*(8947575*c*x + 6774075*d*x**2 + 5311735*e*x**3 + 4279275*f*x**4 + 3522519*g*x**5)/185910725 + 2*x**3*(a + b*x**3)**(sympy.S(3)/2)*(229425*c*x + 205275*d*x**2 + 185725*e*x**3 + 169575*f*x**4 + 156009*g*x**5)/3900225
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_459():
    f = x**2*(a + b*x**3)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3 + g*x**4)
    F = 108*3**(sympy.S(1)/4)*a**(sympy.S(10)/3)*e*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1729*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 216*a**3*e*sqrt(a + b*x**3)/(1729*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 36*3**(sympy.S(3)/4)*a**3*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*e*(43010 - 43010*sqrt(3)) + 13832*a*g - 39767*b*d)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(37182145*b**(sympy.S(7)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 54*a**2*e*x**2*sqrt(a + b*x**3)/(1729*b) + 2*a**2*f*x**3*sqrt(a + b*x**3)/(105*b) + 54*a**2*g*x**4*sqrt(a + b*x**3)/(4301*b) + 54*a**2*x*sqrt(a + b*x**3)*(-8*a*g + 23*b*d)/(21505*b**2) + 2*a**2*sqrt(a + b*x**3)*(-2*a*f + 7*b*c)/(105*b**2) + 2*a*x**2*sqrt(a + b*x**3)*(7436429*c*x + 5368545*d*x**2 + 4064445*e*x**3 + 3187041*f*x**4 + 2567565*g*x**5)/111546435 + 2*x**2*(a + b*x**3)**(sympy.S(3)/2)*(52003*c*x + 45885*d*x**2 + 41055*e*x**3 + 37145*f*x**4 + 33915*g*x**5)/780045
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_460():
    f = x*(a + b*x**3)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3 + g*x**4)
    F = -27*3**(sympy.S(1)/4)*a**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-4*a*f + 19*b*c)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1729*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 18*3**(sympy.S(3)/4)*a**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(3458*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*e + (935 - 935*sqrt(3))*(-4*a*f + 19*b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1616615*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 54*a**2*e*x*sqrt(a + b*x**3)/(935*b) + 54*a**2*f*x**2*sqrt(a + b*x**3)/(1729*b) + 2*a**2*g*x**3*sqrt(a + b*x**3)/(105*b) + 2*a**2*sqrt(a + b*x**3)*(-2*a*g + 7*b*d)/(105*b**2) + 54*a**2*sqrt(a + b*x**3)*(-4*a*f + 19*b*c)/(1729*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*a*x*sqrt(a + b*x**3)*(479655*c*x + 323323*d*x**2 + 233415*e*x**3 + 176715*f*x**4 + 138567*g*x**5)/4849845 + 2*x*(a + b*x**3)**(sympy.S(3)/2)*(33915*c*x + 29393*d*x**2 + 25935*e*x**3 + 23205*f*x**4 + 20995*g*x**5)/440895
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_461():
    f = (a + b*x**3)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3 + g*x**4)
    F = -27*3**(sympy.S(1)/4)*a**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-4*a*g + 19*b*d)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1729*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*a**2*e*sqrt(a + b*x**3)/(15*b) + 54*a**2*f*x*sqrt(a + b*x**3)/(935*b) + 54*a**2*g*x**2*sqrt(a + b*x**3)/(1729*b) + 18*3**(sympy.S(3)/4)*a**2*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*(935 - 935*sqrt(3))*(-4*a*g + 19*b*d) + 1729*b**(sympy.S(1)/3)*(-2*a*f + 17*b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(1616615*b**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 54*a**2*sqrt(a + b*x**3)*(-4*a*g + 19*b*d)/(1729*b**(sympy.S(5)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*a*sqrt(a + b*x**3)*(793611*c*x + 479655*d*x**2 + 323323*e*x**3 + 233415*f*x**4 + 176715*g*x**5)/4849845 + 2*(a + b*x**3)**(sympy.S(3)/2)*(62985*c*x + 53295*d*x**2 + 46189*e*x**3 + 40755*f*x**4 + 36465*g*x**5)/692835
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_462():
    f = (a + b*x**3)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3 + g*x**4)/x
    F = -27*3**(sympy.S(1)/4)*a**(sympy.S(7)/3)*e*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(91*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 2*a**(sympy.S(3)/2)*c*atanh(sqrt(a + b*x**3)/sqrt(a))/3 + 2*a**2*f*sqrt(a + b*x**3)/(15*b) + 54*a**2*g*x*sqrt(a + b*x**3)/(935*b) + 54*a**2*e*sqrt(a + b*x**3)/(91*b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 18*3**(sympy.S(3)/4)*a**2*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*e*(935 - 935*sqrt(3)) - 182*a*g + 1547*b*d)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(85085*b**(sympy.S(4)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*a*sqrt(a + b*x**3)*(85085*c*x + 41769*d*x**2 + 25245*e*x**3 + 17017*f*x**4 + 12285*g*x**5)/(255255*x) + 2*(a + b*x**3)**(sympy.S(3)/2)*(12155*c*x + 9945*d*x**2 + 8415*e*x**3 + 7293*f*x**4 + 6435*g*x**5)/(109395*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_463():
    f = (a + b*x**3)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3 + g*x**4)/x**2
    F = -27*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(2*a*f + 13*b*c)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(182*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 9*3**(sympy.S(3)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(182*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*e - (55 - 55*sqrt(3))*(2*a*f + 13*b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(5005*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 2*a**(sympy.S(3)/2)*d*atanh(sqrt(a + b*x**3)/sqrt(a))/3 + 2*a**2*g*sqrt(a + b*x**3)/(15*b) - 27*a*c*sqrt(a + b*x**3)/(7*x) + 2*a*sqrt(a + b*x**3)*(19305*c*x + 5005*d*x**2 + 2457*e*x**3 + 1485*f*x**4 + 1001*g*x**5)/(15015*x**2) + 27*a*sqrt(a + b*x**3)*(2*a*f + 13*b*c)/(91*b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*(a + b*x**3)**(sympy.S(3)/2)*(6435*c*x + 5005*d*x**2 + 4095*e*x**3 + 3465*f*x**4 + 3003*g*x**5)/(45045*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_464():
    f = (a + b*x**3)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3 + g*x**4)/x**3
    F = -27*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(2*a*g + 13*b*d)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(182*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 2*a**(sympy.S(3)/2)*e*atanh(sqrt(a + b*x**3)/sqrt(a))/3 + 27*a*c*sqrt(a + b*x**3)/(10*x**2) - 27*a*d*sqrt(a + b*x**3)/(7*x) - 2*a*sqrt(a + b*x**3)*(27027*c*x - 19305*d*x**2 - 5005*e*x**3 - 2457*f*x**4 - 1485*g*x**5)/(15015*x**3) + 9*3**(sympy.S(3)/4)*a*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*(110 - 110*sqrt(3))*(2*a*g + 13*b*d) + 91*b**(sympy.S(1)/3)*(4*a*f + 11*b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(10010*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 27*a*sqrt(a + b*x**3)*(2*a*g + 13*b*d)/(91*b**(sympy.S(2)/3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 2*(a + b*x**3)**(sympy.S(3)/2)*(9009*c*x + 6435*d*x**2 + 5005*e*x**3 + 4095*f*x**4 + 3465*g*x**5)/(45045*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_465():
    f = (a + b*x**3)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3 + g*x**4)/x**4
    F = -27*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*b**(sympy.S(1)/3)*e*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(14*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - sqrt(a)*(2*a*f + 3*b*c)*atanh(sqrt(a + b*x**3)/sqrt(a))/3 + 27*a*b**(sympy.S(1)/3)*e*sqrt(a + b*x**3)/(7*a**(sympy.S(1)/3)*(1 + sqrt(3)) + 7*b**(sympy.S(1)/3)*x) + a*c*sqrt(a + b*x**3)/x**3 + 27*a*d*sqrt(a + b*x**3)/(10*x**2) - 27*a*e*sqrt(a + b*x**3)/(7*x) - 2*a*sqrt(a + b*x**3)*(1155*c*x + 2079*d*x**2 - 1485*e*x**3 - 385*f*x**4 - 189*g*x**5)/(1155*x**4) + 9*3**(sympy.S(3)/4)*a*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*e*(110 - 110*sqrt(3)) + 28*a*g + 77*b*d)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(770*b**(sympy.S(1)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 2*(a + b*x**3)**(sympy.S(3)/2)*(1155*c*x + 693*d*x**2 + 495*e*x**3 + 385*f*x**4 + 315*g*x**5)/(3465*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_466():
    f = (a + b*x**3)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3 + g*x**4)/x**5
    F = -27*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(8*a*f + 7*b*c)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(112*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 9*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(28*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*e - (5 - 5*sqrt(3))*(8*a*f + 7*b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(280*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - sqrt(a)*(2*a*g + 3*b*d)*atanh(sqrt(a + b*x**3)/sqrt(a))/3 + 27*a*c*sqrt(a + b*x**3)/(20*x**4) + a*d*sqrt(a + b*x**3)/x**3 + 27*a*e*sqrt(a + b*x**3)/(10*x**2) - 2*a*sqrt(a + b*x**3)*(189*c*x + 105*d*x**2 + 189*e*x**3 - 135*f*x**4 - 35*g*x**5)/(105*x**5) + 27*b**(sympy.S(1)/3)*sqrt(a + b*x**3)*(8*a*f + 7*b*c)/(56*a**(sympy.S(1)/3)*(1 + sqrt(3)) + 56*b**(sympy.S(1)/3)*x) - sqrt(a + b*x**3)*(216*a*f + 189*b*c)/(56*x) + 2*(a + b*x**3)**(sympy.S(3)/2)*(315*c*x + 105*d*x**2 + 63*e*x**3 + 45*f*x**4 + 35*g*x**5)/(315*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_467():
    f = (a + b*x**3)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3 + g*x**4)/x**6
    F = -27*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(8*a*g + 7*b*d)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(112*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - sqrt(a)*b*e*atanh(sqrt(a + b*x**3)/sqrt(a)) + 9*3**(sympy.S(3)/4)*b**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*(5 - 5*sqrt(3))*(8*a*g + 7*b*d) + 14*b**(sympy.S(1)/3)*(2*a*f + b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(280*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 27*b**(sympy.S(1)/3)*sqrt(a + b*x**3)*(8*a*g + 7*b*d)/(56*a**(sympy.S(1)/3)*(1 + sqrt(3)) + 56*b**(sympy.S(1)/3)*x) + 27*b*c*sqrt(a + b*x**3)/(20*x**2) - 27*b*d*sqrt(a + b*x**3)/(8*x) - b*sqrt(a + b*x**3)*(252*c*x - 315*d*x**2 - 140*e*x**3 - 126*f*x**4 - 180*g*x**5)/(140*x**3) - (a + b*x**3)**(sympy.S(3)/2)*(c/(5*x**5) + d/(4*x**4) + e/(3*x**3) + f/(2*x**2) + g/x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_468():
    f = (a + b*x**3)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3 + g*x**4)/x**7
    F = -27*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*b**(sympy.S(4)/3)*e*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(16*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 27*b**(sympy.S(4)/3)*e*sqrt(a + b*x**3)/(8*a**(sympy.S(1)/3)*(1 + sqrt(3)) + 8*b**(sympy.S(1)/3)*x) + 9*3**(sympy.S(3)/4)*b**(sympy.S(2)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*e*(5 - 5*sqrt(3)) + 4*a*g + 2*b*d)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(40*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + b*c*sqrt(a + b*x**3)/(4*x**3) + 27*b*d*sqrt(a + b*x**3)/(20*x**2) - 27*b*e*sqrt(a + b*x**3)/(8*x) - b*sqrt(a + b*x**3)*(10*c*x + 36*d*x**2 - 45*e*x**3 - 20*f*x**4 - 18*g*x**5)/(20*x**4) - (a + b*x**3)**(sympy.S(3)/2)*(c/(6*x**6) + d/(5*x**5) + e/(4*x**4) + f/(3*x**3) + g/(2*x**2)) - b*(4*a*f + b*c)*atanh(sqrt(a + b*x**3)/sqrt(a))/(4*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_469():
    f = (a + b*x**3)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3 + g*x**4)/x**8
    F = 27*b*c*sqrt(a + b*x**3)/(280*x**4) + b*d*sqrt(a + b*x**3)/(4*x**3) + 27*b*e*sqrt(a + b*x**3)/(20*x**2) - b*sqrt(a + b*x**3)*(36*c*x + 70*d*x**2 + 252*e*x**3 - 315*f*x**4 - 140*g*x**5)/(140*x**5) - (a + b*x**3)**(sympy.S(3)/2)*(c/(7*x**7) + d/(6*x**6) + e/(5*x**5) + f/(4*x**4) + g/(3*x**3)) + 27*b**(sympy.S(4)/3)*sqrt(a + b*x**3)*(14*a*f + b*c)/(112*a*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) - 27*b*sqrt(a + b*x**3)*(14*a*f + b*c)/(112*a*x) - b*(4*a*g + b*d)*atanh(sqrt(a + b*x**3)/sqrt(a))/(4*sqrt(a)) - 27*3**(sympy.S(1)/4)*b**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(14*a*f + b*c)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(224*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 9*3**(sympy.S(3)/4)*b**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(28*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*e - (5 - 5*sqrt(3))*(14*a*f + b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(560*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_470():
    f = (a + b*x**3)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3 + g*x**4)/x**9
    F = -b*sqrt(a + b*x**3)*(63*c/x**5 + 90*d/x**4 + 140*e/x**3 + 252*f/x**2 + 630*g/x)/560 - (a + b*x**3)**(sympy.S(3)/2)*(c/(8*x**8) + d/(7*x**7) + e/(6*x**6) + f/(5*x**5) + g/(4*x**4)) - 9*3**(sympy.S(3)/4)*b**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*(20 - 20*sqrt(3))*(14*a*g + b*d) + 7*b**(sympy.S(1)/3)*(-16*a*f + b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(2240*a*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 27*b**(sympy.S(4)/3)*sqrt(a + b*x**3)*(14*a*g + b*d)/(112*a*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) - 27*b**2*c*sqrt(a + b*x**3)/(320*a*x**2) - 27*b**2*d*sqrt(a + b*x**3)/(112*a*x) - b**2*e*atanh(sqrt(a + b*x**3)/sqrt(a))/(4*sqrt(a)) - 27*3**(sympy.S(1)/4)*b**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(14*a*g + b*d)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(224*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_471():
    f = (a + b*x**3)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3 + g*x**4)/x**10
    F = -b*sqrt(a + b*x**3)*(140*c/x**6 + 189*d/x**5 + 270*e/x**4 + 420*f/x**3 + 756*g/x**2)/1680 - (a + b*x**3)**(sympy.S(3)/2)*(280*c/x**9 + 315*d/x**8 + 360*e/x**7 + 420*f/x**6 + 504*g/x**5)/2520 + 27*b**(sympy.S(7)/3)*e*sqrt(a + b*x**3)/(112*a*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) - 9*3**(sympy.S(3)/4)*b**(sympy.S(5)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*e*(20 - 20*sqrt(3)) - 112*a*g + 7*b*d)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(2240*a*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - b**2*c*sqrt(a + b*x**3)/(24*a*x**3) - 27*b**2*d*sqrt(a + b*x**3)/(320*a*x**2) - 27*b**2*e*sqrt(a + b*x**3)/(112*a*x) + b**2*(-6*a*f + b*c)*atanh(sqrt(a + b*x**3)/sqrt(a))/(24*a**(sympy.S(3)/2)) - 27*3**(sympy.S(1)/4)*b**(sympy.S(7)/3)*e*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(224*a**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_472():
    f = (a + b*x**3)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3 + g*x**4)/x**11
    F = -b*sqrt(a + b*x**3)*(108*c/x**7 + 140*d/x**6 + 189*e/x**5 + 270*f/x**4 + 420*g/x**3)/1680 - (a + b*x**3)**(sympy.S(3)/2)*(252*c/x**10 + 280*d/x**9 + 315*e/x**8 + 360*f/x**7 + 420*g/x**6)/2520 - 27*b**2*c*sqrt(a + b*x**3)/(1120*a*x**4) - b**2*d*sqrt(a + b*x**3)/(24*a*x**3) - 27*b**2*e*sqrt(a + b*x**3)/(320*a*x**2) - 27*b**(sympy.S(7)/3)*sqrt(a + b*x**3)*(-4*a*f + b*c)/(448*a**2*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 27*b**2*sqrt(a + b*x**3)*(-4*a*f + b*c)/(448*a**2*x) + b**2*(-6*a*g + b*d)*atanh(sqrt(a + b*x**3)/sqrt(a))/(24*a**(sympy.S(3)/2)) + 27*3**(sympy.S(1)/4)*b**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-4*a*f + b*c)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(896*a**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 9*3**(sympy.S(3)/4)*b**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(7*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*e - (5 - 5*sqrt(3))*(-4*a*f + b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(2240*a**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_473():
    f = (a + b*x**3)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3 + g*x**4)/x**12
    F = -b*sqrt(a + b*x**3)*(945*c/x**8 + 1188*d/x**7 + 1540*e/x**6 + 2079*f/x**5 + 2970*g/x**4)/18480 - (a + b*x**3)**(sympy.S(3)/2)*(2520*c/x**11 + 2772*d/x**10 + 3080*e/x**9 + 3465*f/x**8 + 3960*g/x**7)/27720 - 27*b**2*c*sqrt(a + b*x**3)/(1760*a*x**5) - 27*b**2*d*sqrt(a + b*x**3)/(1120*a*x**4) - b**2*e*sqrt(a + b*x**3)/(24*a*x**3) + 9*3**(sympy.S(3)/4)*b**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*(110 - 110*sqrt(3))*(-4*a*g + b*d) + 7*b**(sympy.S(1)/3)*(-22*a*f + 7*b*c))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(49280*a**2*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 27*b**(sympy.S(7)/3)*sqrt(a + b*x**3)*(-4*a*g + b*d)/(448*a**2*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)) + 27*b**2*sqrt(a + b*x**3)*(-4*a*g + b*d)/(448*a**2*x) + 27*b**2*sqrt(a + b*x**3)*(-22*a*f + 7*b*c)/(7040*a**2*x**2) + b**3*e*atanh(sqrt(a + b*x**3)/sqrt(a))/(24*a**(sympy.S(3)/2)) + 27*3**(sympy.S(1)/4)*b**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-4*a*g + b*d)*elliptic_e(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(896*a**(sympy.S(5)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_474():
    f = (a + b*x**4)*(c + d*x + e*x**2 + f*x**3)
    F = a*c*x + a*d*x**2/2 + a*e*x**3/3 + a*f*x**4/4 + b*c*x**5/5 + b*d*x**6/6 + b*e*x**7/7 + b*f*x**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_475():
    f = x**3*(a + b*x**4)*(c + d*x + e*x**2 + f*x**3)
    F = a*c*x**4/4 + a*d*x**5/5 + a*e*x**6/6 + a*f*x**7/7 + b*c*x**8/8 + b*d*x**9/9 + b*e*x**10/10 + b*f*x**11/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_476():
    f = (a + b*x**4)**2*(c + d*x + e*x**2 + f*x**3)
    F = a**2*c*x + a**2*d*x**2/2 + a**2*e*x**3/3 + 2*a*b*c*x**5/5 + a*b*d*x**6/3 + 2*a*b*e*x**7/7 + b**2*c*x**9/9 + b**2*d*x**10/10 + b**2*e*x**11/11 + f*(a + b*x**4)**3/(12*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_477():
    f = x**3*(a + b*x**4)**2*(c + d*x + e*x**2 + f*x**3)
    F = a**2*d*x**5/5 + a**2*e*x**6/6 + a**2*f*x**7/7 + 2*a*b*d*x**9/9 + a*b*e*x**10/5 + 2*a*b*f*x**11/11 + b**2*d*x**13/13 + b**2*e*x**14/14 + b**2*f*x**15/15 + c*(a + b*x**4)**3/(12*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_478():
    f = (a + b*x**4)**3*(c + d*x + e*x**2 + f*x**3)
    F = a**3*c*x + a**3*d*x**2/2 + a**3*e*x**3/3 + 3*a**2*b*c*x**5/5 + a**2*b*d*x**6/2 + 3*a**2*b*e*x**7/7 + a*b**2*c*x**9/3 + 3*a*b**2*d*x**10/10 + 3*a*b**2*e*x**11/11 + b**3*c*x**13/13 + b**3*d*x**14/14 + b**3*e*x**15/15 + f*(a + b*x**4)**4/(16*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_479():
    f = x**3*(a + b*x**4)**3*(c + d*x + e*x**2 + f*x**3)
    F = a**3*d*x**5/5 + a**3*e*x**6/6 + a**3*f*x**7/7 + a**2*b*d*x**9/3 + 3*a**2*b*e*x**10/10 + 3*a**2*b*f*x**11/11 + 3*a*b**2*d*x**13/13 + 3*a*b**2*e*x**14/14 + a*b**2*f*x**15/5 + b**3*d*x**17/17 + b**3*e*x**18/18 + b**3*f*x**19/19 + c*(a + b*x**4)**4/(16*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_480():
    f = (a + b*x**4)**4*(c + d*x + e*x**2 + f*x**3)
    F = a**4*c*x + a**4*d*x**2/2 + a**4*e*x**3/3 + 4*a**3*b*c*x**5/5 + 2*a**3*b*d*x**6/3 + 4*a**3*b*e*x**7/7 + 2*a**2*b**2*c*x**9/3 + 3*a**2*b**2*d*x**10/5 + 6*a**2*b**2*e*x**11/11 + 4*a*b**3*c*x**13/13 + 2*a*b**3*d*x**14/7 + 4*a*b**3*e*x**15/15 + b**4*c*x**17/17 + b**4*d*x**18/18 + b**4*e*x**19/19 + f*(a + b*x**4)**5/(20*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_481():
    f = x**3*(a + b*x**4)**4*(c + d*x + e*x**2 + f*x**3)
    F = a**4*d*x**5/5 + a**4*e*x**6/6 + a**4*f*x**7/7 + 4*a**3*b*d*x**9/9 + 2*a**3*b*e*x**10/5 + 4*a**3*b*f*x**11/11 + 6*a**2*b**2*d*x**13/13 + 3*a**2*b**2*e*x**14/7 + 2*a**2*b**2*f*x**15/5 + 4*a*b**3*d*x**17/17 + 2*a*b**3*e*x**18/9 + 4*a*b**3*f*x**19/19 + b**4*d*x**21/21 + b**4*e*x**22/22 + b**4*f*x**23/23 + c*(a + b*x**4)**5/(20*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_482():
    f = (c + d*x + e*x**2 + f*x**3)/(a - b*x**4)
    F = -f*log(a - b*x**4)/(4*b) + d*atanh(sqrt(b)*x**2/sqrt(a))/(2*sqrt(a)*sqrt(b)) + (-sqrt(a)*e + sqrt(b)*c)*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)) + (sqrt(a)*e + sqrt(b)*c)*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_483():
    f = x**3*(c + d*x + e*x**2 + f*x**3)/(a - b*x**4)
    F = a**(sympy.S(1)/4)*(-sqrt(a)*f + sqrt(b)*d)*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*b**(sympy.S(7)/4)) + a**(sympy.S(1)/4)*(sqrt(a)*f + sqrt(b)*d)*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*b**(sympy.S(7)/4)) + sqrt(a)*e*atanh(sqrt(b)*x**2/sqrt(a))/(2*b**(sympy.S(3)/2)) - c*log(a - b*x**4)/(4*b) - d*x/b - e*x**2/(2*b) - f*x**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_484():
    f = (c + d*x + e*x**2 + f*x**3)/(a + b*x**4)
    F = f*log(a + b*x**4)/(4*b) + d*atan(sqrt(b)*x**2/sqrt(a))/(2*sqrt(a)*sqrt(b)) - sqrt(2)*(-sqrt(a)*e + sqrt(b)*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(-sqrt(a)*e + sqrt(b)*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)) - sqrt(2)*(sqrt(a)*e + sqrt(b)*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(sqrt(a)*e + sqrt(b)*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_485():
    f = x**3*(c + d*x + e*x**2 + f*x**3)/(a + b*x**4)
    F = sqrt(2)*a**(sympy.S(1)/4)*(-sqrt(a)*f + sqrt(b)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*b**(sympy.S(7)/4)) - sqrt(2)*a**(sympy.S(1)/4)*(-sqrt(a)*f + sqrt(b)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*b**(sympy.S(7)/4)) + sqrt(2)*a**(sympy.S(1)/4)*(sqrt(a)*f + sqrt(b)*d)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*b**(sympy.S(7)/4)) - sqrt(2)*a**(sympy.S(1)/4)*(sqrt(a)*f + sqrt(b)*d)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*b**(sympy.S(7)/4)) - sqrt(a)*e*atan(sqrt(b)*x**2/sqrt(a))/(2*b**(sympy.S(3)/2)) + c*log(a + b*x**4)/(4*b) + d*x/b + e*x**2/(2*b) + f*x**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_486():
    f = (c + d*x + e*x**2 + f*x**3)/(a + b*x**4)**2
    F = -(a*f - b*x*(c + d*x + e*x**2))/(4*a*b*(a + b*x**4)) + d*atan(sqrt(b)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*sqrt(b)) - sqrt(2)*(-sqrt(a)*e + 3*sqrt(b)*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(-sqrt(a)*e + 3*sqrt(b)*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(3)/4)) - sqrt(2)*(sqrt(a)*e + 3*sqrt(b)*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(sqrt(a)*e + 3*sqrt(b)*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_487():
    f = x**3*(c + d*x + e*x**2 + f*x**3)/(a + b*x**4)**2
    F = -(c + d*x + e*x**2 + f*x**3)/(4*b*(a + b*x**4)) + e*atan(sqrt(b)*x**2/sqrt(a))/(4*sqrt(a)*b**(sympy.S(3)/2)) - sqrt(2)*(-3*sqrt(a)*f + sqrt(b)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(3)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(-3*sqrt(a)*f + sqrt(b)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(3)/4)*b**(sympy.S(7)/4)) - sqrt(2)*(3*sqrt(a)*f + sqrt(b)*d)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(3)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(3*sqrt(a)*f + sqrt(b)*d)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(3)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_488():
    f = (c + d*x + e*x**2 + f*x**3)/(a + b*x**4)**3
    F = -(a*f - b*x*(c + d*x + e*x**2))/(8*a*b*(a + b*x**4)**2) + x*(7*c + 6*d*x + 5*e*x**2)/(32*a**2*(a + b*x**4)) + 3*d*atan(sqrt(b)*x**2/sqrt(a))/(16*a**(sympy.S(5)/2)*sqrt(b)) - sqrt(2)*(-5*sqrt(a)*e + 21*sqrt(b)*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(256*a**(sympy.S(11)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(-5*sqrt(a)*e + 21*sqrt(b)*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(256*a**(sympy.S(11)/4)*b**(sympy.S(3)/4)) - sqrt(2)*(5*sqrt(a)*e + 21*sqrt(b)*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(5*sqrt(a)*e + 21*sqrt(b)*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_489():
    f = x**3*(c + d*x + e*x**2 + f*x**3)/(a + b*x**4)**3
    F = -(c + d*x + e*x**2 + f*x**3)/(8*b*(a + b*x**4)**2) + x*(d + 2*e*x + 3*f*x**2)/(32*a*b*(a + b*x**4)) + e*atan(sqrt(b)*x**2/sqrt(a))/(16*a**(sympy.S(3)/2)*b**(sympy.S(3)/2)) - sqrt(2)*(-3*sqrt(a)*f + 3*sqrt(b)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(256*a**(sympy.S(7)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(-3*sqrt(a)*f + 3*sqrt(b)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(256*a**(sympy.S(7)/4)*b**(sympy.S(7)/4)) - sqrt(2)*(3*sqrt(a)*f + 3*sqrt(b)*d)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(7)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(3*sqrt(a)*f + 3*sqrt(b)*d)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(7)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_490():
    f = (c + d*x + e*x**2 + f*x**3)/(a + b*x**4)**4
    F = -(a*f - b*x*(c + d*x + e*x**2))/(12*a*b*(a + b*x**4)**3) + x*(11*c + 10*d*x + 9*e*x**2)/(96*a**2*(a + b*x**4)**2) + x*(77*c + 60*d*x + 45*e*x**2)/(384*a**3*(a + b*x**4)) + 5*d*atan(sqrt(b)*x**2/sqrt(a))/(32*a**(sympy.S(7)/2)*sqrt(b)) - sqrt(2)*(-15*sqrt(a)*e + 77*sqrt(b)*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(1024*a**(sympy.S(15)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(-15*sqrt(a)*e + 77*sqrt(b)*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(1024*a**(sympy.S(15)/4)*b**(sympy.S(3)/4)) - sqrt(2)*(15*sqrt(a)*e + 77*sqrt(b)*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(512*a**(sympy.S(15)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(15*sqrt(a)*e + 77*sqrt(b)*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(512*a**(sympy.S(15)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_491():
    f = x**3*(c + d*x + e*x**2 + f*x**3)/(a + b*x**4)**4
    F = -(c + d*x + e*x**2 + f*x**3)/(12*b*(a + b*x**4)**3) + x*(d + 2*e*x + 3*f*x**2)/(96*a*b*(a + b*x**4)**2) + x*(7*d + 12*e*x + 15*f*x**2)/(384*a**2*b*(a + b*x**4)) + e*atan(sqrt(b)*x**2/sqrt(a))/(32*a**(sympy.S(5)/2)*b**(sympy.S(3)/2)) - sqrt(2)*(-5*sqrt(a)*f + 7*sqrt(b)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(1024*a**(sympy.S(11)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(-5*sqrt(a)*f + 7*sqrt(b)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(1024*a**(sympy.S(11)/4)*b**(sympy.S(7)/4)) - sqrt(2)*(5*sqrt(a)*f + 7*sqrt(b)*d)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(512*a**(sympy.S(11)/4)*b**(sympy.S(7)/4)) + sqrt(2)*(5*sqrt(a)*f + 7*sqrt(b)*d)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(512*a**(sympy.S(11)/4)*b**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_492():
    f = x**4*sqrt(a + b*x**4)*(c + d*x + e*x**2 + f*x**3)
    F = 2*a**(sympy.S(9)/4)*e*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(7)/4)*sqrt(a + b*x**4)) - a**(sympy.S(7)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(7*sqrt(a)*e + 5*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(105*b**(sympy.S(7)/4)*sqrt(a + b*x**4)) - a**2*d*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(16*b**(sympy.S(3)/2)) - 2*a**2*e*x*sqrt(a + b*x**4)/(15*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x**2)) + 2*a*c*x*sqrt(a + b*x**4)/(21*b) - a*d*x**2*sqrt(a + b*x**4)/(16*b) + 2*a*e*x**3*sqrt(a + b*x**4)/(45*b) + x**5*sqrt(a + b*x**4)*(9*c + 7*e*x**2)/63 + f*x**4*(a + b*x**4)**(sympy.S(3)/2)/(10*b) - (a + b*x**4)**(sympy.S(3)/2)*(8*a*f - 15*b*d*x**2)/(120*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_493():
    f = x**3*sqrt(a + b*x**4)*(c + d*x + e*x**2 + f*x**3)
    F = 2*a**(sympy.S(9)/4)*f*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(7)/4)*sqrt(a + b*x**4)) - a**(sympy.S(7)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(7*sqrt(a)*f + 5*sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(105*b**(sympy.S(7)/4)*sqrt(a + b*x**4)) - a**2*e*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(16*b**(sympy.S(3)/2)) - 2*a**2*f*x*sqrt(a + b*x**4)/(15*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x**2)) + 2*a*d*x*sqrt(a + b*x**4)/(21*b) - a*e*x**2*sqrt(a + b*x**4)/(16*b) + 2*a*f*x**3*sqrt(a + b*x**4)/(45*b) + x**5*sqrt(a + b*x**4)*(9*d + 7*f*x**2)/63 + (a + b*x**4)**(sympy.S(3)/2)*(4*c + 3*e*x**2)/(24*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_494():
    f = x**2*sqrt(a + b*x**4)*(c + d*x + e*x**2 + f*x**3)
    F = -2*a**(sympy.S(5)/4)*c*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(3)/4)*sqrt(a + b*x**4)) + a**(sympy.S(5)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(-5*sqrt(a)*e + 21*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(105*b**(sympy.S(5)/4)*sqrt(a + b*x**4)) - a**2*f*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(16*b**(sympy.S(3)/2)) + 2*a*e*x*sqrt(a + b*x**4)/(21*b) - a*f*x**2*sqrt(a + b*x**4)/(16*b) + 2*a*c*x*sqrt(a + b*x**4)/(5*sqrt(b)*(sqrt(a) + sqrt(b)*x**2)) + x**3*sqrt(a + b*x**4)*(7*c + 5*e*x**2)/35 + (a + b*x**4)**(sympy.S(3)/2)*(4*d + 3*f*x**2)/(24*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_495():
    f = x*sqrt(a + b*x**4)*(c + d*x + e*x**2 + f*x**3)
    F = -2*a**(sympy.S(5)/4)*d*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(3)/4)*sqrt(a + b*x**4)) + a**(sympy.S(5)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(-5*sqrt(a)*f + 21*sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(105*b**(sympy.S(5)/4)*sqrt(a + b*x**4)) + 2*a*f*x*sqrt(a + b*x**4)/(21*b) + a*c*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(4*sqrt(b)) + 2*a*d*x*sqrt(a + b*x**4)/(5*sqrt(b)*(sqrt(a) + sqrt(b)*x**2)) + c*x**2*sqrt(a + b*x**4)/4 + x**3*sqrt(a + b*x**4)*(7*d + 5*f*x**2)/35 + e*(a + b*x**4)**(sympy.S(3)/2)/(6*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_496():
    f = sqrt(a + b*x**4)*(c + d*x + e*x**2 + f*x**3)
    F = -2*a**(sympy.S(5)/4)*e*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(3)/4)*sqrt(a + b*x**4)) + a**(sympy.S(3)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(3*sqrt(a)*e + 5*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(3)/4)*sqrt(a + b*x**4)) + a*d*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(4*sqrt(b)) + 2*a*e*x*sqrt(a + b*x**4)/(5*sqrt(b)*(sqrt(a) + sqrt(b)*x**2)) + d*x**2*sqrt(a + b*x**4)/4 + x*sqrt(a + b*x**4)*(5*c + 3*e*x**2)/15 + f*(a + b*x**4)**(sympy.S(3)/2)/(6*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_497():
    f = sqrt(a + b*x**4)*(c + d*x + e*x**2 + f*x**3)/x
    F = -2*a**(sympy.S(5)/4)*f*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(3)/4)*sqrt(a + b*x**4)) + a**(sympy.S(3)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(3*sqrt(a)*f + 5*sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(3)/4)*sqrt(a + b*x**4)) - sqrt(a)*c*atanh(sqrt(a + b*x**4)/sqrt(a))/2 + a*e*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(4*sqrt(b)) + 2*a*f*x*sqrt(a + b*x**4)/(5*sqrt(b)*(sqrt(a) + sqrt(b)*x**2)) + x*sqrt(a + b*x**4)*(5*d + 3*f*x**2)/15 + sqrt(a + b*x**4)*(2*c + e*x**2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_498():
    f = sqrt(a + b*x**4)*(c + d*x + e*x**2 + f*x**3)/x**2
    F = -2*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*c*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/sqrt(a + b*x**4) + a**(sympy.S(1)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(sqrt(a)*e + 3*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(3*b**(sympy.S(1)/4)*sqrt(a + b*x**4)) - sqrt(a)*d*atanh(sqrt(a + b*x**4)/sqrt(a))/2 + a*f*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(4*sqrt(b)) + 2*sqrt(b)*c*x*sqrt(a + b*x**4)/(sqrt(a) + sqrt(b)*x**2) + sqrt(a + b*x**4)*(2*d + f*x**2)/4 - sqrt(a + b*x**4)*(3*c - e*x**2)/(3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_499():
    f = sqrt(a + b*x**4)*(c + d*x + e*x**2 + f*x**3)/x**3
    F = -2*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*d*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/sqrt(a + b*x**4) + a**(sympy.S(1)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(sqrt(a)*f + 3*sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(3*b**(sympy.S(1)/4)*sqrt(a + b*x**4)) - sqrt(a)*e*atanh(sqrt(a + b*x**4)/sqrt(a))/2 + sqrt(b)*c*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/2 + 2*sqrt(b)*d*x*sqrt(a + b*x**4)/(sqrt(a) + sqrt(b)*x**2) - sqrt(a + b*x**4)*(3*d - f*x**2)/(3*x) - sqrt(a + b*x**4)*(c - e*x**2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_500():
    f = sqrt(a + b*x**4)*(c + d*x + e*x**2 + f*x**3)/x**4
    F = -2*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*e*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/sqrt(a + b*x**4) - sqrt(a)*f*atanh(sqrt(a + b*x**4)/sqrt(a))/2 + sqrt(b)*d*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/2 + 2*sqrt(b)*e*x*sqrt(a + b*x**4)/(sqrt(a) + sqrt(b)*x**2) - 2*e*sqrt(a + b*x**4)/x - sqrt(a + b*x**4)*(d - f*x**2)/(2*x**2) - sqrt(a + b*x**4)*(c - 3*e*x**2)/(3*x**3) + b**(sympy.S(1)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(3*sqrt(a)*e + sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(3*a**(sympy.S(1)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_501():
    f = sqrt(a + b*x**4)*(c + d*x + e*x**2 + f*x**3)/x**5
    F = -2*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*f*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/sqrt(a + b*x**4) + sqrt(b)*e*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/2 + 2*sqrt(b)*f*x*sqrt(a + b*x**4)/(sqrt(a) + sqrt(b)*x**2) - sqrt(a + b*x**4)*(3*c/x**4 + 4*d/x**3 + 6*e/x**2 + 12*f/x)/12 - b*c*atanh(sqrt(a + b*x**4)/sqrt(a))/(4*sqrt(a)) + b**(sympy.S(1)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(3*sqrt(a)*f + sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(3*a**(sympy.S(1)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_502():
    f = sqrt(a + b*x**4)*(c + d*x + e*x**2 + f*x**3)/x**6
    F = sqrt(b)*f*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/2 - sqrt(a + b*x**4)*(12*c/x**5 + 15*d/x**4 + 20*e/x**3 + 30*f/x**2)/60 + 2*b**(sympy.S(3)/2)*c*x*sqrt(a + b*x**4)/(5*a*(sqrt(a) + sqrt(b)*x**2)) - 2*b*c*sqrt(a + b*x**4)/(5*a*x) - b*d*atanh(sqrt(a + b*x**4)/sqrt(a))/(4*sqrt(a)) - 2*b**(sympy.S(5)/4)*c*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*a**(sympy.S(3)/4)*sqrt(a + b*x**4)) + b**(sympy.S(3)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(5*sqrt(a)*e + 3*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(15*a**(sympy.S(3)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_503():
    f = sqrt(a + b*x**4)*(c + d*x + e*x**2 + f*x**3)/x**7
    F = -sqrt(a + b*x**4)*(10*c/x**6 + 12*d/x**5 + 15*e/x**4 + 20*f/x**3)/60 + 2*b**(sympy.S(3)/2)*d*x*sqrt(a + b*x**4)/(5*a*(sqrt(a) + sqrt(b)*x**2)) - b*c*sqrt(a + b*x**4)/(6*a*x**2) - 2*b*d*sqrt(a + b*x**4)/(5*a*x) - b*e*atanh(sqrt(a + b*x**4)/sqrt(a))/(4*sqrt(a)) - 2*b**(sympy.S(5)/4)*d*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*a**(sympy.S(3)/4)*sqrt(a + b*x**4)) + b**(sympy.S(3)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(5*sqrt(a)*f + 3*sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(15*a**(sympy.S(3)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_504():
    f = sqrt(a + b*x**4)*(c + d*x + e*x**2 + f*x**3)/x**8
    F = -sqrt(a + b*x**4)*(60*c/x**7 + 70*d/x**6 + 84*e/x**5 + 105*f/x**4)/420 + 2*b**(sympy.S(3)/2)*e*x*sqrt(a + b*x**4)/(5*a*(sqrt(a) + sqrt(b)*x**2)) - 2*b*c*sqrt(a + b*x**4)/(21*a*x**3) - b*d*sqrt(a + b*x**4)/(6*a*x**2) - 2*b*e*sqrt(a + b*x**4)/(5*a*x) - b*f*atanh(sqrt(a + b*x**4)/sqrt(a))/(4*sqrt(a)) - 2*b**(sympy.S(5)/4)*e*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*a**(sympy.S(3)/4)*sqrt(a + b*x**4)) - b**(sympy.S(5)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(-21*sqrt(a)*e + 5*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(105*a**(sympy.S(5)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_505():
    f = sqrt(a + b*x**4)*(c + d*x + e*x**2 + f*x**3)/x**9
    F = -sqrt(a + b*x**4)*(105*c/x**8 + 120*d/x**7 + 140*e/x**6 + 168*f/x**5)/840 + 2*b**(sympy.S(3)/2)*f*x*sqrt(a + b*x**4)/(5*a*(sqrt(a) + sqrt(b)*x**2)) - b*c*sqrt(a + b*x**4)/(16*a*x**4) - 2*b*d*sqrt(a + b*x**4)/(21*a*x**3) - b*e*sqrt(a + b*x**4)/(6*a*x**2) - 2*b*f*sqrt(a + b*x**4)/(5*a*x) + b**2*c*atanh(sqrt(a + b*x**4)/sqrt(a))/(16*a**(sympy.S(3)/2)) - 2*b**(sympy.S(5)/4)*f*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*a**(sympy.S(3)/4)*sqrt(a + b*x**4)) - b**(sympy.S(5)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(-21*sqrt(a)*f + 5*sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(105*a**(sympy.S(5)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_506():
    f = sqrt(a + b*x**4)*(c + d*x + e*x**2 + f*x**3)/x**10
    F = -sqrt(a + b*x**4)*(56*c/x**9 + 63*d/x**8 + 72*e/x**7 + 84*f/x**6)/504 - 2*b*c*sqrt(a + b*x**4)/(45*a*x**5) - b*d*sqrt(a + b*x**4)/(16*a*x**4) - 2*b*e*sqrt(a + b*x**4)/(21*a*x**3) - b*f*sqrt(a + b*x**4)/(6*a*x**2) - 2*b**(sympy.S(5)/2)*c*x*sqrt(a + b*x**4)/(15*a**2*(sqrt(a) + sqrt(b)*x**2)) + 2*b**2*c*sqrt(a + b*x**4)/(15*a**2*x) + b**2*d*atanh(sqrt(a + b*x**4)/sqrt(a))/(16*a**(sympy.S(3)/2)) + 2*b**(sympy.S(9)/4)*c*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(15*a**(sympy.S(7)/4)*sqrt(a + b*x**4)) - b**(sympy.S(7)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(5*sqrt(a)*e + 7*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(105*a**(sympy.S(7)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_507():
    f = x**4*(a + b*x**4)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3)
    F = 4*a**(sympy.S(13)/4)*e*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(65*b**(sympy.S(7)/4)*sqrt(a + b*x**4)) - 2*a**(sympy.S(11)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(77*sqrt(a)*e + 65*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5005*b**(sympy.S(7)/4)*sqrt(a + b*x**4)) - a**3*d*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(32*b**(sympy.S(3)/2)) - 4*a**3*e*x*sqrt(a + b*x**4)/(65*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x**2)) + 4*a**2*c*x*sqrt(a + b*x**4)/(77*b) - a**2*d*x**2*sqrt(a + b*x**4)/(32*b) + 4*a**2*e*x**3*sqrt(a + b*x**4)/(195*b) + 2*a*x**5*sqrt(a + b*x**4)*(117*c + 77*e*x**2)/3003 - a*d*x**2*(a + b*x**4)**(sympy.S(3)/2)/(48*b) + x**5*(a + b*x**4)**(sympy.S(3)/2)*(13*c + 11*e*x**2)/143 + f*x**4*(a + b*x**4)**(sympy.S(5)/2)/(14*b) - (a + b*x**4)**(sympy.S(5)/2)*(12*a*f - 35*b*d*x**2)/(420*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_508():
    f = x**3*(a + b*x**4)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3)
    F = 4*a**(sympy.S(13)/4)*f*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(65*b**(sympy.S(7)/4)*sqrt(a + b*x**4)) - 2*a**(sympy.S(11)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(77*sqrt(a)*f + 65*sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5005*b**(sympy.S(7)/4)*sqrt(a + b*x**4)) - a**3*e*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(32*b**(sympy.S(3)/2)) - 4*a**3*f*x*sqrt(a + b*x**4)/(65*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x**2)) + 4*a**2*d*x*sqrt(a + b*x**4)/(77*b) - a**2*e*x**2*sqrt(a + b*x**4)/(32*b) + 4*a**2*f*x**3*sqrt(a + b*x**4)/(195*b) + 2*a*x**5*sqrt(a + b*x**4)*(117*d + 77*f*x**2)/3003 - a*e*x**2*(a + b*x**4)**(sympy.S(3)/2)/(48*b) + x**5*(a + b*x**4)**(sympy.S(3)/2)*(13*d + 11*f*x**2)/143 + (a + b*x**4)**(sympy.S(5)/2)*(6*c + 5*e*x**2)/(60*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_509():
    f = x**2*(a + b*x**4)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3)
    F = -4*a**(sympy.S(9)/4)*c*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(3)/4)*sqrt(a + b*x**4)) + 2*a**(sympy.S(9)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(-15*sqrt(a)*e + 77*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(1155*b**(sympy.S(5)/4)*sqrt(a + b*x**4)) - a**3*f*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(32*b**(sympy.S(3)/2)) + 4*a**2*e*x*sqrt(a + b*x**4)/(77*b) - a**2*f*x**2*sqrt(a + b*x**4)/(32*b) + 4*a**2*c*x*sqrt(a + b*x**4)/(15*sqrt(b)*(sqrt(a) + sqrt(b)*x**2)) + 2*a*x**3*sqrt(a + b*x**4)*(77*c + 45*e*x**2)/1155 - a*f*x**2*(a + b*x**4)**(sympy.S(3)/2)/(48*b) + x**3*(a + b*x**4)**(sympy.S(3)/2)*(11*c + 9*e*x**2)/99 + (a + b*x**4)**(sympy.S(5)/2)*(6*d + 5*f*x**2)/(60*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_510():
    f = x*(a + b*x**4)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3)
    F = -4*a**(sympy.S(9)/4)*d*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(3)/4)*sqrt(a + b*x**4)) + 2*a**(sympy.S(9)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(-15*sqrt(a)*f + 77*sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(1155*b**(sympy.S(5)/4)*sqrt(a + b*x**4)) + 4*a**2*f*x*sqrt(a + b*x**4)/(77*b) + 3*a**2*c*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(16*sqrt(b)) + 4*a**2*d*x*sqrt(a + b*x**4)/(15*sqrt(b)*(sqrt(a) + sqrt(b)*x**2)) + 3*a*c*x**2*sqrt(a + b*x**4)/16 + 2*a*x**3*sqrt(a + b*x**4)*(77*d + 45*f*x**2)/1155 + c*x**2*(a + b*x**4)**(sympy.S(3)/2)/8 + x**3*(a + b*x**4)**(sympy.S(3)/2)*(11*d + 9*f*x**2)/99 + e*(a + b*x**4)**(sympy.S(5)/2)/(10*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_511():
    f = (a + b*x**4)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3)
    F = -4*a**(sympy.S(9)/4)*e*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(3)/4)*sqrt(a + b*x**4)) + 2*a**(sympy.S(7)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(7*sqrt(a)*e + 15*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(105*b**(sympy.S(3)/4)*sqrt(a + b*x**4)) + 3*a**2*d*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(16*sqrt(b)) + 4*a**2*e*x*sqrt(a + b*x**4)/(15*sqrt(b)*(sqrt(a) + sqrt(b)*x**2)) + 3*a*d*x**2*sqrt(a + b*x**4)/16 + 2*a*x*sqrt(a + b*x**4)*(15*c + 7*e*x**2)/105 + d*x**2*(a + b*x**4)**(sympy.S(3)/2)/8 + x*(a + b*x**4)**(sympy.S(3)/2)*(9*c + 7*e*x**2)/63 + f*(a + b*x**4)**(sympy.S(5)/2)/(10*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_512():
    f = (a + b*x**4)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3)/x
    F = -4*a**(sympy.S(9)/4)*f*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(15*b**(sympy.S(3)/4)*sqrt(a + b*x**4)) + 2*a**(sympy.S(7)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(7*sqrt(a)*f + 15*sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(105*b**(sympy.S(3)/4)*sqrt(a + b*x**4)) - a**(sympy.S(3)/2)*c*atanh(sqrt(a + b*x**4)/sqrt(a))/2 + 3*a**2*e*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(16*sqrt(b)) + 4*a**2*f*x*sqrt(a + b*x**4)/(15*sqrt(b)*(sqrt(a) + sqrt(b)*x**2)) + 2*a*x*sqrt(a + b*x**4)*(15*d + 7*f*x**2)/105 + a*sqrt(a + b*x**4)*(8*c + 3*e*x**2)/16 + x*(a + b*x**4)**(sympy.S(3)/2)*(9*d + 7*f*x**2)/63 + (a + b*x**4)**(sympy.S(3)/2)*(4*c + 3*e*x**2)/24
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_513():
    f = (a + b*x**4)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3)/x**2
    F = -12*a**(sympy.S(5)/4)*b**(sympy.S(1)/4)*c*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*sqrt(a + b*x**4)) + 2*a**(sympy.S(5)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(5*sqrt(a)*e + 21*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(35*b**(sympy.S(1)/4)*sqrt(a + b*x**4)) - a**(sympy.S(3)/2)*d*atanh(sqrt(a + b*x**4)/sqrt(a))/2 + 3*a**2*f*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(16*sqrt(b)) + 12*a*sqrt(b)*c*x*sqrt(a + b*x**4)/(5*sqrt(a) + 5*sqrt(b)*x**2) + a*sqrt(a + b*x**4)*(8*d + 3*f*x**2)/16 + 2*x*sqrt(a + b*x**4)*(5*a*e + 21*b*c*x**2)/35 + (a + b*x**4)**(sympy.S(3)/2)*(4*d + 3*f*x**2)/24 - (a + b*x**4)**(sympy.S(3)/2)*(7*c - e*x**2)/(7*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_514():
    f = (a + b*x**4)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3)/x**3
    F = -12*a**(sympy.S(5)/4)*b**(sympy.S(1)/4)*d*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*sqrt(a + b*x**4)) + 2*a**(sympy.S(5)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(5*sqrt(a)*f + 21*sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(35*b**(sympy.S(1)/4)*sqrt(a + b*x**4)) - a**(sympy.S(3)/2)*e*atanh(sqrt(a + b*x**4)/sqrt(a))/2 + 3*a*sqrt(b)*c*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/4 + 12*a*sqrt(b)*d*x*sqrt(a + b*x**4)/(5*sqrt(a) + 5*sqrt(b)*x**2) + 2*x*sqrt(a + b*x**4)*(5*a*f + 21*b*d*x**2)/35 + sqrt(a + b*x**4)*(2*a*e + 3*b*c*x**2)/4 - (a + b*x**4)**(sympy.S(3)/2)*(7*d - f*x**2)/(7*x) - (a + b*x**4)**(sympy.S(3)/2)*(3*c - e*x**2)/(6*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_515():
    f = (a + b*x**4)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3)/x**4
    F = -12*a**(sympy.S(5)/4)*b**(sympy.S(1)/4)*e*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*sqrt(a + b*x**4)) + 2*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(9*sqrt(a)*e + 5*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(15*sqrt(a + b*x**4)) - a**(sympy.S(3)/2)*f*atanh(sqrt(a + b*x**4)/sqrt(a))/2 + 3*a*sqrt(b)*d*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/4 + 12*a*sqrt(b)*e*x*sqrt(a + b*x**4)/(5*sqrt(a) + 5*sqrt(b)*x**2) + sqrt(a + b*x**4)*(2*a*f + 3*b*d*x**2)/4 - sqrt(a + b*x**4)*(18*a*e - 10*b*c*x**2)/(15*x) - (a + b*x**4)**(sympy.S(3)/2)*(3*d - f*x**2)/(6*x**2) - (a + b*x**4)**(sympy.S(3)/2)*(5*c - 3*e*x**2)/(15*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_516():
    f = (a + b*x**4)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3)/x**5
    F = -12*a**(sympy.S(5)/4)*b**(sympy.S(1)/4)*f*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*sqrt(a + b*x**4)) + 2*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(9*sqrt(a)*f + 5*sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(15*sqrt(a + b*x**4)) - 3*sqrt(a)*b*c*atanh(sqrt(a + b*x**4)/sqrt(a))/4 + 3*a*sqrt(b)*e*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/4 + 12*a*sqrt(b)*f*x*sqrt(a + b*x**4)/(5*sqrt(a) + 5*sqrt(b)*x**2) + 2*b*x*sqrt(a + b*x**4)*(5*d + 9*f*x**2)/15 + 3*b*sqrt(a + b*x**4)*(c + e*x**2)/4 - (a + b*x**4)**(sympy.S(3)/2)*(3*c/x**4 + 4*d/x**3 + 6*e/x**2 + 12*f/x)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_517():
    f = (a + b*x**4)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3)/x**6
    F = -12*a**(sympy.S(1)/4)*b**(sympy.S(5)/4)*c*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*sqrt(a + b*x**4)) + 2*a**(sympy.S(1)/4)*b**(sympy.S(3)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(5*sqrt(a)*e + 9*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(15*sqrt(a + b*x**4)) - 3*sqrt(a)*b*d*atanh(sqrt(a + b*x**4)/sqrt(a))/4 + 3*a*sqrt(b)*f*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/4 + 12*b**(sympy.S(3)/2)*c*x*sqrt(a + b*x**4)/(5*sqrt(a) + 5*sqrt(b)*x**2) + 3*b*sqrt(a + b*x**4)*(d + f*x**2)/4 - 2*b*sqrt(a + b*x**4)*(9*c - 5*e*x**2)/(15*x) - (a + b*x**4)**(sympy.S(3)/2)*(12*c/x**5 + 15*d/x**4 + 20*e/x**3 + 30*f/x**2)/60
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_518():
    f = (a + b*x**4)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3)/x**7
    F = -12*a**(sympy.S(1)/4)*b**(sympy.S(5)/4)*d*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*sqrt(a + b*x**4)) + 2*a**(sympy.S(1)/4)*b**(sympy.S(3)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(5*sqrt(a)*f + 9*sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(15*sqrt(a + b*x**4)) - 3*sqrt(a)*b*e*atanh(sqrt(a + b*x**4)/sqrt(a))/4 + b**(sympy.S(3)/2)*c*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/2 + 12*b**(sympy.S(3)/2)*d*x*sqrt(a + b*x**4)/(5*sqrt(a) + 5*sqrt(b)*x**2) - 2*b*sqrt(a + b*x**4)*(9*d - 5*f*x**2)/(15*x) - b*sqrt(a + b*x**4)*(2*c - 3*e*x**2)/(4*x**2) - (a + b*x**4)**(sympy.S(3)/2)*(10*c/x**6 + 12*d/x**5 + 15*e/x**4 + 20*f/x**3)/60
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_519():
    f = (a + b*x**4)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3)/x**8
    F = -12*a**(sympy.S(1)/4)*b**(sympy.S(5)/4)*e*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*sqrt(a + b*x**4)) - 3*sqrt(a)*b*f*atanh(sqrt(a + b*x**4)/sqrt(a))/4 + b**(sympy.S(3)/2)*d*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/2 + 12*b**(sympy.S(3)/2)*e*x*sqrt(a + b*x**4)/(5*sqrt(a) + 5*sqrt(b)*x**2) - 12*b*e*sqrt(a + b*x**4)/(5*x) - b*sqrt(a + b*x**4)*(2*d - 3*f*x**2)/(4*x**2) - 2*b*sqrt(a + b*x**4)*(5*c - 21*e*x**2)/(35*x**3) - (a + b*x**4)**(sympy.S(3)/2)*(60*c/x**7 + 70*d/x**6 + 84*e/x**5 + 105*f/x**4)/420 + 2*b**(sympy.S(5)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(21*sqrt(a)*e + 5*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(35*a**(sympy.S(1)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_520():
    f = (a + b*x**4)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3)/x**9
    F = -12*a**(sympy.S(1)/4)*b**(sympy.S(5)/4)*f*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*sqrt(a + b*x**4)) + b**(sympy.S(3)/2)*e*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/2 + 12*b**(sympy.S(3)/2)*f*x*sqrt(a + b*x**4)/(5*sqrt(a) + 5*sqrt(b)*x**2) - b*sqrt(a + b*x**4)*(105*c/x**4 + 160*d/x**3 + 280*e/x**2 + 672*f/x)/560 - (a + b*x**4)**(sympy.S(3)/2)*(105*c/x**8 + 120*d/x**7 + 140*e/x**6 + 168*f/x**5)/840 - 3*b**2*c*atanh(sqrt(a + b*x**4)/sqrt(a))/(16*sqrt(a)) + 2*b**(sympy.S(5)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(21*sqrt(a)*f + 5*sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(35*a**(sympy.S(1)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_521():
    f = (a + b*x**4)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3)/x**10
    F = b**(sympy.S(3)/2)*f*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/2 - b*sqrt(a + b*x**4)*(224*c/x**5 + 315*d/x**4 + 480*e/x**3 + 840*f/x**2)/1680 - (a + b*x**4)**(sympy.S(3)/2)*(56*c/x**9 + 63*d/x**8 + 72*e/x**7 + 84*f/x**6)/504 + 4*b**(sympy.S(5)/2)*c*x*sqrt(a + b*x**4)/(15*a*(sqrt(a) + sqrt(b)*x**2)) - 4*b**2*c*sqrt(a + b*x**4)/(15*a*x) - 3*b**2*d*atanh(sqrt(a + b*x**4)/sqrt(a))/(16*sqrt(a)) - 4*b**(sympy.S(9)/4)*c*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(15*a**(sympy.S(3)/4)*sqrt(a + b*x**4)) + 2*b**(sympy.S(7)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(15*sqrt(a)*e + 7*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(105*a**(sympy.S(3)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_522():
    f = (a + b*x**4)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3)/x**11
    F = -b*sqrt(a + b*x**4)*(168*c/x**6 + 224*d/x**5 + 315*e/x**4 + 480*f/x**3)/1680 - (a + b*x**4)**(sympy.S(3)/2)*(252*c/x**10 + 280*d/x**9 + 315*e/x**8 + 360*f/x**7)/2520 + 4*b**(sympy.S(5)/2)*d*x*sqrt(a + b*x**4)/(15*a*(sqrt(a) + sqrt(b)*x**2)) - b**2*c*sqrt(a + b*x**4)/(10*a*x**2) - 4*b**2*d*sqrt(a + b*x**4)/(15*a*x) - 3*b**2*e*atanh(sqrt(a + b*x**4)/sqrt(a))/(16*sqrt(a)) - 4*b**(sympy.S(9)/4)*d*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(15*a**(sympy.S(3)/4)*sqrt(a + b*x**4)) + 2*b**(sympy.S(7)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(15*sqrt(a)*f + 7*sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(105*a**(sympy.S(3)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_523():
    f = (a + b*x**4)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3)/x**12
    F = -b*sqrt(a + b*x**4)*(1440*c/x**7 + 1848*d/x**6 + 2464*e/x**5 + 3465*f/x**4)/18480 - (a + b*x**4)**(sympy.S(3)/2)*(360*c/x**11 + 396*d/x**10 + 440*e/x**9 + 495*f/x**8)/3960 + 4*b**(sympy.S(5)/2)*e*x*sqrt(a + b*x**4)/(15*a*(sqrt(a) + sqrt(b)*x**2)) - 4*b**2*c*sqrt(a + b*x**4)/(77*a*x**3) - b**2*d*sqrt(a + b*x**4)/(10*a*x**2) - 4*b**2*e*sqrt(a + b*x**4)/(15*a*x) - 3*b**2*f*atanh(sqrt(a + b*x**4)/sqrt(a))/(16*sqrt(a)) - 4*b**(sympy.S(9)/4)*e*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(15*a**(sympy.S(3)/4)*sqrt(a + b*x**4)) - 2*b**(sympy.S(9)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(-77*sqrt(a)*e + 15*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(1155*a**(sympy.S(5)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_524():
    f = (a + b*x**4)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3)/x**13
    F = -b*sqrt(a + b*x**4)*(1155*c/x**8 + 1440*d/x**7 + 1848*e/x**6 + 2464*f/x**5)/18480 - (a + b*x**4)**(sympy.S(3)/2)*(165*c/x**12 + 180*d/x**11 + 198*e/x**10 + 220*f/x**9)/1980 + 4*b**(sympy.S(5)/2)*f*x*sqrt(a + b*x**4)/(15*a*(sqrt(a) + sqrt(b)*x**2)) - b**2*c*sqrt(a + b*x**4)/(32*a*x**4) - 4*b**2*d*sqrt(a + b*x**4)/(77*a*x**3) - b**2*e*sqrt(a + b*x**4)/(10*a*x**2) - 4*b**2*f*sqrt(a + b*x**4)/(15*a*x) + b**3*c*atanh(sqrt(a + b*x**4)/sqrt(a))/(32*a**(sympy.S(3)/2)) - 4*b**(sympy.S(9)/4)*f*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(15*a**(sympy.S(3)/4)*sqrt(a + b*x**4)) - 2*b**(sympy.S(9)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(-77*sqrt(a)*f + 15*sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(1155*a**(sympy.S(5)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_525():
    f = (a + b*x**4)**(sympy.S(3)/2)*(c + d*x + e*x**2 + f*x**3)/x**14
    F = -b*sqrt(a + b*x**4)*(12320*c/x**9 + 15015*d/x**8 + 18720*e/x**7 + 24024*f/x**6)/240240 - (a + b*x**4)**(sympy.S(3)/2)*(660*c/x**13 + 715*d/x**12 + 780*e/x**11 + 858*f/x**10)/8580 - 4*b**2*c*sqrt(a + b*x**4)/(195*a*x**5) - b**2*d*sqrt(a + b*x**4)/(32*a*x**4) - 4*b**2*e*sqrt(a + b*x**4)/(77*a*x**3) - b**2*f*sqrt(a + b*x**4)/(10*a*x**2) - 4*b**(sympy.S(7)/2)*c*x*sqrt(a + b*x**4)/(65*a**2*(sqrt(a) + sqrt(b)*x**2)) + 4*b**3*c*sqrt(a + b*x**4)/(65*a**2*x) + b**3*d*atanh(sqrt(a + b*x**4)/sqrt(a))/(32*a**(sympy.S(3)/2)) + 4*b**(sympy.S(13)/4)*c*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(65*a**(sympy.S(7)/4)*sqrt(a + b*x**4)) - 2*b**(sympy.S(11)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(65*sqrt(a)*e + 77*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5005*a**(sympy.S(7)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_526():
    f = x**4*(c + d*x + e*x**2 + f*x**3)/sqrt(a + b*x**4)
    F = 3*a**(sympy.S(5)/4)*e*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(7)/4)*sqrt(a + b*x**4)) - a**(sympy.S(3)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(9*sqrt(a)*e + 5*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(30*b**(sympy.S(7)/4)*sqrt(a + b*x**4)) - a*d*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(4*b**(sympy.S(3)/2)) - 3*a*e*x*sqrt(a + b*x**4)/(5*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x**2)) + c*x*sqrt(a + b*x**4)/(3*b) + e*x**3*sqrt(a + b*x**4)/(5*b) + f*x**4*sqrt(a + b*x**4)/(6*b) - sqrt(a + b*x**4)*(4*a*f - 3*b*d*x**2)/(12*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_527():
    f = x**3*(c + d*x + e*x**2 + f*x**3)/sqrt(a + b*x**4)
    F = 3*a**(sympy.S(5)/4)*f*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*b**(sympy.S(7)/4)*sqrt(a + b*x**4)) - a**(sympy.S(3)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(9*sqrt(a)*f + 5*sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(30*b**(sympy.S(7)/4)*sqrt(a + b*x**4)) - a*e*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(4*b**(sympy.S(3)/2)) - 3*a*f*x*sqrt(a + b*x**4)/(5*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x**2)) + d*x*sqrt(a + b*x**4)/(3*b) + f*x**3*sqrt(a + b*x**4)/(5*b) + sqrt(a + b*x**4)*(2*c + e*x**2)/(4*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_528():
    f = x**2*(c + d*x + e*x**2 + f*x**3)/sqrt(a + b*x**4)
    F = -a**(sympy.S(1)/4)*c*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(b**(sympy.S(3)/4)*sqrt(a + b*x**4)) + a**(sympy.S(1)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(-sqrt(a)*e + 3*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(6*b**(sympy.S(5)/4)*sqrt(a + b*x**4)) - a*f*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(4*b**(sympy.S(3)/2)) + e*x*sqrt(a + b*x**4)/(3*b) + sqrt(a + b*x**4)*(2*d + f*x**2)/(4*b) + c*x*sqrt(a + b*x**4)/(sqrt(b)*(sqrt(a) + sqrt(b)*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_529():
    f = x*(c + d*x + e*x**2 + f*x**3)/sqrt(a + b*x**4)
    F = -a**(sympy.S(1)/4)*d*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(b**(sympy.S(3)/4)*sqrt(a + b*x**4)) + a**(sympy.S(1)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(-sqrt(a)*f + 3*sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(6*b**(sympy.S(5)/4)*sqrt(a + b*x**4)) + e*sqrt(a + b*x**4)/(2*b) + f*x*sqrt(a + b*x**4)/(3*b) + c*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(2*sqrt(b)) + d*x*sqrt(a + b*x**4)/(sqrt(b)*(sqrt(a) + sqrt(b)*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_530():
    f = (c + d*x + e*x**2 + f*x**3)/sqrt(a + b*x**4)
    F = -a**(sympy.S(1)/4)*e*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(b**(sympy.S(3)/4)*sqrt(a + b*x**4)) + a**(sympy.S(1)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(e + sqrt(b)*c/sqrt(a))*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*b**(sympy.S(3)/4)*sqrt(a + b*x**4)) + f*sqrt(a + b*x**4)/(2*b) + d*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(2*sqrt(b)) + e*x*sqrt(a + b*x**4)/(sqrt(b)*(sqrt(a) + sqrt(b)*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_531():
    f = (c + d*x + e*x**2 + f*x**3)/(x*sqrt(a + b*x**4))
    F = -a**(sympy.S(1)/4)*f*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(b**(sympy.S(3)/4)*sqrt(a + b*x**4)) + a**(sympy.S(1)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(f + sqrt(b)*d/sqrt(a))*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*b**(sympy.S(3)/4)*sqrt(a + b*x**4)) + e*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(2*sqrt(b)) + f*x*sqrt(a + b*x**4)/(sqrt(b)*(sqrt(a) + sqrt(b)*x**2)) - c*atanh(sqrt(a + b*x**4)/sqrt(a))/(2*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_532():
    f = (c + d*x + e*x**2 + f*x**3)/(x**2*sqrt(a + b*x**4))
    F = f*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(2*sqrt(b)) + sqrt(b)*c*x*sqrt(a + b*x**4)/(a*(sqrt(a) + sqrt(b)*x**2)) - c*sqrt(a + b*x**4)/(a*x) - d*atanh(sqrt(a + b*x**4)/sqrt(a))/(2*sqrt(a)) - b**(sympy.S(1)/4)*c*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(a**(sympy.S(3)/4)*sqrt(a + b*x**4)) + sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(sqrt(a)*e + sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_533():
    f = (c + d*x + e*x**2 + f*x**3)/(x**3*sqrt(a + b*x**4))
    F = sqrt(b)*d*x*sqrt(a + b*x**4)/(a*(sqrt(a) + sqrt(b)*x**2)) - c*sqrt(a + b*x**4)/(2*a*x**2) - d*sqrt(a + b*x**4)/(a*x) - e*atanh(sqrt(a + b*x**4)/sqrt(a))/(2*sqrt(a)) - b**(sympy.S(1)/4)*d*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(a**(sympy.S(3)/4)*sqrt(a + b*x**4)) + sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(sqrt(a)*f + sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_534():
    f = (c + d*x + e*x**2 + f*x**3)/(x**4*sqrt(a + b*x**4))
    F = sqrt(b)*e*x*sqrt(a + b*x**4)/(a*(sqrt(a) + sqrt(b)*x**2)) - c*sqrt(a + b*x**4)/(3*a*x**3) - d*sqrt(a + b*x**4)/(2*a*x**2) - e*sqrt(a + b*x**4)/(a*x) - f*atanh(sqrt(a + b*x**4)/sqrt(a))/(2*sqrt(a)) - b**(sympy.S(1)/4)*e*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(a**(sympy.S(3)/4)*sqrt(a + b*x**4)) - b**(sympy.S(1)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(-3*sqrt(a)*e + sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(6*a**(sympy.S(5)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_535():
    f = (c + d*x + e*x**2 + f*x**3)/(x**5*sqrt(a + b*x**4))
    F = sqrt(b)*f*x*sqrt(a + b*x**4)/(a*(sqrt(a) + sqrt(b)*x**2)) - c*sqrt(a + b*x**4)/(4*a*x**4) - d*sqrt(a + b*x**4)/(3*a*x**3) - e*sqrt(a + b*x**4)/(2*a*x**2) - f*sqrt(a + b*x**4)/(a*x) + b*c*atanh(sqrt(a + b*x**4)/sqrt(a))/(4*a**(sympy.S(3)/2)) - b**(sympy.S(1)/4)*f*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(a**(sympy.S(3)/4)*sqrt(a + b*x**4)) - b**(sympy.S(1)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(-3*sqrt(a)*f + sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(6*a**(sympy.S(5)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_536():
    f = (c + d*x + e*x**2 + f*x**3)/(x**6*sqrt(a + b*x**4))
    F = -c*sqrt(a + b*x**4)/(5*a*x**5) - d*sqrt(a + b*x**4)/(4*a*x**4) - e*sqrt(a + b*x**4)/(3*a*x**3) - f*sqrt(a + b*x**4)/(2*a*x**2) - 3*b**(sympy.S(3)/2)*c*x*sqrt(a + b*x**4)/(5*a**2*(sqrt(a) + sqrt(b)*x**2)) + 3*b*c*sqrt(a + b*x**4)/(5*a**2*x) + b*d*atanh(sqrt(a + b*x**4)/sqrt(a))/(4*a**(sympy.S(3)/2)) + 3*b**(sympy.S(5)/4)*c*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*a**(sympy.S(7)/4)*sqrt(a + b*x**4)) - b**(sympy.S(3)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(5*sqrt(a)*e + 9*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(30*a**(sympy.S(7)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_537():
    f = x**6*(c + d*x + e*x**2 + f*x**3)/(a + b*x**4)**(sympy.S(3)/2)
    F = -3*a**(sympy.S(1)/4)*c*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*b**(sympy.S(7)/4)*sqrt(a + b*x**4)) + a**(sympy.S(1)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(-5*sqrt(a)*e + 9*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(12*b**(sympy.S(9)/4)*sqrt(a + b*x**4)) - 3*a*f*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(4*b**(sympy.S(5)/2)) + d*sqrt(a + b*x**4)/b**2 + e*x*sqrt(a + b*x**4)/(3*b**2) + f*x**2*sqrt(a + b*x**4)/(4*b**2) + x*(a*e + a*f*x - b*c*x**2 - b*d*x**3)/(2*b**2*sqrt(a + b*x**4)) + 3*c*x*sqrt(a + b*x**4)/(2*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_538():
    f = x**5*(c + d*x + e*x**2 + f*x**3)/(a + b*x**4)**(sympy.S(3)/2)
    F = -3*a**(sympy.S(1)/4)*d*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*b**(sympy.S(7)/4)*sqrt(a + b*x**4)) + a**(sympy.S(1)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(-5*sqrt(a)*f + 9*sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(12*b**(sympy.S(9)/4)*sqrt(a + b*x**4)) + e*sqrt(a + b*x**4)/b**2 + f*x*sqrt(a + b*x**4)/(3*b**2) + x*(a*f - b*c*x - b*d*x**2 - b*e*x**3)/(2*b**2*sqrt(a + b*x**4)) + c*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(2*b**(sympy.S(3)/2)) + 3*d*x*sqrt(a + b*x**4)/(2*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_539():
    f = x**4*(c + d*x + e*x**2 + f*x**3)/(a + b*x**4)**(sympy.S(3)/2)
    F = -3*a**(sympy.S(1)/4)*e*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*b**(sympy.S(7)/4)*sqrt(a + b*x**4)) - x*(c + d*x + e*x**2 + f*x**3)/(2*b*sqrt(a + b*x**4)) + f*sqrt(a + b*x**4)/b**2 + d*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(2*b**(sympy.S(3)/2)) + 3*e*x*sqrt(a + b*x**4)/(2*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x**2)) + sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(3*sqrt(a)*e + sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(4*a**(sympy.S(1)/4)*b**(sympy.S(7)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_540():
    f = x**3*(c + d*x + e*x**2 + f*x**3)/(a + b*x**4)**(sympy.S(3)/2)
    F = -3*a**(sympy.S(1)/4)*f*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*b**(sympy.S(7)/4)*sqrt(a + b*x**4)) + (-c - d*x - e*x**2 - f*x**3)/(2*b*sqrt(a + b*x**4)) + e*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(2*b**(sympy.S(3)/2)) + 3*f*x*sqrt(a + b*x**4)/(2*b**(sympy.S(3)/2)*(sqrt(a) + sqrt(b)*x**2)) + sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(3*sqrt(a)*f + sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(4*a**(sympy.S(1)/4)*b**(sympy.S(7)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_541():
    f = x**2*(c + d*x + e*x**2 + f*x**3)/(a + b*x**4)**(sympy.S(3)/2)
    F = f*atanh(sqrt(b)*x**2/sqrt(a + b*x**4))/(2*b**(sympy.S(3)/2)) - d*sqrt(a + b*x**4)/(2*a*b) - x*(a*e + a*f*x - b*c*x**2 - b*d*x**3)/(2*a*b*sqrt(a + b*x**4)) - c*x*sqrt(a + b*x**4)/(2*a*sqrt(b)*(sqrt(a) + sqrt(b)*x**2)) + c*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*sqrt(a + b*x**4)) - sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(-sqrt(a)*e + sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(4*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_542():
    f = x*(c + d*x + e*x**2 + f*x**3)/(a + b*x**4)**(sympy.S(3)/2)
    F = -e*sqrt(a + b*x**4)/(2*a*b) - x*(a*f - b*c*x - b*d*x**2 - b*e*x**3)/(2*a*b*sqrt(a + b*x**4)) - d*x*sqrt(a + b*x**4)/(2*a*sqrt(b)*(sqrt(a) + sqrt(b)*x**2)) + d*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*sqrt(a + b*x**4)) - sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(-sqrt(a)*f + sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(4*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_543():
    f = (c + d*x + e*x**2 + f*x**3)/(a + b*x**4)**(sympy.S(3)/2)
    F = -(a*f - b*x*(c + d*x + e*x**2))/(2*a*b*sqrt(a + b*x**4)) - e*x*sqrt(a + b*x**4)/(2*a*sqrt(b)*(sqrt(a) + sqrt(b)*x**2)) + e*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*sqrt(a + b*x**4)) + sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(-sqrt(a)*e + sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(4*a**(sympy.S(5)/4)*b**(sympy.S(3)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_544():
    f = (c + d*x + e*x**2 + f*x**3)/(x*(a + b*x**4)**(sympy.S(3)/2))
    F = -f*x*sqrt(a + b*x**4)/(2*a*sqrt(b)*(sqrt(a) + sqrt(b)*x**2)) + c*sqrt(a + b*x**4)/(2*a**2) + x*(a*d + a*e*x + a*f*x**2 - b*c*x**3)/(2*a**2*sqrt(a + b*x**4)) - c*atanh(sqrt(a + b*x**4)/sqrt(a))/(2*a**(sympy.S(3)/2)) + f*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*sqrt(a + b*x**4)) + sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(-sqrt(a)*f + sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(4*a**(sympy.S(5)/4)*b**(sympy.S(3)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_545():
    f = (c + d*x + e*x**2 + f*x**3)/(x**2*(a + b*x**4)**(sympy.S(3)/2))
    F = 3*sqrt(b)*c*x*sqrt(a + b*x**4)/(2*a**2*(sqrt(a) + sqrt(b)*x**2)) - c*sqrt(a + b*x**4)/(a**2*x) + d*sqrt(a + b*x**4)/(2*a**2) + x*(a*e + a*f*x - b*c*x**2 - b*d*x**3)/(2*a**2*sqrt(a + b*x**4)) - d*atanh(sqrt(a + b*x**4)/sqrt(a))/(2*a**(sympy.S(3)/2)) - 3*b**(sympy.S(1)/4)*c*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(7)/4)*sqrt(a + b*x**4)) + sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(sqrt(a)*e + 3*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(4*a**(sympy.S(7)/4)*b**(sympy.S(1)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_546():
    f = (c + d*x + e*x**2 + f*x**3)/(x**3*(a + b*x**4)**(sympy.S(3)/2))
    F = 3*sqrt(b)*d*x*sqrt(a + b*x**4)/(2*a**2*(sqrt(a) + sqrt(b)*x**2)) - c*sqrt(a + b*x**4)/(2*a**2*x**2) - d*sqrt(a + b*x**4)/(a**2*x) + e*sqrt(a + b*x**4)/(2*a**2) + x*(a*f - b*c*x - b*d*x**2 - b*e*x**3)/(2*a**2*sqrt(a + b*x**4)) - e*atanh(sqrt(a + b*x**4)/sqrt(a))/(2*a**(sympy.S(3)/2)) - 3*b**(sympy.S(1)/4)*d*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(7)/4)*sqrt(a + b*x**4)) + sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(sqrt(a)*f + 3*sqrt(b)*d)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(4*a**(sympy.S(7)/4)*b**(sympy.S(1)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_547():
    f = (c + d*x + e*x**2 + f*x**3)/(x**4*(a + b*x**4)**(sympy.S(3)/2))
    F = 3*sqrt(b)*e*x*sqrt(a + b*x**4)/(2*a**2*(sqrt(a) + sqrt(b)*x**2)) - c*sqrt(a + b*x**4)/(3*a**2*x**3) - d*sqrt(a + b*x**4)/(2*a**2*x**2) - e*sqrt(a + b*x**4)/(a**2*x) + f*sqrt(a + b*x**4)/(2*a**2) - x*(b*c + b*d*x + b*e*x**2 + b*f*x**3)/(2*a**2*sqrt(a + b*x**4)) - f*atanh(sqrt(a + b*x**4)/sqrt(a))/(2*a**(sympy.S(3)/2)) - 3*b**(sympy.S(1)/4)*e*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*elliptic_e(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(7)/4)*sqrt(a + b*x**4)) - b**(sympy.S(1)/4)*sqrt((a + b*x**4)/(sqrt(a) + sqrt(b)*x**2)**2)*(sqrt(a) + sqrt(b)*x**2)*(-9*sqrt(a)*e + 5*sqrt(b)*c)*elliptic_f(2*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(12*a**(sympy.S(9)/4)*sqrt(a + b*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_548():
    f = (g*x)**m*(a + b*x**4)**p*(c + d*x + e*x**2 + f*x**3)
    F = c*(g*x)**(m + 1)*(a + b*x**4)**p*hyper((-p, m/4 + sympy.S(1)/4), (m/4 + sympy.S(5)/4,), -b*x**4/a)/(g*(1 + b*x**4/a)**p*(m + 1)) + d*(g*x)**(m + 2)*(a + b*x**4)**p*hyper((-p, m/4 + sympy.S.Half), (m/4 + sympy.S(3)/2,), -b*x**4/a)/(g**2*(1 + b*x**4/a)**p*(m + 2)) + e*(g*x)**(m + 3)*(a + b*x**4)**p*hyper((-p, m/4 + sympy.S(3)/4), (m/4 + sympy.S(7)/4,), -b*x**4/a)/(g**3*(1 + b*x**4/a)**p*(m + 3)) + f*(g*x)**(m + 4)*(a + b*x**4)**p*hyper((-p, m/4 + 1), (m/4 + 2,), -b*x**4/a)/(g**4*(1 + b*x**4/a)**p*(m + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_549():
    f = x**3*(a + b*x**4)**p*(c + d*x + e*x**2 + f*x**3)
    F = d*x**5*(a + b*x**4)**p*hyper((sympy.S(5)/4, -p), (sympy.S(9)/4,), -b*x**4/a)/(5*(1 + b*x**4/a)**p) + e*x**6*(a + b*x**4)**p*hyper((sympy.S(3)/2, -p), (sympy.S(5)/2,), -b*x**4/a)/(6*(1 + b*x**4/a)**p) + f*x**7*(a + b*x**4)**p*hyper((sympy.S(7)/4, -p), (sympy.S(11)/4,), -b*x**4/a)/(7*(1 + b*x**4/a)**p) + c*(a + b*x**4)**(p + 1)/(4*b*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_550():
    f = (x**4 + x**3 + x**2 + x + 1)/(1 - x**5)
    F = -log(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_551():
    f = (-32*x**5 + 48*x**4 - 72*x**3 + 108*x**2 - 162*x + 243)/(729 - 64*x**6)
    F = log(2*x + 3)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_552():
    f = (32*x**5 + 48*x**4 + 72*x**3 + 108*x**2 + 162*x + 243)/(729 - 64*x**6)
    F = -log(3 - 2*x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_553():
    f = (16*x**4 + 36*x**2 + 81)/(729 - 64*x**6)
    F = atanh(2*x/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_554():
    f = (-16*x**4 - 24*x**3 + 54*x + 81)/(729 - 64*x**6)
    F = -sqrt(3)*atan(sqrt(3)*(3 - 4*x)/9)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_555():
    f = (3 - 2*x)/(729 - 64*x**6)
    F = log(2*x + 3)/486 - log(4*x**2 - 6*x + 9)/972 + sqrt(3)*atan(sqrt(3)*(4*x + 3)/9)/486
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_556():
    f = (2*x + 3)/(729 - 64*x**6)
    F = -log(3 - 2*x)/486 + log(4*x**2 + 6*x + 9)/972 - sqrt(3)*atan(sqrt(3)*(3 - 4*x)/9)/486
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_557():
    f = (4*x**2 - 6*x + 9)/(729 - 64*x**6)
    F = -log(3 - 2*x)/324 + log(2*x + 3)/108 - log(4*x**2 + 6*x + 9)/324 + sqrt(3)*atan(sqrt(3)*(4*x + 3)/9)/162
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_558():
    f = (4*x**2 + 6*x + 9)/(729 - 64*x**6)
    F = -log(3 - 2*x)/108 + log(2*x + 3)/324 + log(4*x**2 - 6*x + 9)/324 - sqrt(3)*atan(sqrt(3)*(3 - 4*x)/9)/162
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_559():
    f = (27 - 8*x**3)/(729 - 64*x**6)
    F = log(2*x + 3)/54 - log(4*x**2 - 6*x + 9)/108 - sqrt(3)*atan(sqrt(3)*(3 - 4*x)/9)/54
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_560():
    f = (8*x**3 + 24*x**2 + 36*x + 27)/(729 - 64*x**6)
    F = -log(3 - 2*x)/18 + log(4*x**2 - 6*x + 9)/36 - sqrt(3)*atan(sqrt(3)*(3 - 4*x)/9)/54
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_561():
    f = (-32*x**5 + 48*x**4 - 72*x**3 + 108*x**2 - 162*x + 243)/(729 - 64*x**6)**2
    F = -log(3 - 2*x)/17496 + 5*log(2*x + 3)/17496 - log(4*x**2 - 6*x + 9)/17496 - log(4*x**2 + 6*x + 9)/17496 - sqrt(3)*atan(sqrt(3)*(3 - 4*x)/9)/26244 + sqrt(3)*atan(sqrt(3)*(4*x + 3)/9)/8748 - 1/(5832*x + 8748)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_562():
    f = (32*x**5 + 48*x**4 + 72*x**3 + 108*x**2 + 162*x + 243)/(729 - 64*x**6)**2
    F = -5*log(3 - 2*x)/17496 + log(2*x + 3)/17496 + log(4*x**2 - 6*x + 9)/17496 + log(4*x**2 + 6*x + 9)/17496 - sqrt(3)*atan(sqrt(3)*(3 - 4*x)/9)/8748 + sqrt(3)*atan(sqrt(3)*(4*x + 3)/9)/26244 + 1/(8748 - 5832*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_563():
    f = (16*x**4 + 36*x**2 + 81)/(729 - 64*x**6)**2
    F = -sqrt(3)*atan(sqrt(3)*(3 - 4*x)/9)/39366 + sqrt(3)*atan(sqrt(3)*(4*x + 3)/9)/39366 + atanh(2*x/3)/8748 - 1/(34992*x + 52488) + 1/(52488 - 34992*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_564():
    f = (-16*x**4 - 24*x**3 + 54*x + 81)/(729 - 64*x**6)**2
    F = x/(17496*x**2 - 26244*x + 39366) - log(3 - 2*x)/26244 + log(2*x + 3)/78732 - log(4*x**2 - 6*x + 9)/157464 + log(4*x**2 + 6*x + 9)/52488 - sqrt(3)*atan(sqrt(3)*(3 - 4*x)/9)/13122
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_565():
    f = (3 - 2*x)/(729 - 64*x**6)**2
    F = x/(944784*x**2 + 1417176*x + 2125764) + (3 - x)/(2834352*x**2 - 4251528*x + 6377292) - log(3 - 2*x)/4251528 + log(2*x + 3)/472392 - log(4*x**2 - 6*x + 9)/944784 + log(4*x**2 + 6*x + 9)/8503056 - sqrt(3)*atan(sqrt(3)*(3 - 4*x)/9)/4251528 + sqrt(3)*atan(sqrt(3)*(4*x + 3)/9)/472392 - 1/(1417176*x + 2125764)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_566():
    f = (2*x + 3)/(729 - 64*x**6)**2
    F = x/(944784*x**2 - 1417176*x + 2125764) - (x + 3)/(2834352*x**2 + 4251528*x + 6377292) - log(3 - 2*x)/472392 + log(2*x + 3)/4251528 - log(4*x**2 - 6*x + 9)/8503056 + log(4*x**2 + 6*x + 9)/944784 - sqrt(3)*atan(sqrt(3)*(3 - 4*x)/9)/472392 + sqrt(3)*atan(sqrt(3)*(4*x + 3)/9)/4251528 + 1/(2125764 - 1417176*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_567():
    f = (4*x**2 - 6*x + 9)/(729 - 64*x**6)**2
    F = (4*x + 3)/(944784*x**2 + 1417176*x + 2125764) - log(3 - 2*x)/354294 + log(2*x + 3)/118098 - log(4*x**2 - 6*x + 9)/944784 - 5*log(4*x**2 + 6*x + 9)/2834352 - sqrt(3)*atan(sqrt(3)*(3 - 4*x)/9)/1417176 + sqrt(3)*atan(sqrt(3)*(4*x + 3)/9)/157464 - 1/(314928*x + 472392) + 1/(1417176 - 944784*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_568():
    f = (4*x**2 + 6*x + 9)/(729 - 64*x**6)**2
    F = -(3 - 4*x)/(944784*x**2 - 1417176*x + 2125764) - log(3 - 2*x)/118098 + log(2*x + 3)/354294 + 5*log(4*x**2 - 6*x + 9)/2834352 + log(4*x**2 + 6*x + 9)/944784 - sqrt(3)*atan(sqrt(3)*(3 - 4*x)/9)/157464 + sqrt(3)*atan(sqrt(3)*(4*x + 3)/9)/1417176 - 1/(944784*x + 1417176) + 1/(472392 - 314928*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_569():
    f = (27 - 8*x**3)/(729 - 64*x**6)**2
    F = x/(34992*x**3 + 118098) - log(3 - 2*x)/157464 + 7*log(2*x + 3)/472392 - 7*log(4*x**2 - 6*x + 9)/944784 + log(4*x**2 + 6*x + 9)/314928 - 7*sqrt(3)*atan(sqrt(3)*(3 - 4*x)/9)/472392 + sqrt(3)*atan(sqrt(3)*(4*x + 3)/9)/157464
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_570():
    f = (8*x**3 + 24*x**2 + 36*x + 27)/(729 - 64*x**6)**2
    F = -(3 - 2*x)/(104976*x**2 - 157464*x + 236196) - 7*log(3 - 2*x)/157464 + log(2*x + 3)/472392 + 17*log(4*x**2 - 6*x + 9)/944784 + log(4*x**2 + 6*x + 9)/314928 - 11*sqrt(3)*atan(sqrt(3)*(3 - 4*x)/9)/472392 - sqrt(3)*atan(sqrt(3)*(4*x + 3)/9)/472392 + 1/(78732 - 52488*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_571():
    f = x*(27 - 2*x**3)/(729 - 64*x**6)
    F = -log(3 - 2*x)/96 - 5*log(2*x + 3)/288 + 5*log(4*x**2 - 6*x + 9)/576 + log(4*x**2 + 6*x + 9)/192 - 5*sqrt(3)*atan(sqrt(3)*(3 - 4*x)/9)/288 - sqrt(3)*atan(sqrt(3)*(4*x + 3)/9)/96
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_572():
    f = (c*x)**m*(d + e*x**n + f*x**(2*n) + g*x**(3*n))/(a + b*x**n)
    F = g*x**(2*n + 1)*(c*x)**m/(b*(m + 2*n + 1)) + x**(n + 1)*(c*x)**m*(-a*g + b*f)/(b**2*(m + n + 1)) + (c*x)**(m + 1)*(a**2*g - a*b*f + b**2*e)/(b**3*c*(m + 1)) + (c*x)**(m + 1)*(-a**3*g + a**2*b*f - a*b**2*e + b**3*d)*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(a*b**3*c*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_573():
    f = (a + b*x**n)**3*(c + d*x**(n - 1))
    F = a**3*c*x + 3*a**2*b*c*x**(n + 1)/(n + 1) + 3*a*b**2*c*x**(2*n + 1)/(2*n + 1) + b**3*c*x**(3*n + 1)/(3*n + 1) + d*(a + b*x**n)**4/(4*b*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_574():
    f = (a + b*x**n)**2*(c + d*x**(n - 1))
    F = a**2*c*x + 2*a*b*c*x**(n + 1)/(n + 1) + b**2*c*x**(2*n + 1)/(2*n + 1) + d*(a + b*x**n)**3/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_575():
    f = (a + b*x**n)*(c + d*x**(n - 1))
    F = a*c*x + a*d*x**n/n + b*c*x**(n + 1)/(n + 1) + b*d*x**(2*n)/(2*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_576():
    f = c + d*x**(n - 1)
    F = c*x + d*x**n/n
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_577():
    f = (c + d*x**(n - 1))/(a + b*x**n)
    F = d*log(a + b*x**n)/(b*n) + c*x*hyper((1, 1/n), (1 + 1/n,), -b*x**n/a)/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_578():
    f = (c + d*x**(n - 1))/(a + b*x**n)**2
    F = -d/(b*n*(a + b*x**n)) + c*x*hyper((2, 1/n), (1 + 1/n,), -b*x**n/a)/a**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_579():
    f = (c + d*x**(n - 1))/(a + b*x**n)**3
    F = -d/(2*b*n*(a + b*x**n)**2) + c*x*hyper((3, 1/n), (1 + 1/n,), -b*x**n/a)/a**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_580():
    f = (c*x)**m*(d + e*x**n + f*x**(2*n) + g*x**(3*n))/sqrt(a + b*x**n)
    F = e*x**(n + 1)*(c*x)**m*sqrt(1 + b*x**n/a)*hyper((sympy.S.Half, (m + n + 1)/n), ((m + 2*n + 1)/n,), -b*x**n/a)/(sqrt(a + b*x**n)*(m + n + 1)) + f*x**(2*n + 1)*(c*x)**m*sqrt(1 + b*x**n/a)*hyper((sympy.S.Half, (m + 2*n + 1)/n), ((m + 3*n + 1)/n,), -b*x**n/a)/(sqrt(a + b*x**n)*(m + 2*n + 1)) + g*x**(3*n + 1)*(c*x)**m*sqrt(1 + b*x**n/a)*hyper((sympy.S.Half, (m + 3*n + 1)/n), ((m + 4*n + 1)/n,), -b*x**n/a)/(sqrt(a + b*x**n)*(m + 3*n + 1)) + d*(c*x)**(m + 1)*sqrt(1 + b*x**n/a)*hyper((sympy.S.Half, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(c*sqrt(a + b*x**n)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_581():
    f = (-a*h*x**(n/4 - 1) + b*f*x**(n/2 - 1) + b*g*x**(n - 1) + b*h*x**(5*n/4 - 1))/(a + b*x**n)**(sympy.S(3)/2)
    F = -(2*a*g + 4*a*h*x**(n/4) - 2*b*f*x**(n/2))/(a*n*sqrt(a + b*x**n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_582():
    f = (c*x)**m*(a + b*x**n)**p*(d + e*x + f*x**2 + g*x**3)
    F = d*(c*x)**(m + 1)*(a + b*x**n)**p*hyper((-p, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(c*(1 + b*x**n/a)**p*(m + 1)) + e*(c*x)**(m + 2)*(a + b*x**n)**p*hyper((-p, (m + 2)/n), ((m + n + 2)/n,), -b*x**n/a)/(c**2*(1 + b*x**n/a)**p*(m + 2)) + f*(c*x)**(m + 3)*(a + b*x**n)**p*hyper((-p, (m + 3)/n), ((m + n + 3)/n,), -b*x**n/a)/(c**3*(1 + b*x**n/a)**p*(m + 3)) + g*(c*x)**(m + 4)*(a + b*x**n)**p*hyper((-p, (m + 4)/n), ((m + n + 4)/n,), -b*x**n/a)/(c**4*(1 + b*x**n/a)**p*(m + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_583():
    f = (c*x)**m*(a + b*x**n)**p*(d + e*x**n + f*x**(2*n) + g*x**(3*n))
    F = e*x**(n + 1)*(c*x)**m*(a + b*x**n)**p*hyper((-p, (m + n + 1)/n), ((m + 2*n + 1)/n,), -b*x**n/a)/((1 + b*x**n/a)**p*(m + n + 1)) + f*x**(2*n + 1)*(c*x)**m*(a + b*x**n)**p*hyper((-p, (m + 2*n + 1)/n), ((m + 3*n + 1)/n,), -b*x**n/a)/((1 + b*x**n/a)**p*(m + 2*n + 1)) + g*x**(3*n + 1)*(c*x)**m*(a + b*x**n)**p*hyper((-p, (m + 3*n + 1)/n), ((m + 4*n + 1)/n,), -b*x**n/a)/((1 + b*x**n/a)**p*(m + 3*n + 1)) + d*(c*x)**(m + 1)*(a + b*x**n)**p*hyper((-p, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(c*(1 + b*x**n/a)**p*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_584():
    f = (c + d*x**(n/2) + e*x**n + f*x**(3*n/2))/(a + b*x**n)**2
    F = x*(-a*e + b*c + x**(n/2)*(-a*f + b*d))/(a*b*n*(a + b*x**n)) + x*(a*e - b*c*(1 - n))*hyper((1, 1/n), (1 + 1/n,), -b*x**n/a)/(a**2*b*n) - x**(n/2 + 1)*(-a*f*(n + 2) + b*d*(2 - n))*hyper((1, sympy.S.Half + 1/n), (sympy.S(3)/2 + 1/n,), -b*x**n/a)/(a**2*b*n*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_585():
    f = (a*c + 3*b*d*x**4 + x**2*(2*a*d + 2*b*c))/(sqrt(a + b*x**2)*sqrt(c + d*x**2))
    F = x*sqrt(a + b*x**2)*sqrt(c + d*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_586():
    f = (x**3 + 1)/((1 - x**4)*(x**4 + 1)**(sympy.S(1)/4))
    F = -2**(sympy.S(3)/4)*atan(2**(sympy.S(3)/4)*(x**4 + 1)**(sympy.S(1)/4)/2)/4 + 2**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*x/(x**4 + 1)**(sympy.S(1)/4))/4 + 2**(sympy.S(3)/4)*atanh(2**(sympy.S(3)/4)*(x**4 + 1)**(sympy.S(1)/4)/2)/4 + 2**(sympy.S(3)/4)*atanh(2**(sympy.S(1)/4)*x/(x**4 + 1)**(sympy.S(1)/4))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_587():
    f = (a + b*x**n)**((-n - 1)/n)*(c + d*x**n)**((-n - 1)/n)*(a*c - b*d*x**(2*n))
    F = x/((a + b*x**n)**(1/n)*(c + d*x**n)**(1/n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_588():
    f = (h*x)**(-n*p - n - 1)*(a + b*x**n)**p*(c + d*x**n)**p*(a*c - b*d*x**(2*n))
    F = -(a + b*x**n)**(p + 1)*(c + d*x**n)**(p + 1)/(h*n*(h*x)**(n*(p + 1))*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_589():
    f = (a + b*x**n)**p*(c + d*x**n)**p*(e + b*d*e*x**(2*n)*(2*n*p + 2*n + 1)/(a*c) + e*x**n*(a*d + b*c)*(n*p + n + 1)/(a*c))
    F = e*x*(a + b*x**n)**(p + 1)*(c + d*x**n)**(p + 1)/(a*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_8_P_x_c_x_pow_m_a_plus_b_x_pow_n_pow_p_590():
    f = (h*x)**m*(a + b*x**n)**p*(c + d*x**n)**p*(e + b*d*e*x**(2*n)*(m + 2*n*p + 2*n + 1)/(a*c*(m + 1)) + e*x**n*(a*d + b*c)*(m + n*p + n + 1)/(a*c*(m + 1)))
    F = e*(h*x)**(m + 1)*(a + b*x**n)**(p + 1)*(c + d*x**n)**(p + 1)/(a*c*h*(m + 1))
    assert integrate(f, x) == F

