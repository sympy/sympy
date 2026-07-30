"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.1 Binomial products/1.1.3 General/1.1.3.3 (a+b x^n)^p (c+d x^n)^q.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, m, n, p, q = symbols('a b c d e m n p q')

def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_1():
    f = (a + b*x**3)*(c + d*x**3)**4
    F = a*c**4*x + b*d**4*x**16/16 + c**3*x**4*(4*a*d + b*c)/4 + 2*c**2*d*x**7*(3*a*d + 2*b*c)/7 + c*d**2*x**10*(2*a*d + 3*b*c)/5 + d**3*x**13*(a*d + 4*b*c)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_2():
    f = (a + b*x**3)*(c + d*x**3)**3
    F = a*c**3*x + b*d**3*x**13/13 + c**2*x**4*(3*a*d + b*c)/4 + 3*c*d*x**7*(a*d + b*c)/7 + d**2*x**10*(a*d + 3*b*c)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_3():
    f = (a + b*x**3)*(c + d*x**3)**2
    F = a*c**2*x + b*d**2*x**10/10 + c*x**4*(2*a*d + b*c)/4 + d*x**7*(a*d + 2*b*c)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_4():
    f = (a + b*x**3)*(c + d*x**3)
    F = a*c*x + b*d*x**7/7 + x**4*(a*d/4 + b*c/4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_5():
    f = (a + b*x**3)/(c + d*x**3)
    F = b*x/d - (-a*d + b*c)*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(2)/3)*d**(sympy.S(4)/3)) + (-a*d + b*c)*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(6*c**(sympy.S(2)/3)*d**(sympy.S(4)/3)) + sqrt(3)*(-a*d + b*c)*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(3*c**(sympy.S(2)/3)*d**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_6():
    f = (a + b*x**3)/(c + d*x**3)**2
    F = -x*(-a*d + b*c)/(3*c*d*(c + d*x**3)) + (2*a*d + b*c)*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(9*c**(sympy.S(5)/3)*d**(sympy.S(4)/3)) - (2*a*d + b*c)*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(18*c**(sympy.S(5)/3)*d**(sympy.S(4)/3)) - sqrt(3)*(2*a*d + b*c)*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(9*c**(sympy.S(5)/3)*d**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_7():
    f = (a + b*x**3)/(c + d*x**3)**3
    F = -x*(-a*d + b*c)/(6*c*d*(c + d*x**3)**2) + x*(5*a*d + b*c)/(18*c**2*d*(c + d*x**3)) + (5*a*d + b*c)*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(27*c**(sympy.S(8)/3)*d**(sympy.S(4)/3)) - (5*a*d + b*c)*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(54*c**(sympy.S(8)/3)*d**(sympy.S(4)/3)) - sqrt(3)*(5*a*d + b*c)*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(27*c**(sympy.S(8)/3)*d**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_8():
    f = (a + b*x**3)**2*(c + d*x**3)**3
    F = a**2*c**3*x + a*c**2*x**4*(3*a*d + 2*b*c)/4 + b**2*d**3*x**16/16 + b*d**2*x**13*(2*a*d + 3*b*c)/13 + c*x**7*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2)/7 + d*x**10*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_9():
    f = (a + b*x**3)**2*(c + d*x**3)**2
    F = a**2*c**2*x + a*c*x**4*(a*d + b*c)/2 + b**2*d**2*x**13/13 + b*d*x**10*(a*d + b*c)/5 + x**7*(a**2*d**2/7 + 4*a*b*c*d/7 + b**2*c**2/7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_10():
    f = (a + b*x**3)**2*(c + d*x**3)
    F = a**2*c*x + a*x**4*(a*d + 2*b*c)/4 + b**2*d*x**10/10 + b*x**7*(2*a*d + b*c)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_11():
    f = (a + b*x**3)**2/(c + d*x**3)
    F = b**2*x**4/(4*d) - b*x*(-2*a*d + b*c)/d**2 + (-a*d + b*c)**2*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(2)/3)*d**(sympy.S(7)/3)) - (-a*d + b*c)**2*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(6*c**(sympy.S(2)/3)*d**(sympy.S(7)/3)) - sqrt(3)*(-a*d + b*c)**2*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(3*c**(sympy.S(2)/3)*d**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_12():
    f = (a + b*x**3)**2/(c + d*x**3)**2
    F = b**2*x/d**2 + x*(-a*d + b*c)**2/(3*c*d**2*(c + d*x**3)) - (-2*a*d + 2*b*c)*(a*d + 2*b*c)*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(9*c**(sympy.S(5)/3)*d**(sympy.S(7)/3)) + sqrt(3)*(-2*a*d + 2*b*c)*(a*d + 2*b*c)*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(9*c**(sympy.S(5)/3)*d**(sympy.S(7)/3)) + (-a*d + b*c)*(a*d + 2*b*c)*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(9*c**(sympy.S(5)/3)*d**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_13():
    f = (a + b*x**3)**2/(c + d*x**3)**3
    F = -x*(a + b*x**3)*(-a*d + b*c)/(6*c*d*(c + d*x**3)**2) - x*(-a*d + b*c)*(5*a*d + 4*b*c)/(18*c**2*d**2*(c + d*x**3)) + (5*a**2*d**2 + 2*a*b*c*d + 2*b**2*c**2)*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(27*c**(sympy.S(8)/3)*d**(sympy.S(7)/3)) - (5*a**2*d**2 + 2*a*b*c*d + 2*b**2*c**2)*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(54*c**(sympy.S(8)/3)*d**(sympy.S(7)/3)) - sqrt(3)*(5*a**2*d**2 + 2*a*b*c*d + 2*b**2*c**2)*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(27*c**(sympy.S(8)/3)*d**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_14():
    f = (c + d*x**3)**4/(a + b*x**3)
    F = d**4*x**10/(10*b) + d**3*x**7*(-a*d + 4*b*c)/(7*b**2) + d**2*x**4*(a**2*d**2 - 4*a*b*c*d + 6*b**2*c**2)/(4*b**3) + d*x*(-a*d + 2*b*c)*(a**2*d**2 - 2*a*b*c*d + 2*b**2*c**2)/b**4 + (-a*d + b*c)**4*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(2)/3)*b**(sympy.S(13)/3)) - (-a*d + b*c)**4*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(13)/3)) - sqrt(3)*(-a*d + b*c)**4*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(13)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_15():
    f = (c + d*x**3)**3/(a + b*x**3)
    F = d**3*x**7/(7*b) + d**2*x**4*(-a*d + 3*b*c)/(4*b**2) + d*x*(a**2*d**2 - 3*a*b*c*d + 3*b**2*c**2)/b**3 + (-a*d + b*c)**3*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(2)/3)*b**(sympy.S(10)/3)) - (-a*d + b*c)**3*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(10)/3)) - sqrt(3)*(-a*d + b*c)**3*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(10)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_16():
    f = (c + d*x**3)**2/(a + b*x**3)
    F = d**2*x**4/(4*b) + d*x*(-a*d + 2*b*c)/b**2 + (-a*d + b*c)**2*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(2)/3)*b**(sympy.S(7)/3)) - (-a*d + b*c)**2*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(7)/3)) - sqrt(3)*(-a*d + b*c)**2*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_17():
    f = (c + d*x**3)/(a + b*x**3)
    F = d*x/b + (-a*d + b*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(2)/3)*b**(sympy.S(4)/3)) - (-a*d + b*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(4)/3)) - sqrt(3)*(-a*d + b*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_18():
    f = 1/((a + b*x**3)*(c + d*x**3))
    F = -d**(sympy.S(2)/3)*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(2)/3)*(-a*d + b*c)) + d**(sympy.S(2)/3)*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(6*c**(sympy.S(2)/3)*(-a*d + b*c)) + sqrt(3)*d**(sympy.S(2)/3)*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(3*c**(sympy.S(2)/3)*(-a*d + b*c)) + b**(sympy.S(2)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(2)/3)*(-a*d + b*c)) - b**(sympy.S(2)/3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(2)/3)*(-a*d + b*c)) - sqrt(3)*b**(sympy.S(2)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_19():
    f = 1/((a + b*x**3)*(c + d*x**3)**2)
    F = -d*x/(3*c*(c + d*x**3)*(-a*d + b*c)) - d**(sympy.S(2)/3)*(-2*a*d + 5*b*c)*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(9*c**(sympy.S(5)/3)*(-a*d + b*c)**2) + d**(sympy.S(2)/3)*(-2*a*d + 5*b*c)*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(18*c**(sympy.S(5)/3)*(-a*d + b*c)**2) + sqrt(3)*d**(sympy.S(2)/3)*(-2*a*d + 5*b*c)*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(9*c**(sympy.S(5)/3)*(-a*d + b*c)**2) + b**(sympy.S(5)/3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(2)/3)*(-a*d + b*c)**2) - b**(sympy.S(5)/3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(2)/3)*(-a*d + b*c)**2) - sqrt(3)*b**(sympy.S(5)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_20():
    f = (c + d*x**3)**5/(a + b*x**3)**2
    F = d**5*x**10/(10*b**2) + d**4*x**7*(-2*a*d + 5*b*c)/(7*b**3) + d**3*x**4*(3*a**2*d**2 - 10*a*b*c*d + 10*b**2*c**2)/(4*b**4) + d**2*x*(-4*a**3*d**3 + 15*a**2*b*c*d**2 - 20*a*b**2*c**2*d + 10*b**3*c**3)/b**5 + x*(-a*d + b*c)**5/(3*a*b**5*(a + b*x**3)) + (-a*d + b*c)**4*(13*a*d + 2*b*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(5)/3)*b**(sympy.S(16)/3)) - (-a*d + b*c)**4*(13*a*d + 2*b*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(5)/3)*b**(sympy.S(16)/3)) - sqrt(3)*(-a*d + b*c)**4*(13*a*d + 2*b*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*b**(sympy.S(16)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_21():
    f = (c + d*x**3)**4/(a + b*x**3)**2
    F = d**4*x**7/(7*b**2) + d**3*x**4*(-a*d + 2*b*c)/(2*b**3) + d**2*x*(3*a**2*d**2 - 8*a*b*c*d + 6*b**2*c**2)/b**4 + x*(-a*d + b*c)**4/(3*a*b**4*(a + b*x**3)) + 2*(-a*d + b*c)**3*(5*a*d + b*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(5)/3)*b**(sympy.S(13)/3)) - (-a*d + b*c)**3*(5*a*d + b*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(9*a**(sympy.S(5)/3)*b**(sympy.S(13)/3)) - 2*sqrt(3)*(-a*d + b*c)**3*(5*a*d + b*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*b**(sympy.S(13)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_22():
    f = (c + d*x**3)**3/(a + b*x**3)**2
    F = d**3*x**4/(4*b**2) + d**2*x*(-2*a*d + 3*b*c)/b**3 + x*(-a*d + b*c)**3/(3*a*b**3*(a + b*x**3)) + (-a*d + b*c)**2*(7*a*d + 2*b*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(5)/3)*b**(sympy.S(10)/3)) - (-a*d + b*c)**2*(7*a*d + 2*b*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(5)/3)*b**(sympy.S(10)/3)) - sqrt(3)*(-a*d + b*c)**2*(7*a*d + 2*b*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*b**(sympy.S(10)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_23():
    f = (c + d*x**3)**2/(a + b*x**3)**2
    F = d**2*x/b**2 + x*(-a*d + b*c)**2/(3*a*b**2*(a + b*x**3)) + (-2*a*d + 2*b*c)*(2*a*d + b*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(5)/3)*b**(sympy.S(7)/3)) - sqrt(3)*(-2*a*d + 2*b*c)*(2*a*d + b*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*b**(sympy.S(7)/3)) - (-a*d + b*c)*(2*a*d + b*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(9*a**(sympy.S(5)/3)*b**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_24():
    f = (c + d*x**3)/(a + b*x**3)**2
    F = x*(-a*d + b*c)/(3*a*b*(a + b*x**3)) + (a*d + 2*b*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(5)/3)*b**(sympy.S(4)/3)) - (a*d + 2*b*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(5)/3)*b**(sympy.S(4)/3)) - sqrt(3)*(a*d + 2*b*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*b**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_25():
    f = 1/((a + b*x**3)**2*(c + d*x**3))
    F = d**(sympy.S(5)/3)*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(2)/3)*(-a*d + b*c)**2) - d**(sympy.S(5)/3)*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(6*c**(sympy.S(2)/3)*(-a*d + b*c)**2) - sqrt(3)*d**(sympy.S(5)/3)*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(3*c**(sympy.S(2)/3)*(-a*d + b*c)**2) + b*x/(3*a*(a + b*x**3)*(-a*d + b*c)) + b**(sympy.S(2)/3)*(-5*a*d + 2*b*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(5)/3)*(-a*d + b*c)**2) - b**(sympy.S(2)/3)*(-5*a*d + 2*b*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(18*a**(sympy.S(5)/3)*(-a*d + b*c)**2) - sqrt(3)*b**(sympy.S(2)/3)*(-5*a*d + 2*b*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_26():
    f = 1/((a + b*x**3)**2*(c + d*x**3)**2)
    F = 2*d**(sympy.S(5)/3)*(-a*d + 4*b*c)*log(c**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(9*c**(sympy.S(5)/3)*(-a*d + b*c)**3) - d**(sympy.S(5)/3)*(-a*d + 4*b*c)*log(c**(sympy.S(2)/3) - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(9*c**(sympy.S(5)/3)*(-a*d + b*c)**3) - 2*sqrt(3)*d**(sympy.S(5)/3)*(-a*d + 4*b*c)*atan(sqrt(3)*(c**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*c**(sympy.S(1)/3)))/(9*c**(sympy.S(5)/3)*(-a*d + b*c)**3) + b*x/(3*a*(a + b*x**3)*(c + d*x**3)*(-a*d + b*c)) + d*x*(a*d + b*c)/(3*a*c*(c + d*x**3)*(-a*d + b*c)**2) + 2*b**(sympy.S(5)/3)*(-4*a*d + b*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(9*a**(sympy.S(5)/3)*(-a*d + b*c)**3) - b**(sympy.S(5)/3)*(-4*a*d + b*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(9*a**(sympy.S(5)/3)*(-a*d + b*c)**3) - 2*sqrt(3)*b**(sympy.S(5)/3)*(-4*a*d + b*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_27():
    f = (a + b*x**3)**3/(c + d*x**3)**(sympy.S(13)/3)
    F = 81*a**3*x/(140*c**4*(c + d*x**3)**(sympy.S(1)/3)) + 27*a**2*x*(a + b*x**3)/(140*c**3*(c + d*x**3)**(sympy.S(4)/3)) + 9*a*x*(a + b*x**3)**2/(70*c**2*(c + d*x**3)**(sympy.S(7)/3)) + x*(a + b*x**3)**3/(10*c*(c + d*x**3)**(sympy.S(10)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_28():
    f = (a + b*x**3)**2/(c + d*x**3)**(sympy.S(10)/3)
    F = 9*a**2*x/(14*c**3*(c + d*x**3)**(sympy.S(1)/3)) + 3*a*x*(a + b*x**3)/(14*c**2*(c + d*x**3)**(sympy.S(4)/3)) + x*(a + b*x**3)**2/(7*c*(c + d*x**3)**(sympy.S(7)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_29():
    f = (a + b*x**3)/(c + d*x**3)**(sympy.S(7)/3)
    F = 3*a*x/(4*c**2*(c + d*x**3)**(sympy.S(1)/3)) + x*(a + b*x**3)/(4*c*(c + d*x**3)**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_30():
    f = (c + d*x**3)**(sympy.S(-4)/3)
    F = x/(c*(c + d*x**3)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_31():
    f = 1/((a + b*x**3)*(c + d*x**3)**(sympy.S(1)/3))
    F = log(a**(sympy.S(1)/3) + x*(-a*d + b*c)**(sympy.S(1)/3)/(c + d*x**3)**(sympy.S(1)/3))/(3*a**(sympy.S(2)/3)*(-a*d + b*c)**(sympy.S(1)/3)) - log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*x*(-a*d + b*c)**(sympy.S(1)/3)/(c + d*x**3)**(sympy.S(1)/3) + x**2*(-a*d + b*c)**(sympy.S(2)/3)/(c + d*x**3)**(sympy.S(2)/3))/(6*a**(sympy.S(2)/3)*(-a*d + b*c)**(sympy.S(1)/3)) - sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*x*(-a*d + b*c)**(sympy.S(1)/3)/(c + d*x**3)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*(-a*d + b*c)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_32():
    f = (c + d*x**3)**(sympy.S(2)/3)/(a + b*x**3)**2
    F = x*(c + d*x**3)**(sympy.S(2)/3)/(3*a*(a + b*x**3)) + 2*c*log(a**(sympy.S(1)/3) + x*(-a*d + b*c)**(sympy.S(1)/3)/(c + d*x**3)**(sympy.S(1)/3))/(9*a**(sympy.S(5)/3)*(-a*d + b*c)**(sympy.S(1)/3)) - c*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*x*(-a*d + b*c)**(sympy.S(1)/3)/(c + d*x**3)**(sympy.S(1)/3) + x**2*(-a*d + b*c)**(sympy.S(2)/3)/(c + d*x**3)**(sympy.S(2)/3))/(9*a**(sympy.S(5)/3)*(-a*d + b*c)**(sympy.S(1)/3)) - 2*sqrt(3)*c*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*x*(-a*d + b*c)**(sympy.S(1)/3)/(c + d*x**3)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(9*a**(sympy.S(5)/3)*(-a*d + b*c)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_33():
    f = (c + d*x**3)**(sympy.S(5)/3)/(a + b*x**3)**3
    F = x*(c + d*x**3)**(sympy.S(5)/3)/(6*a*(a + b*x**3)**2) + 5*c*x*(c + d*x**3)**(sympy.S(2)/3)/(18*a**2*(a + b*x**3)) + 5*c**2*log(a**(sympy.S(1)/3) + x*(-a*d + b*c)**(sympy.S(1)/3)/(c + d*x**3)**(sympy.S(1)/3))/(27*a**(sympy.S(8)/3)*(-a*d + b*c)**(sympy.S(1)/3)) - 5*c**2*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*x*(-a*d + b*c)**(sympy.S(1)/3)/(c + d*x**3)**(sympy.S(1)/3) + x**2*(-a*d + b*c)**(sympy.S(2)/3)/(c + d*x**3)**(sympy.S(2)/3))/(54*a**(sympy.S(8)/3)*(-a*d + b*c)**(sympy.S(1)/3)) - 5*sqrt(3)*c**2*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*x*(-a*d + b*c)**(sympy.S(1)/3)/(c + d*x**3)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(8)/3)*(-a*d + b*c)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_34():
    f = (a + b*x**3)**m*(c + d*x**3)**p
    F = x*(a + b*x**3)**m*(c + d*x**3)**p*appellf1(sympy.S(1)/3, -m, -p, sympy.S(4)/3, -b*x**3/a, -d*x**3/c)/((1 + b*x**3/a)**m*(1 + d*x**3/c)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_35():
    f = (c + d*x**3)**q/(a + b*x**3)
    F = x*(c + d*x**3)**q*appellf1(sympy.S(1)/3, 1, -q, sympy.S(4)/3, -b*x**3/a, -d*x**3/c)/(a*(1 + d*x**3/c)**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_36():
    f = (c + d*x**3)**q/(a + b*x**3)**2
    F = x*(c + d*x**3)**q*appellf1(sympy.S(1)/3, 2, -q, sympy.S(4)/3, -b*x**3/a, -d*x**3/c)/(a**2*(1 + d*x**3/c)**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_37():
    f = (a + b*x**3)**m*(c + d*x**3)**3
    F = d*x*(a + b*x**3)**(m + 1)*(c + d*x**3)**2/(b*(3*m + 10)) - d*x*(a + b*x**3)**(m + 1)*(c + d*x**3)*(7*a*d - b*c*(3*m + 16))/(b**2*(3*m + 7)*(3*m + 10)) + d*x*(a + b*x**3)**(m + 1)*(28*a**2*d**2 - a*b*c*d*(15*m + 92) + b**2*c**2*(9*m**2 + 60*m + 118))/(b**3*(3*m + 4)*(3*m + 7)*(3*m + 10)) - x*(a + b*x**3)**m*(28*a**3*d**3 - 12*a**2*b*c*d**2*(3*m + 10) + 3*a*b**2*c**2*d*(9*m**2 + 51*m + 70) - b**3*c**3*(27*m**3 + 189*m**2 + 414*m + 280))*hyper((sympy.S(1)/3, -m), (sympy.S(4)/3,), -b*x**3/a)/(b**3*(1 + b*x**3/a)**m*(3*m + 4)*(3*m + 7)*(3*m + 10))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_38():
    f = (a + b*x**3)**m*(c + d*x**3)**2
    F = d*x*(a + b*x**3)**(m + 1)*(c + d*x**3)/(b*(3*m + 7)) - d*x*(a + b*x**3)**(m + 1)*(4*a*d - b*c*(3*m + 10))/(b**2*(3*m + 4)*(3*m + 7)) + x*(a + b*x**3)**m*(4*a**2*d**2 - 2*a*b*c*d*(3*m + 7) + b**2*c**2*(9*m**2 + 33*m + 28))*hyper((sympy.S(1)/3, -m), (sympy.S(4)/3,), -b*x**3/a)/(b**2*(1 + b*x**3/a)**m*(3*m + 4)*(3*m + 7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_39():
    f = (a + b*x**3)**m*(c + d*x**3)
    F = d*x*(a + b*x**3)**(m + 1)/(b*(3*m + 4)) - x*(a + b*x**3)**m*(a*d - b*c*(3*m + 4))*hyper((sympy.S(1)/3, -m), (sympy.S(4)/3,), -b*x**3/a)/(b*(1 + b*x**3/a)**m*(3*m + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_40():
    f = (a + b*x**3)**m
    F = x*(a + b*x**3)**m*hyper((sympy.S(1)/3, -m), (sympy.S(4)/3,), -b*x**3/a)/(1 + b*x**3/a)**m
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_41():
    f = (a + b*x**3)**m/(c + d*x**3)
    F = x*(a + b*x**3)**m*appellf1(sympy.S(1)/3, 1, -m, sympy.S(4)/3, -d*x**3/c, -b*x**3/a)/(c*(1 + b*x**3/a)**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_42():
    f = (a + b*x**3)**m/(c + d*x**3)**2
    F = x*(a + b*x**3)**m*appellf1(sympy.S(1)/3, 2, -m, sympy.S(4)/3, -d*x**3/c, -b*x**3/a)/(c**2*(1 + b*x**3/a)**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_43():
    f = (a + b*x**3)**m/(c + d*x**3)**3
    F = x*(a + b*x**3)**m*appellf1(sympy.S(1)/3, 3, -m, sympy.S(4)/3, -d*x**3/c, -b*x**3/a)/(c**3*(1 + b*x**3/a)**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_44():
    f = (a + b*x**3)**(-b*c/(-3*a*d + 3*b*c) - 1)*(c + d*x**3)**(a*d/(-3*a*d + 3*b*c) - 1)
    F = x*(c + d*x**3)**(a*d/(-3*a*d + 3*b*c))/(a*c*(a + b*x**3)**(b*c/(-3*a*d + 3*b*c)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_45():
    f = (a + b*x**4)*(c + d*x**4)**4
    F = a*c**4*x + b*d**4*x**21/21 + c**3*x**5*(4*a*d + b*c)/5 + 2*c**2*d*x**9*(3*a*d + 2*b*c)/9 + 2*c*d**2*x**13*(2*a*d + 3*b*c)/13 + d**3*x**17*(a*d + 4*b*c)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_46():
    f = (a + b*x**4)*(c + d*x**4)**3
    F = a*c**3*x + b*d**3*x**17/17 + c**2*x**5*(3*a*d + b*c)/5 + c*d*x**9*(a*d + b*c)/3 + d**2*x**13*(a*d + 3*b*c)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_47():
    f = (a + b*x**4)*(c + d*x**4)**2
    F = a*c**2*x + b*d**2*x**13/13 + c*x**5*(2*a*d + b*c)/5 + d*x**9*(a*d + 2*b*c)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_48():
    f = (a + b*x**4)*(c + d*x**4)
    F = a*c*x + b*d*x**9/9 + x**5*(a*d/5 + b*c/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_49():
    f = (a + b*x**4)/(c + d*x**4)
    F = b*x/d + sqrt(2)*(-a*d + b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*c**(sympy.S(3)/4)*d**(sympy.S(5)/4)) - sqrt(2)*(-a*d + b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*c**(sympy.S(3)/4)*d**(sympy.S(5)/4)) + sqrt(2)*(-a*d + b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*d**(sympy.S(5)/4)) - sqrt(2)*(-a*d + b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*d**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_50():
    f = (a + b*x**4)/(c + d*x**4)**2
    F = -x*(-a*d + b*c)/(4*c*d*(c + d*x**4)) - sqrt(2)*(3*a*d + b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(32*c**(sympy.S(7)/4)*d**(sympy.S(5)/4)) + sqrt(2)*(3*a*d + b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(32*c**(sympy.S(7)/4)*d**(sympy.S(5)/4)) - sqrt(2)*(3*a*d + b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(16*c**(sympy.S(7)/4)*d**(sympy.S(5)/4)) + sqrt(2)*(3*a*d + b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(16*c**(sympy.S(7)/4)*d**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_51():
    f = (a + b*x**4)/(c + d*x**4)**3
    F = -x*(-a*d + b*c)/(8*c*d*(c + d*x**4)**2) + x*(7*a*d + b*c)/(32*c**2*d*(c + d*x**4)) - sqrt(2)*(21*a*d + 3*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(256*c**(sympy.S(11)/4)*d**(sympy.S(5)/4)) + sqrt(2)*(21*a*d + 3*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(256*c**(sympy.S(11)/4)*d**(sympy.S(5)/4)) - sqrt(2)*(21*a*d + 3*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(128*c**(sympy.S(11)/4)*d**(sympy.S(5)/4)) + sqrt(2)*(21*a*d + 3*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(128*c**(sympy.S(11)/4)*d**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_52():
    f = (a + b*x**4)**2*(c + d*x**4)**4
    F = a**2*c**4*x + 2*a*c**3*x**5*(2*a*d + b*c)/5 + b**2*d**4*x**25/25 + 2*b*d**3*x**21*(a*d + 2*b*c)/21 + c**2*x**9*(6*a**2*d**2 + 8*a*b*c*d + b**2*c**2)/9 + 4*c*d*x**13*(a**2*d**2 + 3*a*b*c*d + b**2*c**2)/13 + d**2*x**17*(a**2*d**2 + 8*a*b*c*d + 6*b**2*c**2)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_53():
    f = (a + b*x**4)**2*(c + d*x**4)**3
    F = a**2*c**3*x + a*c**2*x**5*(3*a*d + 2*b*c)/5 + b**2*d**3*x**21/21 + b*d**2*x**17*(2*a*d + 3*b*c)/17 + c*x**9*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2)/9 + d*x**13*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_54():
    f = (a + b*x**4)**2*(c + d*x**4)**2
    F = a**2*c**2*x + 2*a*c*x**5*(a*d + b*c)/5 + b**2*d**2*x**17/17 + 2*b*d*x**13*(a*d + b*c)/13 + x**9*(a**2*d**2/9 + 4*a*b*c*d/9 + b**2*c**2/9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_55():
    f = (a + b*x**4)**2*(c + d*x**4)
    F = a**2*c*x + a*x**5*(a*d + 2*b*c)/5 + b**2*d*x**13/13 + b*x**9*(2*a*d + b*c)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_56():
    f = (a + b*x**4)**2/(c + d*x**4)
    F = b**2*x**5/(5*d) - b*x*(-2*a*d + b*c)/d**2 - sqrt(2)*(-a*d + b*c)**2*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*c**(sympy.S(3)/4)*d**(sympy.S(9)/4)) + sqrt(2)*(-a*d + b*c)**2*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*c**(sympy.S(3)/4)*d**(sympy.S(9)/4)) - sqrt(2)*(-a*d + b*c)**2*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*d**(sympy.S(9)/4)) + sqrt(2)*(-a*d + b*c)**2*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*d**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_57():
    f = (a + b*x**4)**2/(c + d*x**4)**2
    F = b**2*x/d**2 + x*(-a*d + b*c)**2/(4*c*d**2*(c + d*x**4)) + sqrt(2)*(-a*d + b*c)*(3*a*d + 5*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(32*c**(sympy.S(7)/4)*d**(sympy.S(9)/4)) - sqrt(2)*(-a*d + b*c)*(3*a*d + 5*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(32*c**(sympy.S(7)/4)*d**(sympy.S(9)/4)) + sqrt(2)*(-a*d + b*c)*(3*a*d + 5*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(16*c**(sympy.S(7)/4)*d**(sympy.S(9)/4)) - sqrt(2)*(-a*d + b*c)*(3*a*d + 5*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(16*c**(sympy.S(7)/4)*d**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_58():
    f = (a + b*x**4)**2/(c + d*x**4)**3
    F = -x*(a + b*x**4)*(-a*d + b*c)/(8*c*d*(c + d*x**4)**2) - x*(-a*d + b*c)*(7*a*d + 5*b*c)/(32*c**2*d**2*(c + d*x**4)) - sqrt(2)*(21*a**2*d**2 + 6*a*b*c*d + 5*b**2*c**2)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(256*c**(sympy.S(11)/4)*d**(sympy.S(9)/4)) + sqrt(2)*(21*a**2*d**2 + 6*a*b*c*d + 5*b**2*c**2)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(256*c**(sympy.S(11)/4)*d**(sympy.S(9)/4)) - sqrt(2)*(21*a**2*d**2 + 6*a*b*c*d + 5*b**2*c**2)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(128*c**(sympy.S(11)/4)*d**(sympy.S(9)/4)) + sqrt(2)*(21*a**2*d**2 + 6*a*b*c*d + 5*b**2*c**2)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(128*c**(sympy.S(11)/4)*d**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_59():
    f = (c + d*x**4)**4/(a + b*x**4)
    F = d**4*x**13/(13*b) + d**3*x**9*(-a*d + 4*b*c)/(9*b**2) + d**2*x**5*(a**2*d**2 - 4*a*b*c*d + 6*b**2*c**2)/(5*b**3) + d*x*(-a*d + 2*b*c)*(a**2*d**2 - 2*a*b*c*d + 2*b**2*c**2)/b**4 - sqrt(2)*(-a*d + b*c)**4*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(17)/4)) + sqrt(2)*(-a*d + b*c)**4*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(17)/4)) - sqrt(2)*(-a*d + b*c)**4*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(17)/4)) + sqrt(2)*(-a*d + b*c)**4*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(17)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_60():
    f = (c + d*x**4)**3/(a + b*x**4)
    F = d**3*x**9/(9*b) + d**2*x**5*(-a*d + 3*b*c)/(5*b**2) + d*x*(a**2*d**2 - 3*a*b*c*d + 3*b**2*c**2)/b**3 - sqrt(2)*(-a*d + b*c)**3*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(13)/4)) + sqrt(2)*(-a*d + b*c)**3*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(13)/4)) - sqrt(2)*(-a*d + b*c)**3*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(13)/4)) + sqrt(2)*(-a*d + b*c)**3*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(13)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_61():
    f = (c + d*x**4)**2/(a + b*x**4)
    F = d**2*x**5/(5*b) + d*x*(-a*d + 2*b*c)/b**2 - sqrt(2)*(-a*d + b*c)**2*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(9)/4)) + sqrt(2)*(-a*d + b*c)**2*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(9)/4)) - sqrt(2)*(-a*d + b*c)**2*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(9)/4)) + sqrt(2)*(-a*d + b*c)**2*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_62():
    f = (c + d*x**4)/(a + b*x**4)
    F = d*x/b - sqrt(2)*(-a*d + b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(-a*d + b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)) - sqrt(2)*(-a*d + b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(-a*d + b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_63():
    f = 1/((a + b*x**4)*(c + d*x**4))
    F = sqrt(2)*d**(sympy.S(3)/4)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*c**(sympy.S(3)/4)*(-a*d + b*c)) - sqrt(2)*d**(sympy.S(3)/4)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*c**(sympy.S(3)/4)*(-a*d + b*c)) + sqrt(2)*d**(sympy.S(3)/4)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*(-a*d + b*c)) - sqrt(2)*d**(sympy.S(3)/4)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(3)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*(-a*d + b*c)) + sqrt(2)*b**(sympy.S(3)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(3)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(-a*d + b*c)) + sqrt(2)*b**(sympy.S(3)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_64():
    f = 1/((a + b*x**4)*(c + d*x**4)**2)
    F = -d*x/(4*c*(c + d*x**4)*(-a*d + b*c)) + sqrt(2)*d**(sympy.S(3)/4)*(-3*a*d + 7*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(32*c**(sympy.S(7)/4)*(-a*d + b*c)**2) - sqrt(2)*d**(sympy.S(3)/4)*(-3*a*d + 7*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(32*c**(sympy.S(7)/4)*(-a*d + b*c)**2) + sqrt(2)*d**(sympy.S(3)/4)*(-3*a*d + 7*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(16*c**(sympy.S(7)/4)*(-a*d + b*c)**2) - sqrt(2)*d**(sympy.S(3)/4)*(-3*a*d + 7*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(16*c**(sympy.S(7)/4)*(-a*d + b*c)**2) - sqrt(2)*b**(sympy.S(7)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*(-a*d + b*c)**2) + sqrt(2)*b**(sympy.S(7)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*(-a*d + b*c)**2) - sqrt(2)*b**(sympy.S(7)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(-a*d + b*c)**2) + sqrt(2)*b**(sympy.S(7)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_65():
    f = (c + d*x**4)**5/(a + b*x**4)**2
    F = d**5*x**13/(13*b**2) + d**4*x**9*(-2*a*d + 5*b*c)/(9*b**3) + d**3*x**5*(3*a**2*d**2 - 10*a*b*c*d + 10*b**2*c**2)/(5*b**4) + d**2*x*(-4*a**3*d**3 + 15*a**2*b*c*d**2 - 20*a*b**2*c**2*d + 10*b**3*c**3)/b**5 + x*(-a*d + b*c)**5/(4*a*b**5*(a + b*x**4)) - sqrt(2)*(-a*d + b*c)**4*(17*a*d + 3*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(21)/4)) + sqrt(2)*(-a*d + b*c)**4*(17*a*d + 3*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(21)/4)) - sqrt(2)*(-a*d + b*c)**4*(17*a*d + 3*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(21)/4)) + sqrt(2)*(-a*d + b*c)**4*(17*a*d + 3*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(21)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_66():
    f = (c + d*x**4)**4/(a + b*x**4)**2
    F = d**4*x**9/(9*b**2) + 2*d**3*x**5*(-a*d + 2*b*c)/(5*b**3) + d**2*x*(3*a**2*d**2 - 8*a*b*c*d + 6*b**2*c**2)/b**4 + x*(-a*d + b*c)**4/(4*a*b**4*(a + b*x**4)) - sqrt(2)*(-a*d + b*c)**3*(13*a*d + 3*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(17)/4)) + sqrt(2)*(-a*d + b*c)**3*(13*a*d + 3*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(17)/4)) - sqrt(2)*(-a*d + b*c)**3*(13*a*d + 3*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(17)/4)) + sqrt(2)*(-a*d + b*c)**3*(13*a*d + 3*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(17)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_67():
    f = (c + d*x**4)**3/(a + b*x**4)**2
    F = d**3*x**5/(5*b**2) + d**2*x*(-2*a*d + 3*b*c)/b**3 + x*(-a*d + b*c)**3/(4*a*b**3*(a + b*x**4)) - 3*sqrt(2)*(-a*d + b*c)**2*(3*a*d + b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(13)/4)) + 3*sqrt(2)*(-a*d + b*c)**2*(3*a*d + b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(13)/4)) - 3*sqrt(2)*(-a*d + b*c)**2*(3*a*d + b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(13)/4)) + 3*sqrt(2)*(-a*d + b*c)**2*(3*a*d + b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(13)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_68():
    f = (c + d*x**4)**2/(a + b*x**4)**2
    F = d**2*x/b**2 + x*(-a*d + b*c)**2/(4*a*b**2*(a + b*x**4)) - sqrt(2)*(-a*d + b*c)*(5*a*d + 3*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(9)/4)) + sqrt(2)*(-a*d + b*c)*(5*a*d + 3*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(9)/4)) - sqrt(2)*(-a*d + b*c)*(5*a*d + 3*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(9)/4)) + sqrt(2)*(-a*d + b*c)*(5*a*d + 3*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_69():
    f = (c + d*x**4)/(a + b*x**4)**2
    F = x*(-a*d + b*c)/(4*a*b*(a + b*x**4)) - sqrt(2)*(a*d + 3*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(a*d + 3*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*b**(sympy.S(5)/4)) - sqrt(2)*(a*d + 3*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(5)/4)) + sqrt(2)*(a*d + 3*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*b**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_70():
    f = 1/((a + b*x**4)**2*(c + d*x**4))
    F = -sqrt(2)*d**(sympy.S(7)/4)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*c**(sympy.S(3)/4)*(-a*d + b*c)**2) + sqrt(2)*d**(sympy.S(7)/4)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(8*c**(sympy.S(3)/4)*(-a*d + b*c)**2) - sqrt(2)*d**(sympy.S(7)/4)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*(-a*d + b*c)**2) + sqrt(2)*d**(sympy.S(7)/4)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*(-a*d + b*c)**2) + b*x/(4*a*(a + b*x**4)*(-a*d + b*c)) - sqrt(2)*b**(sympy.S(3)/4)*(-7*a*d + 3*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*(-a*d + b*c)**2) + sqrt(2)*b**(sympy.S(3)/4)*(-7*a*d + 3*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*(-a*d + b*c)**2) - sqrt(2)*b**(sympy.S(3)/4)*(-7*a*d + 3*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(-a*d + b*c)**2) + sqrt(2)*b**(sympy.S(3)/4)*(-7*a*d + 3*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_71():
    f = 1/((a + b*x**4)**2*(c + d*x**4)**2)
    F = -sqrt(2)*d**(sympy.S(7)/4)*(-3*a*d + 11*b*c)*log(-sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(32*c**(sympy.S(7)/4)*(-a*d + b*c)**3) + sqrt(2)*d**(sympy.S(7)/4)*(-3*a*d + 11*b*c)*log(sqrt(2)*c**(sympy.S(1)/4)*d**(sympy.S(1)/4)*x + sqrt(c) + sqrt(d)*x**2)/(32*c**(sympy.S(7)/4)*(-a*d + b*c)**3) - sqrt(2)*d**(sympy.S(7)/4)*(-3*a*d + 11*b*c)*atan(1 - sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(16*c**(sympy.S(7)/4)*(-a*d + b*c)**3) + sqrt(2)*d**(sympy.S(7)/4)*(-3*a*d + 11*b*c)*atan(1 + sqrt(2)*d**(sympy.S(1)/4)*x/c**(sympy.S(1)/4))/(16*c**(sympy.S(7)/4)*(-a*d + b*c)**3) + b*x/(4*a*(a + b*x**4)*(c + d*x**4)*(-a*d + b*c)) + d*x*(a*d + b*c)/(4*a*c*(c + d*x**4)*(-a*d + b*c)**2) - sqrt(2)*b**(sympy.S(7)/4)*(-11*a*d + 3*b*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(7)/4)*(-11*a*d + 3*b*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(32*a**(sympy.S(7)/4)*(-a*d + b*c)**3) - sqrt(2)*b**(sympy.S(7)/4)*(-11*a*d + 3*b*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(-a*d + b*c)**3) + sqrt(2)*b**(sympy.S(7)/4)*(-11*a*d + 3*b*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_72():
    f = (a - b*x**4)**(sympy.S(5)/2)/(c - d*x**4)
    F = (Integer(-1) * ((Symbol('b') * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(13) * Symbol('a') * Symbol('d')))) * x * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))) * ((Integer(21) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Symbol('b') * x * ((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(7) * Symbol('d')))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(21) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(56) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(47) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(21) * (Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * (Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * (Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_73():
    f = (a - b*x**4)**(sympy.S(3)/2)/(c - d*x**4)
    F = ((Symbol('b') * x * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1)))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_74():
    f = sqrt(a - b*x**4)/(c - d*x**4)
    F = (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * Symbol('d') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * Symbol('d') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_75():
    f = 1/(sqrt(a - b*x**4)*(c - d*x**4))
    F = (((Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_76():
    f = 1/((a - b*x**4)**(sympy.S(3)/2)*(c - d*x**4))
    F = ((Symbol('b') * x) * ((Integer(2) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))) + (((Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * Symbol('d') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * Symbol('d') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_77():
    f = 1/((a - b*x**4)**(sympy.S(5)/2)*(c - d*x**4))
    F = ((Symbol('b') * x) * ((Integer(6) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(11) * Symbol('a') * Symbol('d')))) * x) * ((Integer(12) * (Symbol('a'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))) + (((Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(11) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(12) * (Symbol('a'))**((Integer(7) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_78():
    f = (a + b*x**4)**(sympy.S(3)/2)/(c + d*x**4)
    F = ((Symbol('b') * x * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * x) * ((((Integer(-1) * Symbol('c')))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(7) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * x) * ((((Integer(-1) * Symbol('c')))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(7) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(6) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1)))) + (((Symbol('b'))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * Symbol('c'))) * (Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))) + (((Symbol('b'))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d')))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * Symbol('c'))) * (Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * (Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * (Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_79():
    f = sqrt(a + b*x**4)/(c + d*x**4)
    F = ((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * x) * ((((Integer(-1) * Symbol('c')))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * x) * ((((Integer(-1) * Symbol('c')))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * Symbol('c'))) * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('b')) + ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((Integer(-1) * Symbol('c'))))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * (((((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * (((((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_80():
    f = 1/(sqrt(a + b*x**4)*(c + d*x**4))
    F = (Integer(-1) * (((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * x) * ((((Integer(-1) * Symbol('c')))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * x) * ((((Integer(-1) * Symbol('c')))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))))**(Integer(-1)))) + (((Symbol('b'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('b')) + ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((Integer(-1) * Symbol('c'))))**(Integer(-1)))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))) + (((Symbol('b'))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('b')) * Symbol('c')) + (sympy.sqrt(Symbol('a')) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d')))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_81():
    f = 1/((a + b*x**4)**(sympy.S(3)/2)*(c + d*x**4))
    F = ((Symbol('b') * x) * ((Integer(2) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))) + (((Symbol('d'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * x) * ((((Integer(-1) * Symbol('c')))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * x) * ((((Integer(-1) * Symbol('c')))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('b')) + ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((Integer(-1) * Symbol('c'))))**(Integer(-1)))) * Symbol('d') * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('b')) * Symbol('c')) + (sympy.sqrt(Symbol('a')) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d')))) * Symbol('d') * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * Symbol('c') * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * (((((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * Symbol('d') * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * (((((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * Symbol('d') * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_82():
    f = 1/((a + b*x**4)**(sympy.S(5)/2)*(c + d*x**4))
    F = ((Symbol('b') * x) * ((Integer(6) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(11) * Symbol('a') * Symbol('d')))) * x) * ((Integer(12) * (Symbol('a'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**((Integer(9) * (Integer(4))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * x) * ((((Integer(-1) * Symbol('c')))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**((Integer(9) * (Integer(4))**(Integer(-1)))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * x) * ((((Integer(-1) * Symbol('c')))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(11) * Symbol('a') * Symbol('d')))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(24) * (Symbol('a'))**((Integer(9) * (Integer(4))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))) + (((Symbol('b'))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('b')) * Symbol('c')) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))) * (Symbol('d'))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))) + (((Symbol('b'))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('b')) * Symbol('c')) + (sympy.sqrt(Symbol('a')) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d')))) * (Symbol('d'))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))) + (((((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_83():
    f = (a - b*x**4)**(sympy.S(7)/2)/(c - d*x**4)**2
    F = (Integer(-1) * ((Symbol('b') * ((Integer(77) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(122) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(21) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) * x * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))) * ((Integer(84) * Symbol('c') * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + ((Symbol('b') * ((Integer(11) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('d')))) * x * ((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(28) * Symbol('c') * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * x * ((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(4) * Symbol('c') * Symbol('d') * (Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + (((Integer(84) * Symbol('c') * (Symbol('d'))**(Integer(4)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1)) * ((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(231) * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3))) + (Integer(-1) * (Integer(553) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * Symbol('d'))) + (Integer(349) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(21) * (Symbol('a'))**(Integer(3)) * (Symbol('d'))**(Integer(3)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * ((Integer(11) * Symbol('b') * Symbol('c')) + (Integer(3) * Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(4)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * ((Integer(11) * Symbol('b') * Symbol('c')) + (Integer(3) * Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(4)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_84():
    f = (a - b*x**4)**(sympy.S(5)/2)/(c - d*x**4)**2
    F = ((Symbol('b') * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * x * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))) * ((Integer(12) * Symbol('c') * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * x * ((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(4) * Symbol('c') * Symbol('d') * (Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(21) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(26) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(12) * Symbol('c') * (Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1)))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(3) * Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(3) * Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_85():
    f = (a - b*x**4)**(sympy.S(3)/2)/(c - d*x**4)**2
    F = (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * x * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))) * ((Integer(4) * Symbol('c') * Symbol('d') * (Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_86():
    f = sqrt(a - b*x**4)/(c - d*x**4)**2
    F = ((x * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))) * ((Integer(4) * Symbol('c') * (Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * Symbol('d') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_87():
    f = 1/(sqrt(a - b*x**4)*(c - d*x**4)**2)
    F = (Integer(-1) * ((Symbol('d') * x * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))) * ((Integer(4) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1)))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_88():
    f = 1/((a - b*x**4)**(sympy.S(3)/2)*(c - d*x**4)**2)
    F = ((Symbol('b') * ((Integer(2) * Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * x) * ((Integer(4) * Symbol('a') * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * x) * ((Integer(4) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4)))))) * (Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + (((Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(2) * Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * (Symbol('a'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * Symbol('d') * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * Symbol('d') * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_89():
    f = 1/((a - b*x**4)**(sympy.S(5)/2)*(c - d*x**4)**2)
    F = ((Symbol('b') * ((Integer(2) * Symbol('b') * Symbol('c')) + (Integer(3) * Symbol('a') * Symbol('d'))) * x) * ((Integer(12) * Symbol('a') * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(5) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(17) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))) * x) * ((Integer(12) * (Symbol('a'))**(Integer(2)) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * x) * ((Integer(4) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c') + (Integer(-1) * (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + (((Symbol('b'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(5) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(17) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(12) * (Symbol('a'))**((Integer(7) * (Integer(4))**(Integer(-1)))) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**(Integer(2)) * ((Integer(13) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**(Integer(2)) * ((Integer(13) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(4))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_90():
    f = sqrt(a + b*x**4)/(a*c - b*c*x**4)
    F = sqrt(2)*atan(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x/sqrt(a + b*x**4))/(4*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*c) + sqrt(2)*atanh(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x/sqrt(a + b*x**4))/(4*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_91():
    f = sqrt(a - b*x**4)/(a*c + b*c*x**4)
    F = atan(b**(sympy.S(1)/4)*x*(sqrt(a) + sqrt(b)*x**2)/(a**(sympy.S(1)/4)*sqrt(a - b*x**4)))/(2*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*c) + atanh(b**(sympy.S(1)/4)*x*(sqrt(a) - sqrt(b)*x**2)/(a**(sympy.S(1)/4)*sqrt(a - b*x**4)))/(2*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_92():
    f = (a + b*x**4)**(sympy.S(7)/4)/(c + d*x**4)
    F = -b**(sympy.S(3)/4)*(-7*a*d + 4*b*c)*atan(b**(sympy.S(1)/4)*x/(a + b*x**4)**(sympy.S(1)/4))/(8*d**2) - b**(sympy.S(3)/4)*(-7*a*d + 4*b*c)*atanh(b**(sympy.S(1)/4)*x/(a + b*x**4)**(sympy.S(1)/4))/(8*d**2) + b*x*(a + b*x**4)**(sympy.S(3)/4)/(4*d) + (-a*d + b*c)**(sympy.S(7)/4)*atan(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(2*c**(sympy.S(3)/4)*d**2) + (-a*d + b*c)**(sympy.S(7)/4)*atanh(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(2*c**(sympy.S(3)/4)*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_93():
    f = (a + b*x**4)**(sympy.S(3)/4)/(c + d*x**4)
    F = b**(sympy.S(3)/4)*atan(b**(sympy.S(1)/4)*x/(a + b*x**4)**(sympy.S(1)/4))/(2*d) + b**(sympy.S(3)/4)*atanh(b**(sympy.S(1)/4)*x/(a + b*x**4)**(sympy.S(1)/4))/(2*d) - (-a*d + b*c)**(sympy.S(3)/4)*atan(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(2*c**(sympy.S(3)/4)*d) - (-a*d + b*c)**(sympy.S(3)/4)*atanh(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(2*c**(sympy.S(3)/4)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_94():
    f = 1/((a + b*x**4)**(sympy.S(1)/4)*(c + d*x**4))
    F = atan(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(2*c**(sympy.S(3)/4)*(-a*d + b*c)**(sympy.S(1)/4)) + atanh(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(2*c**(sympy.S(3)/4)*(-a*d + b*c)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_95():
    f = 1/((a + b*x**4)**(sympy.S(5)/4)*(c + d*x**4))
    F = -d*atan(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(2*c**(sympy.S(3)/4)*(-a*d + b*c)**(sympy.S(5)/4)) - d*atanh(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(2*c**(sympy.S(3)/4)*(-a*d + b*c)**(sympy.S(5)/4)) + b*x/(a*(a + b*x**4)**(sympy.S(1)/4)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_96():
    f = 1/((a + b*x**4)**(sympy.S(9)/4)*(c + d*x**4))
    F = d**2*atan(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(2*c**(sympy.S(3)/4)*(-a*d + b*c)**(sympy.S(9)/4)) + d**2*atanh(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(2*c**(sympy.S(3)/4)*(-a*d + b*c)**(sympy.S(9)/4)) + b*x/(5*a*(a + b*x**4)**(sympy.S(5)/4)*(-a*d + b*c)) + b*x*(-9*a*d + 4*b*c)/(5*a**2*(a + b*x**4)**(sympy.S(1)/4)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_97():
    f = 1/((a + b*x**4)**(sympy.S(13)/4)*(c + d*x**4))
    F = -d**3*atan(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(2*c**(sympy.S(3)/4)*(-a*d + b*c)**(sympy.S(13)/4)) - d**3*atanh(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(2*c**(sympy.S(3)/4)*(-a*d + b*c)**(sympy.S(13)/4)) + b*x/(9*a*(a + b*x**4)**(sympy.S(9)/4)*(-a*d + b*c)) + b*x*(-17*a*d + 8*b*c)/(45*a**2*(a + b*x**4)**(sympy.S(5)/4)*(-a*d + b*c)**2) + b*x*(113*a**2*d**2 - 100*a*b*c*d + 32*b**2*c**2)/(45*a**3*(a + b*x**4)**(sympy.S(1)/4)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_98():
    f = (a + b*x**4)**(sympy.S(9)/4)/(c + d*x**4)
    F = (Integer(-1) * ((Symbol('b') * ((Integer(6) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(11) * Symbol('a') * Symbol('d')))) * x * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1)))) * ((Integer(12) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Symbol('b') * x * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(5) * (Integer(4))**(Integer(-1))))) * ((Integer(6) * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt(Symbol('a')) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(6) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(11) * Symbol('a') * Symbol('d')))) * ((Integer(1) + (Symbol('a') * ((Symbol('b') * (x)**(Integer(4))))**(Integer(-1)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (x)**(Integer(3)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.acot(((sympy.sqrt(Symbol('b')) * (x)**(Integer(2))) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(12) * (Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((Integer(-1) * (sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * (Symbol('d'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_99():
    f = (a + b*x**4)**(sympy.S(5)/4)/(c + d*x**4)
    F = ((Symbol('b') * x * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('a')) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(1) + (Symbol('a') * ((Symbol('b') * (x)**(Integer(4))))**(Integer(-1)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (x)**(Integer(3)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.acot(((sympy.sqrt(Symbol('b')) * (x)**(Integer(2))) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(2) * Symbol('d') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((Integer(-1) * (sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_100():
    f = (a + b*x**4)**(sympy.S(1)/4)/(c + d*x**4)
    F = ((sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((Integer(-1) * (sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_101():
    f = 1/((a + b*x**4)**(sympy.S(3)/4)*(c + d*x**4))
    F = (Integer(-1) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(1) + (Symbol('a') * ((Symbol('b') * (x)**(Integer(4))))**(Integer(-1)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (x)**(Integer(3)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.acot(((sympy.sqrt(Symbol('b')) * (x)**(Integer(2))) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((sympy.sqrt(Symbol('a')) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((Integer(-1) * (sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_102():
    f = 1/((a + b*x**4)**(sympy.S(7)/4)*(c + d*x**4))
    F = ((Symbol('b') * x) * ((Integer(3) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(2) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('d')))) * ((Integer(1) + (Symbol('a') * ((Symbol('b') * (x)**(Integer(4))))**(Integer(-1)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (x)**(Integer(3)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.acot(((sympy.sqrt(Symbol('b')) * (x)**(Integer(2))) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(3) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((Integer(-1) * (sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2))))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_103():
    f = 1/((a + b*x**4)**(sympy.S(11)/4)*(c + d*x**4))
    F = ((Symbol('b') * x) * ((Integer(7) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(7) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(6) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(13) * Symbol('a') * Symbol('d')))) * x) * ((Integer(21) * (Symbol('a'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(12) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(38) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(47) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) * ((Integer(1) + (Symbol('a') * ((Symbol('b') * (x)**(Integer(4))))**(Integer(-1)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (x)**(Integer(3)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.acot(((sympy.sqrt(Symbol('b')) * (x)**(Integer(2))) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(21) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((Integer(-1) * (sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(3)) * sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_104():
    f = (a + b*x**4)**(sympy.S(11)/4)/(c + d*x**4)**2
    F = -b**(sympy.S(7)/4)*(-11*a*d + 8*b*c)*atan(b**(sympy.S(1)/4)*x/(a + b*x**4)**(sympy.S(1)/4))/(8*d**3) - b**(sympy.S(7)/4)*(-11*a*d + 8*b*c)*atanh(b**(sympy.S(1)/4)*x/(a + b*x**4)**(sympy.S(1)/4))/(8*d**3) + b*x*(a + b*x**4)**(sympy.S(3)/4)*(-a*d + 2*b*c)/(4*c*d**2) - x*(a + b*x**4)**(sympy.S(7)/4)*(-a*d + b*c)/(4*c*d*(c + d*x**4)) + (-a*d + b*c)**(sympy.S(7)/4)*(3*a*d + 8*b*c)*atan(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(8*c**(sympy.S(7)/4)*d**3) + (-a*d + b*c)**(sympy.S(7)/4)*(3*a*d + 8*b*c)*atanh(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(8*c**(sympy.S(7)/4)*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_105():
    f = (a + b*x**4)**(sympy.S(7)/4)/(c + d*x**4)**2
    F = b**(sympy.S(7)/4)*atan(b**(sympy.S(1)/4)*x/(a + b*x**4)**(sympy.S(1)/4))/(2*d**2) + b**(sympy.S(7)/4)*atanh(b**(sympy.S(1)/4)*x/(a + b*x**4)**(sympy.S(1)/4))/(2*d**2) - x*(a + b*x**4)**(sympy.S(3)/4)*(-a*d + b*c)/(4*c*d*(c + d*x**4)) - (-a*d + b*c)**(sympy.S(3)/4)*(3*a*d + 4*b*c)*atan(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(8*c**(sympy.S(7)/4)*d**2) - (-a*d + b*c)**(sympy.S(3)/4)*(3*a*d + 4*b*c)*atanh(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(8*c**(sympy.S(7)/4)*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_106():
    f = (a + b*x**4)**(sympy.S(3)/4)/(c + d*x**4)**2
    F = 3*a*atan(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(8*c**(sympy.S(7)/4)*(-a*d + b*c)**(sympy.S(1)/4)) + 3*a*atanh(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(8*c**(sympy.S(7)/4)*(-a*d + b*c)**(sympy.S(1)/4)) + x*(a + b*x**4)**(sympy.S(3)/4)/(4*c*(c + d*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_107():
    f = 1/((a + b*x**4)**(sympy.S(1)/4)*(c + d*x**4)**2)
    F = -d*x*(a + b*x**4)**(sympy.S(3)/4)/(4*c*(c + d*x**4)*(-a*d + b*c)) + (-3*a*d + 4*b*c)*atan(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(8*c**(sympy.S(7)/4)*(-a*d + b*c)**(sympy.S(5)/4)) + (-3*a*d + 4*b*c)*atanh(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(8*c**(sympy.S(7)/4)*(-a*d + b*c)**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_108():
    f = 1/((a + b*x**4)**(sympy.S(5)/4)*(c + d*x**4)**2)
    F = -d*x/(4*c*(a + b*x**4)**(sympy.S(1)/4)*(c + d*x**4)*(-a*d + b*c)) - d*(-3*a*d + 8*b*c)*atan(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(8*c**(sympy.S(7)/4)*(-a*d + b*c)**(sympy.S(9)/4)) - d*(-3*a*d + 8*b*c)*atanh(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(8*c**(sympy.S(7)/4)*(-a*d + b*c)**(sympy.S(9)/4)) + b*x*(a*d + 4*b*c)/(4*a*c*(a + b*x**4)**(sympy.S(1)/4)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_109():
    f = 1/((a + b*x**4)**(sympy.S(9)/4)*(c + d*x**4)**2)
    F = -d*x/(4*c*(a + b*x**4)**(sympy.S(5)/4)*(c + d*x**4)*(-a*d + b*c)) + 3*d**2*(-a*d + 4*b*c)*atan(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(8*c**(sympy.S(7)/4)*(-a*d + b*c)**(sympy.S(13)/4)) + 3*d**2*(-a*d + 4*b*c)*atanh(x*(-a*d + b*c)**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*(a + b*x**4)**(sympy.S(1)/4)))/(8*c**(sympy.S(7)/4)*(-a*d + b*c)**(sympy.S(13)/4)) + b*x*(5*a*d + 4*b*c)/(20*a*c*(a + b*x**4)**(sympy.S(5)/4)*(-a*d + b*c)**2) + b*x*(-5*a**2*d**2 - 56*a*b*c*d + 16*b**2*c**2)/(20*a**2*c*(a + b*x**4)**(sympy.S(1)/4)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_110():
    f = (a + b*x**4)**(sympy.S(9)/4)/(c + d*x**4)**2
    F = ((Symbol('b') * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * x * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1)))) * ((Integer(4) * Symbol('c') * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * x * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(5) * (Integer(4))**(Integer(-1))))) * ((Integer(4) * Symbol('c') * Symbol('d') * (Symbol('c') + (Symbol('d') * (x)**(Integer(4))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('a')) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Integer(1) + (Symbol('a') * ((Symbol('b') * (x)**(Integer(4))))**(Integer(-1)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (x)**(Integer(3)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.acot(((sympy.sqrt(Symbol('b')) * (x)**(Integer(2))) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(4) * Symbol('c') * (Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Integer(2) * Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((Integer(-1) * (sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Integer(2) * Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_111():
    f = (a + b*x**4)**(sympy.S(5)/4)/(c + d*x**4)**2
    F = (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * x * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1)))) * ((Integer(4) * Symbol('c') * Symbol('d') * (Symbol('c') + (Symbol('d') * (x)**(Integer(4))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('a')) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(1) + (Symbol('a') * ((Symbol('b') * (x)**(Integer(4))))**(Integer(-1)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (x)**(Integer(3)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.acot(((sympy.sqrt(Symbol('b')) * (x)**(Integer(2))) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(4) * Symbol('c') * Symbol('d') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(2) * Symbol('b') * Symbol('c')) + (Integer(3) * Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((Integer(-1) * (sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Integer(2) * Symbol('b') * Symbol('c')) + (Integer(3) * Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_112():
    f = (a + b*x**4)**(sympy.S(1)/4)/(c + d*x**4)**2
    F = ((x * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1)))) * ((Integer(4) * Symbol('c') * (Symbol('c') + (Symbol('d') * (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('a')) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(1) + (Symbol('a') * ((Symbol('b') * (x)**(Integer(4))))**(Integer(-1)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (x)**(Integer(3)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.acot(((sympy.sqrt(Symbol('b')) * (x)**(Integer(2))) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(4) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + ((((Integer(2) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((Integer(-1) * (sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1))) + ((((Integer(2) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_113():
    f = 1/((a + b*x**4)**(sympy.S(3)/4)*(c + d*x**4)**2)
    F = (Integer(-1) * ((Symbol('d') * x * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1)))) * ((Integer(4) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('c') + (Symbol('d') * (x)**(Integer(4))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(4) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Integer(1) + (Symbol('a') * ((Symbol('b') * (x)**(Integer(4))))**(Integer(-1)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (x)**(Integer(3)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.acot(((sympy.sqrt(Symbol('b')) * (x)**(Integer(2))) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('d') * ((Integer(2) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((Integer(-1) * (sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('d') * ((Integer(2) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_114():
    f = 1/((a + b*x**4)**(sympy.S(7)/4)*(c + d*x**4)**2)
    F = ((Symbol('b') * ((Integer(4) * Symbol('b') * Symbol('c')) + (Integer(3) * Symbol('a') * Symbol('d'))) * x) * ((Integer(12) * Symbol('a') * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * x) * ((Integer(4) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * (x)**(Integer(4))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(32) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) * ((Integer(1) + (Symbol('a') * ((Symbol('b') * (x)**(Integer(4))))**(Integer(-1)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (x)**(Integer(3)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.acot(((sympy.sqrt(Symbol('b')) * (x)**(Integer(2))) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(12) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('d'))**(Integer(2)) * ((Integer(10) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((Integer(-1) * (sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3))))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * ((Integer(10) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('a') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4))))) * sympy.Function('EllipticPi')((sympy.sqrt(((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))), sympy.asin((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_115():
    f = 1/((x**4 + 1)**(sympy.S(1)/4)*(x**4 + 2))
    F = 2**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(2*(x**4 + 1)**(sympy.S(1)/4)))/4 + 2**(sympy.S(1)/4)*atanh(2**(sympy.S(3)/4)*x/(2*(x**4 + 1)**(sympy.S(1)/4)))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_116():
    f = 1/((a + b*x**4)**(sympy.S(1)/4)*(a - x**4*(a - b)))
    F = atan(a**(sympy.S(1)/4)*x/(a + b*x**4)**(sympy.S(1)/4))/(2*a**(sympy.S(5)/4)) + atanh(a**(sympy.S(1)/4)*x/(a + b*x**4)**(sympy.S(1)/4))/(2*a**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_117():
    f = (a + b*x**4)**p*(c + d*x**4)**q
    F = x*(a + b*x**4)**p*(c + d*x**4)**q*appellf1(sympy.S(1)/4, -p, -q, sympy.S(5)/4, -b*x**4/a, -d*x**4/c)/((1 + b*x**4/a)**p*(1 + d*x**4/c)**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_118():
    f = (a + b*x**4)**2*(c + d*x**4)**q
    F = b*x*(a + b*x**4)*(c + d*x**4)**(q + 1)/(d*(4*q + 9)) - b*x*(c + d*x**4)**(q + 1)*(-a*d*(4*q + 13) + 5*b*c)/(d**2*(4*q + 5)*(4*q + 9)) + x*(c + d*x**4)**q*(a**2*d**2*(16*q**2 + 56*q + 45) - 2*a*b*c*d*(4*q + 9) + 5*b**2*c**2)*hyper((sympy.S(1)/4, -q), (sympy.S(5)/4,), -d*x**4/c)/(d**2*(1 + d*x**4/c)**q*(4*q + 5)*(4*q + 9))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_119():
    f = (a + b*x**4)*(c + d*x**4)**q
    F = b*x*(c + d*x**4)**(q + 1)/(d*(4*q + 5)) - x*(c + d*x**4)**q*(-a*d*(4*q + 5) + b*c)*hyper((sympy.S(1)/4, -q), (sympy.S(5)/4,), -d*x**4/c)/(d*(1 + d*x**4/c)**q*(4*q + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_120():
    f = (c + d*x**4)**q/(a + b*x**4)
    F = x*(c + d*x**4)**q*appellf1(sympy.S(1)/4, 1, -q, sympy.S(5)/4, -b*x**4/a, -d*x**4/c)/(a*(1 + d*x**4/c)**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_121():
    f = (c + d*x**4)**q/(a + b*x**4)**2
    F = x*(c + d*x**4)**q*appellf1(sympy.S(1)/4, 2, -q, sympy.S(5)/4, -b*x**4/a, -d*x**4/c)/(a**2*(1 + d*x**4/c)**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_122():
    f = 1/((a + b*x**5)**(sympy.S(1)/5)*(c + d*x**5))
    F = (1 + sqrt(5))*log((2*c**(sympy.S(2)/5)*(a + b*x**5)**(sympy.S(2)/5) + c**(sympy.S(1)/5)*x*(a + b*x**5)**(sympy.S(1)/5)*(-a*d + b*c)**(sympy.S(1)/5) + sqrt(5)*c**(sympy.S(1)/5)*x*(a + b*x**5)**(sympy.S(1)/5)*(-a*d + b*c)**(sympy.S(1)/5) + 2*x**2*(-a*d + b*c)**(sympy.S(2)/5))/(a + b*x**5)**(sympy.S(2)/5))/(20*c**(sympy.S(4)/5)*(-a*d + b*c)**(sympy.S(1)/5)) + (1 - sqrt(5))*log((2*c**(sympy.S(2)/5)*(a + b*x**5)**(sympy.S(2)/5) - sqrt(5)*c**(sympy.S(1)/5)*x*(a + b*x**5)**(sympy.S(1)/5)*(-a*d + b*c)**(sympy.S(1)/5) + c**(sympy.S(1)/5)*x*(a + b*x**5)**(sympy.S(1)/5)*(-a*d + b*c)**(sympy.S(1)/5) + 2*x**2*(-a*d + b*c)**(sympy.S(2)/5))/(a + b*x**5)**(sympy.S(2)/5))/(20*c**(sympy.S(4)/5)*(-a*d + b*c)**(sympy.S(1)/5)) - log(c**(sympy.S(1)/5) - x*(-a*d + b*c)**(sympy.S(1)/5)/(a + b*x**5)**(sympy.S(1)/5))/(5*c**(sympy.S(4)/5)*(-a*d + b*c)**(sympy.S(1)/5)) - sqrt(sqrt(5)/2 + sympy.S(5)/2)*atan(sqrt(1 - 2*sqrt(5)/5) - 2*sqrt(2)*x*(-a*d + b*c)**(sympy.S(1)/5)/(c**(sympy.S(1)/5)*sqrt(sqrt(5) + 5)*(a + b*x**5)**(sympy.S(1)/5)))/(5*c**(sympy.S(4)/5)*(-a*d + b*c)**(sympy.S(1)/5)) + sqrt(sympy.S(5)/2 - sqrt(5)/2)*atan(sqrt(2*sqrt(5)/5 + 1) + x*sqrt(2*sqrt(5)/5 + 2)*(-a*d + b*c)**(sympy.S(1)/5)/(c**(sympy.S(1)/5)*(a + b*x**5)**(sympy.S(1)/5)))/(5*c**(sympy.S(4)/5)*(-a*d + b*c)**(sympy.S(1)/5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_123():
    f = sqrt(a + b/x)*(c + d/x)**3
    F = -7*d*sqrt(a + b/x)*(c + d/x)**2/5 + x*sqrt(a + b/x)*(c + d/x)**3 - d*sqrt(a + b/x)*(-4*a**2*d**2 + 30*a*b*c*d + 114*b**2*c**2 + b*d*(2*a*d + 33*b*c)/x)/(15*b**2) + c**2*(6*a*d + b*c)*atanh(sqrt(a + b/x)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_124():
    f = sqrt(a + b/x)*(c + d/x)**2
    F = -2*d**2*(a + b/x)**(sympy.S(3)/2)/(3*b) + c**2*x*(a + b/x)**(sympy.S(3)/2)/a - c*sqrt(a + b/x)*(4*a*d + b*c)/a + c*(4*a*d + b*c)*atanh(sqrt(a + b/x)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_125():
    f = sqrt(a + b/x)*(c + d/x)
    F = c*x*(a + b/x)**(sympy.S(3)/2)/a - sqrt(a + b/x)*(2*a*d + b*c)/a + (2*a*d + b*c)*atanh(sqrt(a + b/x)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_126():
    f = sqrt(a + b/x)
    F = x*sqrt(a + b/x) + b*atanh(sqrt(a + b/x)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_127():
    f = sqrt(a + b/x)/(c + d/x)
    F = x*sqrt(a + b/x)/c + 2*sqrt(d)*sqrt(-a*d + b*c)*atan(sqrt(d)*sqrt(a + b/x)/sqrt(-a*d + b*c))/c**2 + (-2*a*d + b*c)*atanh(sqrt(a + b/x)/sqrt(a))/(sqrt(a)*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_128():
    f = sqrt(a + b/x)/(c + d/x)**2
    F = x*sqrt(a + b/x)/(c*(c + d/x)) + 2*d*sqrt(a + b/x)/(c**2*(c + d/x)) + sqrt(d)*(-4*a*d + 3*b*c)*atan(sqrt(d)*sqrt(a + b/x)/sqrt(-a*d + b*c))/(c**3*sqrt(-a*d + b*c)) + (-4*a*d + b*c)*atanh(sqrt(a + b/x)/sqrt(a))/(sqrt(a)*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_129():
    f = sqrt(a + b/x)/(c + d/x)**3
    F = x*sqrt(a + b/x)/(c*(c + d/x)**2) + 3*d*sqrt(a + b/x)/(2*c**2*(c + d/x)**2) + d*sqrt(a + b/x)*(-12*a*d + 11*b*c)/(4*c**3*(c + d/x)*(-a*d + b*c)) + sqrt(d)*(24*a**2*d**2 - 40*a*b*c*d + 15*b**2*c**2)*atan(sqrt(d)*sqrt(a + b/x)/sqrt(-a*d + b*c))/(4*c**4*(-a*d + b*c)**(sympy.S(3)/2)) + (-6*a*d + b*c)*atanh(sqrt(a + b/x)/sqrt(a))/(sqrt(a)*c**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_130():
    f = (a + b/x)**(sympy.S(3)/2)*(c + d/x)**3
    F = 3*sqrt(a)*c**2*(2*a*d + b*c)*atanh(sqrt(a + b/x)/sqrt(a)) - 3*c**2*sqrt(a + b/x)*(2*a*d + b*c) - 9*d*(a + b/x)**(sympy.S(3)/2)*(c + d/x)**2/7 + x*(a + b/x)**(sympy.S(3)/2)*(c + d/x)**3 - d*(a + b/x)**(sympy.S(3)/2)*(3*b*d*(2*a*d + 19*b*c)/x + (-2*a*d + 26*b*c)*(2*a*d + 5*b*c))/(35*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_131():
    f = (a + b/x)**(sympy.S(3)/2)*(c + d/x)**2
    F = sqrt(a)*c*(4*a*d + 3*b*c)*atanh(sqrt(a + b/x)/sqrt(a)) - c*sqrt(a + b/x)*(4*a*d + 3*b*c) - 2*d**2*(a + b/x)**(sympy.S(5)/2)/(5*b) + c**2*x*(a + b/x)**(sympy.S(5)/2)/a - c*(a + b/x)**(sympy.S(3)/2)*(4*a*d + 3*b*c)/(3*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_132():
    f = (a + b/x)**(sympy.S(3)/2)*(c + d/x)
    F = sqrt(a)*(2*a*d + 3*b*c)*atanh(sqrt(a + b/x)/sqrt(a)) - sqrt(a + b/x)*(2*a*d + 3*b*c) + c*x*(a + b/x)**(sympy.S(5)/2)/a - (a + b/x)**(sympy.S(3)/2)*(2*a*d + 3*b*c)/(3*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_133():
    f = (a + b/x)**(sympy.S(3)/2)
    F = 3*sqrt(a)*b*atanh(sqrt(a + b/x)/sqrt(a)) - 3*b*sqrt(a + b/x) + x*(a + b/x)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_134():
    f = (a + b/x)**(sympy.S(3)/2)/(c + d/x)
    F = sqrt(a)*(-2*a*d + 3*b*c)*atanh(sqrt(a + b/x)/sqrt(a))/c**2 + a*x*sqrt(a + b/x)/c - 2*(-a*d + b*c)**(sympy.S(3)/2)*atan(sqrt(d)*sqrt(a + b/x)/sqrt(-a*d + b*c))/(c**2*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_135():
    f = (a + b/x)**(sympy.S(3)/2)/(c + d/x)**2
    F = sqrt(a)*(-4*a*d + 3*b*c)*atanh(sqrt(a + b/x)/sqrt(a))/c**3 + a*x*sqrt(a + b/x)/(c*(c + d/x)) - sqrt(a + b/x)*(-2*a*d + b*c)/(c**2*(c + d/x)) - (-4*a*d + b*c)*sqrt(-a*d + b*c)*atan(sqrt(d)*sqrt(a + b/x)/sqrt(-a*d + b*c))/(c**3*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_136():
    f = (a + b/x)**(sympy.S(3)/2)/(c + d/x)**3
    F = 3*sqrt(a)*(-2*a*d + b*c)*atanh(sqrt(a + b/x)/sqrt(a))/c**4 + a*x*sqrt(a + b/x)/(c*(c + d/x)**2) - sqrt(a + b/x)*(-3*a*d + b*c)/(2*c**2*(c + d/x)**2) - sqrt(a + b/x)*(-12*a*d + 3*b*c)/(4*c**3*(c + d/x)) - (24*a**2*d**2 - 24*a*b*c*d + 3*b**2*c**2)*atan(sqrt(d)*sqrt(a + b/x)/sqrt(-a*d + b*c))/(4*c**4*sqrt(d)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_137():
    f = (a + b/x)**(sympy.S(5)/2)*(c + d/x)**3
    F = a**(sympy.S(3)/2)*c**2*(6*a*d + 5*b*c)*atanh(sqrt(a + b/x)/sqrt(a)) - a*c**2*sqrt(a + b/x)*(6*a*d + 5*b*c) - c**2*(a + b/x)**(sympy.S(3)/2)*(6*a*d + 5*b*c)/3 - 11*d*(a + b/x)**(sympy.S(5)/2)*(c + d/x)**2/9 + x*(a + b/x)**(sympy.S(5)/2)*(c + d/x)**3 - d*(a + b/x)**(sympy.S(5)/2)*(-20*a**2*d**2 + 270*a*b*c*d + 938*b**2*c**2 + 5*b*d*(10*a*d + 89*b*c)/x)/(315*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_138():
    f = (a + b/x)**(sympy.S(5)/2)*(c + d/x)**2
    F = a**(sympy.S(3)/2)*c*(4*a*d + 5*b*c)*atanh(sqrt(a + b/x)/sqrt(a)) - a*c*sqrt(a + b/x)*(4*a*d + 5*b*c) - c*(a + b/x)**(sympy.S(3)/2)*(4*a*d + 5*b*c)/3 - 2*d**2*(a + b/x)**(sympy.S(7)/2)/(7*b) + c**2*x*(a + b/x)**(sympy.S(7)/2)/a - c*(a + b/x)**(sympy.S(5)/2)*(4*a*d + 5*b*c)/(5*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_139():
    f = (a + b/x)**(sympy.S(5)/2)*(c + d/x)
    F = a**(sympy.S(3)/2)*(2*a*d + 5*b*c)*atanh(sqrt(a + b/x)/sqrt(a)) - a*sqrt(a + b/x)*(2*a*d + 5*b*c) - (a + b/x)**(sympy.S(3)/2)*(2*a*d/3 + 5*b*c/3) + c*x*(a + b/x)**(sympy.S(7)/2)/a - (a + b/x)**(sympy.S(5)/2)*(2*a*d + 5*b*c)/(5*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_140():
    f = (a + b/x)**(sympy.S(5)/2)
    F = 5*a**(sympy.S(3)/2)*b*atanh(sqrt(a + b/x)/sqrt(a)) - 5*a*b*sqrt(a + b/x) - 5*b*(a + b/x)**(sympy.S(3)/2)/3 + x*(a + b/x)**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_141():
    f = (a + b/x)**(sympy.S(5)/2)/(c + d/x)
    F = a**(sympy.S(3)/2)*(-2*a*d + 5*b*c)*atanh(sqrt(a + b/x)/sqrt(a))/c**2 + a*x*(a + b/x)**(sympy.S(3)/2)/c - b*sqrt(a + b/x)*(a*d + 2*b*c)/(c*d) + 2*(-a*d + b*c)**(sympy.S(5)/2)*atan(sqrt(d)*sqrt(a + b/x)/sqrt(-a*d + b*c))/(c**2*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_142():
    f = (a + b/x)**(sympy.S(5)/2)/(c + d/x)**2
    F = a**(sympy.S(3)/2)*(-4*a*d + 5*b*c)*atanh(sqrt(a + b/x)/sqrt(a))/c**3 + a*x*(a + b/x)**(sympy.S(3)/2)/(c*(c + d/x)) + sqrt(a + b/x)*(-2*a*d + b*c)*(-a*d + b*c)/(c**2*d*(c + d/x)) - (-a*d + b*c)**(sympy.S(3)/2)*(4*a*d + b*c)*atan(sqrt(d)*sqrt(a + b/x)/sqrt(-a*d + b*c))/(c**3*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_143():
    f = (a + b/x)**(sympy.S(5)/2)/(c + d/x)**3
    F = a**(sympy.S(3)/2)*(-6*a*d + 5*b*c)*atanh(sqrt(a + b/x)/sqrt(a))/c**4 + a*x*(a + b/x)**(sympy.S(3)/2)/(c*(c + d/x)**2) + sqrt(a + b/x)*(-3*a*d + b*c)*(-a*d + b*c)/(2*c**2*d*(c + d/x)**2) - sqrt(a + b/x)*(-12*a**2*d**2 + 7*a*b*c*d + b**2*c**2)/(4*c**3*d*(c + d/x)) - sqrt(-a*d + b*c)*(-24*a**2*d**2 + 8*a*b*c*d + b**2*c**2)*atan(sqrt(d)*sqrt(a + b/x)/sqrt(-a*d + b*c))/(4*c**4*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_144():
    f = (c + d/x)**3/sqrt(a + b/x)
    F = c*x*sqrt(a + b/x)*(c + d/x)**2/a - d*sqrt(a + b/x)*(-4*a**2*d**2 + 18*a*b*c*d + 6*b**2*c**2 + b*d*(2*a*d + 3*b*c)/x)/(3*a*b**2) - c**2*(-6*a*d + b*c)*atanh(sqrt(a + b/x)/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_145():
    f = (c + d/x)**2/sqrt(a + b/x)
    F = -2*d**2*sqrt(a + b/x)/b + c**2*x*sqrt(a + b/x)/a - c*(-4*a*d + b*c)*atanh(sqrt(a + b/x)/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_146():
    f = (c + d/x)/sqrt(a + b/x)
    F = c*x*sqrt(a + b/x)/a - (-2*a*d + b*c)*atanh(sqrt(a + b/x)/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_147():
    f = 1/sqrt(a + b/x)
    F = x*sqrt(a + b/x)/a - b*atanh(sqrt(a + b/x)/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_148():
    f = 1/(sqrt(a + b/x)*(c + d/x))
    F = -2*d**(sympy.S(3)/2)*atan(sqrt(d)*sqrt(a + b/x)/sqrt(-a*d + b*c))/(c**2*sqrt(-a*d + b*c)) + x*sqrt(a + b/x)/(a*c) - (2*a*d + b*c)*atanh(sqrt(a + b/x)/sqrt(a))/(a**(sympy.S(3)/2)*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_149():
    f = 1/(sqrt(a + b/x)*(c + d/x)**2)
    F = -d**(sympy.S(3)/2)*(-4*a*d + 5*b*c)*atan(sqrt(d)*sqrt(a + b/x)/sqrt(-a*d + b*c))/(c**3*(-a*d + b*c)**(sympy.S(3)/2)) + x*sqrt(a + b/x)/(a*c*(c + d/x)) + d*sqrt(a + b/x)*(-2*a*d + b*c)/(a*c**2*(c + d/x)*(-a*d + b*c)) - (4*a*d + b*c)*atanh(sqrt(a + b/x)/sqrt(a))/(a**(sympy.S(3)/2)*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_150():
    f = 1/(sqrt(a + b/x)*(c + d/x)**3)
    F = -d**(sympy.S(3)/2)*(24*a**2*d**2 - 56*a*b*c*d + 35*b**2*c**2)*atan(sqrt(d)*sqrt(a + b/x)/sqrt(-a*d + b*c))/(4*c**4*(-a*d + b*c)**(sympy.S(5)/2)) + x*sqrt(a + b/x)/(a*c*(c + d/x)**2) + d*sqrt(a + b/x)*(-3*a*d + 2*b*c)/(2*a*c**2*(c + d/x)**2*(-a*d + b*c)) + d*sqrt(a + b/x)*(-4*a*d + b*c)*(-3*a*d + 4*b*c)/(4*a*c**3*(c + d/x)*(-a*d + b*c)**2) - (6*a*d + b*c)*atanh(sqrt(a + b/x)/sqrt(a))/(a**(sympy.S(3)/2)*c**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_151():
    f = (c + d/x)**3/(a + b/x)**(sympy.S(3)/2)
    F = c*x*(c + d/x)**2/(a*sqrt(a + b/x)) + (-a*b*d**2*(2*a*d + b*c)/x + (-2*a*d + b*c)*(2*a**2*d**2 - 2*a*b*c*d + 3*b**2*c**2))/(a**2*b**2*sqrt(a + b/x)) - 3*c**2*(-2*a*d + b*c)*atanh(sqrt(a + b/x)/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_152():
    f = (c + d/x)**2/(a + b/x)**(sympy.S(3)/2)
    F = c**2*x/(a*sqrt(a + b/x)) + (2*a**2*d**2 + b*c*(-4*a*d + 3*b*c))/(a**2*b*sqrt(a + b/x)) - c*(-4*a*d + 3*b*c)*atanh(sqrt(a + b/x)/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_153():
    f = (c + d/x)/(a + b/x)**(sympy.S(3)/2)
    F = c*x/(a*sqrt(a + b/x)) + (-2*a*d + 3*b*c)/(a**2*sqrt(a + b/x)) - (-2*a*d + 3*b*c)*atanh(sqrt(a + b/x)/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_154():
    f = (a + b/x)**(sympy.S(-3)/2)
    F = x/(a*sqrt(a + b/x)) + 3*b/(a**2*sqrt(a + b/x)) - 3*b*atanh(sqrt(a + b/x)/sqrt(a))/a**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_155():
    f = 1/((a + b/x)**(sympy.S(3)/2)*(c + d/x))
    F = 2*d**(sympy.S(5)/2)*atan(sqrt(d)*sqrt(a + b/x)/sqrt(-a*d + b*c))/(c**2*(-a*d + b*c)**(sympy.S(3)/2)) + x/(a*c*sqrt(a + b/x)) + b*(-a*d + 3*b*c)/(a**2*c*sqrt(a + b/x)*(-a*d + b*c)) - (2*a*d + 3*b*c)*atanh(sqrt(a + b/x)/sqrt(a))/(a**(sympy.S(5)/2)*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_156():
    f = 1/((a + b/x)**(sympy.S(3)/2)*(c + d/x)**2)
    F = d**(sympy.S(5)/2)*(-4*a*d + 7*b*c)*atan(sqrt(d)*sqrt(a + b/x)/sqrt(-a*d + b*c))/(c**3*(-a*d + b*c)**(sympy.S(5)/2)) + x/(a*c*sqrt(a + b/x)*(c + d/x)) + d*(-2*a*d + b*c)/(a*c**2*sqrt(a + b/x)*(c + d/x)*(-a*d + b*c)) + b*(2*a**2*d**2 - 2*a*b*c*d + 3*b**2*c**2)/(a**2*c**2*sqrt(a + b/x)*(-a*d + b*c)**2) - (4*a*d + 3*b*c)*atanh(sqrt(a + b/x)/sqrt(a))/(a**(sympy.S(5)/2)*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_157():
    f = 1/((a + b/x)**(sympy.S(3)/2)*(c + d/x)**3)
    F = 3*d**(sympy.S(5)/2)*(8*a**2*d**2 - 24*a*b*c*d + 21*b**2*c**2)*atan(sqrt(d)*sqrt(a + b/x)/sqrt(-a*d + b*c))/(4*c**4*(-a*d + b*c)**(sympy.S(7)/2)) + x/(a*c*sqrt(a + b/x)*(c + d/x)**2) + d*(-3*a*d + 2*b*c)/(2*a*c**2*sqrt(a + b/x)*(c + d/x)**2*(-a*d + b*c)) + d*(12*a**2*d**2 - 21*a*b*c*d + 4*b**2*c**2)/(4*a*c**3*sqrt(a + b/x)*(c + d/x)*(-a*d + b*c)**2) + 3*b*(-a*d + 2*b*c)*(4*a**2*d**2 - a*b*c*d + 2*b**2*c**2)/(4*a**2*c**3*sqrt(a + b/x)*(-a*d + b*c)**3) - (6*a*d + 3*b*c)*atanh(sqrt(a + b/x)/sqrt(a))/(a**(sympy.S(5)/2)*c**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_158():
    f = (c + d/x)**3/(a + b/x)**(sympy.S(5)/2)
    F = c*x*(c + d/x)**2/(a*(a + b/x)**(sympy.S(3)/2)) + (-a*d + b*c)*(-4*a**3*d**2*x - 2*a**2*b*d*(5*c*x + 3*d) + a*b**2*c*(20*c*x - 3*d) + 15*b**3*c**2)/(3*a**3*b**2*x*(a + b/x)**(sympy.S(3)/2)) - c**2*(-6*a*d + 5*b*c)*atanh(sqrt(a + b/x)/sqrt(a))/a**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_159():
    f = (c + d/x)**2/(a + b/x)**(sympy.S(5)/2)
    F = c**2*x/(a*(a + b/x)**(sympy.S(3)/2)) + (2*a**2*d**2 + b*c*(-4*a*d + 5*b*c))/(3*a**2*b*(a + b/x)**(sympy.S(3)/2)) + c*(-4*a*d + 5*b*c)/(a**3*sqrt(a + b/x)) - c*(-4*a*d + 5*b*c)*atanh(sqrt(a + b/x)/sqrt(a))/a**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_160():
    f = (c + d/x)/(a + b/x)**(sympy.S(5)/2)
    F = c*x/(a*(a + b/x)**(sympy.S(3)/2)) + (-2*a*d + 5*b*c)/(3*a**2*(a + b/x)**(sympy.S(3)/2)) + (-2*a*d + 5*b*c)/(a**3*sqrt(a + b/x)) - (-2*a*d + 5*b*c)*atanh(sqrt(a + b/x)/sqrt(a))/a**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_161():
    f = (a + b/x)**(sympy.S(-5)/2)
    F = x/(a*(a + b/x)**(sympy.S(3)/2)) + 5*b/(3*a**2*(a + b/x)**(sympy.S(3)/2)) + 5*b/(a**3*sqrt(a + b/x)) - 5*b*atanh(sqrt(a + b/x)/sqrt(a))/a**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_162():
    f = 1/((a + b/x)**(sympy.S(5)/2)*(c + d/x))
    F = -2*d**(sympy.S(7)/2)*atan(sqrt(d)*sqrt(a + b/x)/sqrt(-a*d + b*c))/(c**2*(-a*d + b*c)**(sympy.S(5)/2)) + x/(a*c*(a + b/x)**(sympy.S(3)/2)) + b*(-3*a*d + 5*b*c)/(3*a**2*c*(a + b/x)**(sympy.S(3)/2)*(-a*d + b*c)) + b*(a**2*d**2 - 8*a*b*c*d + 5*b**2*c**2)/(a**3*c*sqrt(a + b/x)*(-a*d + b*c)**2) - (2*a*d + 5*b*c)*atanh(sqrt(a + b/x)/sqrt(a))/(a**(sympy.S(7)/2)*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_163():
    f = 1/((a + b/x)**(sympy.S(5)/2)*(c + d/x)**2)
    F = -d**(sympy.S(7)/2)*(-4*a*d + 9*b*c)*atan(sqrt(d)*sqrt(a + b/x)/sqrt(-a*d + b*c))/(c**3*(-a*d + b*c)**(sympy.S(7)/2)) + x/(a*c*(a + b/x)**(sympy.S(3)/2)*(c + d/x)) + d*(-2*a*d + b*c)/(a*c**2*(a + b/x)**(sympy.S(3)/2)*(c + d/x)*(-a*d + b*c)) + b*(6*a**2*d**2 - 6*a*b*c*d + 5*b**2*c**2)/(3*a**2*c**2*(a + b/x)**(sympy.S(3)/2)*(-a*d + b*c)**2) + b*(-2*a*d + b*c)*(a**2*d**2 - a*b*c*d + 5*b**2*c**2)/(a**3*c**2*sqrt(a + b/x)*(-a*d + b*c)**3) - (4*a*d + 5*b*c)*atanh(sqrt(a + b/x)/sqrt(a))/(a**(sympy.S(7)/2)*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_164():
    f = 1/((a + b/x)**(sympy.S(5)/2)*(c + d/x)**3)
    F = -d**(sympy.S(7)/2)*(24*a**2*d**2 - 88*a*b*c*d + 99*b**2*c**2)*atan(sqrt(d)*sqrt(a + b/x)/sqrt(-a*d + b*c))/(4*c**4*(-a*d + b*c)**(sympy.S(9)/2)) + x/(a*c*(a + b/x)**(sympy.S(3)/2)*(c + d/x)**2) + d*(-3*a*d + 2*b*c)/(2*a*c**2*(a + b/x)**(sympy.S(3)/2)*(c + d/x)**2*(-a*d + b*c)) + d*(12*a**2*d**2 - 23*a*b*c*d + 4*b**2*c**2)/(4*a*c**3*(a + b/x)**(sympy.S(3)/2)*(c + d/x)*(-a*d + b*c)**2) + b*(-36*a**3*d**3 + 87*a**2*b*c*d**2 - 36*a*b**2*c**2*d + 20*b**3*c**3)/(12*a**2*c**3*(a + b/x)**(sympy.S(3)/2)*(-a*d + b*c)**3) + b*(12*a**4*d**4 - 35*a**3*b*c*d**3 + 24*a**2*b**2*c**2*d**2 - 56*a*b**3*c**3*d + 20*b**4*c**4)/(4*a**3*c**3*sqrt(a + b/x)*(-a*d + b*c)**4) - (6*a*d + 5*b*c)*atanh(sqrt(a + b/x)/sqrt(a))/(a**(sympy.S(7)/2)*c**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_165():
    f = sqrt(a + b/x)*sqrt(c + d/x)
    F = -2*sqrt(b)*sqrt(d)*atanh(sqrt(d)*sqrt(a + b/x)/(sqrt(b)*sqrt(c + d/x))) + x*sqrt(a + b/x)*sqrt(c + d/x) + (a*d + b*c)*atanh(sqrt(c)*sqrt(a + b/x)/(sqrt(a)*sqrt(c + d/x)))/(sqrt(a)*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_166():
    f = sqrt(a + b/x)/sqrt(c + d/x)
    F = x*sqrt(a + b/x)*sqrt(c + d/x)/c + (-a*d + b*c)*atanh(sqrt(c)*sqrt(a + b/x)/(sqrt(a)*sqrt(c + d/x)))/(sqrt(a)*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_167():
    f = sqrt(a + b/x)/(c + d/x)**(sympy.S(3)/2)
    F = x*(a + b/x)**(sympy.S(3)/2)/(a*c*sqrt(c + d/x)) - sqrt(a + b/x)*(-3*a*d + b*c)/(a*c**2*sqrt(c + d/x)) + (-3*a*d + b*c)*atanh(sqrt(c)*sqrt(a + b/x)/(sqrt(a)*sqrt(c + d/x)))/(sqrt(a)*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_168():
    f = (a + b/x)**p*(c + d/x)**q
    F = -b*(a + b/x)**(p + 1)*(c + d/x)**q*appellf1(p + 1, 2, -q, p + 2, (a + b/x)/a, -d*(a + b/x)/(-a*d + b*c))/(a**2*(b*(c + d/x)/(-a*d + b*c))**q*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_169():
    f = (a + b/x**2)/(c + d/x**2)
    F = a*x/c + (-a*d + b*c)*atan(sqrt(c)*x/sqrt(d))/(c**(sympy.S(3)/2)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_170():
    f = sqrt(a + b/x**2)*sqrt(c + d/x**2)
    F = 2*sqrt(c)*sqrt(d)*sqrt(a + b/x**2)*elliptic_e(acot(sqrt(c)*x/sqrt(d)), 1 - b*c/(a*d))/(sqrt(c*(a + b/x**2)/(a*(c + d/x**2)))*sqrt(c + d/x**2)) - 2*d*sqrt(a + b/x**2)/(x*sqrt(c + d/x**2)) + x*sqrt(a + b/x**2)*sqrt(c + d/x**2) - sqrt(c)*sqrt(a + b/x**2)*(a*d + b*c)*elliptic_f(acot(sqrt(c)*x/sqrt(d)), 1 - b*c/(a*d))/(a*sqrt(d)*sqrt(c*(a + b/x**2)/(a*(c + d/x**2)))*sqrt(c + d/x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_171():
    f = sqrt(a + b/x**2)/sqrt(c + d/x**2)
    F = -d*sqrt(a + b/x**2)/(c*x*sqrt(c + d/x**2)) + x*sqrt(a + b/x**2)*sqrt(c + d/x**2)/c + sqrt(d)*sqrt(a + b/x**2)*elliptic_e(acot(sqrt(c)*x/sqrt(d)), 1 - b*c/(a*d))/(sqrt(c)*sqrt(c*(a + b/x**2)/(a*(c + d/x**2)))*sqrt(c + d/x**2)) - b*sqrt(c)*sqrt(a + b/x**2)*elliptic_f(acot(sqrt(c)*x/sqrt(d)), 1 - b*c/(a*d))/(a*sqrt(d)*sqrt(c*(a + b/x**2)/(a*(c + d/x**2)))*sqrt(c + d/x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_172():
    f = sqrt(a + b/x**2)/(c + d/x**2)**(sympy.S(3)/2)
    F = -x*sqrt(a + b/x**2)/(c*sqrt(c + d/x**2)) - 2*d*sqrt(a + b/x**2)/(c**2*x*sqrt(c + d/x**2)) + 2*x*sqrt(a + b/x**2)*sqrt(c + d/x**2)/c**2 + 2*sqrt(d)*sqrt(a + b/x**2)*elliptic_e(acot(sqrt(c)*x/sqrt(d)), 1 - b*c/(a*d))/(c**(sympy.S(3)/2)*sqrt(c*(a + b/x**2)/(a*(c + d/x**2)))*sqrt(c + d/x**2)) - b*sqrt(a + b/x**2)*elliptic_f(acot(sqrt(c)*x/sqrt(d)), 1 - b*c/(a*d))/(a*sqrt(c)*sqrt(d)*sqrt(c*(a + b/x**2)/(a*(c + d/x**2)))*sqrt(c + d/x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_173():
    f = (a + b/x**2)**p*(c + d/x**2)**q
    F = x*(a + b/x**2)**p*(c + d/x**2)**q*appellf1(sympy.S(-1)/2, -p, -q, sympy.S.Half, -b/(a*x**2), -d/(c*x**2))/((1 + b/(a*x**2))**p*(1 + d/(c*x**2))**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_174():
    f = (a + b/x**3)/(c + d/x**3)
    F = a*x/c + (-a*d + b*c)*log(c**(sympy.S(1)/3)*x + d**(sympy.S(1)/3))/(3*c**(sympy.S(4)/3)*d**(sympy.S(2)/3)) - (-a*d + b*c)*log(c**(sympy.S(2)/3)*x**2 - c**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3))/(6*c**(sympy.S(4)/3)*d**(sympy.S(2)/3)) - sqrt(3)*(-a*d + b*c)*atan(sqrt(3)*(-2*c**(sympy.S(1)/3)*x + d**(sympy.S(1)/3))/(3*d**(sympy.S(1)/3)))/(3*c**(sympy.S(4)/3)*d**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_175():
    f = (a + b*sqrt(x))/(c + d*sqrt(x))
    F = b*x/d + 2*c*(-a*d + b*c)*log(c + d*sqrt(x))/d**3 - sqrt(x)*(-2*a*d + 2*b*c)/d**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_176():
    f = (x**(sympy.S(1)/3) - 1)/(x**(sympy.S(1)/3) + 1)
    F = -3*x**(sympy.S(2)/3) + 6*x**(sympy.S(1)/3) + x - 6*log(x**(sympy.S(1)/3) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_177():
    f = (x**(sympy.S(2)/3) + 1)/(x**(sympy.S(2)/3) - 1)
    F = 6*x**(sympy.S(1)/3) + x - 6*atanh(x**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_178():
    f = (x**(sympy.S(3)/4) - 16)/(x**(sympy.S(3)/4) + 16)
    F = -128*x**(sympy.S(1)/4) + x + 256*2**(sympy.S(1)/3)*log(x**(sympy.S(1)/4) + 2*2**(sympy.S(1)/3))/3 - 128*2**(sympy.S(1)/3)*log(-2*2**(sympy.S(1)/3)*x**(sympy.S(1)/4) + sqrt(x) + 4*2**(sympy.S(2)/3))/3 - 256*2**(sympy.S(1)/3)*sqrt(3)*atan(2**(sympy.S(2)/3)*sqrt(3)*(-x**(sympy.S(1)/4) + 2**(sympy.S(1)/3))/6)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_179():
    f = (1 + x**(sympy.S(-1)/3))/(-1 + x**(sympy.S(-1)/3))
    F = -3*x**(sympy.S(2)/3) - 6*x**(sympy.S(1)/3) - x - 6*log(1 - x**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_180():
    f = (a - b*x**n)**(sympy.S(3)/2)*(a + b*x**n)**(sympy.S(3)/2)
    F = a**2*x*sqrt(a - b*x**n)*sqrt(a + b*x**n)*hyper((sympy.S(-3)/2, 1/(2*n)), (1 + 1/(2*n),), b**2*x**(2*n)/a**2)/sqrt(1 - b**2*x**(2*n)/a**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_181():
    f = sqrt(a - b*x**n)*sqrt(a + b*x**n)
    F = x*sqrt(a - b*x**n)*sqrt(a + b*x**n)*hyper((sympy.S(-1)/2, 1/(2*n)), (1 + 1/(2*n),), b**2*x**(2*n)/a**2)/sqrt(1 - b**2*x**(2*n)/a**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_182():
    f = (a - b*x**n)**p*(a + b*x**n)**p
    F = x*(a - b*x**n)**p*(a + b*x**n)**p*hyper((-p, 1/(2*n)), (1 + 1/(2*n),), b**2*x**(2*n)/a**2)/(1 - b**2*x**(2*n)/a**2)**p
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_183():
    f = (a + b*x**n)*(c + d*x**n)**4
    F = a*c**4*x + b*d**4*x**(5*n + 1)/(5*n + 1) + c**3*x**(n + 1)*(4*a*d + b*c)/(n + 1) + 2*c**2*d*x**(2*n + 1)*(3*a*d + 2*b*c)/(2*n + 1) + 2*c*d**2*x**(3*n + 1)*(2*a*d + 3*b*c)/(3*n + 1) + d**3*x**(4*n + 1)*(a*d + 4*b*c)/(4*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_184():
    f = (a + b*x**n)*(c + d*x**n)**3
    F = a*c**3*x + b*d**3*x**(4*n + 1)/(4*n + 1) + c**2*x**(n + 1)*(3*a*d + b*c)/(n + 1) + 3*c*d*x**(2*n + 1)*(a*d + b*c)/(2*n + 1) + d**2*x**(3*n + 1)*(a*d + 3*b*c)/(3*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_185():
    f = (a + b*x**n)*(c + d*x**n)**2
    F = a*c**2*x + b*d**2*x**(3*n + 1)/(3*n + 1) + c*x**(n + 1)*(2*a*d + b*c)/(n + 1) + d*x**(2*n + 1)*(a*d + 2*b*c)/(2*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_186():
    f = (a + b*x**n)*(c + d*x**n)
    F = a*c*x + b*d*x**(2*n + 1)/(2*n + 1) + x**(n + 1)*(a*d + b*c)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_187():
    f = (a + b*x**n)/(c + d*x**n)
    F = b*x/d - x*(-a*d + b*c)*hyper((1, 1/n), (1 + 1/n,), -d*x**n/c)/(c*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_188():
    f = (a + b*x**n)/(c + d*x**n)**2
    F = -x*(-a*d + b*c)/(c*d*n*(c + d*x**n)) + x*(-a*d*(1 - n) + b*c)*hyper((1, 1/n), (1 + 1/n,), -d*x**n/c)/(c**2*d*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_189():
    f = (a + b*x**n)/(c + d*x**n)**3
    F = -x*(-a*d + b*c)/(2*c*d*n*(c + d*x**n)**2) + x*(-a*d*(1 - 2*n) + b*c)*hyper((2, 1/n), (1 + 1/n,), -d*x**n/c)/(2*c**3*d*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_190():
    f = (a + b*x**n)/(c + d*x**n)**4
    F = -x*(-a*d + b*c)/(3*c*d*n*(c + d*x**n)**3) + x*(-a*d*(1 - 3*n) + b*c)*hyper((3, 1/n), (1 + 1/n,), -d*x**n/c)/(3*c**4*d*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_191():
    f = (a + b*x**n)**2*(d + e*x**n)**3
    F = a**2*d**3*x + a*d**2*x**(n + 1)*(3*a*e + 2*b*d)/(n + 1) + b**2*e**3*x**(5*n + 1)/(5*n + 1) + b*e**2*x**(4*n + 1)*(2*a*e + 3*b*d)/(4*n + 1) + d*x**(2*n + 1)*(3*a**2*e**2 + 6*a*b*d*e + b**2*d**2)/(2*n + 1) + e*x**(3*n + 1)*(a**2*e**2 + 6*a*b*d*e + 3*b**2*d**2)/(3*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_192():
    f = (a + b*x**n)**2*(d + e*x**n)**2
    F = a**2*d**2*x + 2*a*d*x**(n + 1)*(a*e + b*d)/(n + 1) + b**2*e**2*x**(4*n + 1)/(4*n + 1) + 2*b*e*x**(3*n + 1)*(a*e + b*d)/(3*n + 1) + x**(2*n + 1)*(a**2*e**2 + 4*a*b*d*e + b**2*d**2)/(2*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_193():
    f = (a + b*x**n)**2*(c + d*x**n)
    F = a**2*c*x + a*x**(n + 1)*(a*d + 2*b*c)/(n + 1) + b**2*d*x**(3*n + 1)/(3*n + 1) + b*x**(2*n + 1)*(2*a*d + b*c)/(2*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_194():
    f = (a + b*x**n)**2/(c + d*x**n)
    F = b*x*(a + b*x**n)/(d*(n + 1)) - b*x*(-a*d*(2*n + 1) + b*c*(n + 1))/(d**2*(n + 1)) + x*(-a*d + b*c)**2*hyper((1, 1/n), (1 + 1/n,), -d*x**n/c)/(c*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_195():
    f = (a + b*x**n)**2/(c + d*x**n)**2
    F = -b*x*(a*d - b*c*(n + 1))/(c*d**2*n) - x*(a + b*x**n)*(-a*d + b*c)/(c*d*n*(c + d*x**n)) + x*(-a*d + b*c)*(a*d*(1 - n) - b*c*(n + 1))*hyper((1, 1/n), (1 + 1/n,), -d*x**n/c)/(c**2*d**2*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_196():
    f = (a + b*x**n)**2/(c + d*x**n)**3
    F = -x*(a + b*x**n)*(-a*d + b*c)/(2*c*d*n*(c + d*x**n)**2) + x*(-a*d + b*c)*(a*d*(1 - 2*n) - b*c*(n + 1))/(2*c**2*d**2*n**2*(c + d*x**n)) - x*(-a**2*d**2*(2*n**2 - 3*n + 1) + 2*a*b*c*d*(1 - n) - b**2*c**2*(n + 1))*hyper((1, 1/n), (1 + 1/n,), -d*x**n/c)/(2*c**3*d**2*n**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_197():
    f = (c + d*x**n)**4/(a + b*x**n)
    F = d*x*(c + d*x**n)**3/(b*(3*n + 1)) - d*x*(c + d*x**n)**2*(a*d*(3*n + 1) - b*(6*c*n + c))/(b**2*(6*n**2 + 5*n + 1)) - d*x*(c + d*x**n)*(-a**2*d**2*(6*n**2 + 5*n + 1) + 2*a*b*c*d*(3*n + 1)**2 - b**2*c**2*(18*n**2 + 7*n + 1))/(b**3*(n + 1)*(2*n + 1)*(3*n + 1)) - d*x*(a**3*d**3*(6*n**3 + 11*n**2 + 6*n + 1) - a**2*b*c*d**2*(24*n**3 + 38*n**2 + 19*n + 3) + a*b**2*c**2*d*(36*n**3 + 45*n**2 + 20*n + 3) - b**3*c**3*(24*n**3 + 18*n**2 + 7*n + 1))/(b**4*(n + 1)*(2*n + 1)*(3*n + 1)) + x*(-a*d + b*c)**4*hyper((1, 1/n), (1 + 1/n,), -b*x**n/a)/(a*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_198():
    f = (c + d*x**n)**3/(a + b*x**n)
    F = d*x*(c + d*x**n)**2/(b*(2*n + 1)) - d*x*(c + d*x**n)*(a*d*(2*n + 1) - b*(4*c*n + c))/(b**2*(n + 1)*(2*n + 1)) + d*x*(a**2*d**2*(2*n**2 + 3*n + 1) - a*b*c*d*(6*n**2 + 7*n + 2) + b**2*c**2*(6*n**2 + 4*n + 1))/(b**3*(n + 1)*(2*n + 1)) + x*(-a*d + b*c)**3*hyper((1, 1/n), (1 + 1/n,), -b*x**n/a)/(a*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_199():
    f = (c + d*x**n)**2/(a + b*x**n)
    F = d*x*(c + d*x**n)/(b*(n + 1)) - d*x*(a*d*(n + 1) - b*(2*c*n + c))/(b**2*(n + 1)) + x*(-a*d + b*c)**2*hyper((1, 1/n), (1 + 1/n,), -b*x**n/a)/(a*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_200():
    f = (c + d*x**n)/(a + b*x**n)
    F = d*x/b + x*(-a*d + b*c)*hyper((1, 1/n), (1 + 1/n,), -b*x**n/a)/(a*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_201():
    f = 1/((a + b*x**n)*(c + d*x**n))
    F = -d*x*hyper((1, 1/n), (1 + 1/n,), -d*x**n/c)/(c*(-a*d + b*c)) + b*x*hyper((1, 1/n), (1 + 1/n,), -b*x**n/a)/(a*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_202():
    f = 1/((a + b*x**n)*(c + d*x**n)**2)
    F = -d*x/(c*n*(c + d*x**n)*(-a*d + b*c)) + d*x*(-a*d*(1 - n) + b*c*(1 - 2*n))*hyper((1, 1/n), (1 + 1/n,), -d*x**n/c)/(c**2*n*(-a*d + b*c)**2) + b**2*x*hyper((1, 1/n), (1 + 1/n,), -b*x**n/a)/(a*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_203():
    f = 1/((a + b*x**n)*(c + d*x**n)**3)
    F = -d*x/(2*c*n*(c + d*x**n)**2*(-a*d + b*c)) - d*x*(a*d*(1 - 2*n) - b*(-4*c*n + c))/(2*c**2*n**2*(c + d*x**n)*(-a*d + b*c)**2) - d*x*(a**2*d**2*(2*n**2 - 3*n + 1) - 2*a*b*c*d*(3*n**2 - 4*n + 1) + b**2*c**2*(6*n**2 - 5*n + 1))*hyper((1, 1/n), (1 + 1/n,), -d*x**n/c)/(2*c**3*n**2*(-a*d + b*c)**3) + b**3*x*hyper((1, 1/n), (1 + 1/n,), -b*x**n/a)/(a*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_204():
    f = (c + d*x**n)**4/(a + b*x**n)**2
    F = x*(c + d*x**n)**3*(-a*d + b*c)/(a*b*n*(a + b*x**n)) + d*x*(c + d*x**n)**2*(a*d*(3*n + 1) - b*(2*c*n + c))/(a*b**2*n*(2*n + 1)) - d*x*(c + d*x**n)*(a**2*d**2*(6*n**2 + 5*n + 1) - 2*a*b*c*d*(5*n**2 + 4*n + 1) + b**2*c**2*(2*n**2 + 3*n + 1))/(a*b**3*n*(n + 1)*(2*n + 1)) - d*x*(-a**3*d**3*(6*n**3 + 11*n**2 + 6*n + 1) + a**2*b*c*d**2*(16*n**3 + 26*n**2 + 15*n + 3) - a*b**2*c**2*d*(12*n**3 + 17*n**2 + 12*n + 3) + b**3*c**3*(2*n**2 + 3*n + 1))/(a*b**4*n*(n + 1)*(2*n + 1)) - x*(-a*d + b*c)**3*(-a*d*(3*n + 1) + b*c*(1 - n))*hyper((1, 1/n), (1 + 1/n,), -b*x**n/a)/(a**2*b**4*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_205():
    f = (c + d*x**n)**3/(a + b*x**n)**2
    F = x*(c + d*x**n)**2*(-a*d + b*c)/(a*b*n*(a + b*x**n)) - d*x*(c + d*x**n)*(-a*d*(2*n + 1) + b*c*(n + 1))/(a*b**2*n*(n + 1)) - d*x*(a**2*d**2*(2*n**2 + 3*n + 1) - a*b*c*d*(3*n**2 + 4*n + 2) + b**2*c**2*(n + 1))/(a*b**3*n*(n + 1)) - x*(-a*d + b*c)**2*(-a*d*(2*n + 1) + b*c*(1 - n))*hyper((1, 1/n), (1 + 1/n,), -b*x**n/a)/(a**2*b**3*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_206():
    f = (c + d*x**n)**2/(a + b*x**n)**2
    F = x*(c + d*x**n)*(-a*d + b*c)/(a*b*n*(a + b*x**n)) - d*x*(-a*d*(n + 1) + b*c)/(a*b**2*n) - x*(-a*d + b*c)*(-a*d*(n + 1) + b*c*(1 - n))*hyper((1, 1/n), (1 + 1/n,), -b*x**n/a)/(a**2*b**2*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_207():
    f = (c + d*x**n)/(a + b*x**n)**2
    F = x*(-a*d + b*c)/(a*b*n*(a + b*x**n)) + x*(a*d - b*c*(1 - n))*hyper((1, 1/n), (1 + 1/n,), -b*x**n/a)/(a**2*b*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_208():
    f = 1/((a + b*x**n)**2*(c + d*x**n))
    F = d**2*x*hyper((1, 1/n), (1 + 1/n,), -d*x**n/c)/(c*(-a*d + b*c)**2) + b*x/(a*n*(a + b*x**n)*(-a*d + b*c)) + b*x*(a*d*(1 - 2*n) - b*c*(1 - n))*hyper((1, 1/n), (1 + 1/n,), -b*x**n/a)/(a**2*n*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_209():
    f = 1/((a + b*x**n)**2*(c + d*x**n)**2)
    F = -d**2*x*(-a*d*(1 - n) + b*c*(1 - 3*n))*hyper((1, 1/n), (1 + 1/n,), -d*x**n/c)/(c**2*n*(-a*d + b*c)**3) + b*x/(a*n*(a + b*x**n)*(c + d*x**n)*(-a*d + b*c)) + d*x*(a*d + b*c)/(a*c*n*(c + d*x**n)*(-a*d + b*c)**2) + b**2*x*(a*d*(1 - 3*n) - b*(-c*n + c))*hyper((1, 1/n), (1 + 1/n,), -b*x**n/a)/(a**2*n*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_210():
    f = 1/((a + b*x**n)**2*(c + d*x**n)**3)
    F = d**2*x*(a**2*d**2*(2*n**2 - 3*n + 1) - 2*a*b*c*d*(4*n**2 - 5*n + 1) + b**2*c**2*(12*n**2 - 7*n + 1))*hyper((1, 1/n), (1 + 1/n,), -d*x**n/c)/(2*c**3*n**2*(-a*d + b*c)**4) + b*x/(a*n*(a + b*x**n)*(c + d*x**n)**2*(-a*d + b*c)) + d*x*(a*d + 2*b*c)/(2*a*c*n*(c + d*x**n)**2*(-a*d + b*c)**2) - d*x*(-a**2*d**2*(1 - 2*n) + a*b*c*d*(1 - 6*n) - 2*b**2*c**2*n)/(2*a*c**2*n**2*(c + d*x**n)*(-a*d + b*c)**3) + b**3*x*(a*d*(1 - 4*n) - b*c*(1 - n))*hyper((1, 1/n), (1 + 1/n,), -b*x**n/a)/(a**2*n*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_211():
    f = (a + b*x**n)**p*(c + d*x**n)**q
    F = x*(a + b*x**n)**p*(c + d*x**n)**q*appellf1(1/n, -p, -q, 1 + 1/n, -b*x**n/a, -d*x**n/c)/((1 + b*x**n/a)**p*(1 + d*x**n/c)**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_212():
    f = (a + b*x**n)**p*(c + d*x**n)**3
    F = d*x*(a + b*x**n)**(p + 1)*(c + d*x**n)**2/(b*(n*p + 3*n + 1)) - d*x*(a + b*x**n)**(p + 1)*(c + d*x**n)*(a*d*(2*n + 1) - b*c*(n*(p + 5) + 1))/(b**2*(n*(p + 2) + 1)*(n*(p + 3) + 1)) + d*x*(a + b*x**n)**(p + 1)*(a**2*d**2*(2*n**2 + 3*n + 1) - a*b*c*d*(n**2*(p + 7) + n*(2*p + 9) + 2) + b**2*c**2*(n**2*(p**2 + 6*p + 11) + 2*n*(p + 3) + 1))/(b**3*(n*(p + 2) + 1)*(n*(p + 3) + 1)*(n*p + n + 1)) - x*(a + b*x**n)**p*(a**3*d**3*(2*n**2 + 3*n + 1) - 3*a**2*b*c*d**2*(n + 1)*(n*(p + 3) + 1) + 3*a*b**2*c**2*d*(n**2*(p**2 + 5*p + 6) + n*(2*p + 5) + 1) - b**3*c**3*(n**3*(p**3 + 6*p**2 + 11*p + 6) + n**2*(3*p**2 + 12*p + 11) + 3*n*(p + 2) + 1))*hyper((1/n, -p), (1 + 1/n,), -b*x**n/a)/(b**3*(1 + b*x**n/a)**p*(n*(p + 2) + 1)*(n*(p + 3) + 1)*(n*p + n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_213():
    f = (a + b*x**n)**p*(c + d*x**n)**2
    F = d*x*(a + b*x**n)**(p + 1)*(c + d*x**n)/(b*(n*p + 2*n + 1)) - d*x*(a + b*x**n)**(p + 1)*(a*d*(n + 1) - b*c*(n*(p + 3) + 1))/(b**2*(n*(p + 2) + 1)*(n*p + n + 1)) - x*(a + b*x**n)**p*(-a*d*(a*d*(n + 1) - b*c*(n*(p + 3) + 1)) + b*c*(a*d - b*c*(n*(p + 2) + 1))*(n*p + n + 1))*hyper((1/n, -p), (1 + 1/n,), -b*x**n/a)/(b**2*(1 + b*x**n/a)**p*(n*(p + 2) + 1)*(n*p + n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_214():
    f = (a + b*x**n)**p*(c + d*x**n)
    F = d*x*(a + b*x**n)**(p + 1)/(b*(n*p + n + 1)) - x*(a + b*x**n)**p*(a*d - b*c*(n*p + n + 1))*hyper((1/n, -p), (1 + 1/n,), -b*x**n/a)/(b*(1 + b*x**n/a)**p*(n*p + n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_215():
    f = (a + b*x**n)**p
    F = x*(a + b*x**n)**p*hyper((1/n, -p), (1 + 1/n,), -b*x**n/a)/(1 + b*x**n/a)**p
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_216():
    f = (a + b*x**n)**p/(c + d*x**n)
    F = x*(a + b*x**n)**p*appellf1(1/n, 1, -p, 1 + 1/n, -d*x**n/c, -b*x**n/a)/(c*(1 + b*x**n/a)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_217():
    f = (a + b*x**n)**p/(c + d*x**n)**2
    F = x*(a + b*x**n)**p*appellf1(1/n, 2, -p, 1 + 1/n, -d*x**n/c, -b*x**n/a)/(c**2*(1 + b*x**n/a)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_218():
    f = (a + b*x**n)**p/(c + d*x**n)**3
    F = x*(a + b*x**n)**p*appellf1(1/n, 3, -p, 1 + 1/n, -d*x**n/c, -b*x**n/a)/(c**3*(1 + b*x**n/a)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_219():
    f = (a + b*x**n)**p*(c + d*x**n)**(-p - 1 - 1/n)
    F = x*(a + b*x**n)**p*(c + d*x**n)**(-p - 1/n)*hyper((1/n, -p), (1 + 1/n,), -x**n*(-a*d + b*c)/(a*(c + d*x**n)))/(c*(c*(a + b*x**n)/(a*(c + d*x**n)))**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_220():
    f = (a + b*x**n)**3*(c + d*x**n)**(-4 - 1/n)
    F = 6*a**3*n**3*x/(c**4*(c + d*x**n)**(1/n)*(n + 1)*(2*n + 1)*(3*n + 1)) + 6*a**2*n**2*x*(a + b*x**n)*(c + d*x**n)**(-1 - 1/n)/(c**3*(n + 1)*(2*n + 1)*(3*n + 1)) + 3*a*n*x*(a + b*x**n)**2*(c + d*x**n)**(-2 - 1/n)/(c**2*(6*n**2 + 5*n + 1)) + x*(a + b*x**n)**3*(c + d*x**n)**(-3 - 1/n)/(c*(3*n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_221():
    f = (a + b*x**n)**2*(c + d*x**n)**(-3 - 1/n)
    F = 2*a**2*n**2*x/(c**3*(c + d*x**n)**(1/n)*(n + 1)*(2*n + 1)) + 2*a*n*x*(a + b*x**n)*(c + d*x**n)**(-1 - 1/n)/(c**2*(n + 1)*(2*n + 1)) + x*(a + b*x**n)**2*(c + d*x**n)**(-2 - 1/n)/(c*(2*n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_222():
    f = (a + b*x**n)*(c + d*x**n)**(-2 - 1/n)
    F = a*n*x/(c**2*(c + d*x**n)**(1/n)*(n + 1)) + x*(a + b*x**n)*(c + d*x**n)**(-1 - 1/n)/(c*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_223():
    f = (c + d*x**n)**(-1 - 1/n)
    F = x/(c*(c + d*x**n)**(1/n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_224():
    f = 1/((a + b*x**n)*(c + d*x**n)**(1/n))
    F = x*hyper((1, 1/n), (1 + 1/n,), -x**n*(-a*d + b*c)/(a*(c + d*x**n)))/(a*(c + d*x**n)**(1/n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_225():
    f = (c + d*x**n)**(1 - 1/n)/(a + b*x**n)**2
    F = c*x*hyper((2, 1/n), (1 + 1/n,), -x**n*(-a*d + b*c)/(a*(c + d*x**n)))/(a**2*(c + d*x**n)**(1/n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_226():
    f = (c + d*x**n)**(2 - 1/n)/(a + b*x**n)**3
    F = c**2*x*hyper((3, 1/n), (1 + 1/n,), -x**n*(-a*d + b*c)/(a*(c + d*x**n)))/(a**3*(c + d*x**n)**(1/n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_227():
    f = (a + b*x**n)**p*(c + d*x**n)**(-p - 2 - 1/n)
    F = -b*x*(a + b*x**n)**(p + 1)*(c + d*x**n)**(-p - 1 - 1/n)/(a*n*(p + 1)*(-a*d + b*c)) + x*(c*(a + b*x**n)/(a*(c + d*x**n)))**(-p - 1)*(a + b*x**n)**(p + 1)*(c + d*x**n)**(-p - 1 - 1/n)*(b*c + n*(p + 1)*(-a*d + b*c))*hyper((1/n, -p - 1), (1 + 1/n,), -x**n*(-a*d + b*c)/(a*(c + d*x**n)))/(a*c*n*(p + 1)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_228():
    f = (a + b*x**n)**((a*d*n - b*c*(n + 1))/(n*(-a*d + b*c)))*(c + d*x**n)**((a*d*n + a*d - b*c*n)/(-a*d*n + b*c*n))
    F = x*(c + d*x**n)**(a*d/(n*(-a*d + b*c)))/(a*c*(a + b*x**n)**(b*c/(n*(-a*d + b*c))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_229():
    f = (a + b*x**n)**2*(c + d*x**n)**(-4 - 1/n)
    F = -2*a**2*n**2*x*(3*a*d*n - b*(3*c*n + c))/(c**4*(c + d*x**n)**(1/n)*(n + 1)*(2*n + 1)*(3*n + 1)*(-a*d + b*c)) - 2*a*n*x*(a + b*x**n)*(c + d*x**n)**(-1 - 1/n)*(3*a*d*n - b*(3*c*n + c))/(c**3*(n + 1)*(2*n + 1)*(3*n + 1)*(-a*d + b*c)) - x*(a + b*x**n)**2*(c + d*x**n)**(-2 - 1/n)*(3*a*d*n - b*(3*c*n + c))/(c**2*(-a*d + b*c)*(6*n**2 + 5*n + 1)) - b*x*(a + b*x**n)**3*(c + d*x**n)**(-3 - 1/n)/(3*a*n*(-a*d + b*c)) - x*(a + b*x**n)**3*(c + d*x**n)**(-3 - 1/n)*(3*a*d*n - b*(3*c*n + c))/(3*a*c*n*(3*n + 1)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_230():
    f = (a + b*x**n)*(c + d*x**n)**(-3 - 1/n)
    F = -x*(c + d*x**n)**(-2 - 1/n)*(-a*d + b*c)/(c*d*(2*n + 1)) + x*(c + d*x**n)**(-1 - 1/n)*(2*a*d*n + b*c)/(c**2*d*(n + 1)*(2*n + 1)) + n*x*(2*a*d*n + b*c)/(c**3*d*(c + d*x**n)**(1/n)*(n + 1)*(2*n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_231():
    f = (c + d*x**n)**(-2 - 1/n)
    F = x*(c + d*x**n)**(-1 - 1/n)/(c*(n + 1)) + n*x/(c**2*(c + d*x**n)**(1/n)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_232():
    f = (c + d*x**n)**(-1 - 1/n)/(a + b*x**n)
    F = -d*x/(c*(c + d*x**n)**(1/n)*(-a*d + b*c)) + b*x*hyper((1, 1/n), (1 + 1/n,), -x**n*(-a*d + b*c)/(a*(c + d*x**n)))/(a*(c + d*x**n)**(1/n)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_233():
    f = 1/((a + b*x**n)**2*(c + d*x**n)**(1/n))
    F = b*x/(a*n*(a + b*x**n)*(c + d*x**n)**((1 - n)/n)*(-a*d + b*c)) - x*(a*d*n + b*c*(1 - n))*hyper((1, 1/n), (1 + 1/n,), -x**n*(-a*d + b*c)/(a*(c + d*x**n)))/(a**2*n*(c + d*x**n)**(1/n)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_234():
    f = (c + d*x**n)**(1 - 1/n)/(a + b*x**n)**3
    F = b*x*(c + d*x**n)**(2 - 1/n)/(2*a*n*(a + b*x**n)**2*(-a*d + b*c)) - c*x*(2*a*d*n + b*c*(1 - 2*n))*hyper((2, 1/n), (1 + 1/n,), -x**n*(-a*d + b*c)/(a*(c + d*x**n)))/(2*a**3*n*(c + d*x**n)**(1/n)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_235():
    f = (c + d*x**n)**(2 - 1/n)/(a + b*x**n)**4
    F = b*x*(c + d*x**n)**(3 - 1/n)/(3*a*n*(a + b*x**n)**3*(-a*d + b*c)) - c**2*x*(3*a*d*n + b*c*(1 - 3*n))*hyper((3, 1/n), (1 + 1/n,), -x**n*(-a*d + b*c)/(a*(c + d*x**n)))/(3*a**4*n*(c + d*x**n)**(1/n)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_236():
    f = (a + b*x**2)*sqrt(-c + d*x)*sqrt(c + d*x)/x
    F = -a*c*atan(sqrt(-c + d*x)*sqrt(c + d*x)/c) + a*sqrt(-c + d*x)*sqrt(c + d*x) + b*(-c + d*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2)/(3*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_237():
    f = x**4*(a + b*x**2)*sqrt(-c + d*x)*sqrt(c + d*x)
    F = b*x**5*(-c + d*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2)/(8*d**2) - c**6*(8*a*d**2 + 5*b*c**2)*atanh(sqrt(-c + d*x)/sqrt(c + d*x))/(64*d**7) + c**4*x*sqrt(-c + d*x)*sqrt(c + d*x)*(8*a*d**2 + 5*b*c**2)/(128*d**6) + c**2*x*(-c + d*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2)*(8*a*d**2 + 5*b*c**2)/(64*d**6) + x**3*(-c + d*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2)*(8*a*d**2 + 5*b*c**2)/(48*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_238():
    f = x**2*(a + b*x**2)*sqrt(-c + d*x)*sqrt(c + d*x)
    F = b*x**3*(-c + d*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2)/(6*d**2) - c**4*(2*a*d**2 + b*c**2)*atanh(sqrt(-c + d*x)/sqrt(c + d*x))/(8*d**5) + c**2*x*sqrt(-c + d*x)*sqrt(c + d*x)*(2*a*d**2 + b*c**2)/(16*d**4) + x*(-c + d*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2)*(2*a*d**2 + b*c**2)/(8*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_239():
    f = (a + b*x**2)*sqrt(-c + d*x)*sqrt(c + d*x)
    F = b*x*(-c + d*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2)/(4*d**2) - c**2*(4*a*d**2 + b*c**2)*atanh(sqrt(-c + d*x)/sqrt(c + d*x))/(4*d**3) + x*sqrt(-c + d*x)*sqrt(c + d*x)*(4*a*d**2 + b*c**2)/(8*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_240():
    f = (a + b*x**2)*sqrt(-c + d*x)*sqrt(c + d*x)/x**2
    F = a*(-c + d*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2)/(c**2*x) + x*sqrt(-c + d*x)*sqrt(c + d*x)*(-a*d**2/c**2 + b/2) - (-2*a*d**2 + b*c**2)*atanh(sqrt(-c + d*x)/sqrt(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_241():
    f = (a + b*x**2)*sqrt(-c + d*x)*sqrt(c + d*x)/x**4
    F = a*(-c + d*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2)/(3*c**2*x**3) + 2*b*d*atanh(sqrt(-c + d*x)/sqrt(c + d*x)) - b*sqrt(-c + d*x)*sqrt(c + d*x)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_242():
    f = x**4*(a + b*x**2)/(sqrt(c*x - 1)*sqrt(c*x + 1))
    F = b*x**5*sqrt(c*x - 1)*sqrt(c*x + 1)/(6*c**2) + x**3*(6*a*c**2 + 5*b)*sqrt(c*x - 1)*sqrt(c*x + 1)/(24*c**4) + x*(6*a*c**2 + 5*b)*sqrt(c*x - 1)*sqrt(c*x + 1)/(16*c**6) + (6*a*c**2 + 5*b)*acosh(c*x)/(16*c**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_243():
    f = x**3*(a + b*x**2)/(sqrt(c*x - 1)*sqrt(c*x + 1))
    F = b*x**4*sqrt(c*x - 1)*sqrt(c*x + 1)/(5*c**2) + x**2*(5*a*c**2 + 4*b)*sqrt(c*x - 1)*sqrt(c*x + 1)/(15*c**4) + (10*a*c**2 + 8*b)*sqrt(c*x - 1)*sqrt(c*x + 1)/(15*c**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_244():
    f = x**2*(a + b*x**2)/(sqrt(c*x - 1)*sqrt(c*x + 1))
    F = b*x**3*sqrt(c*x - 1)*sqrt(c*x + 1)/(4*c**2) + x*(4*a*c**2 + 3*b)*sqrt(c*x - 1)*sqrt(c*x + 1)/(8*c**4) + (4*a*c**2 + 3*b)*acosh(c*x)/(8*c**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_245():
    f = x*(a + b*x**2)/(sqrt(c*x - 1)*sqrt(c*x + 1))
    F = b*x**2*sqrt(c*x - 1)*sqrt(c*x + 1)/(3*c**2) + (3*a*c**2 + 2*b)*sqrt(c*x - 1)*sqrt(c*x + 1)/(3*c**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_246():
    f = (a + b*x**2)/(sqrt(c*x - 1)*sqrt(c*x + 1))
    F = b*x*sqrt(c*x - 1)*sqrt(c*x + 1)/(2*c**2) + (2*a*c**2 + b)*acosh(c*x)/(2*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_247():
    f = (a + b*x**2)/(x*sqrt(c*x - 1)*sqrt(c*x + 1))
    F = a*atan(sqrt(c*x - 1)*sqrt(c*x + 1)) + b*sqrt(c*x - 1)*sqrt(c*x + 1)/c**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_248():
    f = (a + b*x**2)/(x**2*sqrt(c*x - 1)*sqrt(c*x + 1))
    F = a*sqrt(c*x - 1)*sqrt(c*x + 1)/x + b*acosh(c*x)/c
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_249():
    f = (a + b*x**2)/(x**3*sqrt(c*x - 1)*sqrt(c*x + 1))
    F = a*sqrt(c*x - 1)*sqrt(c*x + 1)/(2*x**2) + (a*c**2/2 + b)*atan(sqrt(c*x - 1)*sqrt(c*x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_250():
    f = (a + b*x**2)/(x**4*sqrt(c*x - 1)*sqrt(c*x + 1))
    F = a*sqrt(c*x - 1)*sqrt(c*x + 1)/(3*x**3) + (2*a*c**2 + 3*b)*sqrt(c*x - 1)*sqrt(c*x + 1)/(3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_251():
    f = (a + b*x**2)/(x**5*sqrt(c*x - 1)*sqrt(c*x + 1))
    F = a*sqrt(c*x - 1)*sqrt(c*x + 1)/(4*x**4) + c**2*(3*a*c**2 + 4*b)*atan(sqrt(c*x - 1)*sqrt(c*x + 1))/8 + (3*a*c**2 + 4*b)*sqrt(c*x - 1)*sqrt(c*x + 1)/(8*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_252():
    f = x**4*(a + b*x**2)/(sqrt(-c + d*x)*sqrt(c + d*x))
    F = b*x**5*sqrt(-c + d*x)*sqrt(c + d*x)/(6*d**2) + c**4*(6*a*d**2 + 5*b*c**2)*atanh(sqrt(-c + d*x)/sqrt(c + d*x))/(8*d**7) + c**2*x*sqrt(-c + d*x)*sqrt(c + d*x)*(6*a*d**2 + 5*b*c**2)/(16*d**6) + x**3*sqrt(-c + d*x)*sqrt(c + d*x)*(6*a*d**2 + 5*b*c**2)/(24*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_253():
    f = x**3*(a + b*x**2)/(sqrt(-c + d*x)*sqrt(c + d*x))
    F = b*x**4*sqrt(-c + d*x)*sqrt(c + d*x)/(5*d**2) + 2*c**2*sqrt(-c + d*x)*sqrt(c + d*x)*(5*a*d**2 + 4*b*c**2)/(15*d**6) + x**2*sqrt(-c + d*x)*sqrt(c + d*x)*(5*a*d**2 + 4*b*c**2)/(15*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_254():
    f = x**2*(a + b*x**2)/(sqrt(-c + d*x)*sqrt(c + d*x))
    F = b*x**3*sqrt(-c + d*x)*sqrt(c + d*x)/(4*d**2) + c**2*(4*a*d**2 + 3*b*c**2)*atanh(sqrt(-c + d*x)/sqrt(c + d*x))/(4*d**5) + x*sqrt(-c + d*x)*sqrt(c + d*x)*(4*a*d**2 + 3*b*c**2)/(8*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_255():
    f = x*(a + b*x**2)/(sqrt(-c + d*x)*sqrt(c + d*x))
    F = b*x**2*sqrt(-c + d*x)*sqrt(c + d*x)/(3*d**2) + sqrt(-c + d*x)*sqrt(c + d*x)*(3*a*d**2 + 2*b*c**2)/(3*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_256():
    f = (a + b*x**2)/(sqrt(-c + d*x)*sqrt(c + d*x))
    F = b*x*sqrt(-c + d*x)*sqrt(c + d*x)/(2*d**2) + (2*a*d**2 + b*c**2)*atanh(sqrt(-c + d*x)/sqrt(c + d*x))/d**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_257():
    f = (a + b*x**2)/(x*sqrt(-c + d*x)*sqrt(c + d*x))
    F = a*atan(sqrt(-c + d*x)*sqrt(c + d*x)/c)/c + b*sqrt(-c + d*x)*sqrt(c + d*x)/d**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_258():
    f = (a + b*x**2)/(x**2*sqrt(-c + d*x)*sqrt(c + d*x))
    F = a*sqrt(-c + d*x)*sqrt(c + d*x)/(c**2*x) + 2*b*atanh(sqrt(-c + d*x)/sqrt(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_259():
    f = (a + b*x**2)/(x**3*sqrt(-c + d*x)*sqrt(c + d*x))
    F = a*sqrt(-c + d*x)*sqrt(c + d*x)/(2*c**2*x**2) + (a*d**2 + 2*b*c**2)*atan(sqrt(-c + d*x)*sqrt(c + d*x)/c)/(2*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_260():
    f = (a + b*x**2)/(x**4*sqrt(-c + d*x)*sqrt(c + d*x))
    F = a*sqrt(-c + d*x)*sqrt(c + d*x)/(3*c**2*x**3) + sqrt(-c + d*x)*sqrt(c + d*x)*(2*a*d**2 + 3*b*c**2)/(3*c**4*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_261():
    f = (a + b*x**2)/(x**5*sqrt(-c + d*x)*sqrt(c + d*x))
    F = a*sqrt(-c + d*x)*sqrt(c + d*x)/(4*c**2*x**4) + sqrt(-c + d*x)*sqrt(c + d*x)*(3*a*d**2 + 4*b*c**2)/(8*c**4*x**2) + d**2*(3*a*d**2 + 4*b*c**2)*atan(sqrt(-c + d*x)*sqrt(c + d*x)/c)/(8*c**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_262():
    f = x**4*(a + b*x**2)/((-c + d*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2))
    F = b*x**5/(4*d**2*sqrt(-c + d*x)*sqrt(c + d*x)) + 3*c**2*(4*a*d**2 + 5*b*c**2)*atanh(sqrt(-c + d*x)/sqrt(c + d*x))/(4*d**7) - x**3*(4*a*d**2 + 5*b*c**2)/(4*d**4*sqrt(-c + d*x)*sqrt(c + d*x)) + x*sqrt(-c + d*x)*sqrt(c + d*x)*(12*a*d**2 + 15*b*c**2)/(8*d**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_263():
    f = x**3*(a + b*x**2)/((-c + d*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2))
    F = b*x**4/(3*d**2*sqrt(-c + d*x)*sqrt(c + d*x)) - x**2*(3*a*d**2 + 4*b*c**2)/(3*d**4*sqrt(-c + d*x)*sqrt(c + d*x)) + sqrt(-c + d*x)*sqrt(c + d*x)*(6*a*d**2 + 8*b*c**2)/(3*d**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_264():
    f = x**2*(a + b*x**2)/((-c + d*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2))
    F = b*x**3/(2*d**2*sqrt(-c + d*x)*sqrt(c + d*x)) - c*(2*a*d**2 + 3*b*c**2)/(2*d**5*sqrt(-c + d*x)*sqrt(c + d*x)) - sqrt(-c + d*x)*(2*a*d**2 + 3*b*c**2)/(2*d**5*sqrt(c + d*x)) + (2*a*d**2 + 3*b*c**2)*atanh(sqrt(-c + d*x)/sqrt(c + d*x))/d**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_265():
    f = x*(a + b*x**2)/((-c + d*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2))
    F = -x**2*(a/c**2 + b/d**2)/(sqrt(-c + d*x)*sqrt(c + d*x)) + sqrt(-c + d*x)*sqrt(c + d*x)*(a*d**2 + 2*b*c**2)/(c**2*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_266():
    f = (a + b*x**2)/((-c + d*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2))
    F = 2*b*atanh(sqrt(-c + d*x)/sqrt(c + d*x))/d**3 - x*(a/c**2 + b/d**2)/(sqrt(-c + d*x)*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_267():
    f = (a + b*x**2)/(x*(-c + d*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2))
    F = -a*atan(sqrt(-c + d*x)*sqrt(c + d*x)/c)/c**3 - (a/c**2 + b/d**2)/(sqrt(-c + d*x)*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_268():
    f = (a + b*x**2)/(x**2*(-c + d*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2))
    F = a/(c**2*x*sqrt(-c + d*x)*sqrt(c + d*x)) - x*(2*a*d**2 + b*c**2)/(c**4*sqrt(-c + d*x)*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_269():
    f = (a + b*x**2)/(x**3*(-c + d*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2))
    F = a/(2*c**2*x**2*sqrt(-c + d*x)*sqrt(c + d*x)) - (3*a*d**2 + 2*b*c**2)/(2*c**4*sqrt(-c + d*x)*sqrt(c + d*x)) - (3*a*d**2 + 2*b*c**2)*atan(sqrt(-c + d*x)*sqrt(c + d*x)/c)/(2*c**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_270():
    f = (a + b*x**2)/(x**4*(-c + d*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2))
    F = a/(3*c**2*x**3*sqrt(-c + d*x)*sqrt(c + d*x)) + (4*a*d**2 + 3*b*c**2)/(3*c**4*x*sqrt(-c + d*x)*sqrt(c + d*x)) - 2*d**2*x*(4*a*d**2 + 3*b*c**2)/(3*c**6*sqrt(-c + d*x)*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_271():
    f = (a + b*x**2)/(x**5*(-c + d*x)**(sympy.S(3)/2)*(c + d*x)**(sympy.S(3)/2))
    F = a/(4*c**2*x**4*sqrt(-c + d*x)*sqrt(c + d*x)) + (5*a*d**2 + 4*b*c**2)/(8*c**4*x**2*sqrt(-c + d*x)*sqrt(c + d*x)) - 3*d**2*(5*a*d**2 + 4*b*c**2)/(8*c**6*sqrt(-c + d*x)*sqrt(c + d*x)) - 3*d**2*(5*a*d**2 + 4*b*c**2)*atan(sqrt(-c + d*x)*sqrt(c + d*x)/c)/(8*c**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_272():
    f = (c**2*x**2 + 1)/(x*sqrt(c*x - 1)*sqrt(c*x + 1))
    F = sqrt(c*x - 1)*sqrt(c*x + 1) + atan(sqrt(c*x - 1)*sqrt(c*x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_273():
    f = (c + d*x**2)/(x**((a**2*d + 2*b**2*c)/(a**2*d + b**2*c))*sqrt(-a + b*x)*sqrt(a + b*x))
    F = sqrt(-a + b*x)*sqrt(a + b*x)*(d/b**2 + c/a**2)/x**(b**2*c/(a**2*d + b**2*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_274():
    f = 1/(sqrt(-sqrt(x) - 1)*sqrt(sqrt(x) - 1)*sqrt(x + 1))
    F = sqrt(1 - x)*asin(x)/(sqrt(-sqrt(x) - 1)*sqrt(sqrt(x) - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_275():
    f = 1/(sqrt(a - b*sqrt(x))*sqrt(a + b*sqrt(x))*sqrt(a**2 + b**2*x))
    F = -2*sqrt(a**2 - b**2*x)*atan(sqrt(a**2 - b**2*x)/sqrt(a**2 + b**2*x))/(b**2*sqrt(a - b*sqrt(x))*sqrt(a + b*sqrt(x)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_276():
    f = (a - b*x**n)**p*(a + b*x**n)**p*(c + d*x**(2*n))**q
    F = x*(a - b*x**n)**p*(a + b*x**n)**p*(c + d*x**(2*n))**q*appellf1(1/(2*n), -p, -q, 1 + 1/(2*n), b**2*x**(2*n)/a**2, -d*x**(2*n)/c)/((1 - b**2*x**(2*n)/a**2)**p*(1 + d*x**(2*n)/c)**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_277():
    f = (a - b*x**n)**p*(a + b*x**n)**p*(a**2 + b**2*x**(2*n))**p
    F = x*(a - b*x**n)**p*(a + b*x**n)**p*(a**2 + b**2*x**(2*n))**p*hyper((-p, 1/(4*n)), (1 + 1/(4*n),), b**4*x**(4*n)/a**4)/(1 - b**4*x**(4*n)/a**4)**p
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_3_General_1_1_3_3_a_plus_b_x_pow_n_pow_p_c_plus_d_x_pow_n_pow_q_278():
    f = (c + d*x**(2*n))**p/((a - b*x**n)*(a + b*x**n))
    F = x*(c + d*x**(2*n))**p*appellf1(1/(2*n), 1, -p, 1 + 1/(2*n), b**2*x**(2*n)/a**2, -d*x**(2*n)/c)/(a**2*(1 + d*x**(2*n)/c)**p)
    assert integrate(f, x) == F

