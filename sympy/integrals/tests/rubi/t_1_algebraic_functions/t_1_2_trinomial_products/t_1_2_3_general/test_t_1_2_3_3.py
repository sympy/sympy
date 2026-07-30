"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.2 Trinomial products/1.2.3 General/1.2.3.3 (d+e x^n)^q (a+b x^n+c x^(2 n))^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, f, n, p, q = symbols('a b c d e f n p q')

def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_1():
    f = (d + e*x**3)/(a + c*x**6)
    F = -e*log(a**(sympy.S(1)/3) + c**(sympy.S(1)/3)*x**2)/(6*a**(sympy.S(1)/3)*c**(sympy.S(2)/3)) - (-sqrt(a)*e + sqrt(3)*sqrt(c)*d)*log(-sqrt(3)*a**(sympy.S(1)/6)*c**(sympy.S(1)/6)*x + a**(sympy.S(1)/3) + c**(sympy.S(1)/3)*x**2)/(12*a**(sympy.S(5)/6)*c**(sympy.S(2)/3)) + (sqrt(a)*e + sqrt(3)*sqrt(c)*d)*log(sqrt(3)*a**(sympy.S(1)/6)*c**(sympy.S(1)/6)*x + a**(sympy.S(1)/3) + c**(sympy.S(1)/3)*x**2)/(12*a**(sympy.S(5)/6)*c**(sympy.S(2)/3)) + (-sqrt(3)*sqrt(a)*e + sqrt(c)*d)*atan(sqrt(3) + 2*c**(sympy.S(1)/6)*x/a**(sympy.S(1)/6))/(6*a**(sympy.S(5)/6)*c**(sympy.S(2)/3)) - (sqrt(3)*sqrt(a)*e + sqrt(c)*d)*atan(sqrt(3) - 2*c**(sympy.S(1)/6)*x/a**(sympy.S(1)/6))/(6*a**(sympy.S(5)/6)*c**(sympy.S(2)/3)) + d*atan(c**(sympy.S(1)/6)*x/a**(sympy.S(1)/6))/(3*a**(sympy.S(5)/6)*c**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_2():
    f = (d + e*x**3)/(a - c*x**6)
    F = -(sqrt(a)*e + sqrt(c)*d)*log(a**(sympy.S(1)/6) - c**(sympy.S(1)/6)*x)/(6*a**(sympy.S(5)/6)*c**(sympy.S(2)/3)) + (sqrt(a)*e + sqrt(c)*d)*log(a**(sympy.S(1)/6)*c**(sympy.S(1)/6)*x + a**(sympy.S(1)/3) + c**(sympy.S(1)/3)*x**2)/(12*a**(sympy.S(5)/6)*c**(sympy.S(2)/3)) + sqrt(3)*(sqrt(a)*e + sqrt(c)*d)*atan(sqrt(3)*(a**(sympy.S(1)/6) + 2*c**(sympy.S(1)/6)*x)/(3*a**(sympy.S(1)/6)))/(6*a**(sympy.S(5)/6)*c**(sympy.S(2)/3)) + (-sqrt(a)*e/sqrt(c) + d)*log(a**(sympy.S(1)/6) + c**(sympy.S(1)/6)*x)/(6*a**(sympy.S(5)/6)*c**(sympy.S(1)/6)) - (-sqrt(a)*e/sqrt(c) + d)*log(-a**(sympy.S(1)/6)*c**(sympy.S(1)/6)*x + a**(sympy.S(1)/3) + c**(sympy.S(1)/3)*x**2)/(12*a**(sympy.S(5)/6)*c**(sympy.S(1)/6)) - sqrt(3)*(-sqrt(a)*e/sqrt(c) + d)*atan(sqrt(3)*(a**(sympy.S(1)/6) - 2*c**(sympy.S(1)/6)*x)/(3*a**(sympy.S(1)/6)))/(6*a**(sympy.S(5)/6)*c**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_3():
    f = (d + e*x**4)/(a + c*x**8)
    F = (-sqrt(a)*e/sqrt(c) + d + sqrt(2)*d)*log(a**(sympy.S(1)/8)*c**(sympy.S(1)/8)*x*sqrt(sqrt(2) + 2) + a**(sympy.S(1)/4) + c**(sympy.S(1)/4)*x**2)/(8*a**(sympy.S(7)/8)*c**(sympy.S(1)/8)*sqrt(2*sqrt(2) + 4)) + (-sqrt(a)*e + sqrt(c)*d*(1 - sqrt(2)))*log(-a**(sympy.S(1)/8)*c**(sympy.S(1)/8)*x*sqrt(2 - sqrt(2)) + a**(sympy.S(1)/4) + c**(sympy.S(1)/4)*x**2)/(8*a**(sympy.S(7)/8)*c**(sympy.S(5)/8)*sqrt(4 - 2*sqrt(2))) - (-sqrt(a)*e + sqrt(c)*d*(1 - sqrt(2)))*log(a**(sympy.S(1)/8)*c**(sympy.S(1)/8)*x*sqrt(2 - sqrt(2)) + a**(sympy.S(1)/4) + c**(sympy.S(1)/4)*x**2)/(8*a**(sympy.S(7)/8)*c**(sympy.S(5)/8)*sqrt(4 - 2*sqrt(2))) + sqrt(sqrt(2) + 2)*(-sqrt(a)*e + sqrt(c)*d*(1 - sqrt(2)))*atan((a**(sympy.S(1)/8)*sqrt(sqrt(2) + 2) - 2*c**(sympy.S(1)/8)*x)/(a**(sympy.S(1)/8)*sqrt(2 - sqrt(2))))/(8*a**(sympy.S(7)/8)*c**(sympy.S(5)/8)) - sqrt(sqrt(2) + 2)*(-sqrt(a)*e + sqrt(c)*d*(1 - sqrt(2)))*atan((a**(sympy.S(1)/8)*sqrt(sqrt(2) + 2) + 2*c**(sympy.S(1)/8)*x)/(a**(sympy.S(1)/8)*sqrt(2 - sqrt(2))))/(8*a**(sympy.S(7)/8)*c**(sympy.S(5)/8)) - (-sqrt(a)*e + sqrt(c)*d*(1 + sqrt(2)))*log(-a**(sympy.S(1)/8)*c**(sympy.S(1)/8)*x*sqrt(sqrt(2) + 2) + a**(sympy.S(1)/4) + c**(sympy.S(1)/4)*x**2)/(8*a**(sympy.S(7)/8)*c**(sympy.S(5)/8)*sqrt(2*sqrt(2) + 4)) - sqrt(2 - sqrt(2))*(-sqrt(a)*e + sqrt(c)*d*(1 + sqrt(2)))*atan((a**(sympy.S(1)/8)*sqrt(2 - sqrt(2)) - 2*c**(sympy.S(1)/8)*x)/(a**(sympy.S(1)/8)*sqrt(sqrt(2) + 2)))/(8*a**(sympy.S(7)/8)*c**(sympy.S(5)/8)) + sqrt(2 - sqrt(2))*(-sqrt(a)*e + sqrt(c)*d*(1 + sqrt(2)))*atan((a**(sympy.S(1)/8)*sqrt(2 - sqrt(2)) + 2*c**(sympy.S(1)/8)*x)/(a**(sympy.S(1)/8)*sqrt(sqrt(2) + 2)))/(8*a**(sympy.S(7)/8)*c**(sympy.S(5)/8))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_4():
    f = (d + e*x**4)/(a - c*x**8)
    F = -sqrt(2)*(-sqrt(a)*e/sqrt(c) + d)*log(-sqrt(2)*a**(sympy.S(1)/8)*c**(sympy.S(1)/8)*x + a**(sympy.S(1)/4) + c**(sympy.S(1)/4)*x**2)/(16*a**(sympy.S(7)/8)*c**(sympy.S(1)/8)) + sqrt(2)*(-sqrt(a)*e/sqrt(c) + d)*log(sqrt(2)*a**(sympy.S(1)/8)*c**(sympy.S(1)/8)*x + a**(sympy.S(1)/4) + c**(sympy.S(1)/4)*x**2)/(16*a**(sympy.S(7)/8)*c**(sympy.S(1)/8)) - sqrt(2)*(-sqrt(a)*e/sqrt(c) + d)*atan(1 - sqrt(2)*c**(sympy.S(1)/8)*x/a**(sympy.S(1)/8))/(8*a**(sympy.S(7)/8)*c**(sympy.S(1)/8)) + sqrt(2)*(-sqrt(a)*e/sqrt(c) + d)*atan(1 + sqrt(2)*c**(sympy.S(1)/8)*x/a**(sympy.S(1)/8))/(8*a**(sympy.S(7)/8)*c**(sympy.S(1)/8)) + (sqrt(a)*e + sqrt(c)*d)*atan(c**(sympy.S(1)/8)*x/a**(sympy.S(1)/8))/(4*a**(sympy.S(7)/8)*c**(sympy.S(5)/8)) + (sqrt(a)*e + sqrt(c)*d)*atanh(c**(sympy.S(1)/8)*x/a**(sympy.S(1)/8))/(4*a**(sympy.S(7)/8)*c**(sympy.S(5)/8))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_5():
    f = (d + e*x**4)/(b*x**4 + d**2 + e**2*x**8)
    F = -log(sqrt(d) + sqrt(e)*x**2 - x*sqrt(2*sqrt(d)*sqrt(e) + sqrt(-b + 2*d*e)))/(8*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) + sqrt(-b + 2*d*e))) + log(sqrt(d) + sqrt(e)*x**2 + x*sqrt(2*sqrt(d)*sqrt(e) + sqrt(-b + 2*d*e)))/(8*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) + sqrt(-b + 2*d*e))) - atan((-2*sqrt(e)*x + sqrt(2*sqrt(d)*sqrt(e) - sqrt(-b + 2*d*e)))/sqrt(2*sqrt(d)*sqrt(e) + sqrt(-b + 2*d*e)))/(4*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) + sqrt(-b + 2*d*e))) + atan((2*sqrt(e)*x + sqrt(2*sqrt(d)*sqrt(e) - sqrt(-b + 2*d*e)))/sqrt(2*sqrt(d)*sqrt(e) + sqrt(-b + 2*d*e)))/(4*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) + sqrt(-b + 2*d*e))) - log(sqrt(d) + sqrt(e)*x**2 - x*sqrt(2*sqrt(d)*sqrt(e) - sqrt(-b + 2*d*e)))/(8*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) - sqrt(-b + 2*d*e))) + log(sqrt(d) + sqrt(e)*x**2 + x*sqrt(2*sqrt(d)*sqrt(e) - sqrt(-b + 2*d*e)))/(8*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) - sqrt(-b + 2*d*e))) - atan((-2*sqrt(e)*x + sqrt(2*sqrt(d)*sqrt(e) + sqrt(-b + 2*d*e)))/sqrt(2*sqrt(d)*sqrt(e) - sqrt(-b + 2*d*e)))/(4*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) - sqrt(-b + 2*d*e))) + atan((2*sqrt(e)*x + sqrt(2*sqrt(d)*sqrt(e) + sqrt(-b + 2*d*e)))/sqrt(2*sqrt(d)*sqrt(e) - sqrt(-b + 2*d*e)))/(4*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) - sqrt(-b + 2*d*e)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_6():
    f = (d + e*x**4)/(d**2 + e**2*x**8 + f*x**4)
    F = -log(sqrt(d) + sqrt(e)*x**2 - x*sqrt(2*sqrt(d)*sqrt(e) + sqrt(2*d*e - f)))/(8*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) + sqrt(2*d*e - f))) + log(sqrt(d) + sqrt(e)*x**2 + x*sqrt(2*sqrt(d)*sqrt(e) + sqrt(2*d*e - f)))/(8*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) + sqrt(2*d*e - f))) - atan((-2*sqrt(e)*x + sqrt(2*sqrt(d)*sqrt(e) - sqrt(2*d*e - f)))/sqrt(2*sqrt(d)*sqrt(e) + sqrt(2*d*e - f)))/(4*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) + sqrt(2*d*e - f))) + atan((2*sqrt(e)*x + sqrt(2*sqrt(d)*sqrt(e) - sqrt(2*d*e - f)))/sqrt(2*sqrt(d)*sqrt(e) + sqrt(2*d*e - f)))/(4*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) + sqrt(2*d*e - f))) - log(sqrt(d) + sqrt(e)*x**2 - x*sqrt(2*sqrt(d)*sqrt(e) - sqrt(2*d*e - f)))/(8*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) - sqrt(2*d*e - f))) + log(sqrt(d) + sqrt(e)*x**2 + x*sqrt(2*sqrt(d)*sqrt(e) - sqrt(2*d*e - f)))/(8*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) - sqrt(2*d*e - f))) - atan((-2*sqrt(e)*x + sqrt(2*sqrt(d)*sqrt(e) + sqrt(2*d*e - f)))/sqrt(2*sqrt(d)*sqrt(e) - sqrt(2*d*e - f)))/(4*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) - sqrt(2*d*e - f))) + atan((2*sqrt(e)*x + sqrt(2*sqrt(d)*sqrt(e) + sqrt(2*d*e - f)))/sqrt(2*sqrt(d)*sqrt(e) - sqrt(2*d*e - f)))/(4*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) - sqrt(2*d*e - f)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_7():
    f = (d + e*x**4)/(-b*x**4 + d**2 + e**2*x**8)
    F = -sqrt(2)*sqrt(e)*atan(sqrt(2)*sqrt(e)*x/sqrt(sqrt(b - 2*d*e) + sqrt(b + 2*d*e)))/(2*sqrt(b - 2*d*e)*sqrt(sqrt(b - 2*d*e) + sqrt(b + 2*d*e))) - sqrt(2)*sqrt(e)*atanh(sqrt(2)*sqrt(e)*x/sqrt(sqrt(b - 2*d*e) + sqrt(b + 2*d*e)))/(2*sqrt(b - 2*d*e)*sqrt(sqrt(b - 2*d*e) + sqrt(b + 2*d*e))) - sqrt(2)*sqrt(e)*atan(sqrt(2)*sqrt(e)*x/sqrt(sqrt(b - 2*d*e) - sqrt(b + 2*d*e)))/(2*sqrt(b - 2*d*e)*sqrt(sqrt(b - 2*d*e) - sqrt(b + 2*d*e))) - sqrt(2)*sqrt(e)*atanh(sqrt(2)*sqrt(e)*x/sqrt(sqrt(b - 2*d*e) - sqrt(b + 2*d*e)))/(2*sqrt(b - 2*d*e)*sqrt(sqrt(b - 2*d*e) - sqrt(b + 2*d*e)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_8():
    f = (d + e*x**4)/(d**2 + e**2*x**8 - f*x**4)
    F = -log(sqrt(d) + sqrt(e)*x**2 - x*sqrt(2*sqrt(d)*sqrt(e) + sqrt(2*d*e + f)))/(8*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) + sqrt(2*d*e + f))) + log(sqrt(d) + sqrt(e)*x**2 + x*sqrt(2*sqrt(d)*sqrt(e) + sqrt(2*d*e + f)))/(8*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) + sqrt(2*d*e + f))) - atan((-2*sqrt(e)*x + sqrt(2*sqrt(d)*sqrt(e) - sqrt(2*d*e + f)))/sqrt(2*sqrt(d)*sqrt(e) + sqrt(2*d*e + f)))/(4*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) + sqrt(2*d*e + f))) + atan((2*sqrt(e)*x + sqrt(2*sqrt(d)*sqrt(e) - sqrt(2*d*e + f)))/sqrt(2*sqrt(d)*sqrt(e) + sqrt(2*d*e + f)))/(4*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) + sqrt(2*d*e + f))) - log(sqrt(d) + sqrt(e)*x**2 - x*sqrt(2*sqrt(d)*sqrt(e) - sqrt(2*d*e + f)))/(8*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) - sqrt(2*d*e + f))) + log(sqrt(d) + sqrt(e)*x**2 + x*sqrt(2*sqrt(d)*sqrt(e) - sqrt(2*d*e + f)))/(8*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) - sqrt(2*d*e + f))) - atan((-2*sqrt(e)*x + sqrt(2*sqrt(d)*sqrt(e) + sqrt(2*d*e + f)))/sqrt(2*sqrt(d)*sqrt(e) - sqrt(2*d*e + f)))/(4*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) - sqrt(2*d*e + f))) + atan((2*sqrt(e)*x + sqrt(2*sqrt(d)*sqrt(e) + sqrt(2*d*e + f)))/sqrt(2*sqrt(d)*sqrt(e) - sqrt(2*d*e + f)))/(4*sqrt(d)*sqrt(2*sqrt(d)*sqrt(e) - sqrt(2*d*e + f)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_9():
    f = (x**4 + 1)/(b*x**4 + x**8 + 1)
    F = -log(x**2 - x*sqrt(sqrt(2 - b) + 2) + 1)/(8*sqrt(sqrt(2 - b) + 2)) + log(x**2 + x*sqrt(sqrt(2 - b) + 2) + 1)/(8*sqrt(sqrt(2 - b) + 2)) - atan((-2*x + sqrt(2 - sqrt(2 - b)))/sqrt(sqrt(2 - b) + 2))/(4*sqrt(sqrt(2 - b) + 2)) + atan((2*x + sqrt(2 - sqrt(2 - b)))/sqrt(sqrt(2 - b) + 2))/(4*sqrt(sqrt(2 - b) + 2)) - log(x**2 - x*sqrt(2 - sqrt(2 - b)) + 1)/(8*sqrt(2 - sqrt(2 - b))) + log(x**2 + x*sqrt(2 - sqrt(2 - b)) + 1)/(8*sqrt(2 - sqrt(2 - b))) - atan((-2*x + sqrt(sqrt(2 - b) + 2))/sqrt(2 - sqrt(2 - b)))/(4*sqrt(2 - sqrt(2 - b))) + atan((2*x + sqrt(sqrt(2 - b) + 2))/sqrt(2 - sqrt(2 - b)))/(4*sqrt(2 - sqrt(2 - b)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_10():
    f = (x**4 + 1)/(x**8 + 3*x**4 + 1)
    F = -2**(sympy.S(1)/4)*sqrt(5)*(sqrt(5) + 3)**(sympy.S(1)/4)*log(2*x**2 - 2*x*(6 - 2*sqrt(5))**(sympy.S(1)/4) + sqrt(6 - 2*sqrt(5)))/40 + 2**(sympy.S(1)/4)*sqrt(5)*(sqrt(5) + 3)**(sympy.S(1)/4)*log(2*x**2 + 2*x*(6 - 2*sqrt(5))**(sympy.S(1)/4) + sqrt(6 - 2*sqrt(5)))/40 - 2**(sympy.S(1)/4)*sqrt(5)*(3 - sqrt(5))**(sympy.S(1)/4)*log(2*x**2 - 2*x*(2*sqrt(5) + 6)**(sympy.S(1)/4) + sqrt(2*sqrt(5) + 6))/40 + 2**(sympy.S(1)/4)*sqrt(5)*(3 - sqrt(5))**(sympy.S(1)/4)*log(2*x**2 + 2*x*(2*sqrt(5) + 6)**(sympy.S(1)/4) + sqrt(2*sqrt(5) + 6))/40 + 2**(sympy.S(1)/4)*sqrt(5)*(sqrt(5) + 3)**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(3 - sqrt(5))**(sympy.S(1)/4) - 1)/20 + 2**(sympy.S(1)/4)*sqrt(5)*(sqrt(5) + 3)**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(3 - sqrt(5))**(sympy.S(1)/4) + 1)/20 + 2**(sympy.S(1)/4)*sqrt(5)*(3 - sqrt(5))**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4) - 1)/20 + 2**(sympy.S(1)/4)*sqrt(5)*(3 - sqrt(5))**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4) + 1)/20
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_11():
    f = (x**4 + 1)/(x**8 + 2*x**4 + 1)
    F = -sqrt(2)*log(x**2 - sqrt(2)*x + 1)/8 + sqrt(2)*log(x**2 + sqrt(2)*x + 1)/8 + sqrt(2)*atan(sqrt(2)*x - 1)/4 + sqrt(2)*atan(sqrt(2)*x + 1)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_12():
    f = (x**4 + 1)/(x**8 + x**4 + 1)
    F = -log(x**2 - x + 1)/8 + log(x**2 + x + 1)/8 - sqrt(3)*log(x**2 - sqrt(3)*x + 1)/24 + sqrt(3)*log(x**2 + sqrt(3)*x + 1)/24 - sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/12 + sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/12 + atan(2*x - sqrt(3))/4 + atan(2*x + sqrt(3))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_13():
    f = (x**4 + 1)/(x**8 + 1)
    F = -log(x**2 - x*sqrt(2 - sqrt(2)) + 1)/(8*sqrt(2 - sqrt(2))) + log(x**2 + x*sqrt(2 - sqrt(2)) + 1)/(8*sqrt(2 - sqrt(2))) - log(x**2 - x*sqrt(sqrt(2) + 2) + 1)/(8*sqrt(sqrt(2) + 2)) + log(x**2 + x*sqrt(sqrt(2) + 2) + 1)/(8*sqrt(sqrt(2) + 2)) - sqrt(sqrt(2)/2 + 1)*atan((-2*x + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/4 + sqrt(sqrt(2)/2 + 1)*atan((2*x + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/4 - sqrt(1 - sqrt(2)/2)*atan((-2*x + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/4 + sqrt(1 - sqrt(2)/2)*atan((2*x + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_14():
    f = (x**4 + 1)/(x**8 - x**4 + 1)
    F = -log(x**2 - x*sqrt(2 - sqrt(3)) + 1)/(8*sqrt(2 - sqrt(3))) + log(x**2 + x*sqrt(2 - sqrt(3)) + 1)/(8*sqrt(2 - sqrt(3))) - log(x**2 - x*sqrt(sqrt(3) + 2) + 1)/(8*sqrt(sqrt(3) + 2)) + log(x**2 + x*sqrt(sqrt(3) + 2) + 1)/(8*sqrt(sqrt(3) + 2)) - sqrt(sqrt(3) + 2)*atan((-2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/4 + sqrt(sqrt(3) + 2)*atan((2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/4 - sqrt(2 - sqrt(3))*atan((-2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/4 + sqrt(2 - sqrt(3))*atan((2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_15():
    f = (x**4 + 1)/(x**8 - 2*x**4 + 1)
    F = x/(2 - 2*x**4) + atan(x)/4 + atanh(x)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_16():
    f = (x**4 + 1)/(x**8 - 3*x**4 + 1)
    F = atan(sqrt(2)*x/sqrt(-1 + sqrt(5)))/sqrt(-2 + 2*sqrt(5)) - atan(sqrt(2)*x/sqrt(1 + sqrt(5)))/sqrt(2 + 2*sqrt(5)) + atanh(sqrt(2)*x/sqrt(-1 + sqrt(5)))/sqrt(-2 + 2*sqrt(5)) - atanh(sqrt(2)*x/sqrt(1 + sqrt(5)))/sqrt(2 + 2*sqrt(5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_17():
    f = (x**4 + 1)/(x**8 - 4*x**4 + 1)
    F = 2**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*x/sqrt(-1 + sqrt(3)))/(4*sqrt(-1 + sqrt(3))) - 2**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*x/sqrt(1 + sqrt(3)))/(4*sqrt(1 + sqrt(3))) + 2**(sympy.S(3)/4)*atanh(2**(sympy.S(1)/4)*x/sqrt(-1 + sqrt(3)))/(4*sqrt(-1 + sqrt(3))) - 2**(sympy.S(3)/4)*atanh(2**(sympy.S(1)/4)*x/sqrt(1 + sqrt(3)))/(4*sqrt(1 + sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_18():
    f = (x**4 + 1)/(x**8 - 5*x**4 + 1)
    F = atan(sqrt(2)*x/sqrt(-sqrt(3) + sqrt(7)))/sqrt(-6*sqrt(3) + 6*sqrt(7)) - atan(sqrt(2)*x/sqrt(sqrt(3) + sqrt(7)))/sqrt(6*sqrt(3) + 6*sqrt(7)) + atanh(sqrt(2)*x/sqrt(-sqrt(3) + sqrt(7)))/sqrt(-6*sqrt(3) + 6*sqrt(7)) - atanh(sqrt(2)*x/sqrt(sqrt(3) + sqrt(7)))/sqrt(6*sqrt(3) + 6*sqrt(7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_19():
    f = (x**4 + 1)/(x**8 - 6*x**4 + 1)
    F = atan(x/sqrt(-1 + sqrt(2)))/(4*sqrt(-1 + sqrt(2))) - atan(x/sqrt(1 + sqrt(2)))/(4*sqrt(1 + sqrt(2))) + atanh(x/sqrt(-1 + sqrt(2)))/(4*sqrt(-1 + sqrt(2))) - atanh(x/sqrt(1 + sqrt(2)))/(4*sqrt(1 + sqrt(2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_20():
    f = (1 - x**4)/(b*x**4 + x**8 + 1)
    F = sqrt(2 - sqrt(2 - b))*log(x**2 - x*sqrt(2 - sqrt(2 - b)) + 1)/(8*sqrt(2 - b)) - sqrt(2 - sqrt(2 - b))*log(x**2 + x*sqrt(2 - sqrt(2 - b)) + 1)/(8*sqrt(2 - b)) + sqrt(b + 2)*atan((-2*x + sqrt(sqrt(2 - b) + 2))/sqrt(2 - sqrt(2 - b)))/(4*sqrt(2 - b)*sqrt(sqrt(2 - b) + 2)) - sqrt(b + 2)*atan((2*x + sqrt(sqrt(2 - b) + 2))/sqrt(2 - sqrt(2 - b)))/(4*sqrt(2 - b)*sqrt(sqrt(2 - b) + 2)) - sqrt(sqrt(2 - b) + 2)*log(x**2 - x*sqrt(sqrt(2 - b) + 2) + 1)/(8*sqrt(2 - b)) + sqrt(sqrt(2 - b) + 2)*log(x**2 + x*sqrt(sqrt(2 - b) + 2) + 1)/(8*sqrt(2 - b)) - sqrt(b + 2)*atan((-2*x + sqrt(2 - sqrt(2 - b)))/sqrt(sqrt(2 - b) + 2))/(4*sqrt(2 - b)*sqrt(2 - sqrt(2 - b))) + sqrt(b + 2)*atan((2*x + sqrt(2 - sqrt(2 - b)))/sqrt(sqrt(2 - b) + 2))/(4*sqrt(2 - b)*sqrt(2 - sqrt(2 - b)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_21():
    f = (1 - x**4)/(x**8 + 3*x**4 + 1)
    F = -2**(sympy.S(1)/4)*(sqrt(5) + 3)**(sympy.S(1)/4)*log(2*x**2 - 2*x*(6 - 2*sqrt(5))**(sympy.S(1)/4) + sqrt(6 - 2*sqrt(5)))/8 + 2**(sympy.S(1)/4)*(sqrt(5) + 3)**(sympy.S(1)/4)*log(2*x**2 + 2*x*(6 - 2*sqrt(5))**(sympy.S(1)/4) + sqrt(6 - 2*sqrt(5)))/8 + 2**(sympy.S(1)/4)*(3 - sqrt(5))**(sympy.S(1)/4)*log(2*x**2 - 2*x*(2*sqrt(5) + 6)**(sympy.S(1)/4) + sqrt(2*sqrt(5) + 6))/8 - 2**(sympy.S(1)/4)*(3 - sqrt(5))**(sympy.S(1)/4)*log(2*x**2 + 2*x*(2*sqrt(5) + 6)**(sympy.S(1)/4) + sqrt(2*sqrt(5) + 6))/8 + 2**(sympy.S(1)/4)*(sqrt(5) + 3)**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(3 - sqrt(5))**(sympy.S(1)/4) - 1)/4 + 2**(sympy.S(1)/4)*(sqrt(5) + 3)**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(3 - sqrt(5))**(sympy.S(1)/4) + 1)/4 - 2**(sympy.S(1)/4)*(3 - sqrt(5))**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4) - 1)/4 - 2**(sympy.S(1)/4)*(3 - sqrt(5))**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4) + 1)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_22():
    f = (1 - x**4)/(x**8 + 2*x**4 + 1)
    F = x/(2*x**4 + 2) - sqrt(2)*log(x**2 - sqrt(2)*x + 1)/16 + sqrt(2)*log(x**2 + sqrt(2)*x + 1)/16 + sqrt(2)*atan(sqrt(2)*x - 1)/8 + sqrt(2)*atan(sqrt(2)*x + 1)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_23():
    f = (1 - x**4)/(x**8 + x**4 + 1)
    F = log(x**2 - x + 1)/8 - log(x**2 + x + 1)/8 - sqrt(3)*log(x**2 - sqrt(3)*x + 1)/8 + sqrt(3)*log(x**2 + sqrt(3)*x + 1)/8 - sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/4 + sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/4 - atan(2*x - sqrt(3))/4 - atan(2*x + sqrt(3))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_24():
    f = (1 - x**4)/(x**8 + 1)
    F = sqrt(1 - sqrt(2)/2)*log(x**2 - x*sqrt(2 - sqrt(2)) + 1)/8 - sqrt(1 - sqrt(2)/2)*log(x**2 + x*sqrt(2 - sqrt(2)) + 1)/8 - sqrt(sqrt(2)/2 + 1)*log(x**2 - x*sqrt(sqrt(2) + 2) + 1)/8 + sqrt(sqrt(2)/2 + 1)*log(x**2 + x*sqrt(sqrt(2) + 2) + 1)/8 + atan((-2*x + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(4*sqrt(sqrt(2) + 2)) - atan((2*x + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(4*sqrt(sqrt(2) + 2)) - atan((-2*x + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(4*sqrt(2 - sqrt(2))) + atan((2*x + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(4*sqrt(2 - sqrt(2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_25():
    f = (1 - x**4)/(x**8 - x**4 + 1)
    F = sqrt(sympy.S(2)/3 - sqrt(3)/3)*log(x**2 - x*sqrt(2 - sqrt(3)) + 1)/8 - sqrt(sympy.S(2)/3 - sqrt(3)/3)*log(x**2 + x*sqrt(2 - sqrt(3)) + 1)/8 - sqrt(sqrt(3)/3 + sympy.S(2)/3)*log(x**2 - x*sqrt(sqrt(3) + 2) + 1)/8 + sqrt(sqrt(3)/3 + sympy.S(2)/3)*log(x**2 + x*sqrt(sqrt(3) + 2) + 1)/8 + atan((-2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/(4*sqrt(3*sqrt(3) + 6)) - atan((2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/(4*sqrt(3*sqrt(3) + 6)) - atan((-2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/(4*sqrt(6 - 3*sqrt(3))) + atan((2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/(4*sqrt(6 - 3*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_26():
    f = (1 - x**4)/(x**8 - 2*x**4 + 1)
    F = atan(x)/2 + atanh(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_27():
    f = (1 - x**4)/(x**8 - 3*x**4 + 1)
    F = atan(sqrt(2)*x/sqrt(-1 + sqrt(5)))/sqrt(-10 + 10*sqrt(5)) + atan(sqrt(2)*x/sqrt(1 + sqrt(5)))/sqrt(10 + 10*sqrt(5)) + atanh(sqrt(2)*x/sqrt(-1 + sqrt(5)))/sqrt(-10 + 10*sqrt(5)) + atanh(sqrt(2)*x/sqrt(1 + sqrt(5)))/sqrt(10 + 10*sqrt(5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_28():
    f = (1 - x**4)/(x**8 - 4*x**4 + 1)
    F = 2**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*x/sqrt(-1 + sqrt(3)))/(4*sqrt(-3 + 3*sqrt(3))) + 2**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*x/sqrt(1 + sqrt(3)))/(4*sqrt(3 + 3*sqrt(3))) + 2**(sympy.S(3)/4)*atanh(2**(sympy.S(1)/4)*x/sqrt(-1 + sqrt(3)))/(4*sqrt(-3 + 3*sqrt(3))) + 2**(sympy.S(3)/4)*atanh(2**(sympy.S(1)/4)*x/sqrt(1 + sqrt(3)))/(4*sqrt(3 + 3*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_29():
    f = (1 - x**4)/(x**8 - 5*x**4 + 1)
    F = atan(sqrt(2)*x/sqrt(-sqrt(3) + sqrt(7)))/sqrt(-14*sqrt(3) + 14*sqrt(7)) + atan(sqrt(2)*x/sqrt(sqrt(3) + sqrt(7)))/sqrt(14*sqrt(3) + 14*sqrt(7)) + atanh(sqrt(2)*x/sqrt(-sqrt(3) + sqrt(7)))/sqrt(-14*sqrt(3) + 14*sqrt(7)) + atanh(sqrt(2)*x/sqrt(sqrt(3) + sqrt(7)))/sqrt(14*sqrt(3) + 14*sqrt(7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_30():
    f = (1 - x**4)/(x**8 - 6*x**4 + 1)
    F = atan(x/sqrt(-1 + sqrt(2)))/(4*sqrt(-2 + 2*sqrt(2))) + atan(x/sqrt(1 + sqrt(2)))/(4*sqrt(2 + 2*sqrt(2))) + atanh(x/sqrt(-1 + sqrt(2)))/(4*sqrt(-2 + 2*sqrt(2))) + atanh(x/sqrt(1 + sqrt(2)))/(4*sqrt(2 + 2*sqrt(2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_31():
    f = (2*x**4 - 1 + sqrt(3))/(x**8 - x**4 + 1)
    F = -sqrt(2)*log(x**2 - x*sqrt(2 - sqrt(3)) + 1)/4 + sqrt(2)*log(x**2 + x*sqrt(2 - sqrt(3)) + 1)/4 - sqrt(2)*atan((-2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/2 + sqrt(2)*atan((2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_32():
    f = (x**4*(1 + sqrt(3)) + 1)/(x**8 - x**4 + 1)
    F = -sqrt(sqrt(3) + 2)*log(x**2 - x*sqrt(2 - sqrt(3)) + 1)/4 + sqrt(sqrt(3) + 2)*log(x**2 + x*sqrt(2 - sqrt(3)) + 1)/4 - sqrt(sqrt(3) + 2)*atan((-2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/2 + sqrt(sqrt(3) + 2)*atan((2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_33():
    f = (x**4*(-3 + sqrt(3)) - 2*sqrt(3) + 3)/(x**8 - x**4 + 1)
    F = sqrt(6 - 3*sqrt(3))*log(x**2 - x*sqrt(2 - sqrt(3)) + 1)/4 - sqrt(6 - 3*sqrt(3))*log(x**2 + x*sqrt(2 - sqrt(3)) + 1)/4 + sqrt(6 - 3*sqrt(3))*atan((-2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/2 - sqrt(6 - 3*sqrt(3))*atan((2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_34():
    f = (d + e/x)/(a/x**2 + c)
    F = -sqrt(a)*d*atan(sqrt(c)*x/sqrt(a))/c**(sympy.S(3)/2) + d*x/c + e*log(a + c*x**2)/(2*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_35():
    f = (d + e/x)/(a/x**2 + b/x + c)
    F = d*x/c - (b*d - c*e)*log(a + b*x + c*x**2)/(2*c**2) - (-2*a*c*d + b**2*d - b*c*e)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_36():
    f = (d + e/x**2)/(a/x**4 + c)
    F = d*x/c + sqrt(2)*(sqrt(a)*d - sqrt(c)*e)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(1)/4)*c**(sympy.S(5)/4)) - sqrt(2)*(sqrt(a)*d - sqrt(c)*e)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(1)/4)*c**(sympy.S(5)/4)) + sqrt(2)*(sqrt(a)*d + sqrt(c)*e)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(1)/4)*c**(sympy.S(5)/4)) - sqrt(2)*(sqrt(a)*d + sqrt(c)*e)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(1)/4)*c**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_37():
    f = (d + e/x**2)/(a/x**4 + b/x**2 + c)
    F = d*x/c - sqrt(2)*(b*d - c*e + (-2*a*c*d + b**2*d - b*c*e)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2))) - sqrt(2)*(b*d - c*e - (-2*a*c*d + b**2*d - b*c*e)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_38():
    f = (d + e/x**3)/(a/x**6 + c)
    F = -a**(sympy.S(1)/6)*d*atan(c**(sympy.S(1)/6)*x/a**(sympy.S(1)/6))/(3*c**(sympy.S(7)/6)) + d*x/c - e*log(a**(sympy.S(1)/3) + c**(sympy.S(1)/3)*x**2)/(6*a**(sympy.S(1)/3)*c**(sympy.S(2)/3)) + (sqrt(a)*d - sqrt(3)*sqrt(c)*e)*atan(sqrt(3) - 2*c**(sympy.S(1)/6)*x/a**(sympy.S(1)/6))/(6*a**(sympy.S(1)/3)*c**(sympy.S(7)/6)) - (sqrt(a)*d + sqrt(3)*sqrt(c)*e)*atan(sqrt(3) + 2*c**(sympy.S(1)/6)*x/a**(sympy.S(1)/6))/(6*a**(sympy.S(1)/3)*c**(sympy.S(7)/6)) - (sqrt(3)*sqrt(a)*d - sqrt(c)*e)*log(sqrt(3)*a**(sympy.S(1)/6)*c**(sympy.S(1)/6)*x + a**(sympy.S(1)/3) + c**(sympy.S(1)/3)*x**2)/(12*a**(sympy.S(1)/3)*c**(sympy.S(7)/6)) + (sqrt(3)*sqrt(a)*d + sqrt(c)*e)*log(-sqrt(3)*a**(sympy.S(1)/6)*c**(sympy.S(1)/6)*x + a**(sympy.S(1)/3) + c**(sympy.S(1)/3)*x**2)/(12*a**(sympy.S(1)/3)*c**(sympy.S(7)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_39():
    f = (d + e/x**3)/(a/x**6 + b/x**3 + c)
    F = d*x/c - 2**(sympy.S(2)/3)*(b*d - c*e + (-2*a*c*d + b**2*d - b*c*e)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(4)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(b*d - c*e + (-2*a*c*d + b**2*d - b*c*e)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(4)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*sqrt(3)*(b*d - c*e + (-2*a*c*d + b**2*d - b*c*e)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(4)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(b*d - c*e - (-2*a*c*d + b**2*d - b*c*e)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(4)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(b*d - c*e - (-2*a*c*d + b**2*d - b*c*e)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(4)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*sqrt(3)*(b*d - c*e - (-2*a*c*d + b**2*d - b*c*e)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(4)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_40():
    f = (d + e/x**4)/(a/x**8 + c)
    F = d*x/c - (sqrt(a)*(-sqrt(2)*d + d) + sqrt(c)*e)*log(-a**(sympy.S(1)/8)*c**(sympy.S(1)/8)*x*sqrt(2 - sqrt(2)) + a**(sympy.S(1)/4) + c**(sympy.S(1)/4)*x**2)/(8*a**(sympy.S(3)/8)*c**(sympy.S(9)/8)*sqrt(4 - 2*sqrt(2))) + (sqrt(a)*(-sqrt(2)*d + d) + sqrt(c)*e)*log(a**(sympy.S(1)/8)*c**(sympy.S(1)/8)*x*sqrt(2 - sqrt(2)) + a**(sympy.S(1)/4) + c**(sympy.S(1)/4)*x**2)/(8*a**(sympy.S(3)/8)*c**(sympy.S(9)/8)*sqrt(4 - 2*sqrt(2))) - sqrt(sqrt(2) + 2)*(sqrt(a)*(-sqrt(2)*d + d) + sqrt(c)*e)*atan((a**(sympy.S(1)/8)*sqrt(sqrt(2) + 2) - 2*c**(sympy.S(1)/8)*x)/(a**(sympy.S(1)/8)*sqrt(2 - sqrt(2))))/(8*a**(sympy.S(3)/8)*c**(sympy.S(9)/8)) + sqrt(sqrt(2) + 2)*(sqrt(a)*(-sqrt(2)*d + d) + sqrt(c)*e)*atan((a**(sympy.S(1)/8)*sqrt(sqrt(2) + 2) + 2*c**(sympy.S(1)/8)*x)/(a**(sympy.S(1)/8)*sqrt(2 - sqrt(2))))/(8*a**(sympy.S(3)/8)*c**(sympy.S(9)/8)) + (sqrt(a)*d*(1 + sqrt(2)) + sqrt(c)*e)*log(-a**(sympy.S(1)/8)*c**(sympy.S(1)/8)*x*sqrt(sqrt(2) + 2) + a**(sympy.S(1)/4) + c**(sympy.S(1)/4)*x**2)/(8*a**(sympy.S(3)/8)*c**(sympy.S(9)/8)*sqrt(2*sqrt(2) + 4)) - (sqrt(a)*d*(1 + sqrt(2)) + sqrt(c)*e)*log(a**(sympy.S(1)/8)*c**(sympy.S(1)/8)*x*sqrt(sqrt(2) + 2) + a**(sympy.S(1)/4) + c**(sympy.S(1)/4)*x**2)/(8*a**(sympy.S(3)/8)*c**(sympy.S(9)/8)*sqrt(2*sqrt(2) + 4)) + sqrt(2 - sqrt(2))*(sqrt(a)*d*(1 + sqrt(2)) + sqrt(c)*e)*atan((a**(sympy.S(1)/8)*sqrt(2 - sqrt(2)) - 2*c**(sympy.S(1)/8)*x)/(a**(sympy.S(1)/8)*sqrt(sqrt(2) + 2)))/(8*a**(sympy.S(3)/8)*c**(sympy.S(9)/8)) - sqrt(2 - sqrt(2))*(sqrt(a)*d*(1 + sqrt(2)) + sqrt(c)*e)*atan((a**(sympy.S(1)/8)*sqrt(2 - sqrt(2)) + 2*c**(sympy.S(1)/8)*x)/(a**(sympy.S(1)/8)*sqrt(sqrt(2) + 2)))/(8*a**(sympy.S(3)/8)*c**(sympy.S(9)/8))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_41():
    f = (d + e/x**4)/(a/x**8 + b/x**4 + c)
    F = d*x/c + 2**(sympy.S(3)/4)*(b*d - c*e - (-2*a*c*d + b**2*d - b*c*e)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(5)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) + 2**(sympy.S(3)/4)*(b*d - c*e - (-2*a*c*d + b**2*d - b*c*e)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(5)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) + 2**(sympy.S(3)/4)*(b*d - c*e + (-2*a*c*d + b**2*d - b*c*e)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(5)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) + 2**(sympy.S(3)/4)*(b*d - c*e + (-2*a*c*d + b**2*d - b*c*e)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(5)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_42():
    f = (d + e*x**n)**3/(a + c*x**(2*n))
    F = 3*d*e**2*x/c + e**3*x**(n + 1)/(c*(n + 1)) + d*x*(-3*a*e**2 + c*d**2)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(a*c) + e*x**(n + 1)*(-a*e**2 + 3*c*d**2)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(a*c*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_43():
    f = (d + e*x**n)**2/(a + c*x**(2*n))
    F = e**2*x/c + 2*d*e*x**(n + 1)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(a*(n + 1)) + x*(-a*e**2 + c*d**2)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(a*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_44():
    f = (d + e*x**n)/(a + c*x**(2*n))
    F = d*x*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/a + e*x**(n + 1)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(a*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_45():
    f = 1/((a + c*x**(2*n))*(d + e*x**n))
    F = e**2*x*hyper((1, 1/n), (1 + 1/n,), -e*x**n/d)/(d*(a*e**2 + c*d**2)) + c*d*x*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(a*(a*e**2 + c*d**2)) - c*e*x**(n + 1)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(a*(n + 1)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_46():
    f = 1/((a + c*x**(2*n))*(d + e*x**n)**2)
    F = 2*c*e**2*x*hyper((1, 1/n), (1 + 1/n,), -e*x**n/d)/(a*e**2 + c*d**2)**2 + e**2*x*hyper((2, 1/n), (1 + 1/n,), -e*x**n/d)/(d**2*(a*e**2 + c*d**2)) - 2*c**2*d*e*x**(n + 1)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(a*(n + 1)*(a*e**2 + c*d**2)**2) + c*x*(-a*e**2 + c*d**2)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(a*(a*e**2 + c*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_47():
    f = (d + e*x**n)/(a - c*x**(2*n))
    F = d*x*hyper((1, 1/(2*n)), (1 + 1/(2*n),), c*x**(2*n)/a)/a + e*x**(n + 1)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), c*x**(2*n)/a)/(a*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_48():
    f = (d + e*x**n)**3/(a + c*x**(2*n))**2
    F = 3*d*e**2*x*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(a*c) + e**3*x**(n + 1)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(a*c*(n + 1)) + x*(d*(-3*a*e**2 + c*d**2) + e*x**n*(-a*e**2 + 3*c*d**2))/(2*a*c*n*(a + c*x**(2*n))) - d*x*(1 - 2*n)*(-3*a*e**2 + c*d**2)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(2*a**2*c*n) - e*x**(n + 1)*(1 - n)*(-a*e**2 + 3*c*d**2)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(2*a**2*c*n*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_49():
    f = (d + e*x**n)**2/(a + c*x**(2*n))**2
    F = e**2*x*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(a*c) + x*(-a*e**2 + c*d**2 + 2*c*d*e*x**n)/(2*a*c*n*(a + c*x**(2*n))) - d*e*x**(n + 1)*(1 - n)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(a**2*n*(n + 1)) - x*(1 - 2*n)*(-a*e**2 + c*d**2)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(2*a**2*c*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_50():
    f = (d + e*x**n)/(a + c*x**(2*n))**2
    F = x*(d + e*x**n)/(2*a*n*(a + c*x**(2*n))) - d*x*(1 - 2*n)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(2*a**2*n) - e*x**(n + 1)*(1 - n)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(2*a**2*n*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_51():
    f = 1/((a + c*x**(2*n))**2*(d + e*x**n))
    F = e**4*x*hyper((1, 1/n), (1 + 1/n,), -e*x**n/d)/(d*(a*e**2 + c*d**2)**2) + c*d*e**2*x*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(a*(a*e**2 + c*d**2)**2) - c*e**3*x**(n + 1)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(a*(n + 1)*(a*e**2 + c*d**2)**2) + c*x*(d - e*x**n)/(2*a*n*(a + c*x**(2*n))*(a*e**2 + c*d**2)) - c*d*x*(1 - 2*n)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(2*a**2*n*(a*e**2 + c*d**2)) + c*e*x**(n + 1)*(1 - n)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(2*a**2*n*(n + 1)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_52():
    f = 1/((a + c*x**(2*n))**2*(d + e*x**n)**2)
    F = 4*c*e**4*x*hyper((1, 1/n), (1 + 1/n,), -e*x**n/d)/(a*e**2 + c*d**2)**3 + e**4*x*hyper((2, 1/n), (1 + 1/n,), -e*x**n/d)/(d**2*(a*e**2 + c*d**2)**2) - 4*c**2*d*e**3*x**(n + 1)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(a*(n + 1)*(a*e**2 + c*d**2)**3) + c*e**2*x*(-a*e**2 + 3*c*d**2)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(a*(a*e**2 + c*d**2)**3) + c*x*(-a*e**2 + c*d**2 - 2*c*d*e*x**n)/(2*a*n*(a + c*x**(2*n))*(a*e**2 + c*d**2)**2) + c**2*d*e*x**(n + 1)*(1 - n)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(a**2*n*(n + 1)*(a*e**2 + c*d**2)**2) - c*x*(1 - 2*n)*(-a*e**2 + c*d**2)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(2*a**2*n*(a*e**2 + c*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_53():
    f = (d + e*x**n)**3/(a + c*x**(2*n))**3
    F = e**2*x*(3*d + e*x**n)/(2*a*c*n*(a + c*x**(2*n))) + x*(d*(-3*a*e**2 + c*d**2) + e*x**n*(-a*e**2 + 3*c*d**2))/(4*a*c*n*(a + c*x**(2*n))**2) - 3*d*e**2*x*(1 - 2*n)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(2*a**2*c*n) - e**3*x**(n + 1)*(1 - n)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(2*a**2*c*n*(n + 1)) - x*(d*(1 - 4*n)*(-3*a*e**2 + c*d**2) + e*x**n*(1 - 3*n)*(-a*e**2 + 3*c*d**2))/(8*a**2*c*n**2*(a + c*x**(2*n))) + d*x*(1 - 4*n)*(1 - 2*n)*(-3*a*e**2 + c*d**2)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(8*a**3*c*n**2) + e*x**(n + 1)*(1 - 3*n)*(1 - n)*(-a*e**2 + 3*c*d**2)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(8*a**3*c*n**2*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_54():
    f = (d + e*x**n)**2/(a + c*x**(2*n))**3
    F = x*(-a*e**2 + c*d**2 + 2*c*d*e*x**n)/(4*a*c*n*(a + c*x**(2*n))**2) + e**2*x*hyper((2, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(a**2*c) - x*(2*c*d*e*x**n*(1 - 3*n) + (1 - 4*n)*(-a*e**2 + c*d**2))/(8*a**2*c*n**2*(a + c*x**(2*n))) + d*e*x**(n + 1)*(1 - 3*n)*(1 - n)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(4*a**3*n**2*(n + 1)) + x*(1 - 4*n)*(1 - 2*n)*(-a*e**2 + c*d**2)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(8*a**3*c*n**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_55():
    f = (d + e*x**n)/(a + c*x**(2*n))**3
    F = x*(d + e*x**n)/(4*a*n*(a + c*x**(2*n))**2) - x*(d*(1 - 4*n) + e*x**n*(1 - 3*n))/(8*a**2*n**2*(a + c*x**(2*n))) + d*x*(1 - 4*n)*(1 - 2*n)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(8*a**3*n**2) + e*x**(n + 1)*(1 - 3*n)*(1 - n)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(8*a**3*n**2*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_56():
    f = 1/((a + c*x**(2*n))**3*(d + e*x**n))
    F = e**6*x*hyper((1, 1/n), (1 + 1/n,), -e*x**n/d)/(d*(a*e**2 + c*d**2)**3) + c*d*e**4*x*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(a*(a*e**2 + c*d**2)**3) - c*e**5*x**(n + 1)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(a*(n + 1)*(a*e**2 + c*d**2)**3) + c*e**2*x*(d - e*x**n)/(2*a*n*(a + c*x**(2*n))*(a*e**2 + c*d**2)**2) + c*x*(d - e*x**n)/(4*a*n*(a + c*x**(2*n))**2*(a*e**2 + c*d**2)) - c*d*e**2*x*(1 - 2*n)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(2*a**2*n*(a*e**2 + c*d**2)**2) + c*e**3*x**(n + 1)*(1 - n)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(2*a**2*n*(n + 1)*(a*e**2 + c*d**2)**2) - c*x*(d*(1 - 4*n) - e*x**n*(1 - 3*n))/(8*a**2*n**2*(a + c*x**(2*n))*(a*e**2 + c*d**2)) + c*d*x*(1 - 4*n)*(1 - 2*n)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(8*a**3*n**2*(a*e**2 + c*d**2)) - c*e*x**(n + 1)*(1 - 3*n)*(1 - n)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(8*a**3*n**2*(n + 1)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_57():
    f = 1/((a + c*x**(2*n))**3*(d + e*x**n)**2)
    F = 6*c*e**6*x*hyper((1, 1/n), (1 + 1/n,), -e*x**n/d)/(a*e**2 + c*d**2)**4 + e**6*x*hyper((2, 1/n), (1 + 1/n,), -e*x**n/d)/(d**2*(a*e**2 + c*d**2)**3) - 6*c**2*d*e**5*x**(n + 1)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(a*(n + 1)*(a*e**2 + c*d**2)**4) + c*e**4*x*(-a*e**2 + 5*c*d**2)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(a*(a*e**2 + c*d**2)**4) + c*e**2*x*(-a*e**2 + 3*c*d**2 - 4*c*d*e*x**n)/(2*a*n*(a + c*x**(2*n))*(a*e**2 + c*d**2)**3) + c*x*(-a*e**2 + c*d**2 - 2*c*d*e*x**n)/(4*a*n*(a + c*x**(2*n))**2*(a*e**2 + c*d**2)**2) + 2*c**2*d*e**3*x**(n + 1)*(1 - n)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(a**2*n*(n + 1)*(a*e**2 + c*d**2)**3) - c*e**2*x*(1 - 2*n)*(-a*e**2 + 3*c*d**2)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(2*a**2*n*(a*e**2 + c*d**2)**3) - c*x*(-2*c*d*e*x**n*(1 - 3*n) + (1 - 4*n)*(-a*e**2 + c*d**2))/(8*a**2*n**2*(a + c*x**(2*n))*(a*e**2 + c*d**2)**2) - c**2*d*e*x**(n + 1)*(1 - 3*n)*(1 - n)*hyper((1, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/(4*a**3*n**2*(n + 1)*(a*e**2 + c*d**2)**2) + c*x*(1 - 4*n)*(1 - 2*n)*(-a*e**2 + c*d**2)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(8*a**3*n**2*(a*e**2 + c*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_58():
    f = 1/(sqrt(a + c*x**(2*n))*(d + e*x**n))
    F = x*sqrt(1 + c*x**(2*n)/a)*appellf1(1/(2*n), sympy.S.Half, 1, 1 + 1/(2*n), -c*x**(2*n)/a, e**2*x**(2*n)/d**2)/(d*sqrt(a + c*x**(2*n))) - e*x**(n + 1)*sqrt(1 + c*x**(2*n)/a)*appellf1((n + 1)/(2*n), sympy.S.Half, 1, sympy.S(3)/2 + 1/(2*n), -c*x**(2*n)/a, e**2*x**(2*n)/d**2)/(d**2*sqrt(a + c*x**(2*n))*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_59():
    f = (a + c*x**(2*n))**p*(d + e*x**n)**q
    F = sympy.Function('Unintegrable')((((Symbol('d') + (Symbol('e') * (x)**(Symbol('n')))))**(Symbol('q')) * ((Symbol('a') + (Symbol('c') * (x)**((Integer(2) * Symbol('n'))))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_60():
    f = (a + c*x**(2*n))**p*(d + e*x**n)**3
    F = d**3*x*(a + c*x**(2*n))**p*hyper((-p, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(1 + c*x**(2*n)/a)**p + 3*d**2*e*x**(n + 1)*(a + c*x**(2*n))**p*hyper((-p, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/((1 + c*x**(2*n)/a)**p*(n + 1)) + 3*d*e**2*x**(2*n + 1)*(a + c*x**(2*n))**p*hyper((-p, 1 + 1/(2*n)), (2 + 1/(2*n),), -c*x**(2*n)/a)/((1 + c*x**(2*n)/a)**p*(2*n + 1)) + e**3*x**(3*n + 1)*(a + c*x**(2*n))**p*hyper((-p, sympy.S(3)/2 + 1/(2*n)), (sympy.S(5)/2 + 1/(2*n),), -c*x**(2*n)/a)/((1 + c*x**(2*n)/a)**p*(3*n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_61():
    f = (a + c*x**(2*n))**p*(d + e*x**n)**2
    F = d**2*x*(a + c*x**(2*n))**p*hyper((-p, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(1 + c*x**(2*n)/a)**p + 2*d*e*x**(n + 1)*(a + c*x**(2*n))**p*hyper((-p, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/((1 + c*x**(2*n)/a)**p*(n + 1)) + e**2*x**(2*n + 1)*(a + c*x**(2*n))**p*hyper((-p, 1 + 1/(2*n)), (2 + 1/(2*n),), -c*x**(2*n)/a)/((1 + c*x**(2*n)/a)**p*(2*n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_62():
    f = (a + c*x**(2*n))**p*(d + e*x**n)
    F = d*x*(a + c*x**(2*n))**p*hyper((-p, 1/(2*n)), (1 + 1/(2*n),), -c*x**(2*n)/a)/(1 + c*x**(2*n)/a)**p + e*x**(n + 1)*(a + c*x**(2*n))**p*hyper((-p, (n + 1)/(2*n)), (sympy.S(3)/2 + 1/(2*n),), -c*x**(2*n)/a)/((1 + c*x**(2*n)/a)**p*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_63():
    f = (a + c*x**(2*n))**p/(d + e*x**n)
    F = x*(a + c*x**(2*n))**p*appellf1(1/(2*n), 1, -p, 1 + 1/(2*n), e**2*x**(2*n)/d**2, -c*x**(2*n)/a)/(d*(1 + c*x**(2*n)/a)**p) - e*x**(n + 1)*(a + c*x**(2*n))**p*appellf1((n + 1)/(2*n), 1, -p, sympy.S(3)/2 + 1/(2*n), e**2*x**(2*n)/d**2, -c*x**(2*n)/a)/(d**2*(1 + c*x**(2*n)/a)**p*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_64():
    f = (a + c*x**(2*n))**p/(d + e*x**n)**2
    F = x*(a + c*x**(2*n))**p*appellf1(1/(2*n), 2, -p, 1 + 1/(2*n), e**2*x**(2*n)/d**2, -c*x**(2*n)/a)/(d**2*(1 + c*x**(2*n)/a)**p) - 2*e*x**(n + 1)*(a + c*x**(2*n))**p*appellf1((n + 1)/(2*n), 2, -p, sympy.S(3)/2 + 1/(2*n), e**2*x**(2*n)/d**2, -c*x**(2*n)/a)/(d**3*(1 + c*x**(2*n)/a)**p*(n + 1)) + e**2*x**(2*n + 1)*(a + c*x**(2*n))**p*appellf1(1 + 1/(2*n), 2, -p, 2 + 1/(2*n), e**2*x**(2*n)/d**2, -c*x**(2*n)/a)/(d**4*(1 + c*x**(2*n)/a)**p*(2*n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_65():
    f = (a + c*x**(2*n))**p/(d + e*x**n)**3
    F = x*(a + c*x**(2*n))**p*appellf1(1/(2*n), 3, -p, 1 + 1/(2*n), e**2*x**(2*n)/d**2, -c*x**(2*n)/a)/(d**3*(1 + c*x**(2*n)/a)**p) - 3*e*x**(n + 1)*(a + c*x**(2*n))**p*appellf1((n + 1)/(2*n), 3, -p, sympy.S(3)/2 + 1/(2*n), e**2*x**(2*n)/d**2, -c*x**(2*n)/a)/(d**4*(1 + c*x**(2*n)/a)**p*(n + 1)) + 3*e**2*x**(2*n + 1)*(a + c*x**(2*n))**p*appellf1(1 + 1/(2*n), 3, -p, 2 + 1/(2*n), e**2*x**(2*n)/d**2, -c*x**(2*n)/a)/(d**5*(1 + c*x**(2*n)/a)**p*(2*n + 1)) - e**3*x**(3*n + 1)*(a + c*x**(2*n))**p*appellf1(sympy.S(3)/2 + 1/(2*n), 3, -p, sympy.S(5)/2 + 1/(2*n), e**2*x**(2*n)/d**2, -c*x**(2*n)/a)/(d**6*(1 + c*x**(2*n)/a)**p*(3*n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_66():
    f = (d + e*x**n)*(a + b*x**n + c*x**(2*n))
    F = a*d*x + c*e*x**(3*n + 1)/(3*n + 1) + x**(n + 1)*(a*e + b*d)/(n + 1) + x**(2*n + 1)*(b*e + c*d)/(2*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_67():
    f = (d + e*x**n)*(a + b*x**n + c*x**(2*n))**2
    F = a**2*d*x + a*x**(n + 1)*(a*e + 2*b*d)/(n + 1) + c**2*e*x**(5*n + 1)/(5*n + 1) + c*x**(4*n + 1)*(2*b*e + c*d)/(4*n + 1) + x**(2*n + 1)*(2*a*b*e + 2*a*c*d + b**2*d)/(2*n + 1) + x**(3*n + 1)*(2*a*c*e + b**2*e + 2*b*c*d)/(3*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_68():
    f = (d + e*x**n)*(a + b*x**n + c*x**(2*n))**3
    F = a**3*d*x + a**2*x**(n + 1)*(a*e + 3*b*d)/(n + 1) + 3*a*x**(2*n + 1)*(a*b*e + a*c*d + b**2*d)/(2*n + 1) + c**3*e*x**(7*n + 1)/(7*n + 1) + c**2*x**(6*n + 1)*(3*b*e + c*d)/(6*n + 1) + 3*c*x**(5*n + 1)*(a*c*e + b**2*e + b*c*d)/(5*n + 1) + x**(3*n + 1)*(3*a**2*c*e + 3*a*b**2*e + 6*a*b*c*d + b**3*d)/(3*n + 1) + x**(4*n + 1)*(6*a*b*c*e + 3*a*c**2*d + b**3*e + 3*b**2*c*d)/(4*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_69():
    f = (d + e*x**n)**3/(a + b*x**n + c*x**(2*n))
    F = e**3*x**(n + 1)/(c*(n + 1)) + e**2*x*(-b*e + 3*c*d)/c**2 + x*(-a*c*e**3 + b**2*e**3 - 3*b*c*d*e**2 + 3*c**2*d**2*e - (-b*e + 2*c*d)*(b**2*e**2 + c**2*d**2 - c*e*(3*a*e + b*d))/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(c**2*(b + sqrt(-4*a*c + b**2))) + x*(-a*c*e**3 + b**2*e**3 - 3*b*c*d*e**2 + 3*c**2*d**2*e + (-b*e + 2*c*d)*(b**2*e**2 + c**2*d**2 - c*e*(3*a*e + b*d))/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(c**2*(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_70():
    f = (d + e*x**n)**2/(a + b*x**n + c*x**(2*n))
    F = e**2*x/c + x*(-b*e**2 + 2*c*d*e - (b**2*e**2 + 2*c**2*d**2 - 2*c*e*(a*e + b*d))/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(c*(b + sqrt(-4*a*c + b**2))) + x*(-b*e**2 + 2*c*d*e + (b**2*e**2 + 2*c**2*d**2 - 2*c*e*(a*e + b*d))/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(c*(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_71():
    f = (d + e*x**n)/(a + b*x**n + c*x**(2*n))
    F = x*(e - (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(b + sqrt(-4*a*c + b**2)) + x*(e + (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(b - sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_72():
    f = 1/((d + e*x**n)*(a + b*x**n + c*x**(2*n)))
    F = -c*x*(2*c*d - e*(b + sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/((-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)) - c*x*(e + (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((b + sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)) + e**2*x*hyper((1, 1/n), (1 + 1/n,), -e*x**n/d)/(d*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_73():
    f = 1/((d + e*x**n)**2*(a + b*x**n + c*x**(2*n)))
    F = -c*x*(b*e**2*(b - sqrt(-4*a*c + b**2)) + 2*c**2*d**2 - 2*c*e*(a*e + b*d - d*sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)**2) - c*x*(b*e**2*(b + sqrt(-4*a*c + b**2)) + 2*c**2*d**2 - 2*c*e*(a*e + b*d + d*sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/((-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)**2) + e**2*x*(-b*e + 2*c*d)*hyper((1, 1/n), (1 + 1/n,), -e*x**n/d)/(d*(a*e**2 - b*d*e + c*d**2)**2) + e**2*x*hyper((2, 1/n), (1 + 1/n,), -e*x**n/d)/(d**2*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_74():
    f = 1/((d + e*x**n)**3*(a + b*x**n + c*x**(2*n)))
    F = -c*x*(-b**2*e**3*(b - sqrt(-4*a*c + b**2)) + 2*c**3*d**3 - 3*c**2*d*e*(2*a*e + b*d - d*sqrt(-4*a*c + b**2)) + c*e**2*(3*a*b*e - a*e*sqrt(-4*a*c + b**2) + 3*b**2*d - 3*b*d*sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)**3) - c*x*(-b**2*e**3*(b + sqrt(-4*a*c + b**2)) + 2*c**3*d**3 - 3*c**2*d*e*(2*a*e + b*d + d*sqrt(-4*a*c + b**2)) + c*e**2*(a*e*sqrt(-4*a*c + b**2) + 3*b**2*d + 3*b*(a*e + d*sqrt(-4*a*c + b**2))))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/((-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)**3) + e**2*x*(b**2*e**2 + 3*c**2*d**2 - c*e*(a*e + 3*b*d))*hyper((1, 1/n), (1 + 1/n,), -e*x**n/d)/(d*(a*e**2 - b*d*e + c*d**2)**3) + e**2*x*(-b*e + 2*c*d)*hyper((2, 1/n), (1 + 1/n,), -e*x**n/d)/(d**2*(a*e**2 - b*d*e + c*d**2)**2) + e**2*x*hyper((3, 1/n), (1 + 1/n,), -e*x**n/d)/(d**3*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_75():
    f = (d + e*x**n)**3/(a + b*x**n + c*x**(2*n))**2
    F = e**2*x*(e - (-3*b*e + 6*c*d)/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(c*(b + sqrt(-4*a*c + b**2))) + e**2*x*(e + (-3*b*e + 6*c*d)/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(c*(b - sqrt(-4*a*c + b**2))) + x*(-a*b*e*(a*e**2 + 3*c*d**2) - 2*a*c*d*(-3*a*e**2 + c*d**2) + b**2*c*d**3 - x**n*(a*b**2*e**3 + 2*a*c*e*(-a*e**2 + 3*c*d**2) - b*c*d*(3*a*e**2 + c*d**2)))/(a*c*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))) + x*((1 - n)*(a*b**2*e**3 + 2*a*c*e*(-a*e**2 + 3*c*d**2) - b*c*d*(3*a*e**2 + c*d**2)) - (-a*b**3*e**3*(1 - 3*n) + 2*a*b*c*e*(a*e**2*(2 - 5*n) + 3*c*d**2*n) + 4*a*c**2*d*(1 - 2*n)*(-3*a*e**2 + c*d**2) + b**2*c*d*(3*a*e**2*(1 - 3*n) - c*d**2*(1 - n)))/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*c*n*(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + x*((1 - n)*(a*b**2*e**3 + 2*a*c*e*(-a*e**2 + 3*c*d**2) - b*c*d*(3*a*e**2 + c*d**2)) + (-a*b**3*e**3*(1 - 3*n) + 2*a*b*c*e*(a*e**2*(2 - 5*n) + 3*c*d**2*n) + 4*a*c**2*d*(1 - 2*n)*(-3*a*e**2 + c*d**2) + b**2*c*d*(3*a*e**2*(1 - 3*n) - c*d**2*(1 - n)))/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*c*n*(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_76():
    f = (d + e*x**n)**2/(a + b*x**n + c*x**(2*n))**2
    F = -2*e**2*x*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2)) - 2*e**2*x*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2)) + x*(-2*a*b*d*e - 2*a*(-a*e**2 + c*d**2) + b**2*d**2 + x**n*(a*b*e**2 - 4*a*c*d*e + b*c*d**2))/(a*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))) - x*((1 - n)*(a*b*e**2 - 4*a*c*d*e + b*c*d**2) + (4*a*b*c*d*e*n + 4*a*c*(1 - 2*n)*(-a*e**2 + c*d**2) + b**2*(a*e**2*(1 - 3*n) - c*d**2*(1 - n)))/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*n*(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) - x*((1 - n)*(a*b*e**2 - 4*a*c*d*e + b*c*d**2) - (4*a*b*c*d*e*n + 4*a*c*(1 - 2*n)*(-a*e**2 + c*d**2) + b**2*(a*e**2*(1 - 3*n) - c*d**2*(1 - n)))/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*n*(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_77():
    f = (d + e*x**n)/(a + b*x**n + c*x**(2*n))**2
    F = -c*x*(2*a*(c*d*(2 - 4*n) - e*(1 - n)*sqrt(-4*a*c + b**2)) - b**2*d*(1 - n) + b*(2*a*e*n + d*(1 - n)*sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))) - c*x*(2*a*(2*c*d*(1 - 2*n) + e*(1 - n)*sqrt(-4*a*c + b**2)) - b**2*(-d*n + d) - b*(-2*a*e*n + d*(1 - n)*sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))) + x*(-a*b*e - 2*a*c*d + b**2*d + c*x**n*(-2*a*e + b*d))/(a*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_78():
    f = 1/((d + e*x**n)*(a + b*x**n + c*x**(2*n))**2)
    F = -c*e**2*x*(2*c*d - e*(b - sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)**2) - c*e**2*x*(2*c*d - e*(b + sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/((-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)**2) + e**4*x*hyper((1, 1/n), (1 + 1/n,), -e*x**n/d)/(d*(a*e**2 - b*d*e + c*d**2)**2) + c*x*(-2*a*c*(2*c*d*(1 - 2*n) + e*(1 - n)*sqrt(-4*a*c + b**2)) - b**3*e*(1 - n) + b**2*(1 - n)*(c*d + e*sqrt(-4*a*c + b**2)) + b*c*(2*a*e*(2 - 3*n) - d*(1 - n)*sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)) - c*x*((1 - n)*(2*a*c*e - b**2*e + b*c*d) + (2*a*b*c*e*(2 - 3*n) - 4*a*c**2*d*(1 - 2*n) - b**3*e*(1 - n) + b**2*c*d*(1 - n))/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*n*(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2)) + x*(3*a*b*c*e - 2*a*c**2*d - b**3*e + b**2*c*d + c*x**n*(2*a*c*e - b**2*e + b*c*d))/(a*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_79():
    f = 1/((d + e*x**n)**2*(a + b*x**n + c*x**(2*n))**2)
    F = -2*c*e**2*x*(b*e**2*(b - sqrt(-4*a*c + b**2)) + 3*c**2*d**2 - c*e*(a*e + 3*b*d - 2*d*sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)**3) - 2*c*e**2*x*(b*e**2*(b + sqrt(-4*a*c + b**2)) + 3*c**2*d**2 - c*e*(a*e + 3*b*d + 2*d*sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/((-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)**3) + 2*e**4*x*(-b*e + 2*c*d)*hyper((1, 1/n), (1 + 1/n,), -e*x**n/d)/(d*(a*e**2 - b*d*e + c*d**2)**3) + e**4*x*hyper((2, 1/n), (1 + 1/n,), -e*x**n/d)/(d**2*(a*e**2 - b*d*e + c*d**2)**2) + c*x*(4*a*c**2*(-c*d**2*(1 - 2*n) + e*(a*e*(1 - 2*n) - d*(1 - n)*sqrt(-4*a*c + b**2))) + b**4*e**2*(1 - n) - b**3*e*(1 - n)*(2*c*d + e*sqrt(-4*a*c + b**2)) - b**2*c*(-c*d**2*(1 - n) + e*(a*e*(5 - 7*n) - 2*d*(1 - n)*sqrt(-4*a*c + b**2))) + b*c*(3*a*e**2*(1 - n)*sqrt(-4*a*c + b**2) + c*d*(4*a*e*(2 - 3*n) - d*(1 - n)*sqrt(-4*a*c + b**2))))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)**2) + c*x*(4*a*c**2*(-c*d**2*(1 - 2*n) + e*(a*e*(1 - 2*n) + d*(1 - n)*sqrt(-4*a*c + b**2))) + b**4*e**2*(1 - n) - b**3*e*(1 - n)*(2*c*d - e*sqrt(-4*a*c + b**2)) - b**2*c*(-c*d**2*(1 - n) + e*(a*e*(5 - 7*n) + 2*d*(1 - n)*sqrt(-4*a*c + b**2))) + b*c*(-3*a*e**2*(1 - n)*sqrt(-4*a*c + b**2) + c*d*(4*a*e*(2 - 3*n) + d*(1 - n)*sqrt(-4*a*c + b**2))))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)**2) - x*(-6*a*b*c**2*d*e + 2*a*c**2*(-a*e**2 + c*d**2) - b**4*e**2 + 2*b**3*c*d*e - b**2*c*(-4*a*e**2 + c*d**2) + c*x**n*(-4*a*c**2*d*e - b**3*e**2 + 2*b**2*c*d*e - b*c*(-3*a*e**2 + c*d**2)))/(a*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))*(a*e**2 - b*d*e + c*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_80():
    f = (d + e*x**n)**3/(a + b*x**n + c*x**(2*n))**3
    F = e**2*x*(-2*a*c*(6*c*d*(1 - 2*n) - e*(1 - n)*sqrt(-4*a*c + b**2)) - b**3*e*(1 - n) + b**2*(1 - n)*(3*c*d + e*sqrt(-4*a*c + b**2)) + b*c*(2*a*e*(2 - 5*n) - 3*d*(1 - n)*sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*c*n*(-4*a*c + b**2)*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))) + e**2*x*(-2*a*c*(6*c*d*(1 - 2*n) + e*(1 - n)*sqrt(-4*a*c + b**2)) - b**3*e*(1 - n) + b**2*(1 - n)*(3*c*d - e*sqrt(-4*a*c + b**2)) + b*c*(2*a*e*(2 - 5*n) + 3*d*(1 - n)*sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*c*n*(-4*a*c + b**2)*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))) + x*(-a*b*e*(a*e**2 + 3*c*d**2) - 2*a*c*d*(-3*a*e**2 + c*d**2) + b**2*c*d**3 - x**n*(a*b**2*e**3 + 2*a*c*e*(-a*e**2 + 3*c*d**2) - b*c*d*(3*a*e**2 + c*d**2)))/(2*a*c*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))**2) + e**2*x*(a*b*c*e - 6*a*c**2*d - b**3*e + 3*b**2*c*d + c*x**n*(-2*a*c*e - b**2*e + 3*b*c*d))/(a*c**2*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))) + x*((1 - n)*(4*a**2*c**2*e*(1 - 3*n)*(-a*e**2 + 3*c*d**2) - 2*a*b**4*e**3*n - a*b**2*c*e*(-a*e**2*(2*n + 1) + 3*c*d**2) - 2*a*b*c**2*d*(3*a*e**2*n + c*d**2*(2 - 7*n)) + b**3*c*d*(6*a*e**2*n + c*d**2*(1 - 2*n))) + (-4*a**2*b*c**2*e*(a*e**2*(19*n**2 - 11*n + 1) + 3*c*d**2*(-3*n**2 - n + 1)) - 8*a**2*c**3*d*(-3*a*e**2 + c*d**2)*(8*n**2 - 6*n + 1) + 2*a*b**5*e**3*n*(1 - n) + a*b**3*c*e*(a*e**2*(30*n**2 - 19*n + 1) + 3*c*d**2*(1 - n)) + 6*a*b**2*c**2*d*(-a*e**2*(15*n**2 - 10*n + 1) + c*d**2*(3*n**2 - 4*n + 1)) - b**4*c*d*(1 - n)*(6*a*e**2*n + c*d**2*(1 - 2*n)))/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(2*a**2*c*n**2*(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + x*((1 - n)*(4*a**2*c**2*e*(1 - 3*n)*(-a*e**2 + 3*c*d**2) - 2*a*b**4*e**3*n - a*b**2*c*e*(-a*e**2*(2*n + 1) + 3*c*d**2) - 2*a*b*c**2*d*(3*a*e**2*n + c*d**2*(2 - 7*n)) + b**3*c*d*(6*a*e**2*n + c*d**2*(1 - 2*n))) - (-4*a**2*b*c**2*e*(a*e**2*(19*n**2 - 11*n + 1) + 3*c*d**2*(-3*n**2 - n + 1)) - 8*a**2*c**3*d*(-3*a*e**2 + c*d**2)*(8*n**2 - 6*n + 1) + 2*a*b**5*e**3*n*(1 - n) + a*b**3*c*e*(a*e**2*(30*n**2 - 19*n + 1) + 3*c*d**2*(1 - n)) + 6*a*b**2*c**2*d*(-a*e**2*(15*n**2 - 10*n + 1) + c*d**2*(3*n**2 - 4*n + 1)) - b**4*c*d*(1 - n)*(6*a*e**2*n + c*d**2*(1 - 2*n)))/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(2*a**2*c*n**2*(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) - x*(2*a**2*b*c**2*e*(-5*a*e**2*n + 3*c*d**2*(2 - 3*n)) + 4*a**2*c**3*d*(1 - 4*n)*(-3*a*e**2 + c*d**2) - 2*a*b**5*e**3*n - 3*a*b**3*c*e*(-3*a*e**2*n + c*d**2) + a*b**2*c**2*d*(3*a*e**2*(1 - 9*n) - 5*c*d**2*(1 - 3*n)) + b**4*c*d*(6*a*e**2*n + c*d**2*(1 - 2*n)) + c*x**n*(4*a**2*c**2*e*(1 - 3*n)*(-a*e**2 + 3*c*d**2) - 2*a*b**4*e**3*n - a*b**2*c*e*(-a*e**2*(2*n + 1) + 3*c*d**2) - 2*a*b*c**2*d*(3*a*e**2*n + c*d**2*(2 - 7*n)) + b**3*c*d*(6*a*e**2*n + c*d**2*(1 - 2*n))))/(2*a**2*c**2*n**2*(-4*a*c + b**2)**2*(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_81():
    f = (d + e*x**n)**2/(a + b*x**n + c*x**(2*n))**3
    F = -e**2*x*(4*a*c*(1 - 2*n) - b**2*(1 - n) + b*(1 - n)*sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))) - e**2*x*(4*a*c*(1 - 2*n) - b**2*(1 - n) - b*(1 - n)*sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))) + x*(-2*a*b*d*e - 2*a*(-a*e**2 + c*d**2) + b**2*d**2 + x**n*(a*b*e**2 - 4*a*c*d*e + b*c*d**2))/(2*a*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))**2) + e**2*x*(-2*a*c + b**2 + b*c*x**n)/(a*c*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))) - x*((1 - n)*(-8*a**2*c**2*d*e*(1 - 3*n) + 2*a*b**2*c*d*e + 2*a*b*c*(a*e**2*n + c*d**2*(2 - 7*n)) - b**3*(2*a*e**2*n + c*d**2*(1 - 2*n))) - (-8*a**2*b*c**2*d*e*(-3*n**2 - n + 1) - 8*a**2*c**2*(-a*e**2 + c*d**2)*(8*n**2 - 6*n + 1) + 2*a*b**3*c*d*e*(1 - n) + 2*a*b**2*c*(-a*e**2*(15*n**2 - 10*n + 1) + 3*c*d**2*(3*n**2 - 4*n + 1)) - b**4*(1 - n)*(2*a*e**2*n + c*d**2*(1 - 2*n)))/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(2*a**2*n**2*(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) - x*((1 - n)*(-8*a**2*c**2*d*e*(1 - 3*n) + 2*a*b**2*c*d*e + 2*a*b*c*(a*e**2*n + c*d**2*(2 - 7*n)) - b**3*(2*a*e**2*n + c*d**2*(1 - 2*n))) + (-8*a**2*b*c**2*d*e*(-3*n**2 - n + 1) - 8*a**2*c**2*(-a*e**2 + c*d**2)*(8*n**2 - 6*n + 1) + 2*a*b**3*c*d*e*(1 - n) + 2*a*b**2*c*(-a*e**2*(15*n**2 - 10*n + 1) + 3*c*d**2*(3*n**2 - 4*n + 1)) - b**4*(1 - n)*(2*a*e**2*n + c*d**2*(1 - 2*n)))/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(2*a**2*n**2*(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + x*(-4*a**2*b*c**2*d*e*(2 - 3*n) - 4*a**2*c**2*(1 - 4*n)*(-a*e**2 + c*d**2) + 2*a*b**3*c*d*e - a*b**2*c*(a*e**2*(1 - 9*n) - 5*c*d**2*(1 - 3*n)) - b**4*(2*a*e**2*n + c*d**2*(1 - 2*n)) + c*x**n*(-8*a**2*c**2*d*e*(1 - 3*n) + 2*a*b**2*c*d*e + 2*a*b*c*(a*e**2*n + c*d**2*(2 - 7*n)) - b**3*(2*a*e**2*n + c*d**2*(1 - 2*n))))/(2*a**2*c*n**2*(-4*a*c + b**2)**2*(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_82():
    f = (d + e*x**n)/(a + b*x**n + c*x**(2*n))**3
    F = x*(-a*b*e - 2*a*c*d + b**2*d + c*x**n*(-2*a*e + b*d))/(2*a*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))**2) - c*x*(-4*a**2*c*(-2*c*d*(8*n**2 - 6*n + 1) + e*sqrt(-4*a*c + b**2)*(3*n**2 - 4*n + 1)) + a*b**2*(1 - n)*(-6*c*d*(1 - 3*n) + e*sqrt(-4*a*c + b**2)) + 2*a*b*c*(2*a*e*(-3*n**2 - n + 1) + d*sqrt(-4*a*c + b**2)*(7*n**2 - 9*n + 2)) + b**4*d*(2*n**2 - 3*n + 1) - b**3*(1 - n)*(a*e + d*(1 - 2*n)*sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(2*a**2*n**2*(-4*a*c + b**2)**2*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))) + c*x*(-4*a**2*c*(2*c*d*(8*n**2 - 6*n + 1) + e*sqrt(-4*a*c + b**2)*(3*n**2 - 4*n + 1)) + a*b**2*(1 - n)*(6*c*d*(1 - 3*n) + e*sqrt(-4*a*c + b**2)) - 2*a*b*c*(2*a*e*(-3*n**2 - n + 1) - d*sqrt(-4*a*c + b**2)*(7*n**2 - 9*n + 2)) - b**4*d*(2*n**2 - 3*n + 1) + b**3*(1 - n)*(a*e - d*(1 - 2*n)*sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(2*a**2*n**2*(-4*a*c + b**2)**2*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))) + x*(-2*a**2*b*c*e*(2 - 3*n) - 4*a**2*c**2*d*(1 - 4*n) + a*b**3*e + 5*a*b**2*c*d*(1 - 3*n) - b**4*d*(1 - 2*n) + c*x**n*(-4*a**2*c*e*(1 - 3*n) + a*b**2*e + 2*a*b*c*d*(2 - 7*n) - b**3*d*(1 - 2*n)))/(2*a**2*n**2*(-4*a*c + b**2)**2*(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_83():
    f = 1/((d + e*x**n)*(a + b*x**n + c*x**(2*n))**3)
    F = -c*e**4*x*(2*c*d - e*(b - sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)**3) - c*e**4*x*(2*c*d - e*(b + sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/((-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)**3) + e**6*x*hyper((1, 1/n), (1 + 1/n,), -e*x**n/d)/(d*(a*e**2 - b*d*e + c*d**2)**3) + c*e**2*x*(-2*a*c*(2*c*d*(1 - 2*n) + e*(1 - n)*sqrt(-4*a*c + b**2)) - b**3*e*(1 - n) + b**2*(1 - n)*(c*d + e*sqrt(-4*a*c + b**2)) + b*c*(2*a*e*(2 - 3*n) - d*(1 - n)*sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)**2) + c*e**2*x*(-2*a*c*(2*c*d*(1 - 2*n) - e*(1 - n)*sqrt(-4*a*c + b**2)) - b**3*e*(1 - n) + b**2*(1 - n)*(c*d - e*sqrt(-4*a*c + b**2)) + b*c*(2*a*e*(2 - 3*n) + d*(1 - n)*sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)**2) + e**2*x*(3*a*b*c*e - 2*a*c**2*d - b**3*e + b**2*c*d + c*x**n*(2*a*c*e - b**2*e + b*c*d))/(a*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))*(a*e**2 - b*d*e + c*d**2)**2) + x*(3*a*b*c*e - 2*a*c**2*d - b**3*e + b**2*c*d + c*x**n*(2*a*c*e - b**2*e + b*c*d))/(2*a*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))**2*(a*e**2 - b*d*e + c*d**2)) + c*x*(-4*a**2*c**2*(2*c*d*(8*n**2 - 6*n + 1) + e*sqrt(-4*a*c + b**2)*(3*n**2 - 4*n + 1)) + a*b**2*c*(1 - n)*(6*c*d*(1 - 3*n) + e*(5 - 14*n)*sqrt(-4*a*c + b**2)) - 2*a*b*c**2*(-2*a*e*(13*n**2 - 13*n + 3) + d*sqrt(-4*a*c + b**2)*(7*n**2 - 9*n + 2)) + b**5*e*(2*n**2 - 3*n + 1) - b**4*(c*d + e*sqrt(-4*a*c + b**2))*(2*n**2 - 3*n + 1) - b**3*c*(1 - n)*(a*e*(7 - 18*n) - d*(1 - 2*n)*sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(2*a**2*n**2*(-4*a*c + b**2)**2*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)) - c*x*(-4*a**2*c**2*(-2*c*d*(8*n**2 - 6*n + 1) + e*sqrt(-4*a*c + b**2)*(3*n**2 - 4*n + 1)) + a*b**2*c*(1 - n)*(-6*c*d*(1 - 3*n) + e*(5 - 14*n)*sqrt(-4*a*c + b**2)) - 2*a*b*c**2*(2*a*e*(13*n**2 - 13*n + 3) + d*sqrt(-4*a*c + b**2)*(7*n**2 - 9*n + 2)) - b**5*e*(2*n**2 - 3*n + 1) + b**4*(c*d - e*sqrt(-4*a*c + b**2))*(2*n**2 - 3*n + 1) + b**3*c*(1 - n)*(a*e*(7 - 18*n) + d*(1 - 2*n)*sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(2*a**2*n**2*(-4*a*c + b**2)**2*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)) + x*(2*a**2*b*c**2*e*(4 - 11*n) - 4*a**2*c**3*d*(1 - 4*n) - 3*a*b**3*c*e*(2 - 5*n) + 5*a*b**2*c**2*d*(1 - 3*n) + b**5*(-2*e*n + e) - b**4*c*d*(1 - 2*n) - c*x**n*(-4*a**2*c**2*e*(1 - 3*n) + a*b**2*c*e*(5 - 14*n) - 2*a*b*c**2*d*(2 - 7*n) - b**4*e*(1 - 2*n) + b**3*c*d*(1 - 2*n)))/(2*a**2*n**2*(-4*a*c + b**2)**2*(a + b*x**n + c*x**(2*n))*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_84():
    f = 1/((d + e*x**n)**2*(a + b*x**n + c*x**(2*n))**3)
    F = -c*e**4*x*(3*b*e**2*(b - sqrt(-4*a*c + b**2)) + 10*c**2*d**2 - 2*c*e*(a*e + 5*b*d - 3*d*sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)**4) - c*e**4*x*(3*b*e**2*(b + sqrt(-4*a*c + b**2)) + 10*c**2*d**2 - 2*c*e*(a*e + 5*b*d + 3*d*sqrt(-4*a*c + b**2)))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/((-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)**4) + 3*e**6*x*(-b*e + 2*c*d)*hyper((1, 1/n), (1 + 1/n,), -e*x**n/d)/(d*(a*e**2 - b*d*e + c*d**2)**4) + e**6*x*hyper((2, 1/n), (1 + 1/n,), -e*x**n/d)/(d**2*(a*e**2 - b*d*e + c*d**2)**3) + c*e**2*x*(4*a*c**2*(-3*c*d**2*(1 - 2*n) + e*(a*e*(1 - 2*n) - 2*d*(1 - n)*sqrt(-4*a*c + b**2))) + 2*b**4*e**2*(1 - n) - b**3*e*(1 - n)*(5*c*d + 2*e*sqrt(-4*a*c + b**2)) - b**2*c*(-3*c*d**2*(1 - n) + e*(a*e*(9 - 13*n) - 5*d*(1 - n)*sqrt(-4*a*c + b**2))) + b*c*(5*a*e**2*(1 - n)*sqrt(-4*a*c + b**2) + c*d*(4*a*e*(5 - 8*n) - 3*d*(1 - n)*sqrt(-4*a*c + b**2))))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)**3) + c*e**2*x*(4*a*c**2*(-3*c*d**2*(1 - 2*n) + e*(a*e*(1 - 2*n) + 2*d*(1 - n)*sqrt(-4*a*c + b**2))) + 2*b**4*e**2*(1 - n) - b**3*e*(1 - n)*(5*c*d - 2*e*sqrt(-4*a*c + b**2)) - b**2*c*(-3*c*d**2*(1 - n) + e*(a*e*(9 - 13*n) + 5*d*(1 - n)*sqrt(-4*a*c + b**2))) + b*c*(-5*a*e**2*(1 - n)*sqrt(-4*a*c + b**2) + c*d*(4*a*e*(5 - 8*n) + 3*d*(1 - n)*sqrt(-4*a*c + b**2))))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*n*(-4*a*c + b**2)*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)**3) - e**2*x*(-14*a*b*c**2*d*e + 2*a*c**2*(-a*e**2 + 3*c*d**2) - 2*b**4*e**2 + 5*b**3*c*d*e - b**2*c*(-7*a*e**2 + 3*c*d**2) + c*x**n*(-8*a*c**2*d*e - 2*b**3*e**2 + 5*b**2*c*d*e - b*c*(-5*a*e**2 + 3*c*d**2)))/(a*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))*(a*e**2 - b*d*e + c*d**2)**3) - x*(-6*a*b*c**2*d*e + 2*a*c**2*(-a*e**2 + c*d**2) - b**4*e**2 + 2*b**3*c*d*e - b**2*c*(-4*a*e**2 + c*d**2) + c*x**n*(-4*a*c**2*d*e - b**3*e**2 + 2*b**2*c*d*e - b*c*(-3*a*e**2 + c*d**2)))/(2*a*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))**2*(a*e**2 - b*d*e + c*d**2)**2) + c*x*((1 - n)*(-8*a**2*c**3*d*e*(1 - 3*n) + 2*a*b**2*c**2*d*e*(5 - 14*n) + 2*a*b*c**2*(a*e**2*(4 - 13*n) - c*d**2*(2 - 7*n)) + b**5*e**2*(1 - 2*n) - 2*b**4*c*d*e*(1 - 2*n) - b**3*c*(2*a*e**2*(3 - 8*n) - c*d**2*(1 - 2*n))) + (8*a**2*b*c**3*d*e*(13*n**2 - 13*n + 3) - 8*a**2*c**3*(-a*e**2 + c*d**2)*(8*n**2 - 6*n + 1) - 2*a*b**3*c**2*d*e*(18*n**2 - 25*n + 7) + 2*a*b**2*c**2*(-a*e**2*(35*n**2 - 38*n + 9) + 3*c*d**2*(3*n**2 - 4*n + 1)) - b**6*e**2*(2*n**2 - 3*n + 1) + 2*b**5*c*d*e*(2*n**2 - 3*n + 1) + b**4*c*(1 - n)*(4*a*e**2*(2 - 5*n) - c*d**2*(1 - 2*n)))/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(2*a**2*n**2*(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2*(a*e**2 - b*d*e + c*d**2)**2) + c*x*((1 - n)*(-8*a**2*c**3*d*e*(1 - 3*n) + 2*a*b**2*c**2*d*e*(5 - 14*n) + 2*a*b*c**2*(a*e**2*(4 - 13*n) - c*d**2*(2 - 7*n)) + b**5*e**2*(1 - 2*n) - 2*b**4*c*d*e*(1 - 2*n) - b**3*c*(2*a*e**2*(3 - 8*n) - c*d**2*(1 - 2*n))) - (8*a**2*b*c**3*d*e*(13*n**2 - 13*n + 3) - 8*a**2*c**3*(-a*e**2 + c*d**2)*(8*n**2 - 6*n + 1) - 2*a*b**3*c**2*d*e*(18*n**2 - 25*n + 7) + 2*a*b**2*c**2*(-a*e**2*(35*n**2 - 38*n + 9) + 3*c*d**2*(3*n**2 - 4*n + 1)) - b**6*e**2*(2*n**2 - 3*n + 1) + 2*b**5*c*d*e*(2*n**2 - 3*n + 1) + b**4*c*(1 - n)*(4*a*e**2*(2 - 5*n) - c*d**2*(1 - 2*n)))/sqrt(-4*a*c + b**2))*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(2*a**2*n**2*(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2*(a*e**2 - b*d*e + c*d**2)**2) - x*(-4*a**2*b*c**3*d*e*(4 - 11*n) + 4*a**2*c**3*(1 - 4*n)*(-a*e**2 + c*d**2) + 6*a*b**3*c**2*d*e*(2 - 5*n) + a*b**2*c**2*(a*e**2*(13 - 37*n) - 5*c*d**2*(1 - 3*n)) + b**6*e**2*(1 - 2*n) - 2*b**5*c*d*e*(1 - 2*n) - b**4*c*(a*e**2*(7 - 17*n) - c*d**2*(1 - 2*n)) + c*x**n*(-8*a**2*c**3*d*e*(1 - 3*n) + 2*a*b**2*c**2*d*e*(5 - 14*n) + 2*a*b*c**2*(a*e**2*(4 - 13*n) - c*d**2*(2 - 7*n)) + b**5*e**2*(1 - 2*n) - 2*b**4*c*d*e*(1 - 2*n) - b**3*c*(2*a*e**2*(3 - 8*n) - c*d**2*(1 - 2*n))))/(2*a**2*n**2*(-4*a*c + b**2)**2*(a + b*x**n + c*x**(2*n))*(a*e**2 - b*d*e + c*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_85():
    f = (d + e*x**n)*sqrt(a + b*x**n + c*x**(2*n))
    F = d*x*sqrt(a + b*x**n + c*x**(2*n))*appellf1(1/n, sympy.S(-1)/2, sympy.S(-1)/2, 1 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)) + e*x**(n + 1)*sqrt(a + b*x**n + c*x**(2*n))*appellf1(1 + 1/n, sympy.S(-1)/2, sympy.S(-1)/2, 2 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((n + 1)*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_86():
    f = (d + e*x**n)*(a + b*x**n + c*x**(2*n))**(sympy.S(3)/2)
    F = a*d*x*sqrt(a + b*x**n + c*x**(2*n))*appellf1(1/n, sympy.S(-3)/2, sympy.S(-3)/2, 1 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)) + a*e*x**(n + 1)*sqrt(a + b*x**n + c*x**(2*n))*appellf1(1 + 1/n, sympy.S(-3)/2, sympy.S(-3)/2, 2 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((n + 1)*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_87():
    f = (d + e*x**n)/sqrt(a + b*x**n + c*x**(2*n))
    F = d*x*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(1/n, sympy.S.Half, sympy.S.Half, 1 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/sqrt(a + b*x**n + c*x**(2*n)) + e*x**(n + 1)*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(1 + 1/n, sympy.S.Half, sympy.S.Half, 2 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((n + 1)*sqrt(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_88():
    f = (d + e*x**n)/(a + b*x**n + c*x**(2*n))**(sympy.S(3)/2)
    F = d*x*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(1/n, sympy.S(3)/2, sympy.S(3)/2, 1 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*sqrt(a + b*x**n + c*x**(2*n))) + e*x**(n + 1)*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(1 + 1/n, sympy.S(3)/2, sympy.S(3)/2, 2 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*(n + 1)*sqrt(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_89():
    f = (d + e*x**n)/(a + b*x**n + c*x**(2*n))**(sympy.S(5)/2)
    F = d*x*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(1/n, sympy.S(5)/2, sympy.S(5)/2, 1 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a**2*sqrt(a + b*x**n + c*x**(2*n))) + e*x**(n + 1)*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(1 + 1/n, sympy.S(5)/2, sympy.S(5)/2, 2 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a**2*(n + 1)*sqrt(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_90():
    f = (d + e*x**n)**q*(a + b*x**n + c*x**(2*n))**p
    F = sympy.Function('Unintegrable')((((Symbol('d') + (Symbol('e') * (x)**(Symbol('n')))))**(Symbol('q')) * ((Symbol('a') + (Symbol('b') * (x)**(Symbol('n'))) + (Symbol('c') * (x)**((Integer(2) * Symbol('n'))))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_91():
    f = (d + e*x**n)**3*(a + b*x**n + c*x**(2*n))**p
    F = d**3*x*(a + b*x**n + c*x**(2*n))**p*appellf1(1/n, -p, -p, 1 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)**p) + 3*d**2*e*x**(n + 1)*(a + b*x**n + c*x**(2*n))**p*appellf1(1 + 1/n, -p, -p, 2 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((n + 1)*(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)**p) + 3*d*e**2*x**(2*n + 1)*(a + b*x**n + c*x**(2*n))**p*appellf1(2 + 1/n, -p, -p, 3 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((2*n + 1)*(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)**p) + e**3*x**(3*n + 1)*(a + b*x**n + c*x**(2*n))**p*appellf1(3 + 1/n, -p, -p, 4 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((3*n + 1)*(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_92():
    f = (d + e*x**n)**2*(a + b*x**n + c*x**(2*n))**p
    F = d**2*x*(a + b*x**n + c*x**(2*n))**p*appellf1(1/n, -p, -p, 1 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)**p) + 2*d*e*x**(n + 1)*(a + b*x**n + c*x**(2*n))**p*appellf1(1 + 1/n, -p, -p, 2 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((n + 1)*(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)**p) + e**2*x**(2*n + 1)*(a + b*x**n + c*x**(2*n))**p*appellf1(2 + 1/n, -p, -p, 3 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((2*n + 1)*(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_93():
    f = (d + e*x**n)*(a + b*x**n + c*x**(2*n))**p
    F = d*x*(a + b*x**n + c*x**(2*n))**p*appellf1(1/n, -p, -p, 1 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)**p) + e*x**(n + 1)*(a + b*x**n + c*x**(2*n))**p*appellf1(1 + 1/n, -p, -p, 2 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/((n + 1)*(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_94():
    f = (a + b*x**n + c*x**(2*n))**p/(d + e*x**n)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * (x)**(Symbol('n'))) + (Symbol('c') * (x)**((Integer(2) * Symbol('n'))))))**(Symbol('p')) * ((Symbol('d') + (Symbol('e') * (x)**(Symbol('n')))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_95():
    f = (a + b*x**n + c*x**(2*n))**p/(d + e*x**n)**2
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * (x)**(Symbol('n'))) + (Symbol('c') * (x)**((Integer(2) * Symbol('n'))))))**(Symbol('p')) * (((Symbol('d') + (Symbol('e') * (x)**(Symbol('n')))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_3_d_plus_e_x_pow_n_pow_q_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_96():
    f = (a + b*x**n + c*x**(2*n))**p/(d + e*x**n)**3
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * (x)**(Symbol('n'))) + (Symbol('c') * (x)**((Integer(2) * Symbol('n'))))))**(Symbol('p')) * (((Symbol('d') + (Symbol('e') * (x)**(Symbol('n')))))**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F

