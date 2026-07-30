"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.2 Trinomial products/1.2.1 Quadratic/1.2.1.6 (g+h x)^m (a+b x+c x^2)^p (d+e x+f x^2)^q.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, a, b, c, d, e, f, g, h = symbols('A B a b c d e f g h')

def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_1():
    f = (A + B*x)*(a + b*x + c*x**2)/(d + f*x**2)
    F = B*c*x**2/(2*f) + x*(A*c + B*b)/f - (-A*b*f - B*a*f + B*c*d)*log(d + f*x**2)/(2*f**2) - (-A*a*f + A*c*d + B*b*d)*atan(sqrt(f)*x/sqrt(d))/(sqrt(d)*f**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_2():
    f = (A + B*x)*(a + b*x + c*x**2)**2/(d + f*x**2)
    F = B*c**2*x**4/(4*f) + c*x**3*(A*c + 2*B*b)/(3*f) + x**2*(2*A*b*c*f - B*(-2*a*c*f - b**2*f + c**2*d))/(2*f**2) + x*(A*b**2*f - A*c*(-2*a*f + c*d) - B*b*(-2*a*f + 2*c*d))/f**2 - (2*A*b*f*(-a*f + c*d) - B*(-2*a*c*d*f + c**2*d**2 - f*(-a**2*f + b**2*d)))*log(d + f*x**2)/(2*f**3) - (A*b**2*d*f - A*(-a*f + c*d)**2 - 2*B*b*d*(-a*f + c*d))*atan(sqrt(f)*x/sqrt(d))/(sqrt(d)*f**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_3():
    f = (A + B*x)*(a + b*x + c*x**2)**3/(d + f*x**2)
    F = B*c**3*x**6/(6*f) + c**2*x**5*(A*c + 3*B*b)/(5*f) + c*x**4*(3*A*b*c*f - B*(-3*a*c*f - 3*b**2*f + c**2*d))/(4*f**2) + x**3*(3*A*b**2*c*f - A*c**2*(-3*a*f + c*d) + B*b**3*f - 3*B*b*c*(-2*a*f + c*d))/(3*f**2) - x**2*(A*b*f*(-6*a*c*f - b**2*f + 3*c**2*d) - B*(3*a*b**2*f**2 - 3*a*c**2*d*f + c**3*d**2 - 3*c*f*(-a**2*f + b**2*d)))/(2*f**3) - x*(3*A*b**2*f*(-a*f + c*d) - A*c*(3*a**2*f**2 - 3*a*c*d*f + c**2*d**2) + B*b**3*d*f - 3*B*b*(-a*f + c*d)**2)/f**3 + (A*b*f*(-6*a*c*d*f + 3*c**2*d**2 - f*(-3*a**2*f + b**2*d)) - B*(-a*f + c*d)*(-2*a*c*d*f + c**2*d**2 - f*(-a**2*f + 3*b**2*d)))*log(d + f*x**2)/(2*f**4) + (3*A*b**2*d*f*(-a*f + c*d) - A*(-a*f + c*d)**3 + B*b**3*d**2*f - 3*B*b*d*(-a*f + c*d)**2)*atan(sqrt(f)*x/sqrt(d))/(sqrt(d)*f**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_4():
    f = (A + B*x)/((d + f*x**2)*(a + b*x + c*x**2))
    F = -(A*b*f - B*a*f + B*c*d)*log(d + f*x**2)/(-4*a*c*d*f + 2*c**2*d**2 + 2*f*(a**2*f + b**2*d)) + (A*b*f - B*a*f + B*c*d)*log(a + b*x + c*x**2)/(-4*a*c*d*f + 2*c**2*d**2 + 2*f*(a**2*f + b**2*d)) - (A*b**2*f + 2*A*c*(-a*f + c*d) - B*b*(a*f + c*d))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(sqrt(-4*a*c + b**2)*(-2*a*c*d*f + c**2*d**2 + f*(a**2*f + b**2*d))) + sqrt(f)*(A*a*f - A*c*d + B*b*d)*atan(sqrt(f)*x/sqrt(d))/(sqrt(d)*(-2*a*c*d*f + c**2*d**2 + f*(a**2*f + b**2*d)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_5():
    f = (A + B*x)/((d + f*x**2)*(a + b*x + c*x**2)**2)
    F = f*(2*A*b*f*(-a*f + c*d) + B*(-2*a*c*d*f + c**2*d**2 - f*(-a**2*f + b**2*d)))*log(d + f*x**2)/(2*(-2*a*c*d*f + c**2*d**2 + f*(a**2*f + b**2*d))**2) - f*(2*A*b*f*(-a*f + c*d) + B*(-2*a*c*d*f + c**2*d**2 - f*(-a**2*f + b**2*d)))*log(a + b*x + c*x**2)/(2*(-2*a*c*d*f + c**2*d**2 + f*(a**2*f + b**2*d))**2) + (A*b*c*(a*f + c*d) - c*x*(A*b**2*f + 2*A*c*(-a*f + c*d) - B*b*(a*f + c*d)) - (A*b - B*a)*(-2*a*c*f + b**2*f + 2*c**2*d))/((-4*a*c + b**2)*(b**2*d*f + (-a*f + c*d)**2)*(a + b*x + c*x**2)) - (-2*A*b**4*f**2*(-a*f + c*d) - 4*A*b**2*c*f*(3*a**2*f**2 - 3*a*c*d*f + 2*c**2*d**2) - 4*A*c**2*(-3*a*f + c*d)*(-a*f + c*d)**2 + B*b**5*d*f**2 + B*b**3*f*(-a**2*f**2 - 4*a*c*d*f + 5*c**2*d**2) + 2*B*b*c*(3*a**3*f**3 + 3*a**2*c*d*f**2 - 7*a*c**2*d**2*f + c**3*d**3))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/((-4*a*c + b**2)**(sympy.S(3)/2)*(-2*a*c*d*f + c**2*d**2 + f*(a**2*f + b**2*d))**2) - f**(sympy.S(3)/2)*(A*b**2*d*f - A*(-a*f + c*d)**2 + 2*B*b*d*(-a*f + c*d))*atan(sqrt(f)*x/sqrt(d))/(sqrt(d)*(-2*a*c*d*f + c**2*d**2 + f*(a**2*f + b**2*d))**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_6():
    f = (A + B*x)*sqrt(a + b*x + c*x**2)/(d - f*x**2)
    F = -B*sqrt(a + b*x + c*x**2)/f - (-A*sqrt(f) + B*sqrt(d))*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*sqrt(d)*f**(sympy.S(3)/2)) + (A*sqrt(f) + B*sqrt(d))*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*sqrt(d)*f**(sympy.S(3)/2)) - (2*A*c + B*b)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_7():
    f = (A + B*x)/((d - f*x**2)*sqrt(a + b*x + c*x**2))
    F = -(-A*sqrt(f)/sqrt(d) + B)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*sqrt(f)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)) + (A*sqrt(f)/sqrt(d) + B)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*sqrt(f)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_8():
    f = (A + B*x)/((d - f*x**2)*(a + b*x + c*x**2)**(sympy.S(3)/2))
    F = -(2*A*(b**3*f - b*c*(3*a*f + c*d)) + 2*B*a*(2*a*c*f - b**2*f + 2*c**2*d) + 2*c*x*(A*b**2*f - 2*A*c*(a*f + c*d) + B*b*(-a*f + c*d)))/((-4*a*c + b**2)*(b**2*d*f - (a*f + c*d)**2)*sqrt(a + b*x + c*x**2)) - sqrt(f)*(-A*sqrt(f) + B*sqrt(d))*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*sqrt(d)*(a*f - b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)) + sqrt(f)*(A*sqrt(f) + B*sqrt(d))*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*sqrt(d)*(a*f + b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_9():
    f = (A + B*x)/((d - f*x**2)*(a + b*x + c*x**2)**(sympy.S(5)/2))
    F = -(-2*A*b**5*f**2*(6*a*f + 7*c*d) + 2*A*b**3*c*f*(43*a**2*f**2 + 46*a*c*d*f + 15*c**2*d**2) - 8*A*b*c**2*(17*a**3*f**3 + 24*a**2*c*d*f**2 + 9*a*c**2*d**2*f + 2*c**3*d**3) + 48*B*a**2*c**2*f*(a*f + c*d)**2 + 6*B*b**6*d*f**2 - 2*B*b**4*f*(-3*a**2*f**2 + 14*a*c*d*f + 7*c**2*d**2) + 4*B*b**2*c*(-11*a**3*f**3 + 4*a**2*c*d*f**2 + 5*a*c**2*d**2*f + 2*c**3*d**3) + 2*c*x*(-2*A*b**4*f**2*(3*a*f + 4*c*d) + 2*A*b**2*c*f*(19*a**2*f**2 + 22*a*c*d*f + 15*c**2*d**2) - 8*A*c**2*(a*f + c*d)**2*(5*a*f + 2*c*d) + 3*B*b**5*d*f**2 - B*b**3*f*(-3*a**2*f**2 + 10*a*c*d*f + 17*c**2*d**2) + 4*B*b*c*(-5*a**3*f**3 + 4*a**2*c*d*f**2 + 11*a*c**2*d**2*f + 2*c**3*d**3)))/(3*(-4*a*c + b**2)**2*sqrt(a + b*x + c*x**2)*(2*a*c*d*f + c**2*d**2 - f*(-a**2*f + b**2*d))**2) - (2*A*(b**3*f - b*c*(3*a*f + c*d)) + 2*B*a*(2*a*c*f - b**2*f + 2*c**2*d) + 2*c*x*(A*b**2*f - 2*A*c*(a*f + c*d) + B*b*(-a*f + c*d)))/((-12*a*c + 3*b**2)*(b**2*d*f - (a*f + c*d)**2)*(a + b*x + c*x**2)**(sympy.S(3)/2)) - f**(sympy.S(3)/2)*(-A*sqrt(f) + B*sqrt(d))*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*sqrt(d)*(a*f - b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(5)/2)) + f**(sympy.S(3)/2)*(A*sqrt(f) + B*sqrt(d))*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*sqrt(d)*(a*f + b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_10():
    f = (2*x + 1)/((x**2 - 1)*sqrt(x**2 + x - 1))
    F = -atan((x + 3)/(2*sqrt(x**2 + x - 1)))/2 + 3*atanh((1 - 3*x)/(2*sqrt(x**2 + x - 1)))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_11():
    f = (2*x + 1)/((x**2 + 1)*sqrt(x**2 + x - 1))
    F = -sqrt(1 + sqrt(5)/2)*atan((-sqrt(5)*x + 2*sqrt(5) + 5)/(sqrt(20 + 10*sqrt(5))*sqrt(x**2 + x - 1))) + sqrt(-1 + sqrt(5)/2)*atanh((sqrt(5)*x - 2*sqrt(5) + 5)/(sqrt(-20 + 10*sqrt(5))*sqrt(x**2 + x - 1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_12():
    f = (a + b*x - c)/((x**2 + 1)*sqrt(a + b*x + c*x**2))
    F = -sqrt(2)*sqrt(a**2 - a*(2*c - sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c - sqrt(a**2 - 2*a*c + b**2 + c**2)))*atan(sqrt(2)*(b*sqrt(a**2 - 2*a*c + b**2 + c**2) - x*(b**2 + (a - c)*(a - c + sqrt(a**2 - 2*a*c + b**2 + c**2))))/(2*sqrt(a + b*x + c*x**2)*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4)*sqrt(a**2 - a*(2*c - sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c - sqrt(a**2 - 2*a*c + b**2 + c**2)))))/(2*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4)) - sqrt(2)*sqrt(a**2 - a*(2*c + sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c + sqrt(a**2 - 2*a*c + b**2 + c**2)))*atanh(sqrt(2)*(b*sqrt(a**2 - 2*a*c + b**2 + c**2) + x*(b**2 + (a - c)*(a - c - sqrt(a**2 - 2*a*c + b**2 + c**2))))/(2*sqrt(a + b*x + c*x**2)*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4)*sqrt(a**2 - a*(2*c + sqrt(a**2 - 2*a*c + b**2 + c**2)) + b**2 + c*(c + sqrt(a**2 - 2*a*c + b**2 + c**2)))))/(2*(a**2 - 2*a*c + b**2 + c**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_13():
    f = (A + B*x)*(a + b*x + c*x**2)/(d + e*x + f*x**2)
    F = B*c*x**2/(2*f) - x*(-A*c*f - B*b*f + B*c*e)/f**2 - (A*f*(-b*f + c*e) - B*(a*f**2 - b*e*f - c*d*f + c*e**2))*log(d + e*x + f*x**2)/(2*f**3) - (A*f*(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2) + B*(-c*(-3*d*e*f + e**3) + f*(-a*e*f - 2*b*d*f + b*e**2)))*atanh((e + 2*f*x)/sqrt(-4*d*f + e**2))/(f**3*sqrt(-4*d*f + e**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_14():
    f = (A + B*x)*(a + b*x + c*x**2)**2/(d + e*x + f*x**2)
    F = B*c**2*x**4/(4*f) - c*x**3*(-A*c*f - 2*B*b*f + B*c*e)/(3*f**2) - x**2*(A*c*f*(-2*b*f + c*e) - B*(b**2*f**2 + c**2*(-d*f + e**2) - 2*c*f*(-a*f + b*e)))/(2*f**3) + x*(A*f*(b**2*f**2 + c**2*(-d*f + e**2) - 2*c*f*(-a*f + b*e)) + B*(-b*f + c*e)*(-c*(-2*d*f + e**2) + f*(-2*a*f + b*e)))/f**4 + (A*f*(-b*f + c*e)*(-c*(-2*d*f + e**2) + f*(-2*a*f + b*e)) + B*(c**2*(d**2*f**2 - 3*d*e**2*f + e**4) + 2*c*f*(a*f*(-d*f + e**2) - b*(-2*d*e*f + e**3)) - f**2*(-a**2*f**2 + 2*a*b*e*f - b**2*(-d*f + e**2))))*log(d + e*x + f*x**2)/(2*f**5) - (A*f*(c**2*(2*d**2*f**2 - 4*d*e**2*f + e**4) + 2*c*f*(a*f*(-2*d*f + e**2) - b*(-3*d*e*f + e**3)) - f**2*(-2*a**2*f**2 + 2*a*b*e*f - b**2*(-2*d*f + e**2))) - B*(c**2*(5*d**2*e*f**2 - 5*d*e**3*f + e**5) + 2*c*f*(a*e*f*(-3*d*f + e**2) - b*(2*d**2*f**2 - 4*d*e**2*f + e**4)) + f**2*(a**2*e*f**2 - 2*a*b*f*(-2*d*f + e**2) + b**2*(-3*d*e*f + e**3))))*atanh((e + 2*f*x)/sqrt(-4*d*f + e**2))/(f**5*sqrt(-4*d*f + e**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_15():
    f = (A + B*x)/((a + b*x + c*x**2)*(d + e*x + f*x**2))
    F = (-A*(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2) + B*(a*e*f - 2*b*d*f + c*d*e))*atanh((e + 2*f*x)/sqrt(-4*d*f + e**2))/(sqrt(-4*d*f + e**2)*(c**2*d**2 - c*(-a*(-2*d*f + e**2) + b*d*e) + f*(a**2*f - a*b*e + b**2*d))) + (A*b*f - A*c*e - B*a*f + B*c*d)*log(a + b*x + c*x**2)/(2*c**2*d**2 - 2*c*(-a*(-2*d*f + e**2) + b*d*e) + 2*f*(a**2*f - a*b*e + b**2*d)) - (A*b*f - A*c*e - B*a*f + B*c*d)*log(d + e*x + f*x**2)/(2*c**2*d**2 - 2*c*(-a*(-2*d*f + e**2) + b*d*e) + 2*f*(a**2*f - a*b*e + b**2*d)) - (A*b**2*f - b*(A*c*e + B*a*f + B*c*d) + 2*c*(-A*a*f + A*c*d + B*a*e))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(sqrt(-4*a*c + b**2)*(c**2*d**2 - c*(-a*(-2*d*f + e**2) + b*d*e) + f*(a**2*f - a*b*e + b**2*d)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_16():
    f = (A + B*x)/((a + b*x + c*x**2)**2*(d + e*x + f*x**2))
    F = (-A*(c**2*(2*d**2*f**2 - 4*d*e**2*f + e**4) + 2*c*f*(a*f*(-2*d*f + e**2) - b*(-3*d*e*f + e**3)) - f**2*(-2*a**2*f**2 + 2*a*b*e*f - b**2*(-2*d*f + e**2))) + B*(c**2*d*e*(-3*d*f + e**2) - 2*c*d*f*(-a*e*f - 2*b*d*f + b*e**2) + f**2*(a**2*e*f - 4*a*b*d*f + b**2*d*e)))*atanh((e + 2*f*x)/sqrt(-4*d*f + e**2))/(sqrt(-4*d*f + e**2)*(c**2*d**2 - c*(-a*(-2*d*f + e**2) + b*d*e) + f*(a**2*f - a*b*e + b**2*d))**2) + (A*(-b*f + c*e)*(-c*(-2*d*f + e**2) + f*(-2*a*f + b*e)) - B*(-c**2*d*(-d*f + e**2) + 2*c*d*f*(-a*f + b*e) - f**2*(-a**2*f + b**2*d)))*log(a + b*x + c*x**2)/(2*(c**2*d**2 - c*(-a*(-2*d*f + e**2) + b*d*e) + f*(a**2*f - a*b*e + b**2*d))**2) - (A*(-b*f + c*e)*(-c*(-2*d*f + e**2) + f*(-2*a*f + b*e)) - B*(-c**2*d*(-d*f + e**2) + 2*c*d*f*(-a*f + b*e) - f**2*(-a**2*f + b**2*d)))*log(d + e*x + f*x**2)/(2*(c**2*d**2 - c*(-a*(-2*d*f + e**2) + b*d*e) + f*(a**2*f - a*b*e + b**2*d))**2) - (A*c*(2*a*c*e - b*(a*f + c*d)) + c*x*(A*b**2*f - b*(A*c*e + B*a*f + B*c*d) + 2*c*(-A*a*f + A*c*d + B*a*e)) + (A*b - B*a)*(b**2*f + 2*c**2*d - c*(2*a*f + b*e)))/((-4*a*c + b**2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)*(a + b*x + c*x**2)) - (b**5*f**2*(-A*e + B*d) - 2*b**4*f*(-A*(a*f**2 - c*d*f + c*e**2) + B*c*d*e) - b**3*(A*c*e*(-4*a*f**2 - 2*c*d*f + c*e**2) + B*(a**2*f**3 + 4*a*c*d*f**2 - c**2*d*(5*d*f + e**2))) - 4*b**2*c*(A*f*(3*a**2*f**2 + 3*a*c*(-d*f + e**2) + 2*c**2*d**2) + B*c**2*d**2*e) + 2*b*c*(A*c*e*(3*a**2*f**2 + a*c*(2*d*f + 3*e**2) + 3*c**2*d**2) + B*(3*a**3*f**3 + 3*a**2*c*f*(d*f + e**2) + a*c**2*d*(-7*d*f + e**2) + c**3*d**3)) - 4*c**2*(A*(-3*a**3*f**3 - a**2*c*f*(-7*d*f + e**2) + a*c**2*d*(-5*d*f + 3*e**2) + c**3*d**3) - B*a*e*(-3*a**2*f**2 - a*c*(-2*d*f + e**2) + c**2*d**2)))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/((-4*a*c + b**2)**(sympy.S(3)/2)*(c**2*d**2 - c*(-a*(-2*d*f + e**2) + b*d*e) + f*(a**2*f - a*b*e + b**2*d))**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_17():
    f = (g + h*x)/((a + b*x + c*x**2)*(a*d + b*d*x + c*d*x**2)**2)
    F = -6*c*(-b*h + 2*c*g)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(d**2*(-4*a*c + b**2)**(sympy.S(5)/2)) + (b + 2*c*x)*(-3*b*h + 6*c*g)/(2*d**2*(-4*a*c + b**2)**2*(a + b*x + c*x**2)) - (-2*a*h + b*g + x*(-b*h + 2*c*g))/(d**2*(-8*a*c + 2*b**2)*(a + b*x + c*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_18():
    f = (g + h*x)/((a + b*x + c*x**2)**2*(a*d + b*d*x + c*d*x**2))
    F = -6*c*(-b*h + 2*c*g)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(d*(-4*a*c + b**2)**(sympy.S(5)/2)) + (b + 2*c*x)*(-3*b*h + 6*c*g)/(2*d*(-4*a*c + b**2)**2*(a + b*x + c*x**2)) - (-2*a*h + b*g + x*(-b*h + 2*c*g))/(d*(-8*a*c + 2*b**2)*(a + b*x + c*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_19():
    f = (A + B*x)*sqrt(a + b*x + c*x**2)/(d + e*x + f*x**2)
    F = B*sqrt(a + b*x + c*x**2)/f + sqrt(2)*(2*f*(A*f*(-a*f + c*d) - B*d*(-b*f + c*e)) - (e - sqrt(-4*d*f + e**2))*(A*f*(-b*f + c*e) + B*(-c*(-d*f + e**2) + f*(-a*f + b*e))))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*f**2*sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))) - sqrt(2)*(2*f*(A*f*(-a*f + c*d) - B*d*(-b*f + c*e)) - (e + sqrt(-4*d*f + e**2))*(A*f*(-b*f + c*e) + B*(-c*(-d*f + e**2) + f*(-a*f + b*e))))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*f**2*sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))) - (-2*A*c*f - B*b*f + 2*B*c*e)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*sqrt(c)*f**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_20():
    f = (A + B*x)*(a + b*x + c*x**2)**(sympy.S(3)/2)/(d + e*x + f*x**2)
    F = B*(a + b*x + c*x**2)**(sympy.S(3)/2)/(3*f) + sqrt(2)*(2*f*(A*f*(-c**2*d*(-d*f + e**2) + 2*c*d*f*(-a*f + b*e) - f**2*(-a**2*f + b**2*d)) + B*d*(-b*f + c*e)*(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2)) - (e + sqrt(-4*d*f + e**2))*(A*f*(-b*f + c*e)*(-c*(-2*d*f + e**2) + f*(-2*a*f + b*e)) + B*(c**2*(d**2*f**2 - 3*d*e**2*f + e**4) + 2*c*f*(a*f*(-d*f + e**2) - b*(-2*d*e*f + e**3)) - f**2*(-a**2*f**2 + 2*a*b*e*f - b**2*(-d*f + e**2)))))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*f**4*sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))) - sqrt(a + b*x + c*x**2)*(2*A*c*f*(-5*b*f + 4*c*e) - B*(b**2*f**2 + 8*c**2*(-d*f + e**2) - 2*c*f*(-4*a*f + 5*b*e)) + 2*c*f*x*(-2*A*c*f - B*b*f + 2*B*c*e))/(8*c*f**3) - sqrt(2)*(2*c*f*(A*f*(-c**2*d*(-d*f + e**2) + 2*c*d*f*(-a*f + b*e) - f**2*(-a**2*f + b**2*d)) + B*d*(-b*f + c*e)*(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2)) - c*(e - sqrt(-4*d*f + e**2))*(A*f*(-b*f + c*e)*(-c*(-2*d*f + e**2) + f*(-2*a*f + b*e)) + B*(c**2*(d**2*f**2 - 3*d*e**2*f + e**4) + 2*c*f*(a*f*(-d*f + e**2) - b*(-2*d*e*f + e**3)) - f**2*(-a**2*f**2 + 2*a*b*e*f - b**2*(-d*f + e**2)))))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*c*f**4*sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))) + (2*A*c*f*(3*b**2*f**2 + 8*c**2*(-d*f + e**2) - 12*c*f*(-a*f + b*e)) - B*(b**3*f**3 + 6*b*c*f**2*(-2*a*f + b*e) + 16*c**3*(-2*d*e*f + e**3) - 24*c**2*f*(-a*e*f - b*d*f + b*e**2)))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(3)/2)*f**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_21():
    f = (A + B*x)/((a + b*x + c*x**2)*sqrt(d + e*x + f*x**2))
    F = sqrt(2)*(2*A*c - B*(b + sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*(4*c*d - e*(b + sqrt(-4*a*c + b**2)) + x*(2*c*e - 2*f*(b + sqrt(-4*a*c + b**2))))/(4*sqrt(d + e*x + f*x**2)*sqrt(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d - sqrt(-4*a*c + b**2)*(-b*f + c*e))))/(2*sqrt(-4*a*c + b**2)*sqrt(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d - sqrt(-4*a*c + b**2)*(-b*f + c*e))) + sqrt(2)*(-2*A*c + B*b - B*sqrt(-4*a*c + b**2))*atanh(sqrt(2)*(4*c*d - e*(b - sqrt(-4*a*c + b**2)) + x*(2*c*e - 2*f*(b - sqrt(-4*a*c + b**2))))/(4*sqrt(d + e*x + f*x**2)*sqrt(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d + sqrt(-4*a*c + b**2)*(-b*f + c*e))))/(2*sqrt(-4*a*c + b**2)*sqrt(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d + sqrt(-4*a*c + b**2)*(-b*f + c*e)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_22():
    f = (A + B*x)/((a + c*x**2)*sqrt(d + e*x + f*x**2))
    F = sqrt(2)*sqrt(A*(-a*f + c*d - sqrt(a**2*f**2 + a*c*(-2*d*f + e**2) + c**2*d**2)) + B*a*e)*sqrt(-A*c*e + B*(-a*f + c*d + sqrt(a**2*f**2 + a*c*(-2*d*f + e**2) + c**2*d**2)))*atanh(sqrt(2)*sqrt(e)*(a*(A*c*e - B*(-a*f + c*d + sqrt(a**2*f**2 + a*c*(-2*d*f + e**2) + c**2*d**2))) - c*x*(A*(-a*f + c*d - sqrt(a**2*f**2 + a*c*(-2*d*f + e**2) + c**2*d**2)) + B*a*e))/(2*sqrt(a)*sqrt(c)*sqrt(A*(-a*f + c*d - sqrt(a**2*f**2 + a*c*(-2*d*f + e**2) + c**2*d**2)) + B*a*e)*sqrt(-A*c*e + B*(-a*f + c*d + sqrt(a**2*f**2 + a*c*(-2*d*f + e**2) + c**2*d**2)))*sqrt(d + e*x + f*x**2)))/(2*sqrt(a)*sqrt(c)*sqrt(e)*sqrt(a**2*f**2 + a*c*(-2*d*f + e**2) + c**2*d**2)) - sqrt(2)*sqrt(A*(-a*f + c*d + sqrt(a**2*f**2 + a*c*(-2*d*f + e**2) + c**2*d**2)) + B*a*e)*sqrt(-A*c*e + B*(-a*f + c*d - sqrt(a**2*f**2 + a*c*(-2*d*f + e**2) + c**2*d**2)))*atanh(sqrt(2)*sqrt(e)*(a*(A*c*e - B*(-a*f + c*d - sqrt(a**2*f**2 + a*c*(-2*d*f + e**2) + c**2*d**2))) - c*x*(A*(-a*f + c*d + sqrt(a**2*f**2 + a*c*(-2*d*f + e**2) + c**2*d**2)) + B*a*e))/(2*sqrt(a)*sqrt(c)*sqrt(A*(-a*f + c*d + sqrt(a**2*f**2 + a*c*(-2*d*f + e**2) + c**2*d**2)) + B*a*e)*sqrt(-A*c*e + B*(-a*f + c*d - sqrt(a**2*f**2 + a*c*(-2*d*f + e**2) + c**2*d**2)))*sqrt(d + e*x + f*x**2)))/(2*sqrt(a)*sqrt(c)*sqrt(e)*sqrt(a**2*f**2 + a*c*(-2*d*f + e**2) + c**2*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_23():
    f = (A + B*x)/(sqrt(d + f*x**2)*(a + b*x + c*x**2))
    F = sqrt(2)*(2*A*c - B*(b + sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*(2*c*d - f*x*(b + sqrt(-4*a*c + b**2)))/(2*sqrt(d + f*x**2)*sqrt(-2*a*c*f + b*f*(b + sqrt(-4*a*c + b**2)) + 2*c**2*d)))/(2*sqrt(-4*a*c + b**2)*sqrt(-2*a*c*f + b*f*(b + sqrt(-4*a*c + b**2)) + 2*c**2*d)) + sqrt(2)*(-2*A*c + B*b - B*sqrt(-4*a*c + b**2))*atanh(sqrt(2)*(2*c*d - f*x*(b - sqrt(-4*a*c + b**2)))/(2*sqrt(d + f*x**2)*sqrt(-2*a*c*f + b*f*(b - sqrt(-4*a*c + b**2)) + 2*c**2*d)))/(2*sqrt(-4*a*c + b**2)*sqrt(-2*a*c*f + b*f*(b - sqrt(-4*a*c + b**2)) + 2*c**2*d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_24():
    f = (A + B*x)/((a + c*x**2)*sqrt(d + f*x**2))
    F = A*atan(x*sqrt(-a*f + c*d)/(sqrt(a)*sqrt(d + f*x**2)))/(sqrt(a)*sqrt(-a*f + c*d)) - B*atanh(sqrt(c)*sqrt(d + f*x**2)/sqrt(-a*f + c*d))/(sqrt(c)*sqrt(-a*f + c*d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_25():
    f = (x + 2)/((-3*x**2 + 4*x + 2)*sqrt(-2*x**2 + 3*x + 1))
    F = sqrt(sympy.S(-13)/5 + sqrt(10))*atan((x*(1 + 4*sqrt(10)) - 3*sqrt(10) + 12)/(2*sqrt(1 + sqrt(10))*sqrt(-2*x**2 + 3*x + 1)))/2 + sqrt(sympy.S(13)/5 + sqrt(10))*atanh((x*(1 - 4*sqrt(10)) + 3*sqrt(10) + 12)/(2*sqrt(-1 + sqrt(10))*sqrt(-2*x**2 + 3*x + 1)))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_26():
    f = (x + 2)/((-3*x**2 + 4*x + 2)*(-2*x**2 + 3*x + 1)**(sympy.S(3)/2))
    F = -(28*x + 30)/(17*sqrt(-2*x**2 + 3*x + 1)) - 9*sqrt(sympy.S(-3)/5 + sqrt(10)/5)*atan((x*(1 + 4*sqrt(10)) - 3*sqrt(10) + 12)/(2*sqrt(1 + sqrt(10))*sqrt(-2*x**2 + 3*x + 1)))/2 + 9*sqrt(sympy.S(3)/5 + sqrt(10)/5)*atanh((x*(1 - 4*sqrt(10)) + 3*sqrt(10) + 12)/(2*sqrt(-1 + sqrt(10))*sqrt(-2*x**2 + 3*x + 1)))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_27():
    f = (x + 2)/((-3*x**2 + 4*x + 2)*(-2*x**2 + 3*x + 1)**(sympy.S(5)/2))
    F = -(28*x + 30)/(51*(-2*x**2 + 3*x + 1)**(sympy.S(3)/2)) - (9628*x + 582)/(867*sqrt(-2*x**2 + 3*x + 1)) + 9*sqrt(sympy.S(-53)/5 + 17*sqrt(10)/5)*atan((x*(1 + 4*sqrt(10)) - 3*sqrt(10) + 12)/(2*sqrt(1 + sqrt(10))*sqrt(-2*x**2 + 3*x + 1)))/2 + 9*sqrt(sympy.S(53)/5 + 17*sqrt(10)/5)*atanh((x*(1 - 4*sqrt(10)) + 3*sqrt(10) + 12)/(2*sqrt(-1 + sqrt(10))*sqrt(-2*x**2 + 3*x + 1)))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_28():
    f = (x + 2)/((-3*x**2 + 4*x + 2)*sqrt(2*x**2 + 3*x + 1))
    F = -sqrt(7*sqrt(10)/25 + 1)*atanh((x*(17 - 4*sqrt(10)) - 3*sqrt(10) + 12)/(2*sqrt(55 - 17*sqrt(10))*sqrt(2*x**2 + 3*x + 1)))/2 + sqrt(1 - 7*sqrt(10)/25)*atanh((x*(4*sqrt(10) + 17) + 3*sqrt(10) + 12)/(2*sqrt(17*sqrt(10) + 55)*sqrt(2*x**2 + 3*x + 1)))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_29():
    f = (x + 2)/((-3*x**2 + 4*x + 2)*(2*x**2 + 3*x + 1)**(sympy.S(3)/2))
    F = (44*x + 42)/(5*sqrt(2*x**2 + 3*x + 1)) - sqrt(1959*sqrt(10)/5 + 1239)*atanh((x*(17 - 4*sqrt(10)) - 3*sqrt(10) + 12)/(2*sqrt(55 - 17*sqrt(10))*sqrt(2*x**2 + 3*x + 1)))/10 + sqrt(1239 - 1959*sqrt(10)/5)*atanh((x*(4*sqrt(10) + 17) + 3*sqrt(10) + 12)/(2*sqrt(17*sqrt(10) + 55)*sqrt(2*x**2 + 3*x + 1)))/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_30():
    f = (x + 2)/((-3*x**2 + 4*x + 2)*(2*x**2 + 3*x + 1)**(sympy.S(5)/2))
    F = (44*x + 42)/(15*(2*x**2 + 3*x + 1)**(sympy.S(3)/2)) + (460*x + 546)/(15*sqrt(2*x**2 + 3*x + 1)) - sqrt(1544809*sqrt(10)/3 + sympy.S(4885115)/3)*atanh((x*(17 - 4*sqrt(10)) - 3*sqrt(10) + 12)/(2*sqrt(55 - 17*sqrt(10))*sqrt(2*x**2 + 3*x + 1)))/50 + sqrt(sympy.S(4885115)/3 - 1544809*sqrt(10)/3)*atanh((x*(4*sqrt(10) + 17) + 3*sqrt(10) + 12)/(2*sqrt(17*sqrt(10) + 55)*sqrt(2*x**2 + 3*x + 1)))/50
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_31():
    f = (x + 1)/((x**2 + 2*x + 4)*sqrt(x**2 + 2*x + 5))
    F = -atanh(sqrt(x**2 + 2*x + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_32():
    f = (x + 4)/((x**2 + 2*x + 4)*sqrt(x**2 + 2*x + 5))
    F = sqrt(3)*atan(sqrt(3)*(x + 1)/(3*sqrt(x**2 + 2*x + 5))) - atanh(sqrt(x**2 + 2*x + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_33():
    f = (2*x + 1)/((x**2 + x + 3)*sqrt(x**2 + x + 5))
    F = -sqrt(2)*atanh(sqrt(2)*sqrt(x**2 + x + 5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_34():
    f = x/((x**2 + x + 3)*sqrt(x**2 + x + 5))
    F = -sqrt(22)*atan(sqrt(22)*(2*x + 1)/(11*sqrt(x**2 + x + 5)))/22 - sqrt(2)*atanh(sqrt(2)*sqrt(x**2 + x + 5)/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_35():
    f = (A + B*x)/(sqrt(d + e*x + f*x**2)*(a*e + b*e*x + b*f*x**2)**2)
    F = B*atanh(sqrt(b)*sqrt(d + e*x + f*x**2)/sqrt(-a*e + b*d))/(2*sqrt(b)*f*(-a*e + b*d)**(sympy.S(3)/2)) - (-b*x*(-2*A*f + B*e) + e*(A*b - 2*B*a))*sqrt(d + e*x + f*x**2)/(e*(-a*e + b*d)*(-4*a*f + b*e)*(a*e + b*e*x + b*f*x**2)) + (-2*A*f + B*e)*(8*a*e*f - b*(4*d*f + e**2))*atanh((e + 2*f*x)*sqrt(-a*e + b*d)/(sqrt(e)*sqrt(-4*a*f + b*e)*sqrt(d + e*x + f*x**2)))/(2*e**(sympy.S(3)/2)*f*(-a*e + b*d)**(sympy.S(3)/2)*(-4*a*f + b*e)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_36():
    f = (g + h*x)*sqrt(a + b*x + c*x**2)/(a*d + b*d*x + c*d*x**2)**2
    F = -(-4*a*h + 2*b*g + 2*x*(-b*h + 2*c*g))/(d**2*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_37():
    f = (2*x + 3)/(sqrt(-x**2 - 4*x - 3)*(2*x**2 + 4*x + 3))
    F = atanh(x/sqrt(-x**2 - 4*x - 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_38():
    f = (4*x + 3)/(sqrt(-x**2 - 4*x - 3)*(2*x**2 + 4*x + 3))
    F = sqrt(2)*atan(sqrt(2)*(-(x + 3)/sqrt(-x**2 - 4*x - 3) + 1)/2) - sqrt(2)*atan(sqrt(2)*((x + 3)/sqrt(-x**2 - 4*x - 3) + 1)/2) + atanh(x/sqrt(-x**2 - 4*x - 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_39():
    f = (g + h*x)*sqrt(a + b*x + c*x**2)/(a*d + b*d*x + c*d*x**2)**(sympy.S(3)/2)
    F = h*sqrt(a + b*x + c*x**2)*log(a + b*x + c*x**2)/(2*c*d*sqrt(a*d + b*d*x + c*d*x**2)) - (-b*h + 2*c*g)*sqrt(a + b*x + c*x**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c*d*sqrt(-4*a*c + b**2)*sqrt(a*d + b*d*x + c*d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_40():
    f = x**2*sqrt(c + d*x**2)*sqrt(a**2 + 2*a*b*x + b**2*x**2)
    F = -a*c**2*sqrt(a**2 + 2*a*b*x + b**2*x**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(8*d**(sympy.S(3)/2)*(a + b*x)) - a*c*x*sqrt(c + d*x**2)*sqrt(a**2 + 2*a*b*x + b**2*x**2)/(8*d*(a + b*x)) + b*x**2*(c + d*x**2)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x + b**2*x**2)/(5*d*(a + b*x)) - (c + d*x**2)**(sympy.S(3)/2)*(-15*a*d*x + 8*b*c)*sqrt(a**2 + 2*a*b*x + b**2*x**2)/(60*d**2*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_41():
    f = x*sqrt(c + d*x**2)*sqrt(a**2 + 2*a*b*x + b**2*x**2)
    F = -b*c**2*sqrt(a**2 + 2*a*b*x + b**2*x**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(8*d**(sympy.S(3)/2)*(a + b*x)) - b*c*x*sqrt(c + d*x**2)*sqrt(a**2 + 2*a*b*x + b**2*x**2)/(8*d*(a + b*x)) + (4*a + 3*b*x)*(c + d*x**2)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x + b**2*x**2)/(12*d*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_42():
    f = sqrt(c + d*x**2)*sqrt(a**2 + 2*a*b*x + b**2*x**2)
    F = a*c*sqrt(a**2 + 2*a*b*x + b**2*x**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(2*sqrt(d)*(a + b*x)) + a*x*sqrt(c + d*x**2)*sqrt(a**2 + 2*a*b*x + b**2*x**2)/(2*a + 2*b*x) + b*(c + d*x**2)**(sympy.S(3)/2)*sqrt(a**2 + 2*a*b*x + b**2*x**2)/(3*d*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_43():
    f = sqrt(c + d*x**2)*sqrt(a**2 + 2*a*b*x + b**2*x**2)/x
    F = -a*sqrt(c)*sqrt(a**2 + 2*a*b*x + b**2*x**2)*atanh(sqrt(c + d*x**2)/sqrt(c))/(a + b*x) + b*c*sqrt(a**2 + 2*a*b*x + b**2*x**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(2*sqrt(d)*(a + b*x)) + (2*a + b*x)*sqrt(c + d*x**2)*sqrt(a**2 + 2*a*b*x + b**2*x**2)/(2*a + 2*b*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_44():
    f = sqrt(c + d*x**2)*sqrt(a**2 + 2*a*b*x + b**2*x**2)/x**2
    F = a*sqrt(d)*sqrt(a**2 + 2*a*b*x + b**2*x**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(a + b*x) - b*sqrt(c)*sqrt(a**2 + 2*a*b*x + b**2*x**2)*atanh(sqrt(c + d*x**2)/sqrt(c))/(a + b*x) - (a - b*x)*sqrt(c + d*x**2)*sqrt(a**2 + 2*a*b*x + b**2*x**2)/(x*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_45():
    f = sqrt(c + d*x**2)*sqrt(a**2 + 2*a*b*x + b**2*x**2)/x**3
    F = -a*d*sqrt(a**2 + 2*a*b*x + b**2*x**2)*atanh(sqrt(c + d*x**2)/sqrt(c))/(2*sqrt(c)*(a + b*x)) + b*sqrt(d)*sqrt(a**2 + 2*a*b*x + b**2*x**2)*atanh(sqrt(d)*x/sqrt(c + d*x**2))/(a + b*x) - (a + 2*b*x)*sqrt(c + d*x**2)*sqrt(a**2 + 2*a*b*x + b**2*x**2)/(2*x**2*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_46():
    f = x**2*sqrt(a**2 + 2*a*b*x + b**2*x**2)*sqrt(c + d*x**2 + e*x)
    F = b*x**2*sqrt(a**2 + 2*a*b*x + b**2*x**2)*(c + d*x**2 + e*x)**(sympy.S(3)/2)/(5*d*(a + b*x)) - sqrt(a**2 + 2*a*b*x + b**2*x**2)*(c + d*x**2 + e*x)**(sympy.S(3)/2)*(50*a*d*e + 32*b*c*d - 35*b*e**2 - 6*d*x*(10*a*d - 7*b*e))/(240*d**3*(a + b*x)) - (2*d*x + e)*(2*a*d*(4*c*d - 5*e**2) - b*(12*c*d*e - 7*e**3))*sqrt(a**2 + 2*a*b*x + b**2*x**2)*sqrt(c + d*x**2 + e*x)/(128*d**4*(a + b*x)) - (4*c*d - e**2)*sqrt(a**2 + 2*a*b*x + b**2*x**2)*(8*a*c*d**2 - 10*a*d*e**2 - 12*b*c*d*e + 7*b*e**3)*atanh((2*d*x + e)/(2*sqrt(d)*sqrt(c + d*x**2 + e*x)))/(256*d**(sympy.S(9)/2)*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_47():
    f = x*sqrt(a**2 + 2*a*b*x + b**2*x**2)*sqrt(c + d*x**2 + e*x)
    F = sqrt(a**2 + 2*a*b*x + b**2*x**2)*(c + d*x**2 + e*x)**(sympy.S(3)/2)*(8*a*d + 6*b*d*x - 5*b*e)/(24*d**2*(a + b*x)) - (2*d*x + e)*sqrt(a**2 + 2*a*b*x + b**2*x**2)*sqrt(c + d*x**2 + e*x)*(8*a*d*e + 4*b*c*d - 5*b*e**2)/(64*d**3*(a + b*x)) - (4*c*d - e**2)*sqrt(a**2 + 2*a*b*x + b**2*x**2)*(8*a*d*e + 4*b*c*d - 5*b*e**2)*atanh((2*d*x + e)/(2*sqrt(d)*sqrt(c + d*x**2 + e*x)))/(128*d**(sympy.S(7)/2)*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_48():
    f = sqrt(a**2 + 2*a*b*x + b**2*x**2)*sqrt(c + d*x**2 + e*x)
    F = b*sqrt(a**2 + 2*a*b*x + b**2*x**2)*(c + d*x**2 + e*x)**(sympy.S(3)/2)/(3*d*(a + b*x)) + (2*a*d - b*e)*(2*d*x + e)*sqrt(a**2 + 2*a*b*x + b**2*x**2)*sqrt(c + d*x**2 + e*x)/(8*d**2*(a + b*x)) + (2*a*d - b*e)*(4*c*d - e**2)*sqrt(a**2 + 2*a*b*x + b**2*x**2)*atanh((2*d*x + e)/(2*sqrt(d)*sqrt(c + d*x**2 + e*x)))/(16*d**(sympy.S(5)/2)*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_49():
    f = sqrt(a**2 + 2*a*b*x + b**2*x**2)*sqrt(c + d*x**2 + e*x)/x
    F = -a*sqrt(c)*sqrt(a**2 + 2*a*b*x + b**2*x**2)*atanh((2*c + e*x)/(2*sqrt(c)*sqrt(c + d*x**2 + e*x)))/(a + b*x) + sqrt(a**2 + 2*a*b*x + b**2*x**2)*sqrt(c + d*x**2 + e*x)*(4*a*d + 2*b*d*x + b*e)/(4*d*(a + b*x)) + sqrt(a**2 + 2*a*b*x + b**2*x**2)*(4*a*d*e + 4*b*c*d - b*e**2)*atanh((2*d*x + e)/(2*sqrt(d)*sqrt(c + d*x**2 + e*x)))/(8*d**(sympy.S(3)/2)*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_50():
    f = sqrt(a**2 + 2*a*b*x + b**2*x**2)*sqrt(c + d*x**2 + e*x)/x**2
    F = -(a - b*x)*sqrt(a**2 + 2*a*b*x + b**2*x**2)*sqrt(c + d*x**2 + e*x)/(x*(a + b*x)) + (2*a*d + b*e)*sqrt(a**2 + 2*a*b*x + b**2*x**2)*atanh((2*d*x + e)/(2*sqrt(d)*sqrt(c + d*x**2 + e*x)))/(2*sqrt(d)*(a + b*x)) - (a*e + 2*b*c)*sqrt(a**2 + 2*a*b*x + b**2*x**2)*atanh((2*c + e*x)/(2*sqrt(c)*sqrt(c + d*x**2 + e*x)))/(2*sqrt(c)*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_51():
    f = sqrt(a**2 + 2*a*b*x + b**2*x**2)*sqrt(c + d*x**2 + e*x)/x**3
    F = b*sqrt(d)*sqrt(a**2 + 2*a*b*x + b**2*x**2)*atanh((2*d*x + e)/(2*sqrt(d)*sqrt(c + d*x**2 + e*x)))/(a + b*x) - (2*a*c + x*(a*e + 4*b*c))*sqrt(a**2 + 2*a*b*x + b**2*x**2)*sqrt(c + d*x**2 + e*x)/(4*c*x**2*(a + b*x)) - sqrt(a**2 + 2*a*b*x + b**2*x**2)*(4*a*c*d - a*e**2 + 4*b*c*e)*atanh((2*c + e*x)/(2*sqrt(c)*sqrt(c + d*x**2 + e*x)))/(8*c**(sympy.S(3)/2)*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_52():
    f = x**2*sqrt(a + c*x**2)/(d + e*x + f*x**2)
    F = -sqrt(a + c*x**2)*(2*e - f*x)/(2*f**2) + sqrt(2)*(-2*d*f*(a*f**2 + c*(-d*f + e**2)) + e*(e + sqrt(-4*d*f + e**2))*(a*f**2 + c*(-2*d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*f**3*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) - sqrt(2)*(-2*d*f*(a*f**2 + c*(-d*f + e**2)) + e*(e - sqrt(-4*d*f + e**2))*(a*f**2 + c*(-2*d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*f**3*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) + (a*f**2 + 2*c*(-d*f + e**2))*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(2*sqrt(c)*f**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_53():
    f = x*sqrt(a + c*x**2)/(d + e*x + f*x**2)
    F = -sqrt(c)*e*atanh(sqrt(c)*x/sqrt(a + c*x**2))/f**2 + sqrt(a + c*x**2)/f + sqrt(2)*(2*c*d*e*f - (e + sqrt(-4*d*f + e**2))*(a*f**2 + c*(-d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*f**2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) - sqrt(2)*(2*c*d*e*f - (e - sqrt(-4*d*f + e**2))*(a*f**2 + c*(-d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*f**2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_54():
    f = sqrt(a + c*x**2)/(d + e*x + f*x**2)
    F = sqrt(c)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/f - sqrt(2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*f*sqrt(-4*d*f + e**2)) + sqrt(2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*f*sqrt(-4*d*f + e**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_55():
    f = sqrt(a + c*x**2)/(x*(d + e*x + f*x**2))
    F = -sqrt(a)*atanh(sqrt(a + c*x**2)/sqrt(a))/d - sqrt(2)*(2*a*e*f + (e + sqrt(-4*d*f + e**2))*(-a*f + c*d))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*d*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) + sqrt(2)*(2*a*e*f + (e - sqrt(-4*d*f + e**2))*(-a*f + c*d))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*d*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_56():
    f = sqrt(a + c*x**2)/(x**2*(d + e*x + f*x**2))
    F = sqrt(a)*e*atanh(sqrt(a + c*x**2)/sqrt(a))/d**2 - sqrt(a + c*x**2)/(d*x) + sqrt(2)*f*(a*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)) + 2*c*d**2)*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*d**2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) - sqrt(2)*f*(a*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)) + 2*c*d**2)*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*d**2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_57():
    f = sqrt(a + c*x**2)/(x**3*(d + e*x + f*x**2))
    F = -sqrt(a)*(-d*f + e**2)*atanh(sqrt(a + c*x**2)/sqrt(a))/d**3 - sqrt(a + c*x**2)/(2*d*x**2) + e*sqrt(a + c*x**2)/(d**2*x) - sqrt(2)*f*(a*(-3*d*e*f + d*f*sqrt(-4*d*f + e**2) + e**3 - e**2*sqrt(-4*d*f + e**2)) + c*d**2*(e - sqrt(-4*d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*d**3*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) + sqrt(2)*f*(a*(-3*d*e*f - d*f*sqrt(-4*d*f + e**2) + e**3 + e**2*sqrt(-4*d*f + e**2)) + c*d**2*(e + sqrt(-4*d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*d**3*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) - c*atanh(sqrt(a + c*x**2)/sqrt(a))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_58():
    f = x**2*(a + c*x**2)**(sympy.S(3)/2)/(d + e*x + f*x**2)
    F = -(a + c*x**2)**(sympy.S(3)/2)*(4*e - 3*f*x)/(12*f**2) - sqrt(a + c*x**2)*(8*e*(a*f**2 + c*(-2*d*f + e**2)) - f*x*(3*a*f**2 + 4*c*(-d*f + e**2)))/(8*f**4) + sqrt(2)*(a**2*f**4*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)) + 2*a*c*f**2*(2*d**2*f**2 - 4*d*e**2*f - 2*d*e*f*sqrt(-4*d*f + e**2) + e**4 + e**3*sqrt(-4*d*f + e**2)) + c**2*(-2*d**3*f**3 + 9*d**2*e**2*f**2 + 3*d**2*e*f**2*sqrt(-4*d*f + e**2) - 6*d*e**4*f - 4*d*e**3*f*sqrt(-4*d*f + e**2) + e**6 + e**5*sqrt(-4*d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*f**5*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) - sqrt(2)*(a**2*f**4*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)) + 2*a*c*f**2*(2*d**2*f**2 - 4*d*e**2*f + 2*d*e*f*sqrt(-4*d*f + e**2) + e**4 - e**3*sqrt(-4*d*f + e**2)) + c**2*(-2*d**3*f**3 + 9*d**2*e**2*f**2 - 3*d**2*e*f**2*sqrt(-4*d*f + e**2) - 6*d*e**4*f + 4*d*e**3*f*sqrt(-4*d*f + e**2) + e**6 - e**5*sqrt(-4*d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*f**5*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) + (3*a**2*f**4 + 12*a*c*f**2*(-d*f + e**2) + 8*c**2*(d**2*f**2 - 3*d*e**2*f + e**4))*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(8*sqrt(c)*f**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_59():
    f = x*(a + c*x**2)**(sympy.S(3)/2)/(d + e*x + f*x**2)
    F = -sqrt(c)*e*(3*a*f**2 + 2*c*(-2*d*f + e**2))*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(2*f**4) + (a + c*x**2)**(sympy.S(3)/2)/(3*f) + sqrt(a + c*x**2)*(2*a*f**2 - c*e*f*x + 2*c*(-d*f + e**2))/(2*f**3) + sqrt(2)*(2*c*d*e*f*(2*a*f**2 + c*(-2*d*f + e**2)) - (e + sqrt(-4*d*f + e**2))*(a**2*f**4 + 2*a*c*f**2*(-d*f + e**2) + c**2*(d**2*f**2 - 3*d*e**2*f + e**4)))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*f**4*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) - sqrt(2)*(2*c*d*e*f*(2*a*f**2 + c*(-2*d*f + e**2)) - (e - sqrt(-4*d*f + e**2))*(a**2*f**4 + 2*a*c*f**2*(-d*f + e**2) + c**2*(d**2*f**2 - 3*d*e**2*f + e**4)))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*f**4*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_60():
    f = (a + c*x**2)**(sympy.S(3)/2)/(d + e*x + f*x**2)
    F = sqrt(c)*(3*a*f**2 + 2*c*(-d*f + e**2))*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(2*f**3) - c*sqrt(a + c*x**2)*(2*e - f*x)/(2*f**2) + sqrt(2)*(c*e*(e + sqrt(-4*d*f + e**2))*(2*a*f**2 + c*(-2*d*f + e**2)) - 2*f*(-a**2*f**3 + 2*a*c*d*f**2 + c**2*d*(-d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*f**3*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) - sqrt(2)*(c*e*(e - sqrt(-4*d*f + e**2))*(2*a*f**2 + c*(-2*d*f + e**2)) - 2*f*(-a**2*f**3 + 2*a*c*d*f**2 + c**2*d*(-d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*f**3*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_61():
    f = (a + c*x**2)**(sympy.S(3)/2)/(x*(d + e*x + f*x**2))
    F = -a**(sympy.S(3)/2)*atanh(sqrt(a + c*x**2)/sqrt(a))/d + a*sqrt(a + c*x**2)/d - c**(sympy.S(3)/2)*e*atanh(sqrt(c)*x/sqrt(a + c*x**2))/f**2 + sqrt(a + c*x**2)*(-a*f + c*d)/(d*f) + sqrt(2)*(2*e*f*(-a**2*f**2 + c**2*d**2) - (e + sqrt(-4*d*f + e**2))*(c**2*d*e**2 - f*(-a*f + c*d)**2))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*d*f**2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) - sqrt(2)*(2*e*f*(-a**2*f**2 + c**2*d**2) - (e - sqrt(-4*d*f + e**2))*(c**2*d*e**2 - f*(-a*f + c*d)**2))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*d*f**2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_62():
    f = (a + c*x**2)**(sympy.S(3)/2)/(x**2*(d + e*x + f*x**2))
    F = a**(sympy.S(3)/2)*e*atanh(sqrt(a + c*x**2)/sqrt(a))/d**2 + 3*a*sqrt(c)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(2*d) - a*e*sqrt(a + c*x**2)/d**2 + sqrt(c)*(-3*a*f + 2*c*d)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(2*d*f) + 3*c*x*sqrt(a + c*x**2)/(2*d) - (a + c*x**2)**(sympy.S(3)/2)/(d*x) + sqrt(a + c*x**2)*(2*a*e - c*d*x)/(2*d**2) + sqrt(2)*(a**2*f**2*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)) + 4*a*c*d**2*f**2 + c**2*d**2*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*d**2*f*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) - sqrt(2)*(a**2*f**2*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)) + 4*a*c*d**2*f**2 + c**2*d**2*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*d**2*f*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_63():
    f = (a + c*x**2)**(sympy.S(3)/2)/(x**3*(d + e*x + f*x**2))
    F = -a**(sympy.S(3)/2)*(-d*f + e**2)*atanh(sqrt(a + c*x**2)/sqrt(a))/d**3 - 3*sqrt(a)*c*atanh(sqrt(a + c*x**2)/sqrt(a))/(2*d) + a*sqrt(a + c*x**2)*(-d*f + e**2)/d**3 + 3*c*sqrt(a + c*x**2)/(2*d) - 3*c*e*x*sqrt(a + c*x**2)/(2*d**2) - (a + c*x**2)**(sympy.S(3)/2)/(2*d*x**2) + e*(a + c*x**2)**(sympy.S(3)/2)/(d**2*x) - sqrt(a + c*x**2)*(2*a*(-d*f + e**2) + 2*c*d**2 - c*d*e*x)/(2*d**3) - sqrt(2)*(a**2*f*(-3*d*e*f + d*f*sqrt(-4*d*f + e**2) + e**3 - e**2*sqrt(-4*d*f + e**2)) + 2*a*c*d**2*f*(e - sqrt(-4*d*f + e**2)) + c**2*d**3*(e + sqrt(-4*d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*d**3*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) + sqrt(2)*(a**2*f*(-3*d*e*f - d*f*sqrt(-4*d*f + e**2) + e**3 + e**2*sqrt(-4*d*f + e**2)) + 2*a*c*d**2*f*(e + sqrt(-4*d*f + e**2)) + c**2*d**3*(e - sqrt(-4*d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*d**3*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_64():
    f = x**3/(sqrt(a + c*x**2)*(d + e*x + f*x**2))
    F = sqrt(2)*(2*d*e*f - (e + sqrt(-4*d*f + e**2))*(-d*f + e**2))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*f**2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) - sqrt(2)*(2*d*e*f - (e - sqrt(-4*d*f + e**2))*(-d*f + e**2))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*f**2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) + sqrt(a + c*x**2)/(c*f) - e*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(sqrt(c)*f**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_65():
    f = x**2/(sqrt(a + c*x**2)*(d + e*x + f*x**2))
    F = -sqrt(2)*(2*d*f - e*(e + sqrt(-4*d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*f*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) - sqrt(2)*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*f*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) + atanh(sqrt(c)*x/sqrt(a + c*x**2))/(sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_66():
    f = x/(sqrt(a + c*x**2)*(d + e*x + f*x**2))
    F = sqrt(2)*(e - sqrt(-4*d*f + e**2))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) - sqrt(2)*(e + sqrt(-4*d*f + e**2))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_67():
    f = 1/(sqrt(a + c*x**2)*(d + e*x + f*x**2))
    F = sqrt(2)*f*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) - sqrt(2)*f*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_68():
    f = 1/(x*sqrt(a + c*x**2)*(d + e*x + f*x**2))
    F = -sqrt(2)*f*(e - sqrt(-4*d*f + e**2))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*d*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) + sqrt(2)*f*(e + sqrt(-4*d*f + e**2))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*d*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) - atanh(sqrt(a + c*x**2)/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_69():
    f = 1/(x**2*sqrt(a + c*x**2)*(d + e*x + f*x**2))
    F = sqrt(2)*f*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*d**2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) - sqrt(2)*f*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*d**2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) - sqrt(a + c*x**2)/(a*d*x) + e*atanh(sqrt(a + c*x**2)/sqrt(a))/(sqrt(a)*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_70():
    f = 1/(x**3*sqrt(a + c*x**2)*(d + e*x + f*x**2))
    F = -sqrt(2)*f*(-4*d*e*f + 2*e**3 - (e + sqrt(-4*d*f + e**2))*(-d*f + e**2))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*d**3*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) + sqrt(2)*f*(-4*d*e*f + 2*e**3 - (e - sqrt(-4*d*f + e**2))*(-d*f + e**2))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*d**3*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)) - sqrt(a + c*x**2)/(2*a*d*x**2) + e*sqrt(a + c*x**2)/(a*d**2*x) - (-d*f + e**2)*atanh(sqrt(a + c*x**2)/sqrt(a))/(sqrt(a)*d**3) + c*atanh(sqrt(a + c*x**2)/sqrt(a))/(2*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_71():
    f = x**3/((a + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2))
    F = sqrt(2)*(2*a*d*e*f - (e + sqrt(-4*d*f + e**2))*(a*(-d*f + e**2) + c*d**2))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)*(a*c*e**2 + (-a*f + c*d)**2)) - sqrt(2)*(2*a*d*e*f - (e - sqrt(-4*d*f + e**2))*(a*(-d*f + e**2) + c*d**2))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)*(a*c*e**2 + (-a*f + c*d)**2)) - 1/(c*f*sqrt(a + c*x**2)) - e*x/(a*f**2*sqrt(a + c*x**2)) + (a*f*(a*(-d*f + e**2) + c*d**2) + c*e*x*(a*(-2*d*f + e**2) + c*d**2))/(a*f**2*sqrt(a + c*x**2)*(a*c*e**2 + (-a*f + c*d)**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_72():
    f = x**2/((a + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2))
    F = sqrt(2)*f*(a*e*(e + sqrt(-4*d*f + e**2)) + 2*d*(-a*f + c*d))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)*(a*c*e**2 + (-a*f + c*d)**2)) - sqrt(2)*f*(a*e*(e - sqrt(-4*d*f + e**2)) + 2*d*(-a*f + c*d))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)*(a*c*e**2 + (-a*f + c*d)**2)) - (a*e + x*(-a*f + c*d))/(sqrt(a + c*x**2)*(a*c*e**2 + (-a*f + c*d)**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_73():
    f = x/((a + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2))
    F = -sqrt(2)*f*(2*c*d*e - (e + sqrt(-4*d*f + e**2))*(-a*f + c*d))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)*(a*c*e**2 + (-a*f + c*d)**2)) + sqrt(2)*f*(2*c*d*e - (e - sqrt(-4*d*f + e**2))*(-a*f + c*d))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)*(a*c*e**2 + (-a*f + c*d)**2)) - (-a*f + c*d - c*e*x)/(sqrt(a + c*x**2)*(a*c*e**2 + (-a*f + c*d)**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_74():
    f = 1/((a + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2))
    F = sqrt(2)*f*(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)*(a*c*e**2 + (-a*f + c*d)**2)) - sqrt(2)*f*(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)*(a*c*e**2 + (-a*f + c*d)**2)) + c*(a*e + x*(-a*f + c*d))/(a*sqrt(a + c*x**2)*(a*c*e**2 + (-a*f + c*d)**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_75():
    f = 1/(x*(a + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2))
    F = -sqrt(2)*f*(2*e*(a*f**2 + c*(-2*d*f + e**2)) - (e + sqrt(-4*d*f + e**2))*(a*f**2 + c*(-d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*d*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)*(a*c*e**2 + (-a*f + c*d)**2)) + sqrt(2)*f*(2*e*(a*f**2 + c*(-2*d*f + e**2)) - (e - sqrt(-4*d*f + e**2))*(a*f**2 + c*(-d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*d*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)*(a*c*e**2 + (-a*f + c*d)**2)) - (a*(a*f**2 + c*(-d*f + e**2)) + c**2*d*e*x)/(a*d*sqrt(a + c*x**2)*(a*c*e**2 + (-a*f + c*d)**2)) + 1/(a*d*sqrt(a + c*x**2)) - atanh(sqrt(a + c*x**2)/sqrt(a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_76():
    f = 1/(x**2*(a + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2))
    F = -sqrt(2)*f*(-2*a*f**2*(-d*f + e**2) - 2*c*(d**2*f**2 - 3*d*e**2*f + e**4) + e*(e + sqrt(-4*d*f + e**2))*(a*f**2 + c*(-2*d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e + sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))))/(2*d**2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)*(a*c*e**2 + (-a*f + c*d)**2)) + sqrt(2)*f*(-2*a*f**2*(-d*f + e**2) - 2*c*(d**2*f**2 - 3*d*e**2*f + e**4) + e*(e - sqrt(-4*d*f + e**2))*(a*f**2 + c*(-2*d*f + e**2)))*atanh(sqrt(2)*(2*a*f - c*x*(e - sqrt(-4*d*f + e**2)))/(2*sqrt(a + c*x**2)*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))))/(2*d**2*sqrt(2*a*f**2 + c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)))*sqrt(-4*d*f + e**2)*(a*c*e**2 + (-a*f + c*d)**2)) - 1/(a*d*x*sqrt(a + c*x**2)) - e/(a*d**2*sqrt(a + c*x**2)) + (a*e*(a*f**2 + c*(-2*d*f + e**2)) + c*d*x*(a*f**2 + c*(-d*f + e**2)))/(a*d**2*sqrt(a + c*x**2)*(a*c*e**2 + (-a*f + c*d)**2)) - 2*c*x/(a**2*d*sqrt(a + c*x**2)) + e*atanh(sqrt(a + c*x**2)/sqrt(a))/(a**(sympy.S(3)/2)*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_77():
    f = x**3*sqrt(a + b*x + c*x**2)/(d - f*x**2)
    F = b*(b + 2*c*x)*sqrt(a + b*x + c*x**2)/(8*c**2*f) - b*d*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*sqrt(c)*f**2) - b*(-4*a*c + b**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(5)/2)*f) - d*sqrt(a + b*x + c*x**2)/f**2 - d*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*f**(sympy.S(5)/2)) + d*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*f**(sympy.S(5)/2)) - (a + b*x + c*x**2)**(sympy.S(3)/2)/(3*c*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_78():
    f = x**2*sqrt(a + b*x + c*x**2)/(d - f*x**2)
    F = sqrt(d)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*f**2) + sqrt(d)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*f**2) - (b + 2*c*x)*sqrt(a + b*x + c*x**2)/(4*c*f) - (4*a*c*f - b**2*f + 8*c**2*d)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*c**(sympy.S(3)/2)*f**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_79():
    f = x*sqrt(a + b*x + c*x**2)/(d - f*x**2)
    F = -b*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*sqrt(c)*f) - sqrt(a + b*x + c*x**2)/f - sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*f**(sympy.S(3)/2)) + sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*f**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_80():
    f = sqrt(a + b*x + c*x**2)/(d - f*x**2)
    F = -sqrt(c)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/f + sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*sqrt(d)*f) + sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_81():
    f = sqrt(a + b*x + c*x**2)/(x*(d - f*x**2))
    F = -sqrt(a)*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/d - sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*d*sqrt(f)) + sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*d*sqrt(f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_82():
    f = sqrt(a + b*x + c*x**2)/(x**2*(d - f*x**2))
    F = -sqrt(a + b*x + c*x**2)/(d*x) + sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*d**(sympy.S(3)/2)) + sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*d**(sympy.S(3)/2)) - b*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_83():
    f = sqrt(a + b*x + c*x**2)/(x**3*(d - f*x**2))
    F = -sqrt(a)*f*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/d**2 - sqrt(f)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*d**2) + sqrt(f)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*d**2) - (2*a + b*x)*sqrt(a + b*x + c*x**2)/(4*a*d*x**2) + (-4*a*c + b**2)*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(8*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_84():
    f = x**3*(a + b*x + c*x**2)**(sympy.S(3)/2)/(d - f*x**2)
    F = b*(b + 2*c*x)*(a + b*x + c*x**2)**(sympy.S(3)/2)/(16*c**2*f) - 3*b*(b + 2*c*x)*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)/(128*c**3*f) - b*d*(12*a*c*f - b**2*f + 24*c**2*d)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(3)/2)*f**3) + 3*b*(-4*a*c + b**2)**2*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(256*c**(sympy.S(7)/2)*f) - d*(a + b*x + c*x**2)**(sympy.S(3)/2)/(3*f**2) - d*(a*f - b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*f**(sympy.S(7)/2)) + d*(a*f + b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*f**(sympy.S(7)/2)) - d*sqrt(a + b*x + c*x**2)*(8*a*c*f + b**2*f + 2*b*c*f*x + 8*c**2*d)/(8*c*f**3) - (a + b*x + c*x**2)**(sympy.S(5)/2)/(5*c*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_85():
    f = x**2*(a + b*x + c*x**2)**(sympy.S(3)/2)/(d - f*x**2)
    F = sqrt(d)*(a*f - b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*f**3) + sqrt(d)*(a*f + b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*f**3) - (b + 2*c*x)*(a + b*x + c*x**2)**(sympy.S(3)/2)/(8*c*f) - (b*(12*a*c*f - 3*b**2*f + 80*c**2*d) + 2*c*x*(12*a*c*f - 3*b**2*f + 16*c**2*d))*sqrt(a + b*x + c*x**2)/(64*c**2*f**2) - (-24*a*b**2*c*f**2 + 192*a*c**3*d*f + 3*b**4*f**2 + 128*c**4*d**2 + 48*c**2*f*(a**2*f + b**2*d))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(128*c**(sympy.S(5)/2)*f**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_86():
    f = x*(a + b*x + c*x**2)**(sympy.S(3)/2)/(d - f*x**2)
    F = -b*(12*a*c*f - b**2*f + 24*c**2*d)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(3)/2)*f**2) - (a + b*x + c*x**2)**(sympy.S(3)/2)/(3*f) - (a*f - b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*f**(sympy.S(5)/2)) + (a*f + b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*f**(sympy.S(5)/2)) - sqrt(a + b*x + c*x**2)*(8*a*c*f + b**2*f + 2*b*c*f*x + 8*c**2*d)/(8*c*f**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_87():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)/(d - f*x**2)
    F = -(5*b + 2*c*x)*sqrt(a + b*x + c*x**2)/(4*f) + (a*f - b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*sqrt(d)*f**2) + (a*f + b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*sqrt(d)*f**2) - (12*a*c*f + 3*b**2*f + 8*c**2*d)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*sqrt(c)*f**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_88():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)/(x*(d - f*x**2))
    F = -a**(sympy.S(3)/2)*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/d - b*(-12*a*c + b**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(3)/2)*d) - b*(12*a*c*f - b**2*f + 24*c**2*d)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(3)/2)*d*f) - (a*f - b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*d*f**(sympy.S(3)/2)) + (a*f + b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*d*f**(sympy.S(3)/2)) + sqrt(a + b*x + c*x**2)*(8*a*c + b**2 + 2*b*c*x)/(8*c*d) - sqrt(a + b*x + c*x**2)*(8*a*c*f + b**2*f + 2*b*c*f*x + 8*c**2*d)/(8*c*d*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_89():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)/(x**2*(d - f*x**2))
    F = -3*sqrt(a)*b*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(2*d) - (5*b + 2*c*x)*sqrt(a + b*x + c*x**2)/(4*d) + (9*b + 6*c*x)*sqrt(a + b*x + c*x**2)/(4*d) - (a + b*x + c*x**2)**(sympy.S(3)/2)/(d*x) + (a*f - b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*d**(sympy.S(3)/2)*f) + (a*f + b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*d**(sympy.S(3)/2)*f) + (12*a*c + 3*b**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*sqrt(c)*d) - (12*a*c*f + 3*b**2*f + 8*c**2*d)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*sqrt(c)*d*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_90():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)/(x**3*(d - f*x**2))
    F = -a**(sympy.S(3)/2)*f*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/d**2 + 3*b*sqrt(c)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*d) - b*f*(-12*a*c + b**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(3)/2)*d**2) - b*(12*a*c*f - b**2*f + 24*c**2*d)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(3)/2)*d**2) - (3*b - 6*c*x)*sqrt(a + b*x + c*x**2)/(4*d*x) - (a + b*x + c*x**2)**(sympy.S(3)/2)/(2*d*x**2) - (a*f - b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*d**2*sqrt(f)) + (a*f + b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*d**2*sqrt(f)) + f*sqrt(a + b*x + c*x**2)*(8*a*c + b**2 + 2*b*c*x)/(8*c*d**2) - sqrt(a + b*x + c*x**2)*(8*a*c*f + b**2*f + 2*b*c*f*x + 8*c**2*d)/(8*c*d**2) - (12*a*c + 3*b**2)*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(8*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_91():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)/(1 - x**2)
    F = (-5*b/4 - c*x/2)*sqrt(a + b*x + c*x**2) - (a - b + c)**(sympy.S(3)/2)*atanh((2*a - b + x*(b - 2*c))/(2*sqrt(a - b + c)*sqrt(a + b*x + c*x**2)))/2 + (a + b + c)**(sympy.S(3)/2)*atanh((2*a + b + x*(b + 2*c))/(2*sqrt(a + b + c)*sqrt(a + b*x + c*x**2)))/2 - (12*a*c + 3*b**2 + 8*c**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_92():
    f = sqrt(x**2 - x - 1)/(1 - x**2)
    F = -atan((3 - x)/(2*sqrt(x**2 - x - 1)))/2 + atanh((1 - 2*x)/(2*sqrt(x**2 - x - 1))) + atanh((3*x + 1)/(2*sqrt(x**2 - x - 1)))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_93():
    f = (x**2 + x)**(sympy.S(3)/2)/(x**2 + 1)
    F = (x/2 + sympy.S(5)/4)*sqrt(x**2 + x) + sqrt(1 + sqrt(2))*atan((-x + 1 + sqrt(2))/(sqrt(2 + 2*sqrt(2))*sqrt(x**2 + x))) - 5*atanh(x/sqrt(x**2 + x))/4 - sqrt(-1 + sqrt(2))*atanh((-x - sqrt(2) + 1)/(sqrt(-2 + 2*sqrt(2))*sqrt(x**2 + x)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_94():
    f = x**4/((d - f*x**2)*sqrt(a + b*x + c*x**2))
    F = 3*b*sqrt(a + b*x + c*x**2)/(4*c**2*f) + d**(sympy.S(3)/2)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*f**2*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)) + d**(sympy.S(3)/2)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*f**2*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)) - x*sqrt(a + b*x + c*x**2)/(2*c*f) - d*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(sqrt(c)*f**2) - (-4*a*c + 3*b**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_95():
    f = x**3/((d - f*x**2)*sqrt(a + b*x + c*x**2))
    F = b*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*c**(sympy.S(3)/2)*f) + d*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*f**(sympy.S(3)/2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)) - d*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*f**(sympy.S(3)/2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)) - sqrt(a + b*x + c*x**2)/(c*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_96():
    f = x**2/((d - f*x**2)*sqrt(a + b*x + c*x**2))
    F = sqrt(d)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*f*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)) + sqrt(d)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*f*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)) - atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_97():
    f = x/((d - f*x**2)*sqrt(a + b*x + c*x**2))
    F = atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*sqrt(f)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)) - atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*sqrt(f)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_98():
    f = 1/((d - f*x**2)*sqrt(a + b*x + c*x**2))
    F = atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*sqrt(d)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)) + atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*sqrt(d)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_99():
    f = 1/(x*(d - f*x**2)*sqrt(a + b*x + c*x**2))
    F = sqrt(f)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*d*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)) - sqrt(f)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*d*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)) - atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_100():
    f = 1/(x**2*(d - f*x**2)*sqrt(a + b*x + c*x**2))
    F = f*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*d**(sympy.S(3)/2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)) + f*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*d**(sympy.S(3)/2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)) - sqrt(a + b*x + c*x**2)/(a*d*x) + b*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(2*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_101():
    f = 1/(x**3*(d - f*x**2)*sqrt(a + b*x + c*x**2))
    F = f**(sympy.S(3)/2)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*d**2*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)) - f**(sympy.S(3)/2)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*d**2*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)) - sqrt(a + b*x + c*x**2)/(2*a*d*x**2) + 3*b*sqrt(a + b*x + c*x**2)/(4*a**2*d*x) - f*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(sqrt(a)*d**2) - (-4*a*c + 3*b**2)*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(8*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_102():
    f = x**4/((d - f*x**2)*(a + b*x + c*x**2)**(sympy.S(3)/2))
    F = 2*b*sqrt(a + b*x + c*x**2)/(c*f*(-4*a*c + b**2)) + d**(sympy.S(3)/2)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*f*(a*f + b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)) + d**(sympy.S(3)/2)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*f*(a*f - b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)) - 2*d**2*(b*(b**2*f - c*(3*a*f + c*d)) - c*x*(2*a*c*f - b**2*f + 2*c**2*d))/(f**2*(-4*a*c + b**2)*(b**2*d*f - (a*f + c*d)**2)*sqrt(a + b*x + c*x**2)) + 2*d*(b + 2*c*x)/(f**2*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)) - 2*x*(2*a + b*x)/(f*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)) - atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_103():
    f = x**3/((d - f*x**2)*(a + b*x + c*x**2)**(sympy.S(3)/2))
    F = -2*d*(a*(2*a*c*f - b**2*f + 2*c**2*d) + b*c*x*(-a*f + c*d))/(f*(-4*a*c + b**2)*(b**2*d*f - (a*f + c*d)**2)*sqrt(a + b*x + c*x**2)) + d*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*sqrt(f)*(a*f + b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)) - d*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*sqrt(f)*(a*f - b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)) - (4*a + 2*b*x)/(f*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_104():
    f = x**2/((d - f*x**2)*(a + b*x + c*x**2)**(sympy.S(3)/2))
    F = sqrt(d)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*(a*f + b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)) + sqrt(d)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*(a*f - b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)) + (2*a*b*(-a*f + c*d) + 2*c*x*(-2*a*(a*f + c*d) + b**2*d))/((-4*a*c + b**2)*(b**2*d*f - (a*f + c*d)**2)*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_105():
    f = x/((d - f*x**2)*(a + b*x + c*x**2)**(sympy.S(3)/2))
    F = sqrt(f)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*(a*f + b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)) - sqrt(f)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*(a*f - b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)) - (2*a*(2*a*c*f - b**2*f + 2*c**2*d) + 2*b*c*x*(-a*f + c*d))/((-4*a*c + b**2)*(b**2*d*f - (a*f + c*d)**2)*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_106():
    f = 1/((d - f*x**2)*(a + b*x + c*x**2)**(sympy.S(3)/2))
    F = -(2*b*(b**2*f - c*(3*a*f + c*d)) - 2*c*x*(2*a*c*f - b**2*f + 2*c**2*d))/((-4*a*c + b**2)*(b**2*d*f - (a*f + c*d)**2)*sqrt(a + b*x + c*x**2)) + f*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*sqrt(d)*(a*f + b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)) + f*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*sqrt(d)*(a*f - b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_107():
    f = 1/(x*(d - f*x**2)*(a + b*x + c*x**2)**(sympy.S(3)/2))
    F = f**(sympy.S(3)/2)*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*d*(a*f + b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)) - f**(sympy.S(3)/2)*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*d*(a*f - b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)) - 2*f*(a*(2*a*c*f - b**2*f + 2*c**2*d) + b*c*x*(-a*f + c*d))/(d*(-4*a*c + b**2)*(b**2*d*f - (a*f + c*d)**2)*sqrt(a + b*x + c*x**2)) + (-4*a*c + 2*b**2 + 2*b*c*x)/(a*d*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)) - atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_108():
    f = 1/(x**2*(d - f*x**2)*(a + b*x + c*x**2)**(sympy.S(3)/2))
    F = -2*f*(b*(b**2*f - c*(3*a*f + c*d)) - c*x*(2*a*c*f - b**2*f + 2*c**2*d))/(d*(-4*a*c + b**2)*(b**2*d*f - (a*f + c*d)**2)*sqrt(a + b*x + c*x**2)) + f**2*atanh((2*a*sqrt(f) + b*sqrt(d) + x*(b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f + b*sqrt(d)*sqrt(f) + c*d)))/(2*d**(sympy.S(3)/2)*(a*f + b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)) + f**2*atanh((-2*a*sqrt(f) + b*sqrt(d) + x*(-b*sqrt(f) + 2*c*sqrt(d)))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*f - b*sqrt(d)*sqrt(f) + c*d)))/(2*d**(sympy.S(3)/2)*(a*f - b*sqrt(d)*sqrt(f) + c*d)**(sympy.S(3)/2)) + (-4*a*c + 2*b**2 + 2*b*c*x)/(a*d*x*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)) - (-8*a*c + 3*b**2)*sqrt(a + b*x + c*x**2)/(a**2*d*x*(-4*a*c + b**2)) + 3*b*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(2*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_109():
    f = x**2*sqrt(a + b*x + c*x**2)/(d + e*x + f*x**2)
    F = sqrt(2)*(c*(2*d**2*f**2 - 4*d*e**2*f - 2*d*e*f*sqrt(-4*d*f + e**2) + e**4 + e**3*sqrt(-4*d*f + e**2)) + f*(a*f*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)) - b*(-3*d*e*f - d*f*sqrt(-4*d*f + e**2) + e**3 + e**2*sqrt(-4*d*f + e**2))))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*f**3*sqrt(c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)) + f*(2*a*f - b*(e + sqrt(-4*d*f + e**2))))*sqrt(-4*d*f + e**2)) - sqrt(2)*(c*(2*d**2*f**2 - 4*d*e**2*f + 2*d*e*f*sqrt(-4*d*f + e**2) + e**4 - e**3*sqrt(-4*d*f + e**2)) + f*(a*f*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)) - b*(-3*d*e*f + d*f*sqrt(-4*d*f + e**2) + e**3 - e**2*sqrt(-4*d*f + e**2))))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*f**3*sqrt(c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)) + f*(2*a*f - b*(e - sqrt(-4*d*f + e**2))))*sqrt(-4*d*f + e**2)) - sqrt(a + b*x + c*x**2)*(-b*f + 4*c*e - 2*c*f*x)/(4*c*f**2) - (b**2*f**2 - 8*c**2*(-d*f + e**2) + 4*c*f*(-a*f + b*e))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*c**(sympy.S(3)/2)*f**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_110():
    f = x*sqrt(a + b*x + c*x**2)/(d + e*x + f*x**2)
    F = sqrt(a + b*x + c*x**2)/f - sqrt(2)*(2*d*f*(-b*f + c*e) + (e - sqrt(-4*d*f + e**2))*(-c*(-d*f + e**2) + f*(-a*f + b*e)))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*f**2*sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))) + sqrt(2)*(2*d*f*(-b*f + c*e) + (e + sqrt(-4*d*f + e**2))*(-c*(-d*f + e**2) + f*(-a*f + b*e)))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*f**2*sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))) - (-b*f + 2*c*e)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*sqrt(c)*f**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_111():
    f = sqrt(a + b*x + c*x**2)/(d + e*x + f*x**2)
    F = sqrt(c)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/f - sqrt(2)*sqrt(c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)) + f*(2*a*f - b*(e - sqrt(-4*d*f + e**2))))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*f*sqrt(-4*d*f + e**2)) + sqrt(2)*sqrt(c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)) + f*(2*a*f - b*(e + sqrt(-4*d*f + e**2))))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*f*sqrt(-4*d*f + e**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_112():
    f = sqrt(a + b*x + c*x**2)/(x*(d + e*x + f*x**2))
    F = -sqrt(a)*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/d - sqrt(2)*(c*d*(e + sqrt(-4*d*f + e**2)) - f*(-a*(e - sqrt(-4*d*f + e**2)) + 2*b*d))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*d*sqrt(c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)) + f*(2*a*f - b*(e + sqrt(-4*d*f + e**2))))*sqrt(-4*d*f + e**2)) + sqrt(2)*(c*d*(e - sqrt(-4*d*f + e**2)) - f*(-a*(e + sqrt(-4*d*f + e**2)) + 2*b*d))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*d*sqrt(c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)) + f*(2*a*f - b*(e - sqrt(-4*d*f + e**2))))*sqrt(-4*d*f + e**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_113():
    f = sqrt(a + b*x + c*x**2)/(x**2*(d + e*x + f*x**2))
    F = sqrt(a)*e*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/d**2 - b*e*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*sqrt(c)*d**2) + sqrt(c)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/d - sqrt(a + b*x + c*x**2)/(d*x) + sqrt(2)*f*(a*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)) - b*d*(e - sqrt(-4*d*f + e**2)) + 2*c*d**2)*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*d**2*sqrt(c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)) + f*(2*a*f - b*(e + sqrt(-4*d*f + e**2))))*sqrt(-4*d*f + e**2)) - sqrt(2)*f*(a*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)) - b*d*(e + sqrt(-4*d*f + e**2)) + 2*c*d**2)*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*d**2*sqrt(c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)) + f*(2*a*f - b*(e - sqrt(-4*d*f + e**2))))*sqrt(-4*d*f + e**2)) - (-b*e + 2*c*d)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*sqrt(c)*d**2) - b*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_114():
    f = x**3/(sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2))
    F = -b*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*c**(sympy.S(3)/2)*f) - sqrt(2)*(2*d*e*f - (e - sqrt(-4*d*f + e**2))*(-d*f + e**2))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*f**2*sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))) + sqrt(2)*(2*d*e*f - (e + sqrt(-4*d*f + e**2))*(-d*f + e**2))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*f**2*sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))) + sqrt(a + b*x + c*x**2)/(c*f) - e*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(sqrt(c)*f**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_115():
    f = x**2/(sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2))
    F = -sqrt(2)*(2*d*f - e*(e + sqrt(-4*d*f + e**2)))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*f*sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))) - sqrt(2)*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*f*sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))) + atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_116():
    f = x/(sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2))
    F = sqrt(2)*(e - sqrt(-4*d*f + e**2))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))) - sqrt(2)*(e + sqrt(-4*d*f + e**2))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_117():
    f = 1/(sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2))
    F = sqrt(2)*f*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))) - sqrt(2)*f*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_118():
    f = 1/(x*sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2))
    F = -sqrt(2)*f*(e - sqrt(-4*d*f + e**2))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*d*sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))) + sqrt(2)*f*(e + sqrt(-4*d*f + e**2))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*d*sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))) - atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_119():
    f = 1/(x**2*sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2))
    F = sqrt(2)*f*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*d**2*sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))) - sqrt(2)*f*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*d**2*sqrt(c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)) + f*(2*a*f - b*(e - sqrt(-4*d*f + e**2))))*sqrt(-4*d*f + e**2)) - sqrt(a + b*x + c*x**2)/(a*d*x) + e*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(sqrt(a)*d**2) + b*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(2*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_120():
    f = 1/(x**3*sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2))
    F = sqrt(2)*f*(-4*d*e*f + 2*e**3 - (e - sqrt(-4*d*f + e**2))*(-d*f + e**2))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*d**3*sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))) - sqrt(2)*f*(-4*d*e*f + 2*e**3 - (e + sqrt(-4*d*f + e**2))*(-d*f + e**2))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*d**3*sqrt(-4*d*f + e**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))) - sqrt(a + b*x + c*x**2)/(2*a*d*x**2) + e*sqrt(a + b*x + c*x**2)/(a*d**2*x) + 3*b*sqrt(a + b*x + c*x**2)/(4*a**2*d*x) - (-d*f + e**2)*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(sqrt(a)*d**3) - b*e*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(2*a**(sympy.S(3)/2)*d**2) - (-4*a*c + 3*b**2)*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(8*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_121():
    f = x**3/((a + b*x + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2))
    F = 2*e*(b + 2*c*x)/(f**2*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)) + sqrt(2)*(2*d*f*(-a*e + b*d) + (e - sqrt(-4*d*f + e**2))*(a*(-d*f + e**2) - b*d*e + c*d**2))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*sqrt(-4*d*f + e**2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))) - sqrt(2)*(2*d*f*(-a*e + b*d) + (e + sqrt(-4*d*f + e**2))*(a*(-d*f + e**2) - b*d*e + c*d**2))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*sqrt(-4*d*f + e**2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))) + (4*a + 2*b*x)/(f*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)) + (2*c*d*e*(a*b*f - 2*a*c*e + b*c*d) + 2*c*x*(-d*e*(b**2*f + 2*c**2*d - c*(2*a*f + b*e)) + (-d*f + e**2)*(a*b*f - 2*a*c*e + b*c*d)) - 2*(b**2*f + 2*c**2*d - c*(2*a*f + b*e))*(a*d*f - a*e**2 + b*d*e))/(f**2*(-4*a*c + b**2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_122():
    f = x**2/((a + b*x + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2))
    F = -sqrt(2)*f*(2*d*(-a*f + c*d) - (e - sqrt(-4*d*f + e**2))*(-a*e + b*d))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*sqrt(-4*d*f + e**2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))) + sqrt(2)*f*(2*d*(-a*f + c*d) - (e + sqrt(-4*d*f + e**2))*(-a*e + b*d))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*sqrt(-4*d*f + e**2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))) - (2*a*(a*b*f - 2*a*c*e + b*c*d) + 2*c*x*(-a*b*e - 2*a*(-a*f + c*d) + b**2*d))/((-4*a*c + b**2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_123():
    f = x/((a + b*x + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2))
    F = sqrt(2)*f*(2*d*(-b*f + c*e) - (e - sqrt(-4*d*f + e**2))*(-a*f + c*d))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*sqrt(-4*d*f + e**2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))) - sqrt(2)*f*(2*d*(-b*f + c*e) - (e + sqrt(-4*d*f + e**2))*(-a*f + c*d))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*sqrt(-4*d*f + e**2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))) + (2*a*(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d) + 2*c*x*(a*b*f - 2*a*c*e + b*c*d))/((-4*a*c + b**2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_124():
    f = 1/((a + b*x + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2))
    F = sqrt(2)*f*(c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)) + f*(2*a*f - b*(e - sqrt(-4*d*f + e**2))))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*sqrt(c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)) + f*(2*a*f - b*(e + sqrt(-4*d*f + e**2))))*sqrt(-4*d*f + e**2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)) - sqrt(2)*f*(c*(-2*d*f + e**2 + e*sqrt(-4*d*f + e**2)) + f*(2*a*f - b*(e + sqrt(-4*d*f + e**2))))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*sqrt(c*(-2*d*f + e**2 - e*sqrt(-4*d*f + e**2)) + f*(2*a*f - b*(e - sqrt(-4*d*f + e**2))))*sqrt(-4*d*f + e**2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)) + (-4*a*c**2*e - 2*b**3*f + 2*b**2*c*e - 2*b*c*(-3*a*f + c*d) - 2*c*x*(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d))/((-4*a*c + b**2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_125():
    f = 1/(x*(a + b*x + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2))
    F = sqrt(2)*f*(2*c*(-2*d*e*f + e**3) - 2*f*(-a*e*f - b*d*f + b*e**2) + (e - sqrt(-4*d*f + e**2))*(-c*(-d*f + e**2) + f*(-a*f + b*e)))*atanh(sqrt(2)*(4*a*f - b*(e - sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e - sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*d*sqrt(-4*d*f + e**2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 - (-b*f + c*e)*sqrt(-4*d*f + e**2))) - sqrt(2)*f*(2*c*(-2*d*e*f + e**3) - 2*f*(-a*e*f - b*d*f + b*e**2) + (e + sqrt(-4*d*f + e**2))*(-c*(-d*f + e**2) + f*(-a*f + b*e)))*atanh(sqrt(2)*(4*a*f - b*(e + sqrt(-4*d*f + e**2)) + x*(2*b*f - 2*c*(e + sqrt(-4*d*f + e**2))))/(4*sqrt(a + b*x + c*x**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))))/(2*d*sqrt(-4*d*f + e**2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)*sqrt(2*a*f**2 - b*e*f - 2*c*d*f + c*e**2 + (-b*f + c*e)*sqrt(-4*d*f + e**2))) + (2*c*e*(2*a*c*e - b*(a*f + c*d)) + 2*c*x*(-b*c*(d*f + e**2) + b*f*(-a*f + b*e) + 2*c**2*d*e) + 2*(-a*f + b*e)*(b**2*f + 2*c**2*d - c*(2*a*f + b*e)))/(d*(-4*a*c + b**2)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2)*sqrt(a + b*x + c*x**2)) + (-4*a*c + 2*b**2 + 2*b*c*x)/(a*d*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)) - atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_126():
    f = x**4/(sqrt(-x**2 - 4*x - 3)*(2*x**2 + 4*x + 3))
    F = -x*sqrt(-x**2 - 4*x - 3)/4 + 5*sqrt(-x**2 - 4*x - 3)/2 + 11*asin(x + 2)/2 + sqrt(2)*atan(sqrt(2)*(-(x + 3)/sqrt(-x**2 - 4*x - 3) + 1)/2)/4 - sqrt(2)*atan(sqrt(2)*((x + 3)/sqrt(-x**2 - 4*x - 3) + 1)/2)/4 - 5*atanh(x/sqrt(-x**2 - 4*x - 3))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_127():
    f = x**3/(sqrt(-x**2 - 4*x - 3)*(2*x**2 + 4*x + 3))
    F = -sqrt(-x**2 - 4*x - 3)/2 - 2*asin(x + 2) + sqrt(2)*atan(sqrt(2)*(-(x + 3)/sqrt(-x**2 - 4*x - 3) + 1)/2)/4 - sqrt(2)*atan(sqrt(2)*((x + 3)/sqrt(-x**2 - 4*x - 3) + 1)/2)/4 + atanh(x/sqrt(-x**2 - 4*x - 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_128():
    f = x**2/(sqrt(-x**2 - 4*x - 3)*(2*x**2 + 4*x + 3))
    F = asin(x + 2)/2 - sqrt(2)*atan(sqrt(2)*(-(x + 3)/sqrt(-x**2 - 4*x - 3) + 1)/2)/2 + sqrt(2)*atan(sqrt(2)*((x + 3)/sqrt(-x**2 - 4*x - 3) + 1)/2)/2 - atanh(x/sqrt(-x**2 - 4*x - 3))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_129():
    f = 1/(sqrt(-x**2 - 4*x - 3)*(2*x**2 + 4*x + 3))
    F = -sqrt(2)*atan(sqrt(2)*(-(x + 3)/sqrt(-x**2 - 4*x - 3) + 1)/2)/3 + sqrt(2)*atan(sqrt(2)*((x + 3)/sqrt(-x**2 - 4*x - 3) + 1)/2)/3 + atanh(x/sqrt(-x**2 - 4*x - 3))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_130():
    f = 1/(x*sqrt(-x**2 - 4*x - 3)*(2*x**2 + 4*x + 3))
    F = sqrt(2)*atan(sqrt(2)*(-(x + 3)/sqrt(-x**2 - 4*x - 3) + 1)/2)/9 - sqrt(2)*atan(sqrt(2)*((x + 3)/sqrt(-x**2 - 4*x - 3) + 1)/2)/9 - sqrt(3)*atan(sqrt(3)*(2*x + 3)/(3*sqrt(-x**2 - 4*x - 3)))/9 - 4*atanh(x/sqrt(-x**2 - 4*x - 3))/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_131():
    f = 1/(x**2*sqrt(-x**2 - 4*x - 3)*(2*x**2 + 4*x + 3))
    F = 2*sqrt(2)*atan(sqrt(2)*(-(x + 3)/sqrt(-x**2 - 4*x - 3) + 1)/2)/27 - 2*sqrt(2)*atan(sqrt(2)*((x + 3)/sqrt(-x**2 - 4*x - 3) + 1)/2)/27 + 2*sqrt(3)*atan(sqrt(3)*(2*x + 3)/(3*sqrt(-x**2 - 4*x - 3)))/9 + 10*atanh(x/sqrt(-x**2 - 4*x - 3))/27 + sqrt(-x**2 - 4*x - 3)/(9*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_132():
    f = (3*x + 2)**2*(-12*x**2 + 31*x + 30)**2*sqrt(12*x**2 + 17*x + 6)
    F = -(sympy.S(5)/16 - 3*x/32)*(12*x**2 + 17*x + 6)**(sympy.S(7)/2) + (602184*x + 426547)*(12*x**2 + 17*x + 6)**(sympy.S(5)/2)/24576 - (3010920*x + 2132735)*(12*x**2 + 17*x + 6)**(sympy.S(3)/2)/4718592 + (3010920*x + 2132735)*sqrt(12*x**2 + 17*x + 6)/150994944 - 873*(12*x**2 + 17*x + 6)**(sympy.S(7)/2)/1792 - 125455*sqrt(3)*atanh(sqrt(3)*(24*x + 17)/(12*sqrt(12*x**2 + 17*x + 6)))/1811939328
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_133():
    f = (3*x + 2)*(-12*x**2 + 31*x + 30)*sqrt(12*x**2 + 17*x + 6)
    F = (97*x/32 + sympy.S(1649)/768)*(12*x**2 + 17*x + 6)**(sympy.S(3)/2) - (2328*x + 1649)*sqrt(12*x**2 + 17*x + 6)/24576 - (12*x**2 + 17*x + 6)**(sympy.S(5)/2)/20 + 97*sqrt(3)*atanh(sqrt(3)*(24*x + 17)/(12*sqrt(12*x**2 + 17*x + 6)))/294912
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_134():
    f = sqrt(12*x**2 + 17*x + 6)/((3*x + 2)*(-12*x**2 + 31*x + 30))
    F = atanh((291*x + 206)/(84*sqrt(12*x**2 + 17*x + 6)))/42
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_135():
    f = sqrt(12*x**2 + 17*x + 6)/((3*x + 2)**2*(-12*x**2 + 31*x + 30)**2)
    F = 97*atanh((291*x + 206)/(84*sqrt(12*x**2 + 17*x + 6)))/3226944 + 3137*sqrt(12*x**2 + 17*x + 6)/(384160 - 115248*x) - (388*x + 275)/((980 - 294*x)*sqrt(12*x**2 + 17*x + 6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_136():
    f = sqrt(12*x**2 + 17*x + 6)/((3*x + 2)**3*(-12*x**2 + 31*x + 30)**3)
    F = 40325*atanh((291*x + 206)/(84*sqrt(12*x**2 + 17*x + 6)))/637540872192 - 1634466587*sqrt(12*x**2 + 17*x + 6)/(75897722880 - 22769316864*x) - (388*x + 275)/(294*(10 - 3*x)**2*(12*x**2 + 17*x + 6)**(sympy.S(3)/2)) + (1042556*x + 738029)/(8232*(10 - 3*x)**2*sqrt(12*x**2 + 17*x + 6)) - 50555899*sqrt(12*x**2 + 17*x + 6)/(19361664*(10 - 3*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_137():
    f = (2*x - 3)*(x**2 - 3*x)**(sympy.S(2)/3)
    F = 3*(x**2 - 3*x)**(sympy.S(5)/3)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_138():
    f = (x*(x - 3))**(sympy.S(2)/3)*(2*x - 3)
    F = 3*(-x*(3 - x))**(sympy.S(5)/3)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_139():
    f = x*(2*x**2 - 9*x + 9)/(x**2 - 3*x)**(sympy.S(1)/3)
    F = 3*(x**2 - 3*x)**(sympy.S(5)/3)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_140():
    f = x*(2*x**2 - 9*x + 9)/(x*(x - 3))**(sympy.S(1)/3)
    F = 3*(x**2 - 3*x)**(sympy.S(5)/3)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_141():
    f = (g + h*x)/((g**2 + 3*h**2*x**2)*(-c*g**2/h**2 + 9*c*x**2)**(sympy.S(1)/3))
    F = 2**(sympy.S(1)/3)*(1 - 9*h**2*x**2/g**2)**(sympy.S(1)/3)*log(g**2 + 3*h**2*x**2)/(12*h*(-c*g**2/h**2 + 9*c*x**2)**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*(1 - 9*h**2*x**2/g**2)**(sympy.S(1)/3)*log((1 - 3*h*x/g)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(1 + 3*h*x/g)**(sympy.S(1)/3))/(4*h*(-c*g**2/h**2 + 9*c*x**2)**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*sqrt(3)*(1 - 9*h**2*x**2/g**2)**(sympy.S(1)/3)*atan(2**(sympy.S(2)/3)*sqrt(3)*(1 - 3*h*x/g)**(sympy.S(2)/3)/(3*(1 + 3*h*x/g)**(sympy.S(1)/3)) - sqrt(3)/3)/(6*h*(-c*g**2/h**2 + 9*c*x**2)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_6_g_plus_h_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_d_plus_e_x_plus_f_x_pow_2_pow_q_142():
    f = (g + h*x)/((b*x + c*x**2 + (2*b**2*h**2 + b*c*g*h - c**2*g**2)/(9*c*h**2))**(sympy.S(1)/3)*(b*f*x/c + f*x**2 + f*(b**2 - (2*b**2*h**2 + b*c*g*h - c**2*g**2)/(3*h**2))/c**2))
    F = -3*3**(sympy.S(2)/3)*h*(c*h**2*(-9*b*x - 9*c*x**2 + (-2*b*h + c*g)*(b*h + c*g)/(c*h**2))/(-b*h + 2*c*g)**2)**(sympy.S(1)/3)*log((-3*h*(b + 2*c*x)/(-b*h + 2*c*g) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(3*h*(b + 2*c*x)/(-b*h + 2*c*g) + 1)**(sympy.S(1)/3))/(2*f*(9*b*x + 9*c*x**2 - (-2*b*h + c*g)*(b*h + c*g)/(c*h**2))**(sympy.S(1)/3)) + 3**(sympy.S(2)/3)*h*(c*h**2*(-9*b*x - 9*c*x**2 + (-2*b*h + c*g)*(b*h + c*g)/(c*h**2))/(-b*h + 2*c*g)**2)**(sympy.S(1)/3)*log(b*f*x/c + f*x**2 + f*(b**2*h**2 - b*c*g*h + c**2*g**2)/(3*c**2*h**2))/(2*f*(9*b*x + 9*c*x**2 - (-2*b*h + c*g)*(b*h + c*g)/(c*h**2))**(sympy.S(1)/3)) - 3*3**(sympy.S(1)/6)*h*(c*h**2*(-9*b*x - 9*c*x**2 + (-2*b*h + c*g)*(b*h + c*g)/(c*h**2))/(-b*h + 2*c*g)**2)**(sympy.S(1)/3)*atan(2**(sympy.S(2)/3)*sqrt(3)*(-3*h*(b + 2*c*x)/(-b*h + 2*c*g) + 1)**(sympy.S(2)/3)/(3*(3*h*(b + 2*c*x)/(-b*h + 2*c*g) + 1)**(sympy.S(1)/3)) - sqrt(3)/3)/(f*(9*b*x + 9*c*x**2 - (-2*b*h + c*g)*(b*h + c*g)/(c*h**2))**(sympy.S(1)/3))
    assert integrate(f, x) == F

