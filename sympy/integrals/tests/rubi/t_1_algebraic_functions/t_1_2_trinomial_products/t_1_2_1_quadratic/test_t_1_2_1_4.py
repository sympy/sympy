"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.2 Trinomial products/1.2.1 Quadratic/1.2.1.4 (d+e x)^m (f+g x)^n (a+b x+c x^2)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, f, g, m, n, p = symbols('a b c d e f g m n p')

def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_1():
    f = x**2*(d + e*x)*sqrt(d**2 - e**2*x**2)
    F = d**5*atan(e*x/sqrt(d**2 - e**2*x**2))/(8*e**3) + d**3*x*sqrt(d**2 - e**2*x**2)/(8*e**2) - d**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(3*e**3) - d*x*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(4*e**2) + (d**2 - e**2*x**2)**(sympy.S(5)/2)/(5*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_2():
    f = x**4*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)
    F = 3*d**9*atan(e*x/sqrt(d**2 - e**2*x**2))/(128*e**5) + 3*d**7*x*sqrt(d**2 - e**2*x**2)/(128*e**4) + d**5*x*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(64*e**4) - d**3*(128*d + 315*e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(5040*e**5) - 4*d**2*x**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(63*e**3) - d*x**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(8*e**2) - x**4*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(9*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_3():
    f = x**3*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)
    F = 3*d**8*atan(e*x/sqrt(d**2 - e**2*x**2))/(128*e**4) + 3*d**6*x*sqrt(d**2 - e**2*x**2)/(128*e**3) + d**4*x*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(64*e**3) - d**2*(32*d + 35*e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(560*e**4) - d*x**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(7*e**2) - x**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(8*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_4():
    f = x**2*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)
    F = d**7*atan(e*x/sqrt(d**2 - e**2*x**2))/(16*e**3) + d**5*x*sqrt(d**2 - e**2*x**2)/(16*e**2) + d**3*x*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(24*e**2) - d**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(5*e**3) - d*x*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(6*e**2) + (d**2 - e**2*x**2)**(sympy.S(7)/2)/(7*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_5():
    f = x*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)
    F = d**6*atan(e*x/sqrt(d**2 - e**2*x**2))/(16*e**2) + d**4*x*sqrt(d**2 - e**2*x**2)/(16*e) + d**2*x*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(24*e) - (6*d + 5*e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(30*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_6():
    f = x*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)
    F = d**6*atan(e*x/sqrt(d**2 - e**2*x**2))/(16*e**2) + d**4*x*sqrt(d**2 - e**2*x**2)/(16*e) + d**2*x*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(24*e) - (6*d + 5*e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(30*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_7():
    f = (d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/x
    F = 3*d**4*atan(e*x/sqrt(d**2 - e**2*x**2))/8 - d**4*atanh(sqrt(d**2 - e**2*x**2)/d) + d**2*(8*d + 3*e*x)*sqrt(d**2 - e**2*x**2)/8 + (d/3 + e*x/4)*(d**2 - e**2*x**2)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_8():
    f = (d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/x**2
    F = -3*d**3*e*atan(e*x/sqrt(d**2 - e**2*x**2))/2 - d**3*e*atanh(sqrt(d**2 - e**2*x**2)/d) + d*e*(2*d - 3*e*x)*sqrt(d**2 - e**2*x**2)/2 - (3*d - e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_9():
    f = (d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/x**3
    F = -3*d**2*e**2*atan(e*x/sqrt(d**2 - e**2*x**2))/2 + 3*d**2*e**2*atanh(sqrt(d**2 - e**2*x**2)/d)/2 - 3*d*e*(d + e*x)*sqrt(d**2 - e**2*x**2)/(2*x) - (d - e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_10():
    f = (d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/x**4
    F = d*e**3*atan(e*x/sqrt(d**2 - e**2*x**2)) + 3*d*e**3*atanh(sqrt(d**2 - e**2*x**2)/d)/2 + e**2*(2*d - 3*e*x)*sqrt(d**2 - e**2*x**2)/(2*x) - (2*d + 3*e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(6*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_11():
    f = (d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/x**5
    F = e**4*atan(e*x/sqrt(d**2 - e**2*x**2)) - 3*e**4*atanh(sqrt(d**2 - e**2*x**2)/d)/8 + e**2*(3*d + 8*e*x)*sqrt(d**2 - e**2*x**2)/(8*x**2) - (3*d + 4*e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(12*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_12():
    f = (d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/x**6
    F = 3*e**3*sqrt(d**2 - e**2*x**2)/(8*x**2) - e*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(4*x**4) - 3*e**5*atanh(sqrt(d**2 - e**2*x**2)/d)/(8*d) - (d**2 - e**2*x**2)**(sympy.S(5)/2)/(5*d*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_13():
    f = (d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/x**7
    F = e**4*sqrt(d**2 - e**2*x**2)/(16*d*x**2) - e**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(24*d*x**4) - (d**2 - e**2*x**2)**(sympy.S(5)/2)/(6*d*x**6) - e**6*atanh(sqrt(d**2 - e**2*x**2)/d)/(16*d**2) - e*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(5*d**2*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_14():
    f = (d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/x**8
    F = -(d**2 - e**2*x**2)**(sympy.S(5)/2)/(7*d*x**7) + e**5*sqrt(d**2 - e**2*x**2)/(16*d**2*x**2) - e**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(24*d**2*x**4) - e*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(6*d**2*x**6) - e**7*atanh(sqrt(d**2 - e**2*x**2)/d)/(16*d**3) - 2*e**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(35*d**3*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_15():
    f = (d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/x**9
    F = -(d**2 - e**2*x**2)**(sympy.S(5)/2)/(8*d*x**8) - e*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(7*d**2*x**7) + 3*e**6*sqrt(d**2 - e**2*x**2)/(128*d**3*x**2) - e**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(64*d**3*x**4) - e**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(16*d**3*x**6) - 3*e**8*atanh(sqrt(d**2 - e**2*x**2)/d)/(128*d**4) - 2*e**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(35*d**4*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_16():
    f = x**2*(d + e*x)/sqrt(d**2 - e**2*x**2)
    F = d**3*atan(e*x/sqrt(d**2 - e**2*x**2))/(2*e**3) - d**2*sqrt(d**2 - e**2*x**2)/e**3 - d*x*sqrt(d**2 - e**2*x**2)/(2*e**2) + (d**2 - e**2*x**2)**(sympy.S(3)/2)/(3*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_17():
    f = x**2*(d + e*x)/(d**2 - e**2*x**2)**(sympy.S(3)/2)
    F = d*(d + e*x)/(e**3*sqrt(d**2 - e**2*x**2)) - d*atan(e*x/sqrt(d**2 - e**2*x**2))/e**3 + sqrt(d**2 - e**2*x**2)/e**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_18():
    f = x**2*(d + e*x)/(d**2 - e**2*x**2)**(sympy.S(5)/2)
    F = -2/(3*e**3*sqrt(d**2 - e**2*x**2)) + x**2*(d + e*x)/(3*d*e*(d**2 - e**2*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_19():
    f = x**7*(d + e*x)/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = -7*d**2*atan(e*x/sqrt(d**2 - e**2*x**2))/(2*e**8) + x**6*(d + e*x)/(5*e**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - x**4*(6*d + 7*e*x)/(15*e**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + x**2*(24*d + 35*e*x)/(15*e**6*sqrt(d**2 - e**2*x**2)) + (32*d + 35*e*x)*sqrt(d**2 - e**2*x**2)/(10*e**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_20():
    f = x**6*(d + e*x)/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = -d*atan(e*x/sqrt(d**2 - e**2*x**2))/e**7 + x**5*(d + e*x)/(5*e**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - x**3*(5*d + 6*e*x)/(15*e**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + x*(5*d + 8*e*x)/(5*e**6*sqrt(d**2 - e**2*x**2)) + 16*sqrt(d**2 - e**2*x**2)/(5*e**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_21():
    f = x**5*(d + e*x)/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = x**4*(d + e*x)/(5*e**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - x**2*(4*d + 5*e*x)/(15*e**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (8*d + 15*e*x)/(15*e**6*sqrt(d**2 - e**2*x**2)) - atan(e*x/sqrt(d**2 - e**2*x**2))/e**6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_22():
    f = x**4*(d + e*x)/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = -4*d**2/(15*e**5*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + 4/(5*e**5*sqrt(d**2 - e**2*x**2)) + x**4*(d + e*x)/(5*d*e*(d**2 - e**2*x**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_23():
    f = x**3*(d + e*x)/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = x**2*(d + e*x)/(5*e**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - (2*d + 3*e*x)/(15*e**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + x/(5*d**2*e**3*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_24():
    f = x**2*(d + e*x)/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = x**2*(d + e*x)/(5*d*e*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - (2*d - 2*e*x)/(15*d*e**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - 2*x/(15*d**3*e**2*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_25():
    f = x*(d + e*x)/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = (d + e*x)/(5*e**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - x/(15*d**2*e*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - 2*x/(15*d**4*e*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_26():
    f = (d + e*x)/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = (d + e*x)/(5*d*e*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + 4*x/(15*d**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + 8*x/(15*d**5*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_27():
    f = (d + e*x)/(x*(d**2 - e**2*x**2)**(sympy.S(7)/2))
    F = (d + e*x)/(5*d**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + (5*d + 4*e*x)/(15*d**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (15*d + 8*e*x)/(15*d**6*sqrt(d**2 - e**2*x**2)) - atanh(sqrt(d**2 - e**2*x**2)/d)/d**6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_28():
    f = (d + e*x)/(x**2*(d**2 - e**2*x**2)**(sympy.S(7)/2))
    F = (d + e*x)/(5*d**2*x*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + (6*d + 5*e*x)/(15*d**4*x*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (8*d + 5*e*x)/(5*d**6*x*sqrt(d**2 - e**2*x**2)) - e*atanh(sqrt(d**2 - e**2*x**2)/d)/d**7 - 16*sqrt(d**2 - e**2*x**2)/(5*d**7*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_29():
    f = (d + e*x)/(x**3*(d**2 - e**2*x**2)**(sympy.S(7)/2))
    F = (d + e*x)/(5*d**2*x**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + (7*d + 6*e*x)/(15*d**4*x**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (35*d + 24*e*x)/(15*d**6*x**2*sqrt(d**2 - e**2*x**2)) - 7*sqrt(d**2 - e**2*x**2)/(2*d**7*x**2) - 7*e**2*atanh(sqrt(d**2 - e**2*x**2)/d)/(2*d**8) - 16*e*sqrt(d**2 - e**2*x**2)/(5*d**8*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_30():
    f = x**2*(d + e*x)/(d**2 - e**2*x**2)**(sympy.S(9)/2)
    F = x**2*(d + e*x)/(7*d*e*(d**2 - e**2*x**2)**(sympy.S(7)/2)) - (2*d - 4*e*x)/(35*d*e**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 4*x/(105*d**3*e**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - 8*x/(105*d**5*e**2*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_31():
    f = x**2*(d + e*x)/(d**2 - e**2*x**2)**(sympy.S(11)/2)
    F = x**2*(d + e*x)/(9*d*e*(d**2 - e**2*x**2)**(sympy.S(9)/2)) - (2*d - 6*e*x)/(63*d*e**3*(d**2 - e**2*x**2)**(sympy.S(7)/2)) - 2*x/(105*d**3*e**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 8*x/(315*d**5*e**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - 16*x/(315*d**7*e**2*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_32():
    f = x**2*(-a*x + 1)/(-a**2*x**2 + 1)**(sympy.S(3)/2)
    F = -(-a*x + 1)/(a**3*sqrt(-a**2*x**2 + 1)) - sqrt(-a**2*x**2 + 1)/a**3 - asin(a*x)/a**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_33():
    f = x**4*(d + e*x)**2/sqrt(d**2 - e**2*x**2)
    F = 11*d**6*atan(e*x/sqrt(d**2 - e**2*x**2))/(16*e**5) - d**4*(256*d + 165*e*x)*sqrt(d**2 - e**2*x**2)/(240*e**5) - 8*d**3*x**2*sqrt(d**2 - e**2*x**2)/(15*e**3) - 11*d**2*x**3*sqrt(d**2 - e**2*x**2)/(24*e**2) - 2*d*x**4*sqrt(d**2 - e**2*x**2)/(5*e) - x**5*sqrt(d**2 - e**2*x**2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_34():
    f = x**3*(d + e*x)**2/sqrt(d**2 - e**2*x**2)
    F = 3*d**5*atan(e*x/sqrt(d**2 - e**2*x**2))/(4*e**4) - 3*d**3*(8*d + 5*e*x)*sqrt(d**2 - e**2*x**2)/(20*e**4) - 3*d**2*x**2*sqrt(d**2 - e**2*x**2)/(5*e**2) - d*x**3*sqrt(d**2 - e**2*x**2)/(2*e) - x**4*sqrt(d**2 - e**2*x**2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_35():
    f = x**2*(d + e*x)**2/sqrt(d**2 - e**2*x**2)
    F = 7*d**4*atan(e*x/sqrt(d**2 - e**2*x**2))/(8*e**3) - d**2*(32*d + 21*e*x)*sqrt(d**2 - e**2*x**2)/(24*e**3) - 2*d*x**2*sqrt(d**2 - e**2*x**2)/(3*e) - x**3*sqrt(d**2 - e**2*x**2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_36():
    f = x*(d + e*x)**2/sqrt(d**2 - e**2*x**2)
    F = d**3*atan(e*x/sqrt(d**2 - e**2*x**2))/e**2 - d*(5*d + 3*e*x)*sqrt(d**2 - e**2*x**2)/(3*e**2) - x**2*sqrt(d**2 - e**2*x**2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_37():
    f = (d + e*x)**2/sqrt(d**2 - e**2*x**2)
    F = 3*d**2*atan(e*x/sqrt(d**2 - e**2*x**2))/(2*e) - 3*d*sqrt(d**2 - e**2*x**2)/(2*e) - (d + e*x)*sqrt(d**2 - e**2*x**2)/(2*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_38():
    f = (d + e*x)**2/(x*sqrt(d**2 - e**2*x**2))
    F = 2*d*atan(e*x/sqrt(d**2 - e**2*x**2)) - d*atanh(sqrt(d**2 - e**2*x**2)/d) - sqrt(d**2 - e**2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_39():
    f = (d + e*x)**2/(x**2*sqrt(d**2 - e**2*x**2))
    F = e*atan(e*x/sqrt(d**2 - e**2*x**2)) - 2*e*atanh(sqrt(d**2 - e**2*x**2)/d) - sqrt(d**2 - e**2*x**2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_40():
    f = (d + e*x)**2/(x**3*sqrt(d**2 - e**2*x**2))
    F = -sqrt(d**2 - e**2*x**2)/(2*x**2) - 3*e**2*atanh(sqrt(d**2 - e**2*x**2)/d)/(2*d) - 2*e*sqrt(d**2 - e**2*x**2)/(d*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_41():
    f = (d + e*x)**2/(x**4*sqrt(d**2 - e**2*x**2))
    F = -sqrt(d**2 - e**2*x**2)/(3*x**3) - e*sqrt(d**2 - e**2*x**2)/(d*x**2) - e**3*atanh(sqrt(d**2 - e**2*x**2)/d)/d**2 - 5*e**2*sqrt(d**2 - e**2*x**2)/(3*d**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_42():
    f = (d + e*x)**2/(x**5*sqrt(d**2 - e**2*x**2))
    F = -sqrt(d**2 - e**2*x**2)/(4*x**4) - 2*e*sqrt(d**2 - e**2*x**2)/(3*d*x**3) - 7*e**2*sqrt(d**2 - e**2*x**2)/(8*d**2*x**2) - 7*e**4*atanh(sqrt(d**2 - e**2*x**2)/d)/(8*d**3) - 4*e**3*sqrt(d**2 - e**2*x**2)/(3*d**3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_43():
    f = (d + e*x)**2/(x**6*sqrt(d**2 - e**2*x**2))
    F = -sqrt(d**2 - e**2*x**2)/(5*x**5) - e*sqrt(d**2 - e**2*x**2)/(2*d*x**4) - 3*e**2*sqrt(d**2 - e**2*x**2)/(5*d**2*x**3) - 3*e**3*sqrt(d**2 - e**2*x**2)/(4*d**3*x**2) - 3*e**5*atanh(sqrt(d**2 - e**2*x**2)/d)/(4*d**4) - 6*e**4*sqrt(d**2 - e**2*x**2)/(5*d**4*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_44():
    f = x**5*(d + e*x)**2/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = d**4*(d + e*x)**2/(5*e**6*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 22*d**3*(d + e*x)/(15*e**6*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + 2*d*(30*d + 23*e*x)/(15*e**6*sqrt(d**2 - e**2*x**2)) - 2*d*atan(e*x/sqrt(d**2 - e**2*x**2))/e**6 + sqrt(d**2 - e**2*x**2)/e**6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_45():
    f = x**4*(d + e*x)**2/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = d**3*(d + e*x)**2/(5*e**5*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 17*d**2*(d + e*x)/(15*e**5*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (30*d + 26*e*x)/(15*e**5*sqrt(d**2 - e**2*x**2)) - atan(e*x/sqrt(d**2 - e**2*x**2))/e**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_46():
    f = x**3*(d + e*x)**2/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = d**2*(d + e*x)**2/(5*e**4*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 4*d*(d + e*x)/(5*e**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (5*d + 2*e*x)/(5*d*e**4*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_47():
    f = x**2*(d + e*x)**2/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = d*(d + e*x)**2/(5*e**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - (7*d + 7*e*x)/(15*e**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + x/(15*d**2*e**2*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_48():
    f = x*(d + e*x)**2/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = (d + e*x)**2/(5*e**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - (2*d + 2*e*x)/(15*d*e**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - 4*x/(15*d**3*e*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_49():
    f = (d + e*x)**2/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = (2*d + 2*e*x)/(5*e*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + x/(5*d**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + 2*x/(5*d**4*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_50():
    f = (d + e*x)**2/(x*(d**2 - e**2*x**2)**(sympy.S(7)/2))
    F = (2*d + 2*e*x)/(5*d*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + (5*d + 8*e*x)/(15*d**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (15*d + 16*e*x)/(15*d**5*sqrt(d**2 - e**2*x**2)) - atanh(sqrt(d**2 - e**2*x**2)/d)/d**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_51():
    f = (d + e*x)**2/(x**2*(d**2 - e**2*x**2)**(sympy.S(7)/2))
    F = 2*e*(d + e*x)/(5*d**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + e*(10*d + 13*e*x)/(15*d**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + e*(30*d + 41*e*x)/(15*d**6*sqrt(d**2 - e**2*x**2)) - 2*e*atanh(sqrt(d**2 - e**2*x**2)/d)/d**6 - sqrt(d**2 - e**2*x**2)/(d**6*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_52():
    f = (d + e*x)**2/(x**3*(d**2 - e**2*x**2)**(sympy.S(7)/2))
    F = 2*e**2*(d + e*x)/(5*d**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + e**2*(5*d + 6*e*x)/(5*d**5*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - sqrt(d**2 - e**2*x**2)/(2*d**6*x**2) + 2*e**2*(10*d + 11*e*x)/(5*d**7*sqrt(d**2 - e**2*x**2)) - 9*e**2*atanh(sqrt(d**2 - e**2*x**2)/d)/(2*d**7) - 2*e*sqrt(d**2 - e**2*x**2)/(d**7*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_53():
    f = (d + e*x)**2/(x**4*(d**2 - e**2*x**2)**(sympy.S(7)/2))
    F = 2*e**3*(d + e*x)/(5*d**4*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + e**3*(20*d + 23*e*x)/(15*d**6*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - sqrt(d**2 - e**2*x**2)/(3*d**6*x**3) - e*sqrt(d**2 - e**2*x**2)/(d**7*x**2) + 2*e**3*(45*d + 53*e*x)/(15*d**8*sqrt(d**2 - e**2*x**2)) - 7*e**3*atanh(sqrt(d**2 - e**2*x**2)/d)/d**8 - 14*e**2*sqrt(d**2 - e**2*x**2)/(3*d**8*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_54():
    f = x**3*(x + 1)**2/sqrt(1 - x**2)
    F = -x**4*sqrt(1 - x**2)/5 - x**3*sqrt(1 - x**2)/2 - 3*x**2*sqrt(1 - x**2)/5 - sqrt(1 - x**2)*(3*x/4 + sympy.S(6)/5) + 3*asin(x)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_55():
    f = x**2*(x + 1)**2/sqrt(1 - x**2)
    F = -x**3*sqrt(1 - x**2)/4 - 2*x**2*sqrt(1 - x**2)/3 - sqrt(1 - x**2)*(7*x/8 + sympy.S(4)/3) + 7*asin(x)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_56():
    f = x*(x + 1)**2/sqrt(1 - x**2)
    F = -x**2*sqrt(1 - x**2)/3 - sqrt(1 - x**2)*(x + sympy.S(5)/3) + asin(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_57():
    f = (x + 1)**2/sqrt(1 - x**2)
    F = -sqrt(1 - x**2)*(x/2 + sympy.S.Half) - 3*sqrt(1 - x**2)/2 + 3*asin(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_58():
    f = (x + 1)**2/(x*sqrt(1 - x**2))
    F = -sqrt(1 - x**2) + 2*asin(x) - atanh(sqrt(1 - x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_59():
    f = (x + 1)**2/(x**2*sqrt(1 - x**2))
    F = asin(x) - 2*atanh(sqrt(1 - x**2)) - sqrt(1 - x**2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_60():
    f = (x + 1)**2/(x**3*sqrt(1 - x**2))
    F = -3*atanh(sqrt(1 - x**2))/2 - 2*sqrt(1 - x**2)/x - sqrt(1 - x**2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_61():
    f = (x + 1)**2/(x**4*sqrt(1 - x**2))
    F = -atanh(sqrt(1 - x**2)) - 5*sqrt(1 - x**2)/(3*x) - sqrt(1 - x**2)/x**2 - sqrt(1 - x**2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_62():
    f = (x + 1)**2/(x**5*sqrt(1 - x**2))
    F = -7*atanh(sqrt(1 - x**2))/8 - 4*sqrt(1 - x**2)/(3*x) - 7*sqrt(1 - x**2)/(8*x**2) - 2*sqrt(1 - x**2)/(3*x**3) - sqrt(1 - x**2)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_63():
    f = (x + 1)**2/(x**6*sqrt(1 - x**2))
    F = -3*atanh(sqrt(1 - x**2))/4 - 6*sqrt(1 - x**2)/(5*x) - 3*sqrt(1 - x**2)/(4*x**2) - 3*sqrt(1 - x**2)/(5*x**3) - sqrt(1 - x**2)/(2*x**4) - sqrt(1 - x**2)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_64():
    f = (d + e*x)**3*sqrt(d**2 - e**2*x**2)/x**5
    F = -d*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(4*x**4) - e**4*atan(e*x/sqrt(d**2 - e**2*x**2)) + 13*e**4*atanh(sqrt(d**2 - e**2*x**2)/d)/8 - e**2*(13*d + 8*e*x)*sqrt(d**2 - e**2*x**2)/(8*x**2) - e*(d**2 - e**2*x**2)**(sympy.S(3)/2)/x**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_65():
    f = x**5*(d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)
    F = 35*d**14*atan(e*x/sqrt(d**2 - e**2*x**2))/(2048*e**6) + 35*d**12*x*sqrt(d**2 - e**2*x**2)/(2048*e**5) + 35*d**10*x*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(3072*e**5) + 7*d**8*x*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(768*e**5) - d**6*(31744*d + 63063*e*x)*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(1153152*e**6) - 124*d**5*x**2*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(1287*e**4) - 7*d**4*x**3*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(48*e**3) - 31*d**3*x**4*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(143*e**2) - 7*d**2*x**5*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(24*e) - 3*d*x**6*(d**2 - e**2*x**2)**(sympy.S(7)/2)/13 - e*x**7*(d**2 - e**2*x**2)**(sympy.S(7)/2)/14
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_66():
    f = x**4*(d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)
    F = 27*d**13*atan(e*x/sqrt(d**2 - e**2*x**2))/(1024*e**5) + 27*d**11*x*sqrt(d**2 - e**2*x**2)/(1024*e**4) + 9*d**9*x*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(512*e**4) + 9*d**7*x*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(640*e**4) - d**5*(12800*d + 27027*e*x)*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(320320*e**5) - 20*d**4*x**2*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(143*e**3) - 9*d**3*x**3*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(40*e**2) - 45*d**2*x**4*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(143*e) - d*x**5*(d**2 - e**2*x**2)**(sympy.S(7)/2)/4 - e*x**6*(d**2 - e**2*x**2)**(sympy.S(7)/2)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_67():
    f = x**3*(d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)
    F = 41*d**12*atan(e*x/sqrt(d**2 - e**2*x**2))/(1024*e**4) + 41*d**10*x*sqrt(d**2 - e**2*x**2)/(1024*e**3) + 41*d**8*x*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(1536*e**3) + 41*d**6*x*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(1920*e**3) - d**4*(14720*d + 28413*e*x)*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(221760*e**4) - 23*d**3*x**2*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(99*e**2) - 41*d**2*x**3*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(120*e) - 3*d*x**4*(d**2 - e**2*x**2)**(sympy.S(7)/2)/11 - e*x**5*(d**2 - e**2*x**2)**(sympy.S(7)/2)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_68():
    f = x**2*(d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)
    F = 19*d**11*atan(e*x/sqrt(d**2 - e**2*x**2))/(256*e**3) + 19*d**9*x*sqrt(d**2 - e**2*x**2)/(256*e**2) + 19*d**7*x*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(384*e**2) + 19*d**5*x*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(480*e**2) - d**3*(5920*d + 13167*e*x)*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(55440*e**3) - 37*d**2*x**2*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(99*e) - 3*d*x**3*(d**2 - e**2*x**2)**(sympy.S(7)/2)/10 - e*x**4*(d**2 - e**2*x**2)**(sympy.S(7)/2)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_69():
    f = x*(d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)
    F = 33*d**10*atan(e*x/sqrt(d**2 - e**2*x**2))/(256*e**2) + 33*d**8*x*sqrt(d**2 - e**2*x**2)/(256*e) + 11*d**6*x*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(128*e) + 11*d**4*x*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(160*e) - 33*d**3*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(560*e**2) - 11*d**2*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(240*e**2) - d*(d + e*x)**2*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(30*e**2) - (d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(10*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_70():
    f = (d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)
    F = 55*d**9*atan(e*x/sqrt(d**2 - e**2*x**2))/(128*e) + 55*d**7*x*sqrt(d**2 - e**2*x**2)/128 + 55*d**5*x*(d**2 - e**2*x**2)**(sympy.S(3)/2)/192 + 11*d**3*x*(d**2 - e**2*x**2)**(sympy.S(5)/2)/48 - 11*d**2*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(56*e) - 11*d*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(72*e) - (d + e*x)**2*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(9*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_71():
    f = (d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/x
    F = 125*d**8*atan(e*x/sqrt(d**2 - e**2*x**2))/128 - d**8*atanh(sqrt(d**2 - e**2*x**2)/d) + d**6*(128*d + 125*e*x)*sqrt(d**2 - e**2*x**2)/128 + d**4*(64*d + 125*e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/192 + d**2*(48*d + 125*e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)/240 - 3*d*(d**2 - e**2*x**2)**(sympy.S(7)/2)/7 - e*x*(d**2 - e**2*x**2)**(sympy.S(7)/2)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_72():
    f = (d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/x**2
    F = -15*d**7*e*atan(e*x/sqrt(d**2 - e**2*x**2))/16 - 3*d**7*e*atanh(sqrt(d**2 - e**2*x**2)/d) + 3*d**5*e*(16*d - 5*e*x)*sqrt(d**2 - e**2*x**2)/16 + d**3*e*(8*d - 5*e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/8 + d*e*(6*d - 5*e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)/10 - d*(d**2 - e**2*x**2)**(sympy.S(7)/2)/x - e*(d**2 - e**2*x**2)**(sympy.S(7)/2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_73():
    f = (d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/x**3
    F = -85*d**6*e**2*atan(e*x/sqrt(d**2 - e**2*x**2))/16 - d**6*e**2*atanh(sqrt(d**2 - e**2*x**2)/d)/2 + d**4*e**2*(8*d - 85*e*x)*sqrt(d**2 - e**2*x**2)/16 + d**2*e**2*(4*d - 85*e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/24 - d*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(2*x**2) + e**2*(3*d - 85*e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)/30 - 3*e*(d**2 - e**2*x**2)**(sympy.S(7)/2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_74():
    f = (d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/x**4
    F = -25*d**5*e**3*atan(e*x/sqrt(d**2 - e**2*x**2))/8 + 13*d**5*e**3*atanh(sqrt(d**2 - e**2*x**2)/d)/2 - d**3*e**3*(52*d + 25*e*x)*sqrt(d**2 - e**2*x**2)/8 - d*e**3*(26*d + 25*e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/12 - d*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(3*x**3) - e**2*(50*d + 39*e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(30*x) - 3*e*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_75():
    f = (d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/x**5
    F = 45*d**4*e**4*atan(e*x/sqrt(d**2 - e**2*x**2))/8 + 45*d**4*e**4*atanh(sqrt(d**2 - e**2*x**2)/d)/8 - 45*d**2*e**4*(d - e*x)*sqrt(d**2 - e**2*x**2)/8 + 15*d*e**3*(2*d - e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(8*x) - d*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(4*x**4) - 3*e**2*(3*d + 2*e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(8*x**2) - e*(d**2 - e**2*x**2)**(sympy.S(7)/2)/x**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_76():
    f = (d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/x**6
    F = 13*d**3*e**5*atan(e*x/sqrt(d**2 - e**2*x**2))/2 - 25*d**3*e**5*atanh(sqrt(d**2 - e**2*x**2)/d)/8 + d**2*e**4*(52*d + 25*e*x)*sqrt(d**2 - e**2*x**2)/(8*x) + d*e**3*(25*d - 52*e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(24*x**2) - d*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(5*x**5) - e**2*(52*d + 25*e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(60*x**3) - 3*e*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_77():
    f = (d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/x**7
    F = -d**2*e**6*atan(e*x/sqrt(d**2 - e**2*x**2))/2 - 85*d**2*e**6*atanh(sqrt(d**2 - e**2*x**2)/d)/16 - d*e**5*(8*d - 85*e*x)*sqrt(d**2 - e**2*x**2)/(16*x) + d*e**3*(8*d + 85*e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(48*x**3) - d*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(6*x**6) - e**2*(85*d + 12*e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(120*x**4) - 3*e*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_78():
    f = (d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/x**8
    F = -3*d*e**7*atan(e*x/sqrt(d**2 - e**2*x**2)) - 15*d*e**7*atanh(sqrt(d**2 - e**2*x**2)/d)/16 - d*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(7*x**7) - 3*e**6*(16*d - 5*e*x)*sqrt(d**2 - e**2*x**2)/(16*x) + e**4*(16*d + 5*e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(16*x**3) - e**2*(24*d + 5*e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(40*x**5) - e*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(2*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_79():
    f = (d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/x**9
    F = -d*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(8*x**8) - e**8*atan(e*x/sqrt(d**2 - e**2*x**2)) + 125*e**8*atanh(sqrt(d**2 - e**2*x**2)/d)/128 - e**6*(125*d + 128*e*x)*sqrt(d**2 - e**2*x**2)/(128*x**2) + e**4*(125*d + 64*e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(192*x**4) - e**2*(125*d + 48*e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(240*x**6) - 3*e*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(7*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_80():
    f = (d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/x**10
    F = -d*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(9*x**9) - 55*e**7*sqrt(d**2 - e**2*x**2)/(128*x**2) + 55*e**5*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(192*x**4) - 11*e**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(48*x**6) - 3*e*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(8*x**8) + 55*e**9*atanh(sqrt(d**2 - e**2*x**2)/d)/(128*d) - 29*e**2*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(63*d*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_81():
    f = (d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/x**11
    F = -d*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(10*x**10) - e*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(3*x**9) - 33*e**8*sqrt(d**2 - e**2*x**2)/(256*d*x**2) + 11*e**6*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(128*d*x**4) - 11*e**4*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(160*d*x**6) - 33*e**2*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(80*d*x**8) + 33*e**10*atanh(sqrt(d**2 - e**2*x**2)/d)/(256*d**2) - 5*e**3*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(21*d**2*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_82():
    f = (d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/x**12
    F = -d*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(11*x**11) - 3*e*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(10*x**10) - 37*e**2*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(99*d*x**9) - 19*e**9*sqrt(d**2 - e**2*x**2)/(256*d**2*x**2) + 19*e**7*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(384*d**2*x**4) - 19*e**5*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(480*d**2*x**6) - 19*e**3*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(80*d**2*x**8) + 19*e**11*atanh(sqrt(d**2 - e**2*x**2)/d)/(256*d**3) - 74*e**4*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(693*d**3*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_83():
    f = x**5*(d + e*x)**3/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = d**4*(d + e*x)**3/(5*e**6*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 23*d**3*(d + e*x)**2/(15*e**6*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + 127*d**2*(d + e*x)/(15*e**6*sqrt(d**2 - e**2*x**2)) - 13*d**2*atan(e*x/sqrt(d**2 - e**2*x**2))/(2*e**6) + 3*d*sqrt(d**2 - e**2*x**2)/e**6 + x*sqrt(d**2 - e**2*x**2)/(2*e**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_84():
    f = x**4*(d + e*x)**3/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = d**3*(d + e*x)**3/(5*e**5*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 6*d**2*(d + e*x)**2/(5*e**5*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + 24*d*(d + e*x)/(5*e**5*sqrt(d**2 - e**2*x**2)) - 3*d*atan(e*x/sqrt(d**2 - e**2*x**2))/e**5 + sqrt(d**2 - e**2*x**2)/e**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_85():
    f = x**3*(d + e*x)**3/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = d**2*(d + e*x)**3/(5*e**4*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 13*d*(d + e*x)**2/(15*e**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (32*d + 32*e*x)/(15*e**4*sqrt(d**2 - e**2*x**2)) - atan(e*x/sqrt(d**2 - e**2*x**2))/e**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_86():
    f = x**2*(d + e*x)**3/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = d*(d + e*x)**3/(5*e**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 8*(d + e*x)**2/(15*e**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (7*d + 7*e*x)/(15*d*e**3*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_87():
    f = x*(d + e*x)**3/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = (d + e*x)**3/(5*e**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - (2*d + 2*e*x)/(5*e**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - x/(5*d**2*e*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_88():
    f = (d + e*x)**3/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = sqrt(d**2 - e**2*x**2)/(5*d*e*(d - e*x)**3) + 2*sqrt(d**2 - e**2*x**2)/(15*d**2*e*(d - e*x)**2) + 2*sqrt(d**2 - e**2*x**2)/(15*d**3*e*(d - e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_89():
    f = (d + e*x)**3/(x*(d**2 - e**2*x**2)**(sympy.S(7)/2))
    F = (4*d + 4*e*x)/(5*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + (5*d + 11*e*x)/(15*d**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (15*d + 22*e*x)/(15*d**4*sqrt(d**2 - e**2*x**2)) - atanh(sqrt(d**2 - e**2*x**2)/d)/d**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_90():
    f = (d + e*x)**3/(x**2*(d**2 - e**2*x**2)**(sympy.S(7)/2))
    F = 4*e*(d + e*x)/(5*d*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + e*(5*d + 7*e*x)/(5*d**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + e*(15*d + 19*e*x)/(5*d**5*sqrt(d**2 - e**2*x**2)) - 3*e*atanh(sqrt(d**2 - e**2*x**2)/d)/d**5 - sqrt(d**2 - e**2*x**2)/(d**5*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_91():
    f = (d + e*x)**3/(x**3*(d**2 - e**2*x**2)**(sympy.S(7)/2))
    F = 4*e**2*(d + e*x)/(5*d**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + e**2*(25*d + 31*e*x)/(15*d**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - sqrt(d**2 - e**2*x**2)/(2*d**5*x**2) + e**2*(90*d + 107*e*x)/(15*d**6*sqrt(d**2 - e**2*x**2)) - 13*e**2*atanh(sqrt(d**2 - e**2*x**2)/d)/(2*d**6) - 3*e*sqrt(d**2 - e**2*x**2)/(d**6*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_92():
    f = x**4*sqrt(d**2 - e**2*x**2)/(d + e*x)
    F = 3*d**5*atan(e*x/sqrt(d**2 - e**2*x**2))/(8*e**5) + d**3*(64*d - 45*e*x)*sqrt(d**2 - e**2*x**2)/(120*e**5) + 4*d**2*x**2*sqrt(d**2 - e**2*x**2)/(15*e**3) - d*x**3*sqrt(d**2 - e**2*x**2)/(4*e**2) + x**4*sqrt(d**2 - e**2*x**2)/(5*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_93():
    f = x**3*sqrt(d**2 - e**2*x**2)/(d + e*x)
    F = -3*d**4*atan(e*x/sqrt(d**2 - e**2*x**2))/(8*e**4) - d**2*(16*d - 9*e*x)*sqrt(d**2 - e**2*x**2)/(24*e**4) - d*x**2*sqrt(d**2 - e**2*x**2)/(3*e**2) + x**3*sqrt(d**2 - e**2*x**2)/(4*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_94():
    f = x**2*sqrt(d**2 - e**2*x**2)/(d + e*x)
    F = d**3*atan(e*x/sqrt(d**2 - e**2*x**2))/(2*e**3) + d*(2*d - e*x)*sqrt(d**2 - e**2*x**2)/(2*e**3) - (d**2 - e**2*x**2)**(sympy.S(3)/2)/(3*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_95():
    f = x*sqrt(d**2 - e**2*x**2)/(d + e*x)
    F = -d**2*atan(e*x/sqrt(d**2 - e**2*x**2))/(2*e**2) - (2*d - e*x)*sqrt(d**2 - e**2*x**2)/(2*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_96():
    f = sqrt(d**2 - e**2*x**2)/(d + e*x)
    F = d*atan(e*x/sqrt(d**2 - e**2*x**2))/e + sqrt(d**2 - e**2*x**2)/e
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_97():
    f = sqrt(d**2 - e**2*x**2)/(x*(d + e*x))
    F = -atan(e*x/sqrt(d**2 - e**2*x**2)) - atanh(sqrt(d**2 - e**2*x**2)/d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_98():
    f = sqrt(d**2 - e**2*x**2)/(x**2*(d + e*x))
    F = e*atanh(sqrt(d**2 - e**2*x**2)/d)/d - sqrt(d**2 - e**2*x**2)/(d*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_99():
    f = sqrt(d**2 - e**2*x**2)/(x**3*(d + e*x))
    F = -sqrt(d**2 - e**2*x**2)/(2*d*x**2) - e**2*atanh(sqrt(d**2 - e**2*x**2)/d)/(2*d**2) + e*sqrt(d**2 - e**2*x**2)/(d**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_100():
    f = sqrt(d**2 - e**2*x**2)/(x**4*(d + e*x))
    F = -sqrt(d**2 - e**2*x**2)/(3*d*x**3) + e*sqrt(d**2 - e**2*x**2)/(2*d**2*x**2) + e**3*atanh(sqrt(d**2 - e**2*x**2)/d)/(2*d**3) - 2*e**2*sqrt(d**2 - e**2*x**2)/(3*d**3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_101():
    f = sqrt(d**2 - e**2*x**2)/(x**5*(d + e*x))
    F = -sqrt(d**2 - e**2*x**2)/(4*d*x**4) + e*sqrt(d**2 - e**2*x**2)/(3*d**2*x**3) - 3*e**2*sqrt(d**2 - e**2*x**2)/(8*d**3*x**2) - 3*e**4*atanh(sqrt(d**2 - e**2*x**2)/d)/(8*d**4) + 2*e**3*sqrt(d**2 - e**2*x**2)/(3*d**4*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_102():
    f = x**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(d + e*x)
    F = d**5*atan(e*x/sqrt(d**2 - e**2*x**2))/(8*e**3) + d**3*x*sqrt(d**2 - e**2*x**2)/(8*e**2) + d*(4*d - 3*e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(12*e**3) - (d**2 - e**2*x**2)**(sympy.S(5)/2)/(5*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_103():
    f = x**4*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(d + e*x)
    F = 3*d**9*atan(e*x/sqrt(d**2 - e**2*x**2))/(128*e**5) + 3*d**7*x*sqrt(d**2 - e**2*x**2)/(128*e**4) + d**5*x*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(64*e**4) + d**3*(128*d - 315*e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(5040*e**5) + 4*d**2*x**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(63*e**3) - d*x**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(8*e**2) + x**4*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(9*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_104():
    f = x**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(d + e*x)
    F = -3*d**8*atan(e*x/sqrt(d**2 - e**2*x**2))/(128*e**4) - 3*d**6*x*sqrt(d**2 - e**2*x**2)/(128*e**3) - d**4*x*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(64*e**3) - d**2*(32*d - 35*e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(560*e**4) - d*x**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(7*e**2) + x**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(8*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_105():
    f = x**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(d + e*x)
    F = d**7*atan(e*x/sqrt(d**2 - e**2*x**2))/(16*e**3) + d**5*x*sqrt(d**2 - e**2*x**2)/(16*e**2) + d**3*x*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(24*e**2) + d*(6*d - 5*e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(30*e**3) - (d**2 - e**2*x**2)**(sympy.S(7)/2)/(7*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_106():
    f = x*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(d + e*x)
    F = -d**6*atan(e*x/sqrt(d**2 - e**2*x**2))/(16*e**2) - d**4*x*sqrt(d**2 - e**2*x**2)/(16*e) - d**2*x*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(24*e) - (6*d - 5*e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(30*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_107():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(d + e*x)
    F = 3*d**5*atan(e*x/sqrt(d**2 - e**2*x**2))/(8*e) + 3*d**3*x*sqrt(d**2 - e**2*x**2)/8 + d*x*(d**2 - e**2*x**2)**(sympy.S(3)/2)/4 + (d**2 - e**2*x**2)**(sympy.S(5)/2)/(5*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_108():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x*(d + e*x))
    F = -3*d**4*atan(e*x/sqrt(d**2 - e**2*x**2))/8 - d**4*atanh(sqrt(d**2 - e**2*x**2)/d) + d**2*(8*d - 3*e*x)*sqrt(d**2 - e**2*x**2)/8 + (d/3 - e*x/4)*(d**2 - e**2*x**2)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_109():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x**2*(d + e*x))
    F = -3*d**3*e*atan(e*x/sqrt(d**2 - e**2*x**2))/2 + d**3*e*atanh(sqrt(d**2 - e**2*x**2)/d) - d*e*(2*d + 3*e*x)*sqrt(d**2 - e**2*x**2)/2 - (3*d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_110():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x**3*(d + e*x))
    F = 3*d**2*e**2*atan(e*x/sqrt(d**2 - e**2*x**2))/2 + 3*d**2*e**2*atanh(sqrt(d**2 - e**2*x**2)/d)/2 + 3*d*e*(d - e*x)*sqrt(d**2 - e**2*x**2)/(2*x) - (d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_111():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x**4*(d + e*x))
    F = d*e**3*atan(e*x/sqrt(d**2 - e**2*x**2)) - 3*d*e**3*atanh(sqrt(d**2 - e**2*x**2)/d)/2 + e**2*(2*d + 3*e*x)*sqrt(d**2 - e**2*x**2)/(2*x) - (2*d - 3*e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(6*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_112():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x**5*(d + e*x))
    F = -e**4*atan(e*x/sqrt(d**2 - e**2*x**2)) - 3*e**4*atanh(sqrt(d**2 - e**2*x**2)/d)/8 + e**2*(3*d - 8*e*x)*sqrt(d**2 - e**2*x**2)/(8*x**2) - (3*d - 4*e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(12*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_113():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x**6*(d + e*x))
    F = -3*e**3*sqrt(d**2 - e**2*x**2)/(8*x**2) + e*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(4*x**4) + 3*e**5*atanh(sqrt(d**2 - e**2*x**2)/d)/(8*d) - (d**2 - e**2*x**2)**(sympy.S(5)/2)/(5*d*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_114():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x**7*(d + e*x))
    F = e**4*sqrt(d**2 - e**2*x**2)/(16*d*x**2) - e**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(24*d*x**4) - (d**2 - e**2*x**2)**(sympy.S(5)/2)/(6*d*x**6) - e**6*atanh(sqrt(d**2 - e**2*x**2)/d)/(16*d**2) + e*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(5*d**2*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_115():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x**8*(d + e*x))
    F = -(d**2 - e**2*x**2)**(sympy.S(5)/2)/(7*d*x**7) - e**5*sqrt(d**2 - e**2*x**2)/(16*d**2*x**2) + e**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(24*d**2*x**4) + e*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(6*d**2*x**6) + e**7*atanh(sqrt(d**2 - e**2*x**2)/d)/(16*d**3) - 2*e**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(35*d**3*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_116():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x**9*(d + e*x))
    F = -(d**2 - e**2*x**2)**(sympy.S(5)/2)/(8*d*x**8) + e*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(7*d**2*x**7) + 3*e**6*sqrt(d**2 - e**2*x**2)/(128*d**3*x**2) - e**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(64*d**3*x**4) - e**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(16*d**3*x**6) - 3*e**8*atanh(sqrt(d**2 - e**2*x**2)/d)/(128*d**4) + 2*e**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(35*d**4*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_117():
    f = x*sqrt(1 - x**2)/(x + 1)
    F = sqrt(1 - x**2)*(x/2 - 1) - asin(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_118():
    f = (-a**2*x**2 + 1)**(sympy.S(3)/2)/(x**2*(-a*x + 1))
    F = -a*asin(a*x) - a*atanh(sqrt(-a**2*x**2 + 1)) - (-a*x + 1)*sqrt(-a**2*x**2 + 1)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_119():
    f = x**4/((d + e*x)*sqrt(d**2 - e**2*x**2))
    F = -3*d**3*atan(e*x/sqrt(d**2 - e**2*x**2))/(2*e**5) - d*(16*d - 9*e*x)*sqrt(d**2 - e**2*x**2)/(6*e**5) + x**3*(d - e*x)/(e**2*sqrt(d**2 - e**2*x**2)) - 4*x**2*sqrt(d**2 - e**2*x**2)/(3*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_120():
    f = x**3/((d + e*x)*sqrt(d**2 - e**2*x**2))
    F = 3*d**2*atan(e*x/sqrt(d**2 - e**2*x**2))/(2*e**4) + x**2*(d - e*x)/(e**2*sqrt(d**2 - e**2*x**2)) + (4*d - 3*e*x)*sqrt(d**2 - e**2*x**2)/(2*e**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_121():
    f = x**2/((d + e*x)*sqrt(d**2 - e**2*x**2))
    F = -d*atan(e*x/sqrt(d**2 - e**2*x**2))/e**3 - d*sqrt(d**2 - e**2*x**2)/(e**3*(d + e*x)) - sqrt(d**2 - e**2*x**2)/e**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_122():
    f = x/((d + e*x)*sqrt(d**2 - e**2*x**2))
    F = atan(e*x/sqrt(d**2 - e**2*x**2))/e**2 + sqrt(d**2 - e**2*x**2)/(e**2*(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_123():
    f = 1/((d + e*x)*sqrt(d**2 - e**2*x**2))
    F = -sqrt(d**2 - e**2*x**2)/(d*e*(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_124():
    f = 1/(x*(d + e*x)*sqrt(d**2 - e**2*x**2))
    F = -atanh(sqrt(d**2 - e**2*x**2)/d)/d**2 + sqrt(d**2 - e**2*x**2)/(d**2*(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_125():
    f = 1/(x**2*(d + e*x)*sqrt(d**2 - e**2*x**2))
    F = sqrt(d**2 - e**2*x**2)/(d**2*x*(d + e*x)) + e*atanh(sqrt(d**2 - e**2*x**2)/d)/d**3 - 2*sqrt(d**2 - e**2*x**2)/(d**3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_126():
    f = 1/(x**3*(d + e*x)*sqrt(d**2 - e**2*x**2))
    F = sqrt(d**2 - e**2*x**2)/(d**2*x**2*(d + e*x)) - 3*sqrt(d**2 - e**2*x**2)/(2*d**3*x**2) - 3*e**2*atanh(sqrt(d**2 - e**2*x**2)/d)/(2*d**4) + 2*e*sqrt(d**2 - e**2*x**2)/(d**4*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_127():
    f = x**5/((d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2))
    F = -5*d**2*atan(e*x/sqrt(d**2 - e**2*x**2))/(2*e**6) + x**4*(d - e*x)/(3*e**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - x**2*(4*d - 5*e*x)/(3*e**4*sqrt(d**2 - e**2*x**2)) - (16*d - 15*e*x)*sqrt(d**2 - e**2*x**2)/(6*e**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_128():
    f = x**4/((d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2))
    F = d*atan(e*x/sqrt(d**2 - e**2*x**2))/e**5 + x**3*(d - e*x)/(3*e**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - x*(3*d - 4*e*x)/(3*e**4*sqrt(d**2 - e**2*x**2)) + 8*sqrt(d**2 - e**2*x**2)/(3*e**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_129():
    f = x**3/((d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2))
    F = x**2*(d - e*x)/(3*e**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - (2*d - 3*e*x)/(3*e**4*sqrt(d**2 - e**2*x**2)) - atan(e*x/sqrt(d**2 - e**2*x**2))/e**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_130():
    f = x**2/((d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2))
    F = 2/(3*e**3*sqrt(d**2 - e**2*x**2)) - x**2/(3*d*e*(d + e*x)*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_131():
    f = x/((d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2))
    F = 1/(3*e**2*(d + e*x)*sqrt(d**2 - e**2*x**2)) + x/(3*d**2*e*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_132():
    f = 1/((d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2))
    F = -1/(3*d*e*(d + e*x)*sqrt(d**2 - e**2*x**2)) + 2*x/(3*d**3*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_133():
    f = 1/(x*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2))
    F = 1/(3*d**2*(d + e*x)*sqrt(d**2 - e**2*x**2)) + (3*d - 2*e*x)/(3*d**4*sqrt(d**2 - e**2*x**2)) - atanh(sqrt(d**2 - e**2*x**2)/d)/d**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_134():
    f = 1/(x**2*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2))
    F = 1/(3*d**2*x*(d + e*x)*sqrt(d**2 - e**2*x**2)) + (4*d - 3*e*x)/(3*d**4*x*sqrt(d**2 - e**2*x**2)) + e*atanh(sqrt(d**2 - e**2*x**2)/d)/d**5 - 8*sqrt(d**2 - e**2*x**2)/(3*d**5*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_135():
    f = 1/(x**3*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2))
    F = 1/(3*d**2*x**2*(d + e*x)*sqrt(d**2 - e**2*x**2)) + (5*d - 4*e*x)/(3*d**4*x**2*sqrt(d**2 - e**2*x**2)) - 5*sqrt(d**2 - e**2*x**2)/(2*d**5*x**2) - 5*e**2*atanh(sqrt(d**2 - e**2*x**2)/d)/(2*d**6) + 8*e*sqrt(d**2 - e**2*x**2)/(3*d**6*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_136():
    f = x**7/((d + e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2))
    F = 7*d**2*atan(e*x/sqrt(d**2 - e**2*x**2))/(2*e**8) + x**6*(d - e*x)/(5*e**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - x**4*(6*d - 7*e*x)/(15*e**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + x**2*(24*d - 35*e*x)/(15*e**6*sqrt(d**2 - e**2*x**2)) + (32*d - 35*e*x)*sqrt(d**2 - e**2*x**2)/(10*e**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_137():
    f = x**6/((d + e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2))
    F = -d*atan(e*x/sqrt(d**2 - e**2*x**2))/e**7 + x**5*(d - e*x)/(5*e**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - x**3*(5*d - 6*e*x)/(15*e**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + x*(5*d - 8*e*x)/(5*e**6*sqrt(d**2 - e**2*x**2)) - 16*sqrt(d**2 - e**2*x**2)/(5*e**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_138():
    f = x**5/((d + e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2))
    F = x**4*(d - e*x)/(5*e**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - x**2*(4*d - 5*e*x)/(15*e**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (8*d - 15*e*x)/(15*e**6*sqrt(d**2 - e**2*x**2)) + atan(e*x/sqrt(d**2 - e**2*x**2))/e**6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_139():
    f = x**4/((d + e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2))
    F = 4*d**2/(15*e**5*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - 4/(5*e**5*sqrt(d**2 - e**2*x**2)) - x**4*(d - e*x)/(5*d*e*(d**2 - e**2*x**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_140():
    f = x**3/((d + e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2))
    F = x**2*(d - e*x)/(5*e**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - (2*d - 3*e*x)/(15*e**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - x/(5*d**2*e**3*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_141():
    f = x**2/((d + e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2))
    F = -x**2/(5*d*e*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (2*d + 2*e*x)/(15*d*e**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - 2*x/(15*d**3*e**2*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_142():
    f = x/((d + e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2))
    F = 1/(5*e**2*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + x/(15*d**2*e*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + 2*x/(15*d**4*e*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_143():
    f = 1/((d + e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2))
    F = -1/(5*d*e*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + 4*x/(15*d**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + 8*x/(15*d**5*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_144():
    f = 1/(x*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2))
    F = 1/(5*d**2*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (5*d - 4*e*x)/(15*d**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (15*d - 8*e*x)/(15*d**6*sqrt(d**2 - e**2*x**2)) - atanh(sqrt(d**2 - e**2*x**2)/d)/d**6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_145():
    f = 1/(x**2*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2))
    F = 1/(5*d**2*x*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (6*d - 5*e*x)/(15*d**4*x*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (8*d - 5*e*x)/(5*d**6*x*sqrt(d**2 - e**2*x**2)) + e*atanh(sqrt(d**2 - e**2*x**2)/d)/d**7 - 16*sqrt(d**2 - e**2*x**2)/(5*d**7*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_146():
    f = 1/(x**3*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2))
    F = 1/(5*d**2*x**2*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (7*d - 6*e*x)/(15*d**4*x**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (35*d - 24*e*x)/(15*d**6*x**2*sqrt(d**2 - e**2*x**2)) - 7*sqrt(d**2 - e**2*x**2)/(2*d**7*x**2) - 7*e**2*atanh(sqrt(d**2 - e**2*x**2)/d)/(2*d**8) + 16*e*sqrt(d**2 - e**2*x**2)/(5*d**8*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_147():
    f = 1/(x**4*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2))
    F = 1/(5*d**2*x**3*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (8*d - 7*e*x)/(15*d**4*x**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (48*d - 35*e*x)/(15*d**6*x**3*sqrt(d**2 - e**2*x**2)) - 64*sqrt(d**2 - e**2*x**2)/(15*d**7*x**3) + 7*e*sqrt(d**2 - e**2*x**2)/(2*d**8*x**2) + 7*e**3*atanh(sqrt(d**2 - e**2*x**2)/d)/(2*d**9) - 128*e**2*sqrt(d**2 - e**2*x**2)/(15*d**9*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_148():
    f = x**3/((d + e*x)*(d**2 - e**2*x**2)**(sympy.S(7)/2))
    F = x**2*(d - e*x)/(7*e**2*(d**2 - e**2*x**2)**(sympy.S(7)/2)) - (2*d - 3*e*x)/(35*e**4*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - x/(35*d**2*e**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - 2*x/(35*d**4*e**3*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_149():
    f = x**2/((d + e*x)*(d**2 - e**2*x**2)**(sympy.S(7)/2))
    F = -x**2/(7*d*e*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + (2*d + 4*e*x)/(35*d*e**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 4*x/(105*d**3*e**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - 8*x/(105*d**5*e**2*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_150():
    f = x**3/((a*x + 1)*sqrt(-a**2*x**2 + 1))
    F = x**2*(-a*x + 1)/(a**2*sqrt(-a**2*x**2 + 1)) + (-3*a*x + 4)*sqrt(-a**2*x**2 + 1)/(2*a**4) + 3*asin(a*x)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_151():
    f = x**2/((a*x + 1)*sqrt(-a**2*x**2 + 1))
    F = -sqrt(-a**2*x**2 + 1)/a**3 - asin(a*x)/a**3 - sqrt(-a**2*x**2 + 1)/(a**3*(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_152():
    f = x/((a*x + 1)*sqrt(-a**2*x**2 + 1))
    F = asin(a*x)/a**2 + sqrt(-a**2*x**2 + 1)/(a**2*(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_153():
    f = 1/((a*x + 1)*sqrt(-a**2*x**2 + 1))
    F = -sqrt(-a**2*x**2 + 1)/(a*(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_154():
    f = 1/(x*(-a*x + 1)*sqrt(-a**2*x**2 + 1))
    F = -atanh(sqrt(-a**2*x**2 + 1)) + sqrt(-a**2*x**2 + 1)/(-a*x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_155():
    f = 1/(x**2*(-a*x + 1)*sqrt(-a**2*x**2 + 1))
    F = -a*atanh(sqrt(-a**2*x**2 + 1)) - 2*sqrt(-a**2*x**2 + 1)/x + sqrt(-a**2*x**2 + 1)/(x*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_156():
    f = 1/(x**3*(-a*x + 1)*sqrt(-a**2*x**2 + 1))
    F = -3*a**2*atanh(sqrt(-a**2*x**2 + 1))/2 - 2*a*sqrt(-a**2*x**2 + 1)/x - 3*sqrt(-a**2*x**2 + 1)/(2*x**2) + sqrt(-a**2*x**2 + 1)/(x**2*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_157():
    f = x**5*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(d + e*x)**2
    F = -5*d**9*atan(e*x/sqrt(d**2 - e**2*x**2))/(64*e**6) - 5*d**7*x*sqrt(d**2 - e**2*x**2)/(64*e**5) - d**5*(256*d - 315*e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(2016*e**6) - 4*d**4*x**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(21*e**4) + 5*d**3*x**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(24*e**3) - 5*d**2*x**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(21*e**2) + d*x**5*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(4*e) - x**6*(d**2 - e**2*x**2)**(sympy.S(3)/2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_158():
    f = x**4*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(d + e*x)**2
    F = 13*d**8*atan(e*x/sqrt(d**2 - e**2*x**2))/(128*e**5) + 13*d**6*x*sqrt(d**2 - e**2*x**2)/(128*e**4) + d**4*(1024*d - 1365*e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(6720*e**5) + 8*d**3*x**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(35*e**3) - 13*d**2*x**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(48*e**2) + 2*d*x**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(7*e) - x**5*(d**2 - e**2*x**2)**(sympy.S(3)/2)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_159():
    f = x**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(d + e*x)**2
    F = -d**7*atan(e*x/sqrt(d**2 - e**2*x**2))/(8*e**4) - d**5*x*sqrt(d**2 - e**2*x**2)/(8*e**3) - d**3*(88*d - 105*e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(420*e**4) - 11*d**2*x**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(35*e**2) + d*x**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(3*e) - x**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_160():
    f = x**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(d + e*x)**2
    F = 3*d**6*atan(e*x/sqrt(d**2 - e**2*x**2))/(16*e**3) + 3*d**4*x*sqrt(d**2 - e**2*x**2)/(16*e**2) + d**2*(32*d - 45*e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(120*e**3) + 2*d*x**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(5*e) - x**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_161():
    f = x*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(d + e*x)**2
    F = -d**5*atan(e*x/sqrt(d**2 - e**2*x**2))/(4*e**2) - d**3*x*sqrt(d**2 - e**2*x**2)/(4*e) - d*x*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(6*e) - 2*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(15*e**2) - (d**2 - e**2*x**2)**(sympy.S(7)/2)/(3*e**2*(d + e*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_162():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(d + e*x)**2
    F = 5*d**4*atan(e*x/sqrt(d**2 - e**2*x**2))/(8*e) + 5*d**2*x*sqrt(d**2 - e**2*x**2)/8 + 5*d*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(12*e) + (d - e*x)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(4*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_163():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x*(d + e*x)**2)
    F = -d**3*atan(e*x/sqrt(d**2 - e**2*x**2)) - d**3*atanh(sqrt(d**2 - e**2*x**2)/d) + d*(d - e*x)*sqrt(d**2 - e**2*x**2) - (d**2 - e**2*x**2)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_164():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x**2*(d + e*x)**2)
    F = -d**2*e*atan(e*x/sqrt(d**2 - e**2*x**2))/2 + 2*d**2*e*atanh(sqrt(d**2 - e**2*x**2)/d) - e*(4*d + e*x)*sqrt(d**2 - e**2*x**2)/2 - (d**2 - e**2*x**2)**(sympy.S(3)/2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_165():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x**3*(d + e*x)**2)
    F = 2*d*e**2*atan(e*x/sqrt(d**2 - e**2*x**2)) - d*e**2*atanh(sqrt(d**2 - e**2*x**2)/d)/2 + e*(4*d + e*x)*sqrt(d**2 - e**2*x**2)/(2*x) - (d**2 - e**2*x**2)**(sympy.S(3)/2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_166():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x**4*(d + e*x)**2)
    F = -e**3*atan(e*x/sqrt(d**2 - e**2*x**2)) - e**3*atanh(sqrt(d**2 - e**2*x**2)/d) + e*(d - e*x)*sqrt(d**2 - e**2*x**2)/x**2 - (d**2 - e**2*x**2)**(sympy.S(3)/2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_167():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x**5*(d + e*x)**2)
    F = -5*e**2*sqrt(d**2 - e**2*x**2)/(8*x**2) - (d**2 - e**2*x**2)**(sympy.S(3)/2)/(4*x**4) + 5*e**4*atanh(sqrt(d**2 - e**2*x**2)/d)/(8*d) + 2*e*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(3*d*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_168():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x**6*(d + e*x)**2)
    F = -(d**2 - e**2*x**2)**(sympy.S(3)/2)/(5*x**5) + e**3*sqrt(d**2 - e**2*x**2)/(4*d*x**2) + e*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(2*d*x**4) - e**5*atanh(sqrt(d**2 - e**2*x**2)/d)/(4*d**2) - 7*e**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(15*d**2*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_169():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x**7*(d + e*x)**2)
    F = -(d**2 - e**2*x**2)**(sympy.S(3)/2)/(6*x**6) + 2*e*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(5*d*x**5) - 3*e**4*sqrt(d**2 - e**2*x**2)/(16*d**2*x**2) - 3*e**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(8*d**2*x**4) + 3*e**6*atanh(sqrt(d**2 - e**2*x**2)/d)/(16*d**3) + 4*e**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(15*d**3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_170():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x**8*(d + e*x)**2)
    F = -(d**2 - e**2*x**2)**(sympy.S(3)/2)/(7*x**7) + e*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(3*d*x**6) - 11*e**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(35*d**2*x**5) + e**5*sqrt(d**2 - e**2*x**2)/(8*d**3*x**2) + e**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(4*d**3*x**4) - e**7*atanh(sqrt(d**2 - e**2*x**2)/d)/(8*d**4) - 22*e**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(105*d**4*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_171():
    f = x**4/((d + e*x)**2*(d**2 - e**2*x**2)**(sympy.S(3)/2))
    F = -d**3*(d - e*x)**2/(5*e**5*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + 17*d**2*(d - e*x)/(15*e**5*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - (30*d - 26*e*x)/(15*e**5*sqrt(d**2 - e**2*x**2)) - atan(e*x/sqrt(d**2 - e**2*x**2))/e**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_172():
    f = x**3/((d + e*x)**2*(d**2 - e**2*x**2)**(sympy.S(3)/2))
    F = d**2*(d - e*x)**2/(5*e**4*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 4*d*(d - e*x)/(5*e**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (5*d - 2*e*x)/(5*d*e**4*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_173():
    f = x/((d + e*x)**2*(d**2 - e**2*x**2)**(sympy.S(3)/2))
    F = 1/(5*e**2*(d + e*x)**2*sqrt(d**2 - e**2*x**2)) - 2/(15*d*e**2*(d + e*x)*sqrt(d**2 - e**2*x**2)) + 4*x/(15*d**3*e*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_174():
    f = 1/((d + e*x)**2*(d**2 - e**2*x**2)**(sympy.S(3)/2))
    F = -1/(5*d*e*(d + e*x)**2*sqrt(d**2 - e**2*x**2)) - 1/(5*d**2*e*(d + e*x)*sqrt(d**2 - e**2*x**2)) + 2*x/(5*d**4*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_175():
    f = 1/(x*(d + e*x)**2*(d**2 - e**2*x**2)**(sympy.S(3)/2))
    F = (2*d - 2*e*x)/(5*d*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + (5*d - 8*e*x)/(15*d**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (15*d - 16*e*x)/(15*d**5*sqrt(d**2 - e**2*x**2)) - atanh(sqrt(d**2 - e**2*x**2)/d)/d**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_176():
    f = 1/(x**2*(d + e*x)**2*(d**2 - e**2*x**2)**(sympy.S(3)/2))
    F = -2*e*(d - e*x)/(5*d**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - e*(10*d - 13*e*x)/(15*d**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - e*(30*d - 41*e*x)/(15*d**6*sqrt(d**2 - e**2*x**2)) + 2*e*atanh(sqrt(d**2 - e**2*x**2)/d)/d**6 - sqrt(d**2 - e**2*x**2)/(d**6*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_177():
    f = 1/(x**3*(d + e*x)**2*(d**2 - e**2*x**2)**(sympy.S(3)/2))
    F = 2*e**2*(d - e*x)/(5*d**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + e**2*(5*d - 6*e*x)/(5*d**5*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - sqrt(d**2 - e**2*x**2)/(2*d**6*x**2) + 2*e**2*(10*d - 11*e*x)/(5*d**7*sqrt(d**2 - e**2*x**2)) - 9*e**2*atanh(sqrt(d**2 - e**2*x**2)/d)/(2*d**7) + 2*e*sqrt(d**2 - e**2*x**2)/(d**7*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_178():
    f = x**5/((d + e*x)**3*sqrt(d**2 - e**2*x**2))
    F = d**4*(d - e*x)**3/(5*e**6*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 23*d**3*(d - e*x)**2/(15*e**6*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + 127*d**2*(d - e*x)/(15*e**6*sqrt(d**2 - e**2*x**2)) + 13*d**2*atan(e*x/sqrt(d**2 - e**2*x**2))/(2*e**6) + 3*d*sqrt(d**2 - e**2*x**2)/e**6 - x*sqrt(d**2 - e**2*x**2)/(2*e**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_179():
    f = x**4/((d + e*x)**3*sqrt(d**2 - e**2*x**2))
    F = -d**3*(d - e*x)**3/(5*e**5*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + 6*d**2*(d - e*x)**2/(5*e**5*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - 24*d*(d - e*x)/(5*e**5*sqrt(d**2 - e**2*x**2)) - 3*d*atan(e*x/sqrt(d**2 - e**2*x**2))/e**5 - sqrt(d**2 - e**2*x**2)/e**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_180():
    f = x**3/((d + e*x)**3*sqrt(d**2 - e**2*x**2))
    F = d**2*(d - e*x)**3/(5*e**4*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 13*d*(d - e*x)**2/(15*e**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (32*d - 32*e*x)/(15*e**4*sqrt(d**2 - e**2*x**2)) + atan(e*x/sqrt(d**2 - e**2*x**2))/e**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_181():
    f = x**2/((d + e*x)**3*sqrt(d**2 - e**2*x**2))
    F = -d*sqrt(d**2 - e**2*x**2)/(5*e**3*(d + e*x)**3) + 8*sqrt(d**2 - e**2*x**2)/(15*e**3*(d + e*x)**2) - 7*sqrt(d**2 - e**2*x**2)/(15*d*e**3*(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_182():
    f = x/((d + e*x)**3*sqrt(d**2 - e**2*x**2))
    F = sqrt(d**2 - e**2*x**2)/(5*e**2*(d + e*x)**3) - sqrt(d**2 - e**2*x**2)/(5*d*e**2*(d + e*x)**2) - sqrt(d**2 - e**2*x**2)/(5*d**2*e**2*(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_183():
    f = 1/((d + e*x)**3*sqrt(d**2 - e**2*x**2))
    F = -sqrt(d**2 - e**2*x**2)/(5*d*e*(d + e*x)**3) - 2*sqrt(d**2 - e**2*x**2)/(15*d**2*e*(d + e*x)**2) - 2*sqrt(d**2 - e**2*x**2)/(15*d**3*e*(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_184():
    f = 1/(x*(d + e*x)**3*sqrt(d**2 - e**2*x**2))
    F = (4*d - 4*e*x)/(5*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + (5*d - 11*e*x)/(15*d**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (15*d - 22*e*x)/(15*d**4*sqrt(d**2 - e**2*x**2)) - atanh(sqrt(d**2 - e**2*x**2)/d)/d**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_185():
    f = 1/(x**2*(d + e*x)**3*sqrt(d**2 - e**2*x**2))
    F = -4*e*(d - e*x)/(5*d*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - e*(5*d - 7*e*x)/(5*d**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - e*(15*d - 19*e*x)/(5*d**5*sqrt(d**2 - e**2*x**2)) + 3*e*atanh(sqrt(d**2 - e**2*x**2)/d)/d**5 - sqrt(d**2 - e**2*x**2)/(d**5*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_186():
    f = 1/(x**3*(d + e*x)**3*sqrt(d**2 - e**2*x**2))
    F = 4*e**2*(d - e*x)/(5*d**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + e**2*(25*d - 31*e*x)/(15*d**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - sqrt(d**2 - e**2*x**2)/(2*d**5*x**2) + e**2*(90*d - 107*e*x)/(15*d**6*sqrt(d**2 - e**2*x**2)) - 13*e**2*atanh(sqrt(d**2 - e**2*x**2)/d)/(2*d**6) + 3*e*sqrt(d**2 - e**2*x**2)/(d**6*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_187():
    f = x**5*sqrt(d**2 - e**2*x**2)/(d + e*x)**4
    F = d**4*(d - e*x)**4/(5*e**6*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 8*d**3*(d - e*x)**3/(5*e**6*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + 18*d**3*atan(e*x/sqrt(d**2 - e**2*x**2))/e**6 + 10*d**2*(d - e*x)**2/(e**6*sqrt(d**2 - e**2*x**2)) + 59*d**2*sqrt(d**2 - e**2*x**2)/(3*e**6) - 2*d*x*sqrt(d**2 - e**2*x**2)/e**5 + x**2*sqrt(d**2 - e**2*x**2)/(3*e**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_188():
    f = x**4*sqrt(d**2 - e**2*x**2)/(d + e*x)**4
    F = -d**3*(d - e*x)**4/(5*e**5*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + 19*d**2*(d - e*x)**3/(15*e**5*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - 19*d**2*atan(e*x/sqrt(d**2 - e**2*x**2))/(2*e**5) - 6*d*(d - e*x)**2/(e**5*sqrt(d**2 - e**2*x**2)) - (20*d - e*x)*sqrt(d**2 - e**2*x**2)/(2*e**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_189():
    f = x**3*sqrt(d**2 - e**2*x**2)/(d + e*x)**4
    F = d**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(5*e**4*(d + e*x)**4) + 4*d*atan(e*x/sqrt(d**2 - e**2*x**2))/e**4 + 8*d*sqrt(d**2 - e**2*x**2)/(e**4*(d + e*x)) - 14*d*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(15*e**4*(d + e*x)**3) - (d**2 - e**2*x**2)**(sympy.S(3)/2)/(e**4*(d + e*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_190():
    f = x**2*sqrt(d**2 - e**2*x**2)/(d + e*x)**4
    F = -d*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(5*e**3*(d + e*x)**4) - atan(e*x/sqrt(d**2 - e**2*x**2))/e**3 - 2*sqrt(d**2 - e**2*x**2)/(e**3*(d + e*x)) + 3*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(5*e**3*(d + e*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_191():
    f = x*sqrt(d**2 - e**2*x**2)/(d + e*x)**4
    F = (d**2 - e**2*x**2)**(sympy.S(3)/2)/(5*e**2*(d + e*x)**4) - 4*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(15*d*e**2*(d + e*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_192():
    f = sqrt(d**2 - e**2*x**2)/(d + e*x)**4
    F = -(d**2 - e**2*x**2)**(sympy.S(3)/2)/(5*d*e*(d + e*x)**4) - (d**2 - e**2*x**2)**(sympy.S(3)/2)/(15*d**2*e*(d + e*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_193():
    f = sqrt(d**2 - e**2*x**2)/(x*(d + e*x)**4)
    F = 8*d*(d - e*x)/(5*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 4*e*x/(5*d*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (5*d - 8*e*x)/(5*d**3*sqrt(d**2 - e**2*x**2)) - atanh(sqrt(d**2 - e**2*x**2)/d)/d**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_194():
    f = sqrt(d**2 - e**2*x**2)/(x**2*(d + e*x)**4)
    F = -8*e*(d - e*x)/(5*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 4*e*(5*d - 8*e*x)/(15*d**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - e*(60*d - 79*e*x)/(15*d**4*sqrt(d**2 - e**2*x**2)) + 4*e*atanh(sqrt(d**2 - e**2*x**2)/d)/d**4 - sqrt(d**2 - e**2*x**2)/(d**4*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_195():
    f = sqrt(d**2 - e**2*x**2)/(x**3*(d + e*x)**4)
    F = 8*e**2*(d - e*x)/(5*d*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + 4*e**2*(10*d - 13*e*x)/(15*d**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - sqrt(d**2 - e**2*x**2)/(2*d**4*x**2) + e**2*(135*d - 164*e*x)/(15*d**5*sqrt(d**2 - e**2*x**2)) - 19*e**2*atanh(sqrt(d**2 - e**2*x**2)/d)/(2*d**5) + 4*e*sqrt(d**2 - e**2*x**2)/(d**5*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_196():
    f = sqrt(d**2 - e**2*x**2)/(x**4*(d + e*x)**4)
    F = -8*e**3*(d - e*x)/(5*d**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 4*e**3*(5*d - 6*e*x)/(5*d**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - sqrt(d**2 - e**2*x**2)/(3*d**4*x**3) + 2*e*sqrt(d**2 - e**2*x**2)/(d**5*x**2) - e**3*(80*d - 93*e*x)/(5*d**6*sqrt(d**2 - e**2*x**2)) + 18*e**3*atanh(sqrt(d**2 - e**2*x**2)/d)/d**6 - 29*e**2*sqrt(d**2 - e**2*x**2)/(3*d**6*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_197():
    f = x**5*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(d + e*x)**4
    F = 65*d**7*atan(e*x/sqrt(d**2 - e**2*x**2))/(4*e**6) + 515*d**6*sqrt(d**2 - e**2*x**2)/(21*e**6) - 49*d**5*x*sqrt(d**2 - e**2*x**2)/(4*e**5) + 121*d**4*x**2*sqrt(d**2 - e**2*x**2)/(21*e**4) + d**4*(d - e*x)**4/(e**6*sqrt(d**2 - e**2*x**2)) - 17*d**3*x**3*sqrt(d**2 - e**2*x**2)/(6*e**3) + 11*d**2*x**4*sqrt(d**2 - e**2*x**2)/(7*e**2) - 2*d*x**5*sqrt(d**2 - e**2*x**2)/(3*e) + x**6*sqrt(d**2 - e**2*x**2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_198():
    f = x**4*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(d + e*x)**4
    F = -239*d**6*atan(e*x/sqrt(d**2 - e**2*x**2))/(16*e**5) - 337*d**5*sqrt(d**2 - e**2*x**2)/(15*e**5) + 175*d**4*x*sqrt(d**2 - e**2*x**2)/(16*e**4) - 71*d**3*x**2*sqrt(d**2 - e**2*x**2)/(15*e**3) - d**3*(d - e*x)**4/(e**5*sqrt(d**2 - e**2*x**2)) + 47*d**2*x**3*sqrt(d**2 - e**2*x**2)/(24*e**2) - 4*d*x**4*sqrt(d**2 - e**2*x**2)/(5*e) + x**5*sqrt(d**2 - e**2*x**2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_199():
    f = x**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(d + e*x)**4
    F = 27*d**5*atan(e*x/sqrt(d**2 - e**2*x**2))/(2*e**4) + 101*d**4*sqrt(d**2 - e**2*x**2)/(5*e**4) - 19*d**3*x*sqrt(d**2 - e**2*x**2)/(2*e**3) + 18*d**2*x**2*sqrt(d**2 - e**2*x**2)/(5*e**2) + d**2*(d - e*x)**4/(e**4*sqrt(d**2 - e**2*x**2)) - d*x**3*sqrt(d**2 - e**2*x**2)/e + x**4*sqrt(d**2 - e**2*x**2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_200():
    f = x**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(d + e*x)**4
    F = -95*d**4*atan(e*x/sqrt(d**2 - e**2*x**2))/(8*e**3) - 95*d**3*sqrt(d**2 - e**2*x**2)/(8*e**3) - 95*d**2*(d - e*x)*sqrt(d**2 - e**2*x**2)/(24*e**3) - d*(d - e*x)**4/(e**3*sqrt(d**2 - e**2*x**2)) - 19*d*(d - e*x)**2*sqrt(d**2 - e**2*x**2)/(12*e**3) - (d - e*x)**3*sqrt(d**2 - e**2*x**2)/(4*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_201():
    f = x*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(d + e*x)**4
    F = 10*d**3*atan(e*x/sqrt(d**2 - e**2*x**2))/e**2 + 10*d*x*sqrt(d**2 - e**2*x**2)/e + 20*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(3*e**2) + 8*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(e**2*(d + e*x)**2) + (d**2 - e**2*x**2)**(sympy.S(7)/2)/(e**2*(d + e*x)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_202():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(d + e*x)**4
    F = -15*d**2*atan(e*x/sqrt(d**2 - e**2*x**2))/(2*e) - 15*d*sqrt(d**2 - e**2*x**2)/(2*e) - 5*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(2*e*(d + e*x)) - 2*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(e*(d + e*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_203():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x*(d + e*x)**4)
    F = 8*d*(d - e*x)/sqrt(d**2 - e**2*x**2) + 4*d*atan(e*x/sqrt(d**2 - e**2*x**2)) - d*atanh(sqrt(d**2 - e**2*x**2)/d) + sqrt(d**2 - e**2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_204():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x**2*(d + e*x)**4)
    F = -8*e*(d - e*x)/sqrt(d**2 - e**2*x**2) - e*atan(e*x/sqrt(d**2 - e**2*x**2)) + 4*e*atanh(sqrt(d**2 - e**2*x**2)/d) - sqrt(d**2 - e**2*x**2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_205():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x**3*(d + e*x)**4)
    F = -sqrt(d**2 - e**2*x**2)/(2*x**2) + 8*e**2*(d - e*x)/(d*sqrt(d**2 - e**2*x**2)) - 15*e**2*atanh(sqrt(d**2 - e**2*x**2)/d)/(2*d) + 4*e*sqrt(d**2 - e**2*x**2)/(d*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_206():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x**4*(d + e*x)**4)
    F = -sqrt(d**2 - e**2*x**2)/(3*x**3) + 2*e*sqrt(d**2 - e**2*x**2)/(d*x**2) - 8*e**3*(d - e*x)/(d**2*sqrt(d**2 - e**2*x**2)) + 10*e**3*atanh(sqrt(d**2 - e**2*x**2)/d)/d**2 - 23*e**2*sqrt(d**2 - e**2*x**2)/(3*d**2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_207():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x**5*(d + e*x)**4)
    F = -sqrt(d**2 - e**2*x**2)/(4*x**4) + 4*e*sqrt(d**2 - e**2*x**2)/(3*d*x**3) - 31*e**2*sqrt(d**2 - e**2*x**2)/(8*d**2*x**2) + 8*e**4*(d - e*x)/(d**3*sqrt(d**2 - e**2*x**2)) - 95*e**4*atanh(sqrt(d**2 - e**2*x**2)/d)/(8*d**3) + 32*e**3*sqrt(d**2 - e**2*x**2)/(3*d**3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_208():
    f = (d**2 - e**2*x**2)**(sympy.S(5)/2)/(x**6*(d + e*x)**4)
    F = -sqrt(d**2 - e**2*x**2)/(5*x**5) + e*sqrt(d**2 - e**2*x**2)/(d*x**4) - 13*e**2*sqrt(d**2 - e**2*x**2)/(5*d**2*x**3) + 11*e**3*sqrt(d**2 - e**2*x**2)/(2*d**3*x**2) - 8*e**5*(d - e*x)/(d**4*sqrt(d**2 - e**2*x**2)) + 27*e**5*atanh(sqrt(d**2 - e**2*x**2)/d)/(2*d**4) - 66*e**4*sqrt(d**2 - e**2*x**2)/(5*d**4*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_209():
    f = x**2*sqrt(-a**2*x**2 + 1)/(-a*x + 1)**4
    F = -asin(a*x)/a**3 + 2*sqrt(-a**2*x**2 + 1)/(a**3*(-a*x + 1)) - 3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*a**3*(-a*x + 1)**3) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*a**3*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_210():
    f = x**2*sqrt(-a**2*x**2 + 1)/(-a*x + 1)**5
    F = 23*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(105*a**3*(-a*x + 1)**3) - 12*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(35*a**3*(-a*x + 1)**4) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(7*a**3*(-a*x + 1)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_211():
    f = x**3/((d + e*x)**4*(d**2 - e**2*x**2)**(sympy.S(7)/2))
    F = d**2/(13*e**4*(d + e*x)**4*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 30*d/(143*e**4*(d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + 21/(143*e**4*(d + e*x)**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + 4/(1001*d*e**4*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 24*x/(5005*d**3*e**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 32*x/(5005*d**5*e**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - 64*x/(5005*d**7*e**3*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_212():
    f = x**2/((d + e*x)**4*(d**2 - e**2*x**2)**(sympy.S(7)/2))
    F = -d/(13*e**3*(d + e*x)**4*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + 17/(143*e**3*(d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 7/(1287*d*e**3*(d + e*x)**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 7/(1287*d**2*e**3*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + 14*x/(2145*d**4*e**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + 56*x/(6435*d**6*e**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + 112*x/(6435*d**8*e**2*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_213():
    f = x/((d + e*x)**4*(d**2 - e**2*x**2)**(sympy.S(7)/2))
    F = 1/(13*e**2*(d + e*x)**4*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 4/(143*d*e**2*(d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 32/(1287*d**2*e**2*(d + e*x)**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 32/(1287*d**3*e**2*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + 64*x/(2145*d**5*e*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + 256*x/(6435*d**7*e*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + 512*x/(6435*d**9*e*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_214():
    f = 1/((d + e*x)**4*(d**2 - e**2*x**2)**(sympy.S(7)/2))
    F = -1/(13*d*e*(d + e*x)**4*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 9/(143*d**2*e*(d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 8/(143*d**3*e*(d + e*x)**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - 8/(143*d**4*e*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + 48*x/(715*d**6*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + 64*x/(715*d**8*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + 128*x/(715*d**10*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_215():
    f = 1/(x*(d + e*x)**4*(d**2 - e**2*x**2)**(sympy.S(7)/2))
    F = 8*d*(d - e*x)/(13*(d**2 - e**2*x**2)**(sympy.S(13)/2)) - 4*e*x/(13*d*(d**2 - e**2*x**2)**(sympy.S(11)/2)) + (13*d - 40*e*x)/(117*d**3*(d**2 - e**2*x**2)**(sympy.S(9)/2)) + (117*d - 320*e*x)/(819*d**5*(d**2 - e**2*x**2)**(sympy.S(7)/2)) + (273*d - 640*e*x)/(1365*d**7*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + (273*d - 512*e*x)/(819*d**9*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (819*d - 1024*e*x)/(819*d**11*sqrt(d**2 - e**2*x**2)) - atanh(sqrt(d**2 - e**2*x**2)/d)/d**11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_216():
    f = 1/(x**2*(d + e*x)**4*(d**2 - e**2*x**2)**(sympy.S(7)/2))
    F = -8*e*(d - e*x)/(13*(d**2 - e**2*x**2)**(sympy.S(13)/2)) - 4*e*(13*d - 24*e*x)/(143*d**2*(d**2 - e**2*x**2)**(sympy.S(11)/2)) - e*(572*d - 1103*e*x)/(1287*d**4*(d**2 - e**2*x**2)**(sympy.S(9)/2)) - e*(5148*d - 10111*e*x)/(9009*d**6*(d**2 - e**2*x**2)**(sympy.S(7)/2)) - e*(12012*d - 23225*e*x)/(15015*d**8*(d**2 - e**2*x**2)**(sympy.S(5)/2)) - e*(12012*d - 21583*e*x)/(9009*d**10*(d**2 - e**2*x**2)**(sympy.S(3)/2)) - e*(36036*d - 52175*e*x)/(9009*d**12*sqrt(d**2 - e**2*x**2)) + 4*e*atanh(sqrt(d**2 - e**2*x**2)/d)/d**12 - sqrt(d**2 - e**2*x**2)/(d**12*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_217():
    f = sqrt(-a**2*x**2 + 1)*sqrt(-a*c*x + c)/x**2
    F = a*sqrt(c)*atanh(sqrt(c)*sqrt(-a**2*x**2 + 1)/sqrt(-a*c*x + c)) - a*c*sqrt(-a**2*x**2 + 1)/sqrt(-a*c*x + c) - c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(x*(-a*c*x + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_218():
    f = sqrt(-a*c*x + c)/(x*sqrt(-a**2*x**2 + 1))
    F = -2*sqrt(c)*atanh(sqrt(c)*sqrt(-a**2*x**2 + 1)/sqrt(-a*c*x + c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_219():
    f = sqrt(-a*x + 1)/sqrt(x)
    F = sqrt(x)*sqrt(-a*x + 1) + asin(sqrt(a)*sqrt(x))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_220():
    f = sqrt(-a**2*x**2 + 1)/(sqrt(x)*sqrt(a*x + 1))
    F = sqrt(x)*sqrt(-a*x + 1) + asin(sqrt(a)*sqrt(x))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_221():
    f = sqrt(a*x + 1)/sqrt(x)
    F = sqrt(x)*sqrt(a*x + 1) + asinh(sqrt(a)*sqrt(x))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_222():
    f = sqrt(-a**2*x**2 + 1)/(sqrt(x)*sqrt(-a*x + 1))
    F = sqrt(x)*sqrt(a*x + 1) + asinh(sqrt(a)*sqrt(x))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_223():
    f = sqrt(x)*sqrt(-a*x + 1)
    F = x**(sympy.S(3)/2)*sqrt(-a*x + 1)/2 - sqrt(x)*sqrt(-a*x + 1)/(4*a) + asin(sqrt(a)*sqrt(x))/(4*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_224():
    f = sqrt(x)*sqrt(-a**2*x**2 + 1)/sqrt(a*x + 1)
    F = x**(sympy.S(3)/2)*sqrt(-a*x + 1)/2 - sqrt(x)*sqrt(-a*x + 1)/(4*a) + asin(sqrt(a)*sqrt(x))/(4*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_225():
    f = (g*x)**m*(d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)
    F = d**7*(g*x)**(m + 1)*sqrt(d**2 - e**2*x**2)*(4*m + 11)*hyper((sympy.S(-5)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), e**2*x**2/d**2)/(g*sqrt(1 - e**2*x**2/d**2)*(m + 1)*(m + 8)) + d**6*e*(g*x)**(m + 2)*sqrt(d**2 - e**2*x**2)*(4*m + 29)*hyper((sympy.S(-5)/2, m/2 + 1), (m/2 + 2,), e**2*x**2/d**2)/(g**2*sqrt(1 - e**2*x**2/d**2)*(m + 2)*(m + 9)) - 3*d*(g*x)**(m + 1)*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(g*(m + 8)) - e*(g*x)**(m + 2)*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(g**2*(m + 9))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_226():
    f = (g*x)**m*(d + e*x)**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)
    F = d**6*(g*x)**(m + 1)*sqrt(d**2 - e**2*x**2)*(2*m + 9)*hyper((sympy.S(-5)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), e**2*x**2/d**2)/(g*sqrt(1 - e**2*x**2/d**2)*(m + 1)*(m + 8)) + 2*d**5*e*(g*x)**(m + 2)*sqrt(d**2 - e**2*x**2)*hyper((sympy.S(-5)/2, m/2 + 1), (m/2 + 2,), e**2*x**2/d**2)/(g**2*sqrt(1 - e**2*x**2/d**2)*(m + 2)) - (g*x)**(m + 1)*(d**2 - e**2*x**2)**(sympy.S(7)/2)/(g*(m + 8))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_227():
    f = (g*x)**m*(d + e*x)*(d**2 - e**2*x**2)**(sympy.S(5)/2)
    F = d**5*(g*x)**(m + 1)*sqrt(d**2 - e**2*x**2)*hyper((sympy.S(-5)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), e**2*x**2/d**2)/(g*sqrt(1 - e**2*x**2/d**2)*(m + 1)) + d**4*e*(g*x)**(m + 2)*sqrt(d**2 - e**2*x**2)*hyper((sympy.S(-5)/2, m/2 + 1), (m/2 + 2,), e**2*x**2/d**2)/(g**2*sqrt(1 - e**2*x**2/d**2)*(m + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_228():
    f = (g*x)**m*(d**2 - e**2*x**2)**(sympy.S(5)/2)
    F = d**4*(g*x)**(m + 1)*sqrt(d**2 - e**2*x**2)*hyper((sympy.S(-5)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), e**2*x**2/d**2)/(g*sqrt(1 - e**2*x**2/d**2)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_229():
    f = (g*x)**m*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(d + e*x)
    F = d**3*(g*x)**(m + 1)*sqrt(d**2 - e**2*x**2)*hyper((sympy.S(-3)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), e**2*x**2/d**2)/(g*sqrt(1 - e**2*x**2/d**2)*(m + 1)) - d**2*e*(g*x)**(m + 2)*sqrt(d**2 - e**2*x**2)*hyper((sympy.S(-3)/2, m/2 + 1), (m/2 + 2,), e**2*x**2/d**2)/(g**2*sqrt(1 - e**2*x**2/d**2)*(m + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_230():
    f = (g*x)**m*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(d + e*x)**2
    F = d**2*(g*x)**(m + 1)*sqrt(d**2 - e**2*x**2)*(2*m + 5)*hyper((sympy.S(-1)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), e**2*x**2/d**2)/(g*sqrt(1 - e**2*x**2/d**2)*(m + 1)*(m + 4)) - 2*d*e*(g*x)**(m + 2)*sqrt(d**2 - e**2*x**2)*hyper((sympy.S(-1)/2, m/2 + 1), (m/2 + 2,), e**2*x**2/d**2)/(g**2*sqrt(1 - e**2*x**2/d**2)*(m + 2)) - (g*x)**(m + 1)*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(g*(m + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_231():
    f = (g*x)**m*(d**2 - e**2*x**2)**(sympy.S(5)/2)/(d + e*x)**3
    F = d**3*(g*x)**(m + 1)*sqrt(1 - e**2*x**2/d**2)*(4*m + 5)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), e**2*x**2/d**2)/(g*sqrt(d**2 - e**2*x**2)*(m + 1)*(m + 2)) - d**2*e*(g*x)**(m + 2)*sqrt(1 - e**2*x**2/d**2)*(4*m + 11)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), e**2*x**2/d**2)/(g**2*sqrt(d**2 - e**2*x**2)*(m + 2)*(m + 3)) - 3*d*(g*x)**(m + 1)*sqrt(d**2 - e**2*x**2)/(g*(m + 2)) + e*(g*x)**(m + 2)*sqrt(d**2 - e**2*x**2)/(g**2*(m + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_232():
    f = (g*x)**m*(d + e*x)**3/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = 4*(g*x)**(m + 1)*(d + e*x)/(5*g*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + (g*x)**(m + 1)*(1 - 4*m)*sqrt(1 - e**2*x**2/d**2)*hyper((sympy.S(5)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), e**2*x**2/d**2)/(5*d**3*g*sqrt(d**2 - e**2*x**2)*(m + 1)) + e*(g*x)**(m + 2)*sqrt(1 - e**2*x**2/d**2)*(7 - 4*m)*hyper((sympy.S(5)/2, m/2 + 1), (m/2 + 2,), e**2*x**2/d**2)/(5*d**4*g**2*sqrt(d**2 - e**2*x**2)*(m + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_233():
    f = (g*x)**m*(d + e*x)**2/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = 2*(g*x)**(m + 1)*(d + e*x)/(5*d*g*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + (g*x)**(m + 1)*sqrt(1 - e**2*x**2/d**2)*(3 - 2*m)*hyper((sympy.S(5)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), e**2*x**2/d**2)/(5*d**4*g*sqrt(d**2 - e**2*x**2)*(m + 1)) + 2*e*(g*x)**(m + 2)*sqrt(1 - e**2*x**2/d**2)*(3 - m)*hyper((sympy.S(5)/2, m/2 + 1), (m/2 + 2,), e**2*x**2/d**2)/(5*d**5*g**2*sqrt(d**2 - e**2*x**2)*(m + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_234():
    f = (g*x)**m/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = (g*x)**(m + 1)*sqrt(1 - e**2*x**2/d**2)*hyper((sympy.S(7)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), e**2*x**2/d**2)/(d**6*g*sqrt(d**2 - e**2*x**2)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_235():
    f = (g*x)**m/((d + e*x)*(d**2 - e**2*x**2)**(sympy.S(7)/2))
    F = (g*x)**(m + 1)*sqrt(1 - e**2*x**2/d**2)*hyper((sympy.S(9)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), e**2*x**2/d**2)/(d**7*g*sqrt(d**2 - e**2*x**2)*(m + 1)) - e*(g*x)**(m + 2)*sqrt(1 - e**2*x**2/d**2)*hyper((sympy.S(9)/2, m/2 + 1), (m/2 + 2,), e**2*x**2/d**2)/(d**8*g**2*sqrt(d**2 - e**2*x**2)*(m + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_236():
    f = (g*x)**m/((d + e*x)**2*(d**2 - e**2*x**2)**(sympy.S(7)/2))
    F = 2*(g*x)**(m + 1)*(d - e*x)/(9*d*g*(d**2 - e**2*x**2)**(sympy.S(9)/2)) + (g*x)**(m + 1)*sqrt(1 - e**2*x**2/d**2)*(7 - 2*m)*hyper((sympy.S(9)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), e**2*x**2/d**2)/(9*d**8*g*sqrt(d**2 - e**2*x**2)*(m + 1)) - 2*e*(g*x)**(m + 2)*sqrt(1 - e**2*x**2/d**2)*(7 - m)*hyper((sympy.S(9)/2, m/2 + 1), (m/2 + 2,), e**2*x**2/d**2)/(9*d**9*g**2*sqrt(d**2 - e**2*x**2)*(m + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_237():
    f = (g*x)**m/((d + e*x)**3*(d**2 - e**2*x**2)**(sympy.S(7)/2))
    F = 4*(g*x)**(m + 1)*(d - e*x)/(11*g*(d**2 - e**2*x**2)**(sympy.S(11)/2)) + (g*x)**(m + 1)*sqrt(1 - e**2*x**2/d**2)*(7 - 4*m)*hyper((sympy.S(11)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), e**2*x**2/d**2)/(11*d**9*g*sqrt(d**2 - e**2*x**2)*(m + 1)) - e*(g*x)**(m + 2)*sqrt(1 - e**2*x**2/d**2)*(25 - 4*m)*hyper((sympy.S(11)/2, m/2 + 1), (m/2 + 2,), e**2*x**2/d**2)/(11*d**10*g**2*sqrt(d**2 - e**2*x**2)*(m + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_238():
    f = x**5*(d + e*x)*(d**2 - e**2*x**2)**p
    F = -d**5*(d**2 - e**2*x**2)**(p + 1)/(2*e**6*(p + 1)) + d**3*(d**2 - e**2*x**2)**(p + 2)/(e**6*(p + 2)) - d*(d**2 - e**2*x**2)**(p + 3)/(2*e**6*(p + 3)) + e*x**7*(d**2 - e**2*x**2)**p*hyper((sympy.S(7)/2, -p), (sympy.S(9)/2,), e**2*x**2/d**2)/(7*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_239():
    f = x**4*(d + e*x)*(d**2 - e**2*x**2)**p
    F = -d**4*(d**2 - e**2*x**2)**(p + 1)/(2*e**5*(p + 1)) + d**2*(d**2 - e**2*x**2)**(p + 2)/(e**5*(p + 2)) + d*x**5*(d**2 - e**2*x**2)**p*hyper((sympy.S(5)/2, -p), (sympy.S(7)/2,), e**2*x**2/d**2)/(5*(1 - e**2*x**2/d**2)**p) - (d**2 - e**2*x**2)**(p + 3)/(2*e**5*(p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_240():
    f = x**3*(d + e*x)*(d**2 - e**2*x**2)**p
    F = -d**3*(d**2 - e**2*x**2)**(p + 1)/(2*e**4*(p + 1)) + d*(d**2 - e**2*x**2)**(p + 2)/(2*e**4*(p + 2)) + e*x**5*(d**2 - e**2*x**2)**p*hyper((sympy.S(5)/2, -p), (sympy.S(7)/2,), e**2*x**2/d**2)/(5*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_241():
    f = x**2*(d + e*x)*(d**2 - e**2*x**2)**p
    F = -d**2*(d**2 - e**2*x**2)**(p + 1)/(2*e**3*(p + 1)) + d*x**3*(d**2 - e**2*x**2)**p*hyper((sympy.S(3)/2, -p), (sympy.S(5)/2,), e**2*x**2/d**2)/(3*(1 - e**2*x**2/d**2)**p) + (d**2 - e**2*x**2)**(p + 2)/(2*e**3*(p + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_242():
    f = x*(d + e*x)*(d**2 - e**2*x**2)**p
    F = -d*(d**2 - e**2*x**2)**(p + 1)/(2*e**2*(p + 1)) + e*x**3*(d**2 - e**2*x**2)**p*hyper((sympy.S(3)/2, -p), (sympy.S(5)/2,), e**2*x**2/d**2)/(3*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_243():
    f = (d + e*x)*(d**2 - e**2*x**2)**p
    F = d*x*(d**2 - e**2*x**2)**p*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), e**2*x**2/d**2)/(1 - e**2*x**2/d**2)**p - (d**2 - e**2*x**2)**(p + 1)/(2*e*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_244():
    f = (d + e*x)*(d**2 - e**2*x**2)**p/x
    F = e*x*(d**2 - e**2*x**2)**p*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), e**2*x**2/d**2)/(1 - e**2*x**2/d**2)**p - (d**2 - e**2*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 - e**2*x**2/d**2)/(2*d*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_245():
    f = (d + e*x)*(d**2 - e**2*x**2)**p/x**2
    F = -d*(d**2 - e**2*x**2)**p*hyper((sympy.S(-1)/2, -p), (sympy.S.Half,), e**2*x**2/d**2)/(x*(1 - e**2*x**2/d**2)**p) - e*(d**2 - e**2*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 - e**2*x**2/d**2)/(2*d**2*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_246():
    f = (d + e*x)*(d**2 - e**2*x**2)**p/x**3
    F = -e*(d**2 - e**2*x**2)**p*hyper((sympy.S(-1)/2, -p), (sympy.S.Half,), e**2*x**2/d**2)/(x*(1 - e**2*x**2/d**2)**p) - e**2*(d**2 - e**2*x**2)**(p + 1)*hyper((2, p + 1), (p + 2,), 1 - e**2*x**2/d**2)/(2*d**3*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_247():
    f = x**5*(d + e*x)**2*(d**2 - e**2*x**2)**p
    F = -d**6*(d**2 - e**2*x**2)**(p + 1)/(e**6*(p + 1)) + 5*d**4*(d**2 - e**2*x**2)**(p + 2)/(2*e**6*(p + 2)) - 2*d**2*(d**2 - e**2*x**2)**(p + 3)/(e**6*(p + 3)) + 2*d*e*x**7*(d**2 - e**2*x**2)**p*hyper((sympy.S(7)/2, -p), (sympy.S(9)/2,), e**2*x**2/d**2)/(7*(1 - e**2*x**2/d**2)**p) + (d**2 - e**2*x**2)**(p + 4)/(2*e**6*(p + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_248():
    f = x**4*(d + e*x)**2*(d**2 - e**2*x**2)**p
    F = -d**5*(d**2 - e**2*x**2)**(p + 1)/(e**5*(p + 1)) + 2*d**3*(d**2 - e**2*x**2)**(p + 2)/(e**5*(p + 2)) + 2*d**2*x**5*(d**2 - e**2*x**2)**p*(p + 6)*hyper((sympy.S(5)/2, -p), (sympy.S(7)/2,), e**2*x**2/d**2)/((1 - e**2*x**2/d**2)**p*(10*p + 35)) - d*(d**2 - e**2*x**2)**(p + 3)/(e**5*(p + 3)) - x**5*(d**2 - e**2*x**2)**(p + 1)/(2*p + 7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_249():
    f = x**3*(d + e*x)**2*(d**2 - e**2*x**2)**p
    F = -d**4*(d**2 - e**2*x**2)**(p + 1)/(e**4*(p + 1)) + 3*d**2*(d**2 - e**2*x**2)**(p + 2)/(2*e**4*(p + 2)) + 2*d*e*x**5*(d**2 - e**2*x**2)**p*hyper((sympy.S(5)/2, -p), (sympy.S(7)/2,), e**2*x**2/d**2)/(5*(1 - e**2*x**2/d**2)**p) - (d**2 - e**2*x**2)**(p + 3)/(2*e**4*(p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_250():
    f = x**2*(d + e*x)**2*(d**2 - e**2*x**2)**p
    F = -d**3*(d**2 - e**2*x**2)**(p + 1)/(e**3*(p + 1)) + 2*d**2*x**3*(d**2 - e**2*x**2)**p*(p + 4)*hyper((sympy.S(3)/2, -p), (sympy.S(5)/2,), e**2*x**2/d**2)/((1 - e**2*x**2/d**2)**p*(6*p + 15)) + d*(d**2 - e**2*x**2)**(p + 2)/(e**3*(p + 2)) - x**3*(d**2 - e**2*x**2)**(p + 1)/(2*p + 5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_251():
    f = x*(d + e*x)**2*(d**2 - e**2*x**2)**p
    F = -d**2*(d**2 - e**2*x**2)**(p + 1)/(e**2*(p + 1)) + 2*d*e*x**3*(d**2 - e**2*x**2)**p*hyper((sympy.S(3)/2, -p), (sympy.S(5)/2,), e**2*x**2/d**2)/(3*(1 - e**2*x**2/d**2)**p) + (d**2 - e**2*x**2)**(p + 2)/(2*e**2*(p + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_252():
    f = (d + e*x)**2*(d**2 - e**2*x**2)**p
    F = -2**(p + 2)*d*(1 + e*x/d)**(-p - 1)*(d**2 - e**2*x**2)**(p + 1)*hyper((p + 1, -p - 2), (p + 2,), (d - e*x)/(2*d))/(e*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_253():
    f = (d + e*x)**2*(d**2 - e**2*x**2)**p/x
    F = 2*d*e*x*(d**2 - e**2*x**2)**p*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), e**2*x**2/d**2)/(1 - e**2*x**2/d**2)**p - (d**2 - e**2*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 - e**2*x**2/d**2)/(2*p + 2) - (d**2 - e**2*x**2)**(p + 1)/(2*p + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_254():
    f = (d + e*x)**2*(d**2 - e**2*x**2)**p/x**2
    F = -2*e**2*p*x*(d**2 - e**2*x**2)**p*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), e**2*x**2/d**2)/(1 - e**2*x**2/d**2)**p - (d**2 - e**2*x**2)**(p + 1)/x - e*(d**2 - e**2*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 - e**2*x**2/d**2)/(d*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_255():
    f = (d + e*x)**2*(d**2 - e**2*x**2)**p/x**3
    F = -2*d*e*(d**2 - e**2*x**2)**p*hyper((sympy.S(-1)/2, -p), (sympy.S.Half,), e**2*x**2/d**2)/(x*(1 - e**2*x**2/d**2)**p) - (d**2 - e**2*x**2)**(p + 1)/(2*x**2) - e**2*(1 - p)*(d**2 - e**2*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 - e**2*x**2/d**2)/(2*d**2*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_256():
    f = x**5*(d + e*x)**3*(d**2 - e**2*x**2)**p
    F = -2*d**7*(d**2 - e**2*x**2)**(p + 1)/(e**6*(p + 1)) + 11*d**5*(d**2 - e**2*x**2)**(p + 2)/(2*e**6*(p + 2)) - 5*d**3*(d**2 - e**2*x**2)**(p + 3)/(e**6*(p + 3)) + 2*d**2*e*x**7*(d**2 - e**2*x**2)**p*(3*p + 17)*hyper((sympy.S(7)/2, -p), (sympy.S(9)/2,), e**2*x**2/d**2)/((1 - e**2*x**2/d**2)**p*(14*p + 63)) + 3*d*(d**2 - e**2*x**2)**(p + 4)/(2*e**6*(p + 4)) - e*x**7*(d**2 - e**2*x**2)**(p + 1)/(2*p + 9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_257():
    f = x**4*(d + e*x)**3*(d**2 - e**2*x**2)**p
    F = -2*d**6*(d**2 - e**2*x**2)**(p + 1)/(e**5*(p + 1)) + 9*d**4*(d**2 - e**2*x**2)**(p + 2)/(2*e**5*(p + 2)) + 2*d**3*x**5*(d**2 - e**2*x**2)**p*(p + 11)*hyper((sympy.S(5)/2, -p), (sympy.S(7)/2,), e**2*x**2/d**2)/((1 - e**2*x**2/d**2)**p*(10*p + 35)) - 3*d**2*(d**2 - e**2*x**2)**(p + 3)/(e**5*(p + 3)) - 3*d*x**5*(d**2 - e**2*x**2)**(p + 1)/(2*p + 7) + (d**2 - e**2*x**2)**(p + 4)/(2*e**5*(p + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_258():
    f = x**3*(d + e*x)**3*(d**2 - e**2*x**2)**p
    F = -2*d**5*(d**2 - e**2*x**2)**(p + 1)/(e**4*(p + 1)) + 7*d**3*(d**2 - e**2*x**2)**(p + 2)/(2*e**4*(p + 2)) + 2*d**2*e*x**5*(d**2 - e**2*x**2)**p*(3*p + 13)*hyper((sympy.S(5)/2, -p), (sympy.S(7)/2,), e**2*x**2/d**2)/((1 - e**2*x**2/d**2)**p*(10*p + 35)) - 3*d*(d**2 - e**2*x**2)**(p + 3)/(2*e**4*(p + 3)) - e*x**5*(d**2 - e**2*x**2)**(p + 1)/(2*p + 7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_259():
    f = x**2*(d + e*x)**3*(d**2 - e**2*x**2)**p
    F = -2*d**4*(d**2 - e**2*x**2)**(p + 1)/(e**3*(p + 1)) + 2*d**3*x**3*(d**2 - e**2*x**2)**p*(p + 7)*hyper((sympy.S(3)/2, -p), (sympy.S(5)/2,), e**2*x**2/d**2)/((1 - e**2*x**2/d**2)**p*(6*p + 15)) + 5*d**2*(d**2 - e**2*x**2)**(p + 2)/(2*e**3*(p + 2)) - 3*d*x**3*(d**2 - e**2*x**2)**(p + 1)/(2*p + 5) - (d**2 - e**2*x**2)**(p + 3)/(2*e**3*(p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_260():
    f = x*(d + e*x)**3*(d**2 - e**2*x**2)**p
    F = -3*2**(p + 3)*d**3*(1 + e*x/d)**(-p - 1)*(d**2 - e**2*x**2)**(p + 1)*hyper((p + 1, -p - 3), (p + 2,), (d - e*x)/(2*d))/(e**2*(p + 1)*(2*p + 5)) - (d + e*x)**3*(d**2 - e**2*x**2)**(p + 1)/(e**2*(2*p + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_261():
    f = (d + e*x)**3*(d**2 - e**2*x**2)**p
    F = -2**(p + 3)*d**2*(1 + e*x/d)**(-p - 1)*(d**2 - e**2*x**2)**(p + 1)*hyper((p + 1, -p - 3), (p + 2,), (d - e*x)/(2*d))/(e*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_262():
    f = (d + e*x)**3*(d**2 - e**2*x**2)**p/x
    F = 2*d**2*e*x*(d**2 - e**2*x**2)**p*(3*p + 5)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), e**2*x**2/d**2)/((1 - e**2*x**2/d**2)**p*(2*p + 3)) - d*(d**2 - e**2*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 - e**2*x**2/d**2)/(2*p + 2) - 3*d*(d**2 - e**2*x**2)**(p + 1)/(2*p + 2) - e*x*(d**2 - e**2*x**2)**(p + 1)/(2*p + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_263():
    f = (d + e*x)**3*(d**2 - e**2*x**2)**p/x**2
    F = 2*d*e**2*x*(1 - p)*(d**2 - e**2*x**2)**p*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), e**2*x**2/d**2)/(1 - e**2*x**2/d**2)**p - d*(d**2 - e**2*x**2)**(p + 1)/x - 3*e*(d**2 - e**2*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 - e**2*x**2/d**2)/(2*p + 2) - e*(d**2 - e**2*x**2)**(p + 1)/(2*p + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_264():
    f = (d + e*x)**3*(d**2 - e**2*x**2)**p/x**3
    F = -d*(d**2 - e**2*x**2)**(p + 1)/(2*x**2) - 2*e**3*x*(d**2 - e**2*x**2)**p*(3*p + 1)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), e**2*x**2/d**2)/(1 - e**2*x**2/d**2)**p - 3*e*(d**2 - e**2*x**2)**(p + 1)/x - e**2*(3 - p)*(d**2 - e**2*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 - e**2*x**2/d**2)/(2*d*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_265():
    f = x**4*(d**2 - e**2*x**2)**p/(d + e*x)
    F = d**4*(d**2 - e**2*x**2)**p/(2*e**5*p) - d**2*(d**2 - e**2*x**2)**(p + 1)/(e**5*(p + 1)) + (d**2 - e**2*x**2)**(p + 2)/(2*e**5*(p + 2)) + x**5*(d**2 - e**2*x**2)**p*hyper((sympy.S(5)/2, 1 - p), (sympy.S(7)/2,), e**2*x**2/d**2)/(5*d*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_266():
    f = x**3*(d**2 - e**2*x**2)**p/(d + e*x)
    F = -d**3*(d**2 - e**2*x**2)**p/(2*e**4*p) + d*(d**2 - e**2*x**2)**(p + 1)/(2*e**4*(p + 1)) - e*x**5*(d**2 - e**2*x**2)**p*hyper((sympy.S(5)/2, 1 - p), (sympy.S(7)/2,), e**2*x**2/d**2)/(5*d**2*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_267():
    f = x**2*(d**2 - e**2*x**2)**p/(d + e*x)
    F = d**2*(d**2 - e**2*x**2)**p/(2*e**3*p) - (d**2 - e**2*x**2)**(p + 1)/(2*e**3*(p + 1)) + x**3*(d**2 - e**2*x**2)**p*hyper((sympy.S(3)/2, 1 - p), (sympy.S(5)/2,), e**2*x**2/d**2)/(3*d*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_268():
    f = x*(d**2 - e**2*x**2)**p/(d + e*x)
    F = -d*(d**2 - e**2*x**2)**p/(2*e**2*p) - e*x**3*(d**2 - e**2*x**2)**p*hyper((sympy.S(3)/2, 1 - p), (sympy.S(5)/2,), e**2*x**2/d**2)/(3*d**2*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_269():
    f = (d**2 - e**2*x**2)**p/(d + e*x)
    F = -2**(p - 1)*(1 + e*x/d)**(-p - 1)*(d**2 - e**2*x**2)**(p + 1)*hyper((p + 1, 1 - p), (p + 2,), (d - e*x)/(2*d))/(d**2*e*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_270():
    f = (d**2 - e**2*x**2)**p/(x*(d + e*x))
    F = -(d**2 - e**2*x**2)**p*hyper((1, p), (p + 1,), 1 - e**2*x**2/d**2)/(2*d*p) - e*x*(d**2 - e**2*x**2)**p*hyper((sympy.S.Half, 1 - p), (sympy.S(3)/2,), e**2*x**2/d**2)/(d**2*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_271():
    f = (d**2 - e**2*x**2)**p/(x**2*(d + e*x))
    F = -(d**2 - e**2*x**2)**p*hyper((sympy.S(-1)/2, 1 - p), (sympy.S.Half,), e**2*x**2/d**2)/(d*x*(1 - e**2*x**2/d**2)**p) + e*(d**2 - e**2*x**2)**p*hyper((1, p), (p + 1,), 1 - e**2*x**2/d**2)/(2*d**2*p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_272():
    f = (d**2 - e**2*x**2)**p/(x**3*(d + e*x))
    F = e*(d**2 - e**2*x**2)**p*hyper((sympy.S(-1)/2, 1 - p), (sympy.S.Half,), e**2*x**2/d**2)/(d**2*x*(1 - e**2*x**2/d**2)**p) - e**2*(d**2 - e**2*x**2)**p*hyper((2, p), (p + 1,), 1 - e**2*x**2/d**2)/(2*d**3*p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_273():
    f = x**5*(d**2 - e**2*x**2)**p/(d + e*x)**2
    F = d**6*(d**2 - e**2*x**2)**(p - 1)/(e**6*(1 - p)) + 5*d**4*(d**2 - e**2*x**2)**p/(2*e**6*p) - 2*d**2*(d**2 - e**2*x**2)**(p + 1)/(e**6*(p + 1)) + (d**2 - e**2*x**2)**(p + 2)/(2*e**6*(p + 2)) - 2*e*x**7*(d**2 - e**2*x**2)**p*hyper((sympy.S(7)/2, 2 - p), (sympy.S(9)/2,), e**2*x**2/d**2)/(7*d**3*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_274():
    f = x**4*(d**2 - e**2*x**2)**p/(d + e*x)**2
    F = -d**5*(d**2 - e**2*x**2)**(p - 1)/(e**5*(1 - p)) - 2*d**3*(d**2 - e**2*x**2)**p/(e**5*p) + d*(d**2 - e**2*x**2)**(p + 1)/(e**5*(p + 1)) - x**5*(d**2 - e**2*x**2)**(p - 1)/(2*p + 3) + x**5*(d**2 - e**2*x**2)**p*(2*p + 8)*hyper((sympy.S(5)/2, 2 - p), (sympy.S(7)/2,), e**2*x**2/d**2)/(5*d**2*(1 - e**2*x**2/d**2)**p*(2*p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_275():
    f = x**3*(d**2 - e**2*x**2)**p/(d + e*x)**2
    F = d**4*(d**2 - e**2*x**2)**(p - 1)/(e**4*(1 - p)) + 3*d**2*(d**2 - e**2*x**2)**p/(2*e**4*p) - (d**2 - e**2*x**2)**(p + 1)/(2*e**4*(p + 1)) - 2*e*x**5*(d**2 - e**2*x**2)**p*hyper((sympy.S(5)/2, 2 - p), (sympy.S(7)/2,), e**2*x**2/d**2)/(5*d**3*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_276():
    f = x**2*(d**2 - e**2*x**2)**p/(d + e*x)**2
    F = -d**3*(d**2 - e**2*x**2)**(p - 1)/(e**3*(1 - p)) - d*(d**2 - e**2*x**2)**p/(e**3*p) - x**3*(d**2 - e**2*x**2)**(p - 1)/(2*p + 1) + x**3*(d**2 - e**2*x**2)**p*(2*p + 4)*hyper((sympy.S(3)/2, 2 - p), (sympy.S(5)/2,), e**2*x**2/d**2)/(3*d**2*(1 - e**2*x**2/d**2)**p*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_277():
    f = x*(d**2 - e**2*x**2)**p/(d + e*x)**2
    F = -2**(p - 1)*(1 + e*x/d)**(-p - 1)*(d**2 - e**2*x**2)**(p + 1)*hyper((p + 1, 1 - p), (p + 2,), (d - e*x)/(2*d))/(d**2*e**2*(1 - p**2)) + (d**2 - e**2*x**2)**(p + 1)/(2*e**2*(1 - p)*(d + e*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_278():
    f = (d**2 - e**2*x**2)**p/(d + e*x)**2
    F = -2**(p - 2)*(1 + e*x/d)**(-p - 1)*(d**2 - e**2*x**2)**(p + 1)*hyper((p + 1, 2 - p), (p + 2,), (d - e*x)/(2*d))/(d**3*e*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_279():
    f = (d**2 - e**2*x**2)**p/(x*(d + e*x)**2)
    F = (d**2 - e**2*x**2)**(p - 1)/(1 - p) - (d**2 - e**2*x**2)**p*hyper((1, p), (p + 1,), 1 - e**2*x**2/d**2)/(2*d**2*p) - 2*e*x*(d**2 - e**2*x**2)**p*hyper((sympy.S.Half, 2 - p), (sympy.S(3)/2,), e**2*x**2/d**2)/(d**3*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_280():
    f = (d**2 - e**2*x**2)**p/(x**2*(d + e*x)**2)
    F = -(d**2 - e**2*x**2)**(p - 1)/x - e*(d**2 - e**2*x**2)**(p - 1)*hyper((1, p - 1), (p,), 1 - e**2*x**2/d**2)/(d*(1 - p)) + 2*e**2*x*(2 - p)*(d**2 - e**2*x**2)**p*hyper((sympy.S.Half, 2 - p), (sympy.S(3)/2,), e**2*x**2/d**2)/(d**4*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_281():
    f = (d**2 - e**2*x**2)**p/(x**3*(d + e*x)**2)
    F = -(d**2 - e**2*x**2)**(p - 1)/(2*x**2) + e**2*(3 - p)*(d**2 - e**2*x**2)**(p - 1)*hyper((1, p - 1), (p,), 1 - e**2*x**2/d**2)/(2*d**2*(1 - p)) + 2*e*(d**2 - e**2*x**2)**p*hyper((sympy.S(-1)/2, 2 - p), (sympy.S.Half,), e**2*x**2/d**2)/(d**3*x*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_282():
    f = (d**2 - e**2*x**2)**p/(x**4*(d + e*x)**2)
    F = -(d**2 - e**2*x**2)**(p - 1)/(3*x**3) - e**3*(d**2 - e**2*x**2)**(p - 1)*hyper((2, p - 1), (p,), 1 - e**2*x**2/d**2)/(d**3*(1 - p)) - 2*e**2*(4 - p)*(d**2 - e**2*x**2)**p*hyper((sympy.S(-1)/2, 2 - p), (sympy.S.Half,), e**2*x**2/d**2)/(3*d**4*x*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_283():
    f = (d**2 - e**2*x**2)**p/(x**5*(d + e*x)**2)
    F = -(d**2 - e**2*x**2)**(p - 1)/(4*x**4) + 2*e*(d**2 - e**2*x**2)**p*hyper((sympy.S(-3)/2, 2 - p), (sympy.S(-1)/2,), e**2*x**2/d**2)/(3*d**3*x**3*(1 - e**2*x**2/d**2)**p) + e**4*(5 - p)*(d**2 - e**2*x**2)**(p - 1)*hyper((2, p - 1), (p,), 1 - e**2*x**2/d**2)/(4*d**4*(1 - p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_284():
    f = x**4*(d**2 - e**2*x**2)**p/(d + e*x)**3
    F = -2*d**6*(d**2 - e**2*x**2)**(p - 2)/(e**5*(2 - p)) + 9*d**4*(d**2 - e**2*x**2)**(p - 1)/(2*e**5*(1 - p)) + 3*d**2*(d**2 - e**2*x**2)**p/(e**5*p) - 3*d*x**5*(d**2 - e**2*x**2)**(p - 2)/(2*p + 1) - (d**2 - e**2*x**2)**(p + 1)/(2*e**5*(p + 1)) + x**5*(d**2 - e**2*x**2)**p*(2*p + 16)*hyper((sympy.S(5)/2, 3 - p), (sympy.S(7)/2,), e**2*x**2/d**2)/(5*d**3*(1 - e**2*x**2/d**2)**p*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_285():
    f = x**3*(d**2 - e**2*x**2)**p/(d + e*x)**3
    F = 2*d**5*(d**2 - e**2*x**2)**(p - 2)/(e**4*(2 - p)) - 7*d**3*(d**2 - e**2*x**2)**(p - 1)/(2*e**4*(1 - p)) - 3*d*(d**2 - e**2*x**2)**p/(2*e**4*p) + e*x**5*(d**2 - e**2*x**2)**(p - 2)/(2*p + 1) - 2*e*x**5*(d**2 - e**2*x**2)**p*(3*p + 4)*hyper((sympy.S(5)/2, 3 - p), (sympy.S(7)/2,), e**2*x**2/d**2)/(5*d**4*(1 - e**2*x**2/d**2)**p*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_286():
    f = x**2*(d**2 - e**2*x**2)**p/(d + e*x)**3
    F = 2**(p - 3)*(1 + e*x/d)**(-p - 1)*(d**2 - e**2*x**2)**(p + 1)*(p + 4)*hyper((p + 1, 2 - p), (p + 2,), (d - e*x)/(2*d))/(d**2*e**3*p*(2 - p)*(p + 1)) - d*(d**2 - e**2*x**2)**(p + 1)/(2*e**3*(2 - p)*(d + e*x)**3) - (d**2 - e**2*x**2)**(p + 1)/(2*e**3*p*(d + e*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_287():
    f = x*(d**2 - e**2*x**2)**p/(d + e*x)**3
    F = -3*2**(p - 3)*(1 + e*x/d)**(-p - 1)*(d**2 - e**2*x**2)**(p + 1)*hyper((p + 1, 2 - p), (p + 2,), (d - e*x)/(2*d))/(d**3*e**2*(2 - p)*(p + 1)) + (d**2 - e**2*x**2)**(p + 1)/(2*e**2*(2 - p)*(d + e*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_288():
    f = (d**2 - e**2*x**2)**p/(d + e*x)**3
    F = -2**(p - 3)*(1 + e*x/d)**(-p - 1)*(d**2 - e**2*x**2)**(p + 1)*hyper((p + 1, 3 - p), (p + 2,), (d - e*x)/(2*d))/(d**4*e*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_289():
    f = (d**2 - e**2*x**2)**p/(x*(d + e*x)**3)
    F = 2*d*(d**2 - e**2*x**2)**(p - 2)/(2 - p) - e*x*(d**2 - e**2*x**2)**(p - 2)/(3 - 2*p) + (d**2 - e**2*x**2)**(p - 1)*hyper((1, p - 1), (p,), 1 - e**2*x**2/d**2)/(2*d*(1 - p)) - 2*e*x*(4 - 3*p)*(d**2 - e**2*x**2)**p*hyper((sympy.S.Half, 3 - p), (sympy.S(3)/2,), e**2*x**2/d**2)/(d**4*(1 - e**2*x**2/d**2)**p*(3 - 2*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_290():
    f = (d**2 - e**2*x**2)**p/(x**2*(d + e*x)**3)
    F = -d*(d**2 - e**2*x**2)**(p - 2)/x - 2*e*(d**2 - e**2*x**2)**(p - 2)/(2 - p) - 3*e*(d**2 - e**2*x**2)**(p - 1)*hyper((1, p - 1), (p,), 1 - e**2*x**2/d**2)/(2*d**2*(1 - p)) + 2*e**2*x*(4 - p)*(d**2 - e**2*x**2)**p*hyper((sympy.S.Half, 3 - p), (sympy.S(3)/2,), e**2*x**2/d**2)/(d**5*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_291():
    f = (d**2 - e**2*x**2)**p/(x**3*(d + e*x)**3)
    F = -d*(d**2 - e**2*x**2)**(p - 2)/(2*x**2) + 3*e*(d**2 - e**2*x**2)**(p - 2)/x + e**2*(6 - p)*(d**2 - e**2*x**2)**(p - 2)*hyper((1, p - 2), (p - 1,), 1 - e**2*x**2/d**2)/(2*d*(2 - p)) - 2*e**3*x*(8 - 3*p)*(d**2 - e**2*x**2)**p*hyper((sympy.S.Half, 3 - p), (sympy.S(3)/2,), e**2*x**2/d**2)/(d**6*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_292():
    f = (d**2 - e**2*x**2)**p/(x**4*(d + e*x)**3)
    F = -d*(d**2 - e**2*x**2)**(p - 2)/(3*x**3) + 3*e*(d**2 - e**2*x**2)**(p - 2)/(2*x**2) - e**3*(10 - 3*p)*(d**2 - e**2*x**2)**(p - 2)*hyper((1, p - 2), (p - 1,), 1 - e**2*x**2/d**2)/(2*d**2*(2 - p)) - 2*e**2*(8 - p)*(d**2 - e**2*x**2)**p*hyper((sympy.S(-1)/2, 3 - p), (sympy.S.Half,), e**2*x**2/d**2)/(3*d**5*x*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_293():
    f = (d**2 - e**2*x**2)**p/(x**5*(d + e*x)**3)
    F = -d*(d**2 - e**2*x**2)**(p - 2)/(4*x**4) + e*(d**2 - e**2*x**2)**(p - 2)/x**3 + e**4*(10 - p)*(d**2 - e**2*x**2)**(p - 2)*hyper((2, p - 2), (p - 1,), 1 - e**2*x**2/d**2)/(4*d**3*(2 - p)) + 2*e**3*(4 - p)*(d**2 - e**2*x**2)**p*hyper((sympy.S(-1)/2, 3 - p), (sympy.S.Half,), e**2*x**2/d**2)/(d**6*x*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_294():
    f = x**4*(d**2 - e**2*x**2)**p/(d + e*x)**4
    F = -4*d**7*(d**2 - e**2*x**2)**(p - 3)/(e**5*(3 - p)) + 10*d**5*(d**2 - e**2*x**2)**(p - 2)/(e**5*(2 - p)) - 8*d**3*(d**2 - e**2*x**2)**(p - 1)/(e**5*(1 - p)) + d**2*x**5*(d**2 - e**2*x**2)**(p - 3)*(12*p + 13)/(1 - 4*p**2) - 2*d*(d**2 - e**2*x**2)**p/(e**5*p) - e**2*x**7*(d**2 - e**2*x**2)**(p - 3)/(2*p + 1) - x**5*(d**2 - e**2*x**2)**p*(4*p**2 + 60*p + 64)*hyper((sympy.S(5)/2, 4 - p), (sympy.S(7)/2,), e**2*x**2/d**2)/(5*d**4*(1 - 4*p**2)*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_295():
    f = x**3*(d**2 - e**2*x**2)**p/(d + e*x)**4
    F = 3*2**(p - 2)*(1 + e*x/d)**(-p - 1)*(d**2 - e**2*x**2)**(p + 1)*(p + 2)*hyper((p + 1, 3 - p), (p + 2,), (d - e*x)/(2*d))/(d**2*e**4*p*(1 - 2*p)*(3 - p)*(p + 1)) + d**2*(d**2 - e**2*x**2)**(p + 1)/(2*e**4*(3 - p)*(d + e*x)**4) - d*(d**2 - e**2*x**2)**(p + 1)*(2*p + 1)/(e**4*p*(1 - 2*p)*(d + e*x)**3) - (d**2 - e**2*x**2)**(p + 1)/(2*e**4*p*(d + e*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_296():
    f = x**2*(d**2 - e**2*x**2)**p/(d + e*x)**4
    F = -2**(p - 3)*(1 + e*x/d)**(-p - 1)*(d**2 - e**2*x**2)**(p + 1)*(p + 7)*hyper((p + 1, 3 - p), (p + 2,), (d - e*x)/(2*d))/(d**3*e**3*(1 - 2*p)*(3 - p)*(p + 1)) - d*(d**2 - e**2*x**2)**(p + 1)/(2*e**3*(3 - p)*(d + e*x)**4) + (d**2 - e**2*x**2)**(p + 1)/(e**3*(1 - 2*p)*(d + e*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_297():
    f = x*(d**2 - e**2*x**2)**p/(d + e*x)**4
    F = -2**(p - 2)*(1 + e*x/d)**(-p - 1)*(d**2 - e**2*x**2)**(p + 1)*hyper((p + 1, 3 - p), (p + 2,), (d - e*x)/(2*d))/(d**4*e**2*(3 - p)*(p + 1)) + (d**2 - e**2*x**2)**(p + 1)/(2*e**2*(3 - p)*(d + e*x)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_298():
    f = (d**2 - e**2*x**2)**p/(d + e*x)**4
    F = -2**(p - 4)*(1 + e*x/d)**(-p - 1)*(d**2 - e**2*x**2)**(p + 1)*hyper((p + 1, 4 - p), (p + 2,), (d - e*x)/(2*d))/(d**5*e*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_299():
    f = (d**2 - e**2*x**2)**p/(x*(d + e*x)**4)
    F = 4*d**2*(d**2 - e**2*x**2)**(p - 3)/(3 - p) - 4*d*e*x*(d**2 - e**2*x**2)**(p - 3)/(5 - 2*p) + (d**2 - e**2*x**2)**(p - 2)*hyper((1, p - 2), (p - 1,), 1 - e**2*x**2/d**2)/(4 - 2*p) - (d**2 - e**2*x**2)**(p - 2)/(4 - 2*p) - 8*e*x*(2 - p)*(d**2 - e**2*x**2)**p*hyper((sympy.S.Half, 4 - p), (sympy.S(3)/2,), e**2*x**2/d**2)/(d**5*(1 - e**2*x**2/d**2)**p*(5 - 2*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_300():
    f = (d**2 - e**2*x**2)**p/(x**2*(d + e*x)**4)
    F = -d**2*(d**2 - e**2*x**2)**(p - 3)/x - 4*d*e*(d**2 - e**2*x**2)**(p - 3)/(3 - p) + e**2*x*(d**2 - e**2*x**2)**(p - 3)/(5 - 2*p) - 2*e*(d**2 - e**2*x**2)**(p - 2)*hyper((1, p - 2), (p - 1,), 1 - e**2*x**2/d**2)/(d*(2 - p)) + 4*e**2*x*(d**2 - e**2*x**2)**p*(p**2 - 9*p + 16)*hyper((sympy.S.Half, 4 - p), (sympy.S(3)/2,), e**2*x**2/d**2)/(d**6*(1 - e**2*x**2/d**2)**p*(5 - 2*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_301():
    f = (d**2 - e**2*x**2)**p/(x**3*(d + e*x)**4)
    F = -d**2*(d**2 - e**2*x**2)**(p - 3)/(2*x**2) + 4*d*e*(d**2 - e**2*x**2)**(p - 3)/x + e**2*(11 - p)*(d**2 - e**2*x**2)**(p - 3)/(6 - 2*p) + e**2*(10 - p)*(d**2 - e**2*x**2)**(p - 2)*hyper((1, p - 2), (p - 1,), 1 - e**2*x**2/d**2)/(2*d**2*(2 - p)) - 8*e**3*x*(4 - p)*(d**2 - e**2*x**2)**p*hyper((sympy.S.Half, 4 - p), (sympy.S(3)/2,), e**2*x**2/d**2)/(d**7*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_302():
    f = (d**2 - e**2*x**2)**p/(x**4*(d + e*x)**4)
    F = -d**2*(d**2 - e**2*x**2)**(p - 3)/(3*x**3) + 2*d*e*(d**2 - e**2*x**2)**(p - 3)/x**2 - e**2*(27 - 2*p)*(d**2 - e**2*x**2)**(p - 3)/(3*x) - 2*e**3*(5 - p)*(d**2 - e**2*x**2)**(p - 3)*hyper((1, p - 3), (p - 2,), 1 - e**2*x**2/d**2)/(d*(3 - p)) + 4*e**4*x*(d**2 - e**2*x**2)**p*(p**2 - 17*p + 48)*hyper((sympy.S.Half, 4 - p), (sympy.S(3)/2,), e**2*x**2/d**2)/(3*d**8*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_303():
    f = (d**2 - e**2*x**2)**p/(x**5*(d + e*x)**4)
    F = -d**2*(d**2 - e**2*x**2)**(p - 3)/(4*x**4) + 4*d*e*(d**2 - e**2*x**2)**(p - 3)/(3*x**3) - e**2*(17 - p)*(d**2 - e**2*x**2)**(p - 3)/(4*x**2) + e**4*(d**2 - e**2*x**2)**(p - 3)*(p**2 - 21*p + 70)*hyper((1, p - 3), (p - 2,), 1 - e**2*x**2/d**2)/(4*d**2*(3 - p)) + 8*e**3*(6 - p)*(d**2 - e**2*x**2)**p*hyper((sympy.S(-1)/2, 4 - p), (sympy.S.Half,), e**2*x**2/d**2)/(3*d**7*x*(1 - e**2*x**2/d**2)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_304():
    f = (g*x)**m*(d + e*x)**3*(d**2 - e**2*x**2)**p
    F = 2*d**3*(g*x)**(m + 1)*(d**2 - e**2*x**2)**p*(2*m + p + 3)*hyper((-p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), e**2*x**2/d**2)/(g*(1 - e**2*x**2/d**2)**p*(m + 1)*(m + 2*p + 3)) + 2*d**2*e*(g*x)**(m + 2)*(d**2 - e**2*x**2)**p*(2*m + 3*p + 7)*hyper((-p, m/2 + 1), (m/2 + 2,), e**2*x**2/d**2)/(g**2*(1 - e**2*x**2/d**2)**p*(m + 2)*(m + 2*p + 4)) - 3*d*(g*x)**(m + 1)*(d**2 - e**2*x**2)**(p + 1)/(g*(m + 2*p + 3)) - e*(g*x)**(m + 2)*(d**2 - e**2*x**2)**(p + 1)/(g**2*(m + 2*p + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_305():
    f = (g*x)**m*(d + e*x)**2*(d**2 - e**2*x**2)**p
    F = 2*d**2*(g*x)**(m + 1)*(d**2 - e**2*x**2)**p*(m + p + 2)*hyper((-p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), e**2*x**2/d**2)/(g*(1 - e**2*x**2/d**2)**p*(m + 1)*(m + 2*p + 3)) + 2*d*e*(g*x)**(m + 2)*(d**2 - e**2*x**2)**p*hyper((-p, m/2 + 1), (m/2 + 2,), e**2*x**2/d**2)/(g**2*(1 - e**2*x**2/d**2)**p*(m + 2)) - (g*x)**(m + 1)*(d**2 - e**2*x**2)**(p + 1)/(g*(m + 2*p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_306():
    f = (g*x)**m*(d + e*x)*(d**2 - e**2*x**2)**p
    F = d*(g*x)**(m + 1)*(d**2 - e**2*x**2)**p*hyper((-p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), e**2*x**2/d**2)/(g*(1 - e**2*x**2/d**2)**p*(m + 1)) + e*(g*x)**(m + 2)*(d**2 - e**2*x**2)**p*hyper((-p, m/2 + 1), (m/2 + 2,), e**2*x**2/d**2)/(g**2*(1 - e**2*x**2/d**2)**p*(m + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_307():
    f = (g*x)**m*(d**2 - e**2*x**2)**p
    F = (g*x)**(m + 1)*(d**2 - e**2*x**2)**p*hyper((-p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), e**2*x**2/d**2)/(g*(1 - e**2*x**2/d**2)**p*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_308():
    f = (g*x)**m*(d**2 - e**2*x**2)**p/(d + e*x)
    F = (g*x)**(m + 1)*(d**2 - e**2*x**2)**p*hyper((1 - p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), e**2*x**2/d**2)/(d*g*(1 - e**2*x**2/d**2)**p*(m + 1)) - e*(g*x)**(m + 2)*(d**2 - e**2*x**2)**p*hyper((1 - p, m/2 + 1), (m/2 + 2,), e**2*x**2/d**2)/(d**2*g**2*(1 - e**2*x**2/d**2)**p*(m + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_309():
    f = (g*x)**m*(d**2 - e**2*x**2)**p/(d + e*x)**2
    F = (g*x)**(m + 1)*(d**2 - e**2*x**2)**(p - 1)/(g*(-m - 2*p + 1)) - (g*x)**(m + 1)*(d**2 - e**2*x**2)**p*(2*m + 2*p)*hyper((2 - p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), e**2*x**2/d**2)/(d**2*g*(1 - e**2*x**2/d**2)**p*(m + 1)*(-m - 2*p + 1)) - 2*e*(g*x)**(m + 2)*(d**2 - e**2*x**2)**p*hyper((2 - p, m/2 + 1), (m/2 + 2,), e**2*x**2/d**2)/(d**3*g**2*(1 - e**2*x**2/d**2)**p*(m + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_310():
    f = (g*x)**m*(d**2 - e**2*x**2)**p/(d + e*x)**3
    F = 3*d*(g*x)**(m + 1)*(d**2 - e**2*x**2)**(p - 2)/(g*(-m - 2*p + 3)) - e*(g*x)**(m + 2)*(d**2 - e**2*x**2)**(p - 2)/(g**2*(-m - 2*p + 2)) - (g*x)**(m + 1)*(d**2 - e**2*x**2)**p*(4*m + 2*p)*hyper((3 - p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), e**2*x**2/d**2)/(d**3*g*(1 - e**2*x**2/d**2)**p*(m + 1)*(-m - 2*p + 3)) - 2*e*(g*x)**(m + 2)*(d**2 - e**2*x**2)**p*(-2*m - 3*p + 2)*hyper((3 - p, m/2 + 1), (m/2 + 2,), e**2*x**2/d**2)/(d**4*g**2*(1 - e**2*x**2/d**2)**p*(m + 2)*(-m - 2*p + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_311():
    f = (g*x)**m*(-a**2*x**2 + 1)**p/(a*x + 1)
    F = -a*(g*x)**(m + 2)*hyper((1 - p, m/2 + 1), (m/2 + 2,), a**2*x**2)/(g**2*(m + 2)) + (g*x)**(m + 1)*hyper((1 - p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/(g*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_312():
    f = (g*x)**m*(d + e*x)**n*(d**2 - e**2*x**2)**p
    F = (g*x)**(m + 1)*(1 + e*x/d)**(-n - p)*(d + e*x)**n*(d**2 - e**2*x**2)**p*appellf1(m + 1, -p, -n - p, m + 2, e*x/d, -e*x/d)/(g*(1 - e*x/d)**p*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_313():
    f = x*sqrt(x + 1)/(x**2 + 1)
    F = 2*sqrt(x + 1) + sqrt(sympy.S.Half + sqrt(2)/2)*log(x - sqrt(2 + 2*sqrt(2))*sqrt(x + 1) + 1 + sqrt(2))/2 - sqrt(sympy.S.Half + sqrt(2)/2)*log(x + sqrt(2 + 2*sqrt(2))*sqrt(x + 1) + 1 + sqrt(2))/2 + atan((-2*sqrt(x + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/sqrt(2 + 2*sqrt(2)) - atan((2*sqrt(x + 1) + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/sqrt(2 + 2*sqrt(2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_314():
    f = x**4*sqrt(a + c*x**2)/(d + e*x)
    F = -d**4*sqrt(a*e**2 + c*d**2)*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/e**6 - 13*d*(a + c*x**2)**(sympy.S(3)/2)*(d + e*x)/(20*c*e**3) + d*sqrt(a + c*x**2)*(8*c*d**3 - e*x*(-a*e**2 + 4*c*d**2))/(8*c*e**5) + (a + c*x**2)**(sympy.S(3)/2)*(d + e*x)**2/(5*c*e**3) + (a + c*x**2)**(sympy.S(3)/2)*(-8*a*e**2 + 47*c*d**2)/(60*c**2*e**3) - d*(-a**2*e**4 + 4*a*c*d**2*e**2 + 8*c**2*d**4)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(8*c**(sympy.S(3)/2)*e**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_315():
    f = x**3*sqrt(a + c*x**2)/(d + e*x)
    F = d**3*sqrt(a*e**2 + c*d**2)*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/e**5 - 7*d*(a + c*x**2)**(sympy.S(3)/2)/(12*c*e**2) + (a + c*x**2)**(sympy.S(3)/2)*(d + e*x)/(4*c*e**2) - sqrt(a + c*x**2)*(8*c*d**3 - e*x*(-a*e**2 + 4*c*d**2))/(8*c*e**4) + (-a**2*e**4 + 4*a*c*d**2*e**2 + 8*c**2*d**4)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(8*c**(sympy.S(3)/2)*e**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_316():
    f = x**2*sqrt(a + c*x**2)/(d + e*x)
    F = -d**2*sqrt(a*e**2 + c*d**2)*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/e**4 + d*sqrt(a + c*x**2)*(2*d - e*x)/(2*e**3) + (a + c*x**2)**(sympy.S(3)/2)/(3*c*e) - d*(a*e**2 + 2*c*d**2)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(2*sqrt(c)*e**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_317():
    f = x*sqrt(a + c*x**2)/(d + e*x)
    F = d*sqrt(a*e**2 + c*d**2)*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/e**3 - sqrt(a + c*x**2)*(2*d - e*x)/(2*e**2) + (a*e**2 + 2*c*d**2)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(2*sqrt(c)*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_318():
    f = sqrt(a + c*x**2)/(d + e*x)
    F = -sqrt(c)*d*atanh(sqrt(c)*x/sqrt(a + c*x**2))/e**2 + sqrt(a + c*x**2)/e - sqrt(a*e**2 + c*d**2)*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/e**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_319():
    f = sqrt(a + c*x**2)/(x*(d + e*x))
    F = -sqrt(a)*atanh(sqrt(a + c*x**2)/sqrt(a))/d + sqrt(c)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/e + sqrt(a*e**2 + c*d**2)*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(d*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_320():
    f = sqrt(a + c*x**2)/(x**2*(d + e*x))
    F = sqrt(a)*e*atanh(sqrt(a + c*x**2)/sqrt(a))/d**2 - sqrt(a + c*x**2)/(d*x) - sqrt(a*e**2 + c*d**2)*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/d**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_321():
    f = sqrt(a + c*x**2)/(x**3*(d + e*x))
    F = -sqrt(a)*e**2*atanh(sqrt(a + c*x**2)/sqrt(a))/d**3 - sqrt(a + c*x**2)/(2*d*x**2) + e*sqrt(a + c*x**2)/(d**2*x) + e*sqrt(a*e**2 + c*d**2)*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/d**3 - c*atanh(sqrt(a + c*x**2)/sqrt(a))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_322():
    f = sqrt(a + c*x**2)/(x**4*(d + e*x))
    F = sqrt(a)*e**3*atanh(sqrt(a + c*x**2)/sqrt(a))/d**4 + e*sqrt(a + c*x**2)/(2*d**2*x**2) - e**2*sqrt(a + c*x**2)/(d**3*x) - e**2*sqrt(a*e**2 + c*d**2)*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/d**4 - (a + c*x**2)**(sympy.S(3)/2)/(3*a*d*x**3) + c*e*atanh(sqrt(a + c*x**2)/sqrt(a))/(2*sqrt(a)*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_323():
    f = sqrt(a + c*x**2)/(x**5*(d + e*x))
    F = -sqrt(a)*e**4*atanh(sqrt(a + c*x**2)/sqrt(a))/d**5 - sqrt(a + c*x**2)/(4*d*x**4) - e**2*sqrt(a + c*x**2)/(2*d**3*x**2) + e**3*sqrt(a + c*x**2)/(d**4*x) + e**3*sqrt(a*e**2 + c*d**2)*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/d**5 - c*sqrt(a + c*x**2)/(8*a*d*x**2) + e*(a + c*x**2)**(sympy.S(3)/2)/(3*a*d**2*x**3) - c*e**2*atanh(sqrt(a + c*x**2)/sqrt(a))/(2*sqrt(a)*d**3) + c**2*atanh(sqrt(a + c*x**2)/sqrt(a))/(8*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_324():
    f = x**4/(sqrt(a + c*x**2)*(d + e*x))
    F = -d**4*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(e**4*sqrt(a*e**2 + c*d**2)) - 7*d*sqrt(a + c*x**2)*(d + e*x)/(6*c*e**3) + sqrt(a + c*x**2)*(d + e*x)**2/(3*c*e**3) + sqrt(a + c*x**2)*(-4*a*e**2 + 11*c*d**2)/(6*c**2*e**3) - d*(-a*e**2 + 2*c*d**2)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(2*c**(sympy.S(3)/2)*e**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_325():
    f = x**3/(sqrt(a + c*x**2)*(d + e*x))
    F = d**3*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(e**3*sqrt(a*e**2 + c*d**2)) - 3*d*sqrt(a + c*x**2)/(2*c*e**2) + sqrt(a + c*x**2)*(d + e*x)/(2*c*e**2) + (-a*e**2 + 2*c*d**2)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(2*c**(sympy.S(3)/2)*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_326():
    f = x**2/(sqrt(a + c*x**2)*(d + e*x))
    F = -d**2*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(e**2*sqrt(a*e**2 + c*d**2)) + sqrt(a + c*x**2)/(c*e) - d*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(sqrt(c)*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_327():
    f = x/(sqrt(a + c*x**2)*(d + e*x))
    F = d*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(e*sqrt(a*e**2 + c*d**2)) + atanh(sqrt(c)*x/sqrt(a + c*x**2))/(sqrt(c)*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_328():
    f = 1/(sqrt(a + c*x**2)*(d + e*x))
    F = -atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/sqrt(a*e**2 + c*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_329():
    f = 1/(x*sqrt(a + c*x**2)*(d + e*x))
    F = e*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(d*sqrt(a*e**2 + c*d**2)) - atanh(sqrt(a + c*x**2)/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_330():
    f = 1/(x**2*sqrt(a + c*x**2)*(d + e*x))
    F = -e**2*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(d**2*sqrt(a*e**2 + c*d**2)) - sqrt(a + c*x**2)/(a*d*x) + e*atanh(sqrt(a + c*x**2)/sqrt(a))/(sqrt(a)*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_331():
    f = 1/(x**3*sqrt(a + c*x**2)*(d + e*x))
    F = e**3*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(d**3*sqrt(a*e**2 + c*d**2)) - sqrt(a + c*x**2)/(2*a*d*x**2) + e*sqrt(a + c*x**2)/(a*d**2*x) - e**2*atanh(sqrt(a + c*x**2)/sqrt(a))/(sqrt(a)*d**3) + c*atanh(sqrt(a + c*x**2)/sqrt(a))/(2*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_332():
    f = x**4/((a + c*x**2)**(sympy.S(3)/2)*(d + e*x))
    F = a*(a*e + c*d*x)/(c**2*sqrt(a + c*x**2)*(a*e**2 + c*d**2)) - d**4*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(e**2*(a*e**2 + c*d**2)**(sympy.S(3)/2)) + sqrt(a + c*x**2)/(c**2*e) - d*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(c**(sympy.S(3)/2)*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_333():
    f = x**3/((a + c*x**2)**(sympy.S(3)/2)*(d + e*x))
    F = a*(d - e*x)/(c*sqrt(a + c*x**2)*(a*e**2 + c*d**2)) + d**3*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(e*(a*e**2 + c*d**2)**(sympy.S(3)/2)) + atanh(sqrt(c)*x/sqrt(a + c*x**2))/(c**(sympy.S(3)/2)*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_334():
    f = x**2/((a + c*x**2)**(sympy.S(3)/2)*(d + e*x))
    F = -d**2*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(a*e**2 + c*d**2)**(sympy.S(3)/2) - (a*e + c*d*x)/(c*sqrt(a + c*x**2)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_335():
    f = x/((a + c*x**2)**(sympy.S(3)/2)*(d + e*x))
    F = d*e*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(a*e**2 + c*d**2)**(sympy.S(3)/2) - (d - e*x)/(sqrt(a + c*x**2)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_336():
    f = 1/((a + c*x**2)**(sympy.S(3)/2)*(d + e*x))
    F = -e**2*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(a*e**2 + c*d**2)**(sympy.S(3)/2) + (a*e + c*d*x)/(a*sqrt(a + c*x**2)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_337():
    f = 1/(x*(a + c*x**2)**(sympy.S(3)/2)*(d + e*x))
    F = e**3*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(d*(a*e**2 + c*d**2)**(sympy.S(3)/2)) - e*(a*e + c*d*x)/(a*d*sqrt(a + c*x**2)*(a*e**2 + c*d**2)) + 1/(a*d*sqrt(a + c*x**2)) - atanh(sqrt(a + c*x**2)/sqrt(a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_338():
    f = 1/(x**2*(a + c*x**2)**(sympy.S(3)/2)*(d + e*x))
    F = -e**4*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(d**2*(a*e**2 + c*d**2)**(sympy.S(3)/2)) - 1/(a*d*x*sqrt(a + c*x**2)) + e**2*(a*e + c*d*x)/(a*d**2*sqrt(a + c*x**2)*(a*e**2 + c*d**2)) - e/(a*d**2*sqrt(a + c*x**2)) - 2*c*x/(a**2*d*sqrt(a + c*x**2)) + e*atanh(sqrt(a + c*x**2)/sqrt(a))/(a**(sympy.S(3)/2)*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_339():
    f = 1/(x**3*(a + c*x**2)**(sympy.S(3)/2)*(d + e*x))
    F = e**5*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(d**3*(a*e**2 + c*d**2)**(sympy.S(3)/2)) - 1/(2*a*d*x**2*sqrt(a + c*x**2)) + e/(a*d**2*x*sqrt(a + c*x**2)) - e**3*(a*e + c*d*x)/(a*d**3*sqrt(a + c*x**2)*(a*e**2 + c*d**2)) + e**2/(a*d**3*sqrt(a + c*x**2)) - 3*c/(2*a**2*d*sqrt(a + c*x**2)) + 2*c*e*x/(a**2*d**2*sqrt(a + c*x**2)) - e**2*atanh(sqrt(a + c*x**2)/sqrt(a))/(a**(sympy.S(3)/2)*d**3) + 3*c*atanh(sqrt(a + c*x**2)/sqrt(a))/(2*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_340():
    f = x**5/(sqrt(a + c*x**2)*(d + e*x)**2)
    F = d**5*sqrt(a + c*x**2)/(e**4*(d + e*x)*(a*e**2 + c*d**2)) - d**4*(5*a*e**2 + 4*c*d**2)*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(e**5*(a*e**2 + c*d**2)**(sympy.S(3)/2)) - 5*d*sqrt(a + c*x**2)*(d + e*x)/(3*c*e**4) + sqrt(a + c*x**2)*(d + e*x)**2/(3*c*e**4) + sqrt(a + c*x**2)*(-2*a*e**2 + 13*c*d**2)/(3*c**2*e**4) - d*(-a*e**2 + 4*c*d**2)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(c**(sympy.S(3)/2)*e**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_341():
    f = x**4/(sqrt(a + c*x**2)*(d + e*x)**2)
    F = -d**4*sqrt(a + c*x**2)/(e**3*(d + e*x)*(a*e**2 + c*d**2)) + d**3*(4*a*e**2 + 3*c*d**2)*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(e**4*(a*e**2 + c*d**2)**(sympy.S(3)/2)) - 5*d*sqrt(a + c*x**2)/(2*c*e**3) + sqrt(a + c*x**2)*(d + e*x)/(2*c*e**3) + (-a*e**2 + 6*c*d**2)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(2*c**(sympy.S(3)/2)*e**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_342():
    f = x**3/(sqrt(a + c*x**2)*(d + e*x)**2)
    F = d**3*sqrt(a + c*x**2)/(e**2*(d + e*x)*(a*e**2 + c*d**2)) - d**2*(3*a*e**2 + 2*c*d**2)*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(e**3*(a*e**2 + c*d**2)**(sympy.S(3)/2)) + sqrt(a + c*x**2)/(c*e**2) - 2*d*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(sqrt(c)*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_343():
    f = x**2/(sqrt(a + c*x**2)*(d + e*x)**2)
    F = -d**2*sqrt(a + c*x**2)/(e*(d + e*x)*(a*e**2 + c*d**2)) + d*(2*a*e**2 + c*d**2)*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(e**2*(a*e**2 + c*d**2)**(sympy.S(3)/2)) + atanh(sqrt(c)*x/sqrt(a + c*x**2))/(sqrt(c)*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_344():
    f = x/(sqrt(a + c*x**2)*(d + e*x)**2)
    F = -a*e*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(a*e**2 + c*d**2)**(sympy.S(3)/2) + d*sqrt(a + c*x**2)/((d + e*x)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_345():
    f = 1/(sqrt(a + c*x**2)*(d + e*x)**2)
    F = -c*d*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(a*e**2 + c*d**2)**(sympy.S(3)/2) - e*sqrt(a + c*x**2)/((d + e*x)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_346():
    f = 1/(x*sqrt(a + c*x**2)*(d + e*x)**2)
    F = c*e*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(a*e**2 + c*d**2)**(sympy.S(3)/2) + e**2*sqrt(a + c*x**2)/(d*(d + e*x)*(a*e**2 + c*d**2)) + e*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(d**2*sqrt(a*e**2 + c*d**2)) - atanh(sqrt(a + c*x**2)/sqrt(a))/(sqrt(a)*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_347():
    f = 1/(x**2*sqrt(a + c*x**2)*(d + e*x)**2)
    F = -c*e**2*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(d*(a*e**2 + c*d**2)**(sympy.S(3)/2)) - e**3*sqrt(a + c*x**2)/(d**2*(d + e*x)*(a*e**2 + c*d**2)) - 2*e**2*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(d**3*sqrt(a*e**2 + c*d**2)) - sqrt(a + c*x**2)/(a*d**2*x) + 2*e*atanh(sqrt(a + c*x**2)/sqrt(a))/(sqrt(a)*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_348():
    f = 1/(x**3*sqrt(a + c*x**2)*(d + e*x)**2)
    F = c*e**3*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(d**2*(a*e**2 + c*d**2)**(sympy.S(3)/2)) + e**4*sqrt(a + c*x**2)/(d**3*(d + e*x)*(a*e**2 + c*d**2)) + 3*e**3*atanh((a*e - c*d*x)/(sqrt(a + c*x**2)*sqrt(a*e**2 + c*d**2)))/(d**4*sqrt(a*e**2 + c*d**2)) - sqrt(a + c*x**2)/(2*a*d**2*x**2) + 2*e*sqrt(a + c*x**2)/(a*d**3*x) - 3*e**2*atanh(sqrt(a + c*x**2)/sqrt(a))/(sqrt(a)*d**4) + c*atanh(sqrt(a + c*x**2)/sqrt(a))/(2*a**(sympy.S(3)/2)*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_349():
    f = x**2*(a + b*x)**n*(c + d*x**2)
    F = a**2*(a + b*x)**(n + 1)*(a**2*d + b**2*c)/(b**5*(n + 1)) - 4*a*d*(a + b*x)**(n + 4)/(b**5*(n + 4)) - 2*a*(a + b*x)**(n + 2)*(2*a**2*d + b**2*c)/(b**5*(n + 2)) + d*(a + b*x)**(n + 5)/(b**5*(n + 5)) + (a + b*x)**(n + 3)*(6*a**2*d + b**2*c)/(b**5*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_350():
    f = x*(a + b*x)**n*(c + d*x**2)
    F = -3*a*d*(a + b*x)**(n + 3)/(b**4*(n + 3)) - a*(a + b*x)**(n + 1)*(a**2*d + b**2*c)/(b**4*(n + 1)) + d*(a + b*x)**(n + 4)/(b**4*(n + 4)) + (a + b*x)**(n + 2)*(3*a**2*d + b**2*c)/(b**4*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_351():
    f = (a + b*x)**n*(c + d*x**2)
    F = -2*a*d*(a + b*x)**(n + 2)/(b**3*(n + 2)) + d*(a + b*x)**(n + 3)/(b**3*(n + 3)) + (a + b*x)**(n + 1)*(a**2*d + b**2*c)/(b**3*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_352():
    f = (a + b*x)**n*(c + d*x**2)/x
    F = -a*d*(a + b*x)**(n + 1)/(b**2*(n + 1)) + d*(a + b*x)**(n + 2)/(b**2*(n + 2)) - c*(a + b*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b*x/a)/(a*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_353():
    f = x**2*(a + b*x)**n*(c + d*x**2)**2
    F = a**2*(a + b*x)**(n + 1)*(a**2*d + b**2*c)**2/(b**7*(n + 1)) - 6*a*d**2*(a + b*x)**(n + 6)/(b**7*(n + 6)) - 4*a*d*(a + b*x)**(n + 4)*(5*a**2*d + 2*b**2*c)/(b**7*(n + 4)) - 2*a*(a + b*x)**(n + 2)*(a**2*d + b**2*c)*(3*a**2*d + b**2*c)/(b**7*(n + 2)) + d**2*(a + b*x)**(n + 7)/(b**7*(n + 7)) + d*(a + b*x)**(n + 5)*(15*a**2*d + 2*b**2*c)/(b**7*(n + 5)) + (a + b*x)**(n + 3)*(15*a**4*d**2 + 12*a**2*b**2*c*d + b**4*c**2)/(b**7*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_354():
    f = x*(a + b*x)**n*(c + d*x**2)**2
    F = -5*a*d**2*(a + b*x)**(n + 5)/(b**6*(n + 5)) - 2*a*d*(a + b*x)**(n + 3)*(5*a**2*d + 3*b**2*c)/(b**6*(n + 3)) - a*(a + b*x)**(n + 1)*(a**2*d + b**2*c)**2/(b**6*(n + 1)) + d**2*(a + b*x)**(n + 6)/(b**6*(n + 6)) + 2*d*(a + b*x)**(n + 4)*(5*a**2*d + b**2*c)/(b**6*(n + 4)) + (a + b*x)**(n + 2)*(a**2*d + b**2*c)*(5*a**2*d + b**2*c)/(b**6*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_355():
    f = (a + b*x)**n*(c + d*x**2)**2
    F = -4*a*d**2*(a + b*x)**(n + 4)/(b**5*(n + 4)) - 4*a*d*(a + b*x)**(n + 2)*(a**2*d + b**2*c)/(b**5*(n + 2)) + d**2*(a + b*x)**(n + 5)/(b**5*(n + 5)) + 2*d*(a + b*x)**(n + 3)*(3*a**2*d + b**2*c)/(b**5*(n + 3)) + (a + b*x)**(n + 1)*(a**2*d + b**2*c)**2/(b**5*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_356():
    f = (a + b*x)**n*(c + d*x**2)**2/x
    F = -3*a*d**2*(a + b*x)**(n + 3)/(b**4*(n + 3)) - a*d*(a + b*x)**(n + 1)*(a**2*d + 2*b**2*c)/(b**4*(n + 1)) + d**2*(a + b*x)**(n + 4)/(b**4*(n + 4)) + d*(a + b*x)**(n + 2)*(3*a**2*d + 2*b**2*c)/(b**4*(n + 2)) - c**2*(a + b*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b*x/a)/(a*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_357():
    f = x**2*(a + b*x)**n*(c + d*x**2)**3
    F = a**2*(a + b*x)**(n + 1)*(a**2*d + b**2*c)**3/(b**9*(n + 1)) - 8*a*d**3*(a + b*x)**(n + 8)/(b**9*(n + 8)) - 2*a*d**2*(a + b*x)**(n + 6)*(28*a**2*d + 9*b**2*c)/(b**9*(n + 6)) - 4*a*d*(a + b*x)**(n + 4)*(14*a**4*d**2 + 15*a**2*b**2*c*d + 3*b**4*c**2)/(b**9*(n + 4)) - 2*a*(a + b*x)**(n + 2)*(a**2*d + b**2*c)**2*(4*a**2*d + b**2*c)/(b**9*(n + 2)) + d**3*(a + b*x)**(n + 9)/(b**9*(n + 9)) + d**2*(a + b*x)**(n + 7)*(28*a**2*d + 3*b**2*c)/(b**9*(n + 7)) + d*(a + b*x)**(n + 5)*(70*a**4*d**2 + 45*a**2*b**2*c*d + 3*b**4*c**2)/(b**9*(n + 5)) + (a + b*x)**(n + 3)*(a**2*d + b**2*c)*(28*a**4*d**2 + 17*a**2*b**2*c*d + b**4*c**2)/(b**9*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_358():
    f = x*(a + b*x)**n*(c + d*x**2)**3
    F = -7*a*d**3*(a + b*x)**(n + 7)/(b**8*(n + 7)) - 5*a*d**2*(a + b*x)**(n + 5)*(7*a**2*d + 3*b**2*c)/(b**8*(n + 5)) - 3*a*d*(a + b*x)**(n + 3)*(a**2*d + b**2*c)*(7*a**2*d + 3*b**2*c)/(b**8*(n + 3)) - a*(a + b*x)**(n + 1)*(a**2*d + b**2*c)**3/(b**8*(n + 1)) + d**3*(a + b*x)**(n + 8)/(b**8*(n + 8)) + 3*d**2*(a + b*x)**(n + 6)*(7*a**2*d + b**2*c)/(b**8*(n + 6)) + d*(a + b*x)**(n + 4)*(35*a**4*d**2 + 30*a**2*b**2*c*d + 3*b**4*c**2)/(b**8*(n + 4)) + (a + b*x)**(n + 2)*(a**2*d + b**2*c)**2*(7*a**2*d + b**2*c)/(b**8*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_359():
    f = (a + b*x)**n*(c + d*x**2)**3
    F = -6*a*d**3*(a + b*x)**(n + 6)/(b**7*(n + 6)) - 4*a*d**2*(a + b*x)**(n + 4)*(5*a**2*d + 3*b**2*c)/(b**7*(n + 4)) - 6*a*d*(a + b*x)**(n + 2)*(a**2*d + b**2*c)**2/(b**7*(n + 2)) + d**3*(a + b*x)**(n + 7)/(b**7*(n + 7)) + 3*d**2*(a + b*x)**(n + 5)*(5*a**2*d + b**2*c)/(b**7*(n + 5)) + 3*d*(a + b*x)**(n + 3)*(a**2*d + b**2*c)*(5*a**2*d + b**2*c)/(b**7*(n + 3)) + (a + b*x)**(n + 1)*(a**2*d + b**2*c)**3/(b**7*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_360():
    f = (a + b*x)**n*(c + d*x**2)**3/x
    F = -5*a*d**3*(a + b*x)**(n + 5)/(b**6*(n + 5)) - a*d**2*(a + b*x)**(n + 3)*(10*a**2*d + 9*b**2*c)/(b**6*(n + 3)) - a*d*(a + b*x)**(n + 1)*(a**4*d**2 + 3*a**2*b**2*c*d + 3*b**4*c**2)/(b**6*(n + 1)) + d**3*(a + b*x)**(n + 6)/(b**6*(n + 6)) + d**2*(a + b*x)**(n + 4)*(10*a**2*d + 3*b**2*c)/(b**6*(n + 4)) + d*(a + b*x)**(n + 2)*(5*a**4*d**2 + 9*a**2*b**2*c*d + 3*b**4*c**2)/(b**6*(n + 2)) - c**3*(a + b*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b*x/a)/(a*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_361():
    f = x**4*(d + e*x)**n/(a + c*x**2)
    F = -2*d*(d + e*x)**(n + 2)/(c*e**3*(n + 2)) + (d + e*x)**(n + 3)/(c*e**3*(n + 3)) - (-a)**(sympy.S(3)/2)*(d + e*x)**(n + 1)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d + e*sqrt(-a)))/(2*c**2*(n + 1)*(sqrt(c)*d + e*sqrt(-a))) + (-a)**(sympy.S(3)/2)*(d + e*x)**(n + 1)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d - e*sqrt(-a)))/(2*c**2*(n + 1)*(sqrt(c)*d - e*sqrt(-a))) + (d + e*x)**(n + 1)*(-a*e**2 + c*d**2)/(c**2*e**3*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_362():
    f = x**3*(d + e*x)**n/(a + c*x**2)
    F = a*(d + e*x)**(n + 1)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d + e*sqrt(-a)))/(2*c**(sympy.S(3)/2)*(n + 1)*(sqrt(c)*d + e*sqrt(-a))) + a*(d + e*x)**(n + 1)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d - e*sqrt(-a)))/(2*c**(sympy.S(3)/2)*(n + 1)*(sqrt(c)*d - e*sqrt(-a))) - d*(d + e*x)**(n + 1)/(c*e**2*(n + 1)) + (d + e*x)**(n + 2)/(c*e**2*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_363():
    f = x**2*(d + e*x)**n/(a + c*x**2)
    F = -sqrt(-a)*(d + e*x)**(n + 1)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d + e*sqrt(-a)))/(2*c*(n + 1)*(sqrt(c)*d + e*sqrt(-a))) + sqrt(-a)*(d + e*x)**(n + 1)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d - e*sqrt(-a)))/(2*c*(n + 1)*(sqrt(c)*d - e*sqrt(-a))) + (d + e*x)**(n + 1)/(c*e*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_364():
    f = x*(d + e*x)**n/(a + c*x**2)
    F = -(d + e*x)**(n + 1)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d + e*sqrt(-a)))/(2*sqrt(c)*(n + 1)*(sqrt(c)*d + e*sqrt(-a))) - (d + e*x)**(n + 1)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d - e*sqrt(-a)))/(2*sqrt(c)*(n + 1)*(sqrt(c)*d - e*sqrt(-a)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_365():
    f = (d + e*x)**n/(a + c*x**2)
    F = -(d + e*x)**(n + 1)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d + e*sqrt(-a)))/(2*sqrt(-a)*(n + 1)*(sqrt(c)*d + e*sqrt(-a))) + (d + e*x)**(n + 1)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d - e*sqrt(-a)))/(2*sqrt(-a)*(n + 1)*(sqrt(c)*d - e*sqrt(-a)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_366():
    f = (d + e*x)**n/(x*(a + c*x**2))
    F = sqrt(c)*(d + e*x)**(n + 1)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d + e*sqrt(-a)))/(2*a*(n + 1)*(sqrt(c)*d + e*sqrt(-a))) + sqrt(c)*(d + e*x)**(n + 1)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d - e*sqrt(-a)))/(2*a*(n + 1)*(sqrt(c)*d - e*sqrt(-a))) - (d + e*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + e*x/d)/(a*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_367():
    f = (d + e*x)**n/(x**2*(a + c*x**2))
    F = -c*(d + e*x)**(n + 1)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d + e*sqrt(-a)))/(2*(-a)**(sympy.S(3)/2)*(n + 1)*(sqrt(c)*d + e*sqrt(-a))) + c*(d + e*x)**(n + 1)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d - e*sqrt(-a)))/(2*(-a)**(sympy.S(3)/2)*(n + 1)*(sqrt(c)*d - e*sqrt(-a))) + e*(d + e*x)**(n + 1)*hyper((2, n + 1), (n + 2,), 1 + e*x/d)/(a*d**2*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_368():
    f = x**4*(d + e*x)**n/(a + c*x**2)**2
    F = a*(d + e*x)**(n + 1)*(a*e + c*d*x)/(2*c**2*(a + c*x**2)*(a*e**2 + c*d**2)) - (d + e*x)**(n + 1)*(-a*sqrt(c)*d*e*n + a*e**2*sqrt(-a)*(n + 3) + 3*c*d**2*sqrt(-a))*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d + e*sqrt(-a)))/(4*c**2*(n + 1)*(a*e**2 + c*d**2)*(sqrt(c)*d + e*sqrt(-a))) + (d + e*x)**(n + 1)*(a*sqrt(c)*d*e*n + a*e**2*sqrt(-a)*(n + 3) + 3*c*d**2*sqrt(-a))*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d - e*sqrt(-a)))/(4*c**2*(n + 1)*(a*e**2 + c*d**2)*(sqrt(c)*d - e*sqrt(-a))) + (d + e*x)**(n + 1)/(c**2*e*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_369():
    f = x**3*(d + e*x)**n/(a + c*x**2)**2
    F = a*(d - e*x)*(d + e*x)**(n + 1)/(2*c*(a + c*x**2)*(a*e**2 + c*d**2)) + (d + e*x)**(n + 1)*(d*e*n*sqrt(-a) - (a*e**2*(n + 2) + 2*c*d**2)/sqrt(c))*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d - e*sqrt(-a)))/(4*c*(n + 1)*(a*e**2 + c*d**2)*(sqrt(c)*d - e*sqrt(-a))) - (d + e*x)**(n + 1)*(a*e**2*(n + 2) + sqrt(c)*d*e*n*sqrt(-a) + 2*c*d**2)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d + e*sqrt(-a)))/(4*c**(sympy.S(3)/2)*(n + 1)*(a*e**2 + c*d**2)*(sqrt(c)*d + e*sqrt(-a)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_370():
    f = x**2*(d + e*x)**n/(a + c*x**2)**2
    F = -(d + e*x)**(n + 1)*(a*e + c*d*x)/(2*c*(a + c*x**2)*(a*e**2 + c*d**2)) - (d + e*x)**(n + 1)*(a*e**2*(n + 1) + sqrt(c)*d*e*n*sqrt(-a) + c*d**2)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d + e*sqrt(-a)))/(4*c*sqrt(-a)*(n + 1)*(a*e**2 + c*d**2)*(sqrt(c)*d + e*sqrt(-a))) + (d + e*x)**(n + 1)*(a*e**2*(n + 1) - sqrt(c)*d*e*n*sqrt(-a) + c*d**2)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d - e*sqrt(-a)))/(4*c*sqrt(-a)*(n + 1)*(a*e**2 + c*d**2)*(sqrt(c)*d - e*sqrt(-a)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_371():
    f = x*(d + e*x)**n/(a + c*x**2)**2
    F = -(d - e*x)*(d + e*x)**(n + 1)/((a + c*x**2)*(2*a*e**2 + 2*c*d**2)) + e*n*(d + e*x)**(n + 1)*(sqrt(c)*d + e*sqrt(-a))*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d - e*sqrt(-a)))/(4*sqrt(c)*sqrt(-a)*(n + 1)*(a*e**2 + c*d**2)*(sqrt(c)*d - e*sqrt(-a))) + e*n*(d + e*x)**(n + 1)*(a*e + sqrt(c)*d*sqrt(-a))*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d + e*sqrt(-a)))/(4*a*sqrt(c)*(n + 1)*(a*e**2 + c*d**2)*(sqrt(c)*d + e*sqrt(-a)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_372():
    f = (d + e*x)**n/(a + c*x**2)**2
    F = (d + e*x)**(n + 1)*(a*e**2*(1 - n) - sqrt(c)*d*e*n*sqrt(-a) + c*d**2)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d + e*sqrt(-a)))/(4*(-a)**(sympy.S(3)/2)*(n + 1)*(a*e**2 + c*d**2)*(sqrt(c)*d + e*sqrt(-a))) - (d + e*x)**(n + 1)*(a*e**2*(1 - n) + sqrt(c)*d*e*n*sqrt(-a) + c*d**2)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d - e*sqrt(-a)))/(4*(-a)**(sympy.S(3)/2)*(n + 1)*(a*e**2 + c*d**2)*(sqrt(c)*d - e*sqrt(-a))) + (d + e*x)**(n + 1)*(a*e + c*d*x)/(2*a*(a + c*x**2)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_373():
    f = (d + e*x)**n/(x*(a + c*x**2)**2)
    F = sqrt(c)*e*n*(d + e*x)**(n + 1)*(sqrt(c)*d + e*sqrt(-a))*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d - e*sqrt(-a)))/(4*(-a)**(sympy.S(3)/2)*(n + 1)*(a*e**2 + c*d**2)*(sqrt(c)*d - e*sqrt(-a))) + c*(d - e*x)*(d + e*x)**(n + 1)/(2*a*(a + c*x**2)*(a*e**2 + c*d**2)) - sqrt(c)*e*n*(d + e*x)**(n + 1)*(a*e + sqrt(c)*d*sqrt(-a))*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d + e*sqrt(-a)))/(4*a**2*(n + 1)*(a*e**2 + c*d**2)*(sqrt(c)*d + e*sqrt(-a))) + sqrt(c)*(d + e*x)**(n + 1)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d + e*sqrt(-a)))/(2*a**2*(n + 1)*(sqrt(c)*d + e*sqrt(-a))) + sqrt(c)*(d + e*x)**(n + 1)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d - e*sqrt(-a)))/(2*a**2*(n + 1)*(sqrt(c)*d - e*sqrt(-a))) - (d + e*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + e*x/d)/(a**2*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_374():
    f = (d + e*x)**n/(x**2*(a + c*x**2)**2)
    F = c*(d + e*x)**(n + 1)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d + e*sqrt(-a)))/(2*(-a)**(sympy.S(5)/2)*(n + 1)*(sqrt(c)*d + e*sqrt(-a))) - c*(d + e*x)**(n + 1)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d - e*sqrt(-a)))/(2*(-a)**(sympy.S(5)/2)*(n + 1)*(sqrt(c)*d - e*sqrt(-a))) + c*(d + e*x)**(n + 1)*(a*e**2*(1 - n) - sqrt(c)*d*e*n*sqrt(-a) + c*d**2)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d + e*sqrt(-a)))/(4*(-a)**(sympy.S(5)/2)*(n + 1)*(a*e**2 + c*d**2)*(sqrt(c)*d + e*sqrt(-a))) - c*(d + e*x)**(n + 1)*(a*e**2*(1 - n) + sqrt(c)*d*e*n*sqrt(-a) + c*d**2)*hyper((1, n + 1), (n + 2,), sqrt(c)*(d + e*x)/(sqrt(c)*d - e*sqrt(-a)))/(4*(-a)**(sympy.S(5)/2)*(n + 1)*(a*e**2 + c*d**2)*(sqrt(c)*d - e*sqrt(-a))) - c*(d + e*x)**(n + 1)*(a*e + c*d*x)/(2*a**2*(a + c*x**2)*(a*e**2 + c*d**2)) + e*(d + e*x)**(n + 1)*hyper((2, n + 1), (n + 2,), 1 + e*x/d)/(a**2*d**2*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_375():
    f = (g*x)**m*(a + c*x**2)**2*(d + e*x)**n
    F = -c**2*d*(g*x)**(m + 3)*(d + e*x)**(n + 1)*(m + 4)/(e**2*g**3*(m + n + 4)*(m + n + 5)) + c**2*(g*x)**(m + 4)*(d + e*x)**(n + 1)/(e*g**4*(m + n + 5)) - c*d*(g*x)**(m + 1)*(d + e*x)**(n + 1)*(m + 2)*(2*a*e**2*(m**2 + m*(2*n + 9) + n**2 + 9*n + 20) + c*d**2*(m**2 + 7*m + 12))/(e**4*g*(m + n + 2)*(m + n + 3)*(m + n + 4)*(m + n + 5)) + c*(g*x)**(m + 2)*(d + e*x)**(n + 1)*(2*a*e**2*(m**2 + m*(2*n + 9) + n**2 + 9*n + 20) + c*d**2*(m**2 + 7*m + 12))/(e**3*g**2*(m + n + 3)*(m + n + 4)*(m + n + 5)) + (g*x)**(m + 1)*(d + e*x)**n*(a**2*e**4*(m + n + 2)*(m + n + 3)*(m + n + 4)*(m + n + 5) + c*d**2*(m + 1)*(m + 2)*(2*a*e**2*(m**2 + m*(2*n + 9) + n**2 + 9*n + 20) + c*d**2*(m**2 + 7*m + 12)))*hyper((-n, m + 1), (m + 2,), -e*x/d)/(e**4*g*(1 + e*x/d)**n*(m + 1)*(m + n + 2)*(m + n + 3)*(m + n + 4)*(m + n + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_376():
    f = (g*x)**m*(a + c*x**2)*(d + e*x)**n
    F = -c*d*(g*x)**(m + 1)*(d + e*x)**(n + 1)*(m + 2)/(e**2*g*(m + n + 2)*(m + n + 3)) + c*(g*x)**(m + 2)*(d + e*x)**(n + 1)/(e*g**2*(m + n + 3)) + (g*x)**(m + 1)*(d + e*x)**n*(a*e**2*(m + n + 2)*(m + n + 3) + c*d**2*(m + 1)*(m + 2))*hyper((-n, m + 1), (m + 2,), -e*x/d)/(e**2*g*(1 + e*x/d)**n*(m + 1)*(m + n + 2)*(m + n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_377():
    f = (g*x)**m*(d + e*x)**n/(a + c*x**2)
    F = (g*x)**(m + 1)*(d + e*x)**n*appellf1(m + 1, 1, -n, m + 2, -sqrt(c)*x/sqrt(-a), -e*x/d)/(2*a*g*(1 + e*x/d)**n*(m + 1)) + (g*x)**(m + 1)*(d + e*x)**n*appellf1(m + 1, 1, -n, m + 2, sqrt(c)*x/sqrt(-a), -e*x/d)/(2*a*g*(1 + e*x/d)**n*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_378():
    f = (g*x)**m*(d + e*x)**n/(a + c*x**2)**2
    F = (g*x)**(m + 1)*(d + e*x)**n*appellf1(m + 1, 1, -n, m + 2, -sqrt(c)*x/sqrt(-a), -e*x/d)/(4*a**2*g*(1 + e*x/d)**n*(m + 1)) + (g*x)**(m + 1)*(d + e*x)**n*appellf1(m + 1, 1, -n, m + 2, sqrt(c)*x/sqrt(-a), -e*x/d)/(4*a**2*g*(1 + e*x/d)**n*(m + 1)) + (g*x)**(m + 1)*(d + e*x)**n*appellf1(m + 1, 2, -n, m + 2, -sqrt(c)*x/sqrt(-a), -e*x/d)/(4*a**2*g*(1 + e*x/d)**n*(m + 1)) + (g*x)**(m + 1)*(d + e*x)**n*appellf1(m + 1, 2, -n, m + 2, sqrt(c)*x/sqrt(-a), -e*x/d)/(4*a**2*g*(1 + e*x/d)**n*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_379():
    f = x**5*(a + b*x**2)**p*(d + e*x)
    F = a**2*d*(a + b*x**2)**(p + 1)/(2*b**3*(p + 1)) - a*d*(a + b*x**2)**(p + 2)/(b**3*(p + 2)) + e*x**7*(a + b*x**2)**p*hyper((sympy.S(7)/2, -p), (sympy.S(9)/2,), -b*x**2/a)/(7*(1 + b*x**2/a)**p) + d*(a + b*x**2)**(p + 3)/(2*b**3*(p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_380():
    f = x**4*(a + b*x**2)**p*(d + e*x)
    F = a**2*e*(a + b*x**2)**(p + 1)/(2*b**3*(p + 1)) - a*e*(a + b*x**2)**(p + 2)/(b**3*(p + 2)) + d*x**5*(a + b*x**2)**p*hyper((sympy.S(5)/2, -p), (sympy.S(7)/2,), -b*x**2/a)/(5*(1 + b*x**2/a)**p) + e*(a + b*x**2)**(p + 3)/(2*b**3*(p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_381():
    f = x**3*(a + b*x**2)**p*(d + e*x)
    F = -a*d*(a + b*x**2)**(p + 1)/(2*b**2*(p + 1)) + e*x**5*(a + b*x**2)**p*hyper((sympy.S(5)/2, -p), (sympy.S(7)/2,), -b*x**2/a)/(5*(1 + b*x**2/a)**p) + d*(a + b*x**2)**(p + 2)/(2*b**2*(p + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_382():
    f = x**2*(a + b*x**2)**p*(d + e*x)
    F = -a*e*(a + b*x**2)**(p + 1)/(2*b**2*(p + 1)) + d*x**3*(a + b*x**2)**p*hyper((sympy.S(3)/2, -p), (sympy.S(5)/2,), -b*x**2/a)/(3*(1 + b*x**2/a)**p) + e*(a + b*x**2)**(p + 2)/(2*b**2*(p + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_383():
    f = x*(a + b*x**2)**p*(d + e*x)
    F = e*x**3*(a + b*x**2)**p*hyper((sympy.S(3)/2, -p), (sympy.S(5)/2,), -b*x**2/a)/(3*(1 + b*x**2/a)**p) + d*(a + b*x**2)**(p + 1)/(2*b*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_384():
    f = (a + b*x**2)**p*(d + e*x)
    F = d*x*(a + b*x**2)**p*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*x**2/a)/(1 + b*x**2/a)**p + e*(a + b*x**2)**(p + 1)/(2*b*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_385():
    f = (a + b*x**2)**p*(d + e*x)/x
    F = e*x*(a + b*x**2)**p*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*x**2/a)/(1 + b*x**2/a)**p - d*(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*x**2/a)/(2*a*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_386():
    f = (a + b*x**2)**p*(d + e*x)/x**2
    F = -d*(a + b*x**2)**p*hyper((sympy.S(-1)/2, -p), (sympy.S.Half,), -b*x**2/a)/(x*(1 + b*x**2/a)**p) - e*(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*x**2/a)/(2*a*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_387():
    f = (a + b*x**2)**p*(d + e*x)/x**3
    F = -e*(a + b*x**2)**p*hyper((sympy.S(-1)/2, -p), (sympy.S.Half,), -b*x**2/a)/(x*(1 + b*x**2/a)**p) + b*d*(a + b*x**2)**(p + 1)*hyper((2, p + 1), (p + 2,), 1 + b*x**2/a)/(2*a**2*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_388():
    f = x**5*(a + b*x**2)**p*(d + e*x)**2
    F = a**2*(a + b*x**2)**(p + 1)*(-a*e**2 + b*d**2)/(2*b**4*(p + 1)) - a*(a + b*x**2)**(p + 2)*(-3*a*e**2 + 2*b*d**2)/(2*b**4*(p + 2)) + 2*d*e*x**7*(a + b*x**2)**p*hyper((sympy.S(7)/2, -p), (sympy.S(9)/2,), -b*x**2/a)/(7*(1 + b*x**2/a)**p) + e**2*(a + b*x**2)**(p + 4)/(2*b**4*(p + 4)) + (a + b*x**2)**(p + 3)*(-3*a*e**2 + b*d**2)/(2*b**4*(p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_389():
    f = x**4*(a + b*x**2)**p*(d + e*x)**2
    F = a**2*d*e*(a + b*x**2)**(p + 1)/(b**3*(p + 1)) - 2*a*d*e*(a + b*x**2)**(p + 2)/(b**3*(p + 2)) + e**2*x**5*(a + b*x**2)**(p + 1)/(b*(2*p + 7)) - x**5*(a + b*x**2)**p*(5*a*e**2 - b*d**2*(2*p + 7))*hyper((sympy.S(5)/2, -p), (sympy.S(7)/2,), -b*x**2/a)/(5*b*(1 + b*x**2/a)**p*(2*p + 7)) + d*e*(a + b*x**2)**(p + 3)/(b**3*(p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_390():
    f = x**3*(a + b*x**2)**p*(d + e*x)**2
    F = -a*(a + b*x**2)**(p + 1)*(-a*e**2 + b*d**2)/(2*b**3*(p + 1)) + 2*d*e*x**5*(a + b*x**2)**p*hyper((sympy.S(5)/2, -p), (sympy.S(7)/2,), -b*x**2/a)/(5*(1 + b*x**2/a)**p) + e**2*(a + b*x**2)**(p + 3)/(2*b**3*(p + 3)) + (a + b*x**2)**(p + 2)*(-2*a*e**2 + b*d**2)/(2*b**3*(p + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_391():
    f = x**2*(a + b*x**2)**p*(d + e*x)**2
    F = -a*d*e*(a + b*x**2)**(p + 1)/(b**2*(p + 1)) + e**2*x**3*(a + b*x**2)**(p + 1)/(b*(2*p + 5)) - x**3*(a + b*x**2)**p*(3*a*e**2 - b*d**2*(2*p + 5))*hyper((sympy.S(3)/2, -p), (sympy.S(5)/2,), -b*x**2/a)/(3*b*(1 + b*x**2/a)**p*(2*p + 5)) + d*e*(a + b*x**2)**(p + 2)/(b**2*(p + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_392():
    f = x*(a + b*x**2)**p*(d + e*x)**2
    F = 2*d*e*x**3*(a + b*x**2)**p*hyper((sympy.S(3)/2, -p), (sympy.S(5)/2,), -b*x**2/a)/(3*(1 + b*x**2/a)**p) + e**2*(a + b*x**2)**(p + 2)/(2*b**2*(p + 2)) + (a + b*x**2)**(p + 1)*(-a*e**2 + b*d**2)/(2*b**2*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_393():
    f = (a + b*x**2)**p*(d + e*x)**2
    F = d*e*(a + b*x**2)**(p + 1)*(p + 2)/(b*(p + 1)*(2*p + 3)) + e*(a + b*x**2)**(p + 1)*(d + e*x)/(b*(2*p + 3)) - x*(a + b*x**2)**p*(a*e**2 - b*d**2*(2*p + 3))*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*x**2/a)/(b*(1 + b*x**2/a)**p*(2*p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_394():
    f = (a + b*x**2)**p*(d + e*x)**2/x
    F = 2*d*e*x*(a + b*x**2)**p*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*x**2/a)/(1 + b*x**2/a)**p + e**2*(a + b*x**2)**(p + 1)/(2*b*(p + 1)) - d**2*(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*x**2/a)/(2*a*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_395():
    f = (a + b*x**2)**p*(d + e*x)**2/x**2
    F = -d**2*(a + b*x**2)**(p + 1)/(a*x) - d*e*(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*x**2/a)/(a*(p + 1)) + x*(a + b*x**2)**p*(a*e**2 + b*d**2*(2*p + 1))*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*x**2/a)/(a*(1 + b*x**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_396():
    f = (a + b*x**2)**p*(d + e*x)**2/x**3
    F = -2*d*e*(a + b*x**2)**p*hyper((sympy.S(-1)/2, -p), (sympy.S.Half,), -b*x**2/a)/(x*(1 + b*x**2/a)**p) - d**2*(a + b*x**2)**(p + 1)/(2*a*x**2) - (a + b*x**2)**(p + 1)*(a*e**2 + b*d**2*p)*hyper((1, p + 1), (p + 2,), 1 + b*x**2/a)/(2*a**2*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_397():
    f = x**5*(a + b*x**2)**p*(d + e*x)**3
    F = a**2*d*(a + b*x**2)**(p + 1)*(-3*a*e**2 + b*d**2)/(2*b**4*(p + 1)) - a*d*(a + b*x**2)**(p + 2)*(-9*a*e**2 + 2*b*d**2)/(2*b**4*(p + 2)) + e**3*x**7*(a + b*x**2)**(p + 1)/(b*(2*p + 9)) - e*x**7*(a + b*x**2)**p*(7*a*e**2 - 3*b*d**2*(2*p + 9))*hyper((sympy.S(7)/2, -p), (sympy.S(9)/2,), -b*x**2/a)/(7*b*(1 + b*x**2/a)**p*(2*p + 9)) + 3*d*e**2*(a + b*x**2)**(p + 4)/(2*b**4*(p + 4)) + d*(a + b*x**2)**(p + 3)*(-9*a*e**2 + b*d**2)/(2*b**4*(p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_398():
    f = x**4*(a + b*x**2)**p*(d + e*x)**3
    F = a**2*e*(a + b*x**2)**(p + 1)*(-a*e**2 + 3*b*d**2)/(2*b**4*(p + 1)) - 3*a*e*(a + b*x**2)**(p + 2)*(-a*e**2 + 2*b*d**2)/(2*b**4*(p + 2)) + 3*d*e**2*x**5*(a + b*x**2)**(p + 1)/(b*(2*p + 7)) - d*x**5*(a + b*x**2)**p*(15*a*e**2 - b*d**2*(2*p + 7))*hyper((sympy.S(5)/2, -p), (sympy.S(7)/2,), -b*x**2/a)/(5*b*(1 + b*x**2/a)**p*(2*p + 7)) + e**3*(a + b*x**2)**(p + 4)/(2*b**4*(p + 4)) + 3*e*(a + b*x**2)**(p + 3)*(-a*e**2 + b*d**2)/(2*b**4*(p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_399():
    f = x**3*(a + b*x**2)**p*(d + e*x)**3
    F = -a*d*(a + b*x**2)**(p + 1)*(-3*a*e**2 + b*d**2)/(2*b**3*(p + 1)) + e**3*x**5*(a + b*x**2)**(p + 1)/(b*(2*p + 7)) - e*x**5*(a + b*x**2)**p*(5*a*e**2 - 3*b*d**2*(2*p + 7))*hyper((sympy.S(5)/2, -p), (sympy.S(7)/2,), -b*x**2/a)/(5*b*(1 + b*x**2/a)**p*(2*p + 7)) + 3*d*e**2*(a + b*x**2)**(p + 3)/(2*b**3*(p + 3)) + d*(a + b*x**2)**(p + 2)*(-6*a*e**2 + b*d**2)/(2*b**3*(p + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_400():
    f = x**2*(a + b*x**2)**p*(d + e*x)**3
    F = -a*e*(a + b*x**2)**(p + 1)*(-a*e**2 + 3*b*d**2)/(2*b**3*(p + 1)) + 3*d*e**2*x**3*(a + b*x**2)**(p + 1)/(b*(2*p + 5)) - d*x**3*(a + b*x**2)**p*(9*a*e**2 - b*d**2*(2*p + 5))*hyper((sympy.S(3)/2, -p), (sympy.S(5)/2,), -b*x**2/a)/(3*b*(1 + b*x**2/a)**p*(2*p + 5)) + e**3*(a + b*x**2)**(p + 3)/(2*b**3*(p + 3)) + e*(a + b*x**2)**(p + 2)*(-2*a*e**2 + 3*b*d**2)/(2*b**3*(p + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_401():
    f = x*(a + b*x**2)**p*(d + e*x)**3
    F = e**3*x**3*(a + b*x**2)**(p + 1)/(b*(2*p + 5)) - e*x**3*(a + b*x**2)**p*(a*e**2 - b*d**2*(2*p + 5))*hyper((sympy.S(3)/2, -p), (sympy.S(5)/2,), -b*x**2/a)/(b*(1 + b*x**2/a)**p*(2*p + 5)) + 3*d*e**2*(a + b*x**2)**(p + 2)/(2*b**2*(p + 2)) + d*(a + b*x**2)**(p + 1)*(-3*a*e**2 + b*d**2)/(2*b**2*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_402():
    f = (a + b*x**2)**p*(d + e*x)**3/x
    F = 3*d*e**2*(a + b*x**2)**(p + 1)/(2*b*(p + 1)) + e**3*x*(a + b*x**2)**(p + 1)/(b*(2*p + 3)) - e*x*(a + b*x**2)**p*(a*e**2 - 3*b*d**2*(2*p + 3))*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*x**2/a)/(b*(1 + b*x**2/a)**p*(2*p + 3)) - d**3*(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*x**2/a)/(2*a*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_403():
    f = (a + b*x**2)**p*(d + e*x)**3/x**2
    F = e**3*(a + b*x**2)**(p + 1)/(2*b*(p + 1)) - d**3*(a + b*x**2)**(p + 1)/(a*x) - 3*d**2*e*(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*x**2/a)/(2*a*(p + 1)) + d*x*(a + b*x**2)**p*(3*a*e**2 + b*d**2*(2*p + 1))*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*x**2/a)/(a*(1 + b*x**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_404():
    f = (a + b*x**2)**p*(d + e*x)**3/x**3
    F = -d**3*(a + b*x**2)**(p + 1)/(2*a*x**2) - 3*d**2*e*(a + b*x**2)**(p + 1)/(a*x) + e*x*(a + b*x**2)**p*(a*e**2 + 3*b*d**2*(2*p + 1))*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*x**2/a)/(a*(1 + b*x**2/a)**p) - d*(a + b*x**2)**(p + 1)*(3*a*e**2 + b*d**2*p)*hyper((1, p + 1), (p + 2,), 1 + b*x**2/a)/(2*a**2*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_405():
    f = x**4*(a + b*x**2)**p/(d + e*x)
    F = -d**4*(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(2*e**3*(p + 1)*(a*e**2 + b*d**2)) + x**5*(a + b*x**2)**p*appellf1(sympy.S(5)/2, 1, -p, sympy.S(7)/2, e**2*x**2/d**2, -b*x**2/a)/(5*d*(1 + b*x**2/a)**p) + (a + b*x**2)**(p + 2)/(2*b**2*e*(p + 2)) + (a + b*x**2)**(p + 1)*(-a*e**2 + b*d**2)/(2*b**2*e**3*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_406():
    f = x**3*(a + b*x**2)**p/(d + e*x)
    F = d**3*(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(2*e**2*(p + 1)*(a*e**2 + b*d**2)) - e*x**5*(a + b*x**2)**p*appellf1(sympy.S(5)/2, 1, -p, sympy.S(7)/2, e**2*x**2/d**2, -b*x**2/a)/(5*d**2*(1 + b*x**2/a)**p) - d*(a + b*x**2)**(p + 1)/(2*b*e**2*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_407():
    f = x**2*(a + b*x**2)**p/(d + e*x)
    F = -d**2*(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(2*e*(p + 1)*(a*e**2 + b*d**2)) + x**3*(a + b*x**2)**p*appellf1(sympy.S(3)/2, 1, -p, sympy.S(5)/2, e**2*x**2/d**2, -b*x**2/a)/(3*d*(1 + b*x**2/a)**p) + (a + b*x**2)**(p + 1)/(2*b*e*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_408():
    f = x*(a + b*x**2)**p/(d + e*x)
    F = d*(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/((p + 1)*(2*a*e**2 + 2*b*d**2)) - x*(a + b*x**2)**p*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(e*(1 + b*x**2/a)**p) + x*(a + b*x**2)**p*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*x**2/a)/(e*(1 + b*x**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_409():
    f = (a + b*x**2)**p/(d + e*x)
    F = -e*(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/((p + 1)*(2*a*e**2 + 2*b*d**2)) + x*(a + b*x**2)**p*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(d*(1 + b*x**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_410():
    f = (a + b*x**2)**p/(x*(d + e*x))
    F = e**2*(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(2*d*(p + 1)*(a*e**2 + b*d**2)) - e*x*(a + b*x**2)**p*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(d**2*(1 + b*x**2/a)**p) - (a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*x**2/a)/(2*a*d*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_411():
    f = (a + b*x**2)**p/(x**2*(d + e*x))
    F = -(a + b*x**2)**p*appellf1(sympy.S(-1)/2, 1, -p, sympy.S.Half, e**2*x**2/d**2, -b*x**2/a)/(d*x*(1 + b*x**2/a)**p) - e**3*(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(2*d**2*(p + 1)*(a*e**2 + b*d**2)) + e*(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*x**2/a)/(2*a*d**2*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_412():
    f = (a + b*x**2)**p/(x**3*(d + e*x))
    F = e*(a + b*x**2)**p*appellf1(sympy.S(-1)/2, 1, -p, sympy.S.Half, e**2*x**2/d**2, -b*x**2/a)/(d**2*x*(1 + b*x**2/a)**p) + e**4*(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(2*d**3*(p + 1)*(a*e**2 + b*d**2)) - (a + b*x**2)**(p + 1)/(2*a*d*x**2) - (a + b*x**2)**(p + 1)*(a*e**2 + b*d**2*p)*hyper((1, p + 1), (p + 2,), 1 + b*x**2/a)/(2*a**2*d**3*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_413():
    f = x**4*(a + b*x**2)**p/(d + e*x)**2
    F = -d**4*(a + b*x**2)**(p + 1)/(e**3*(d + e*x)*(a*e**2 + b*d**2)) + d**3*(a + b*x**2)**(p + 1)*(2*a*e**2 + b*d**2*(p + 2))*hyper((1, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(e**3*(p + 1)*(a*e**2 + b*d**2)**2) - 2*d**2*x*(a + b*x**2)**p*(2*a*e**2 + b*d**2*(p + 2))*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(e**4*(1 + b*x**2/a)**p*(a*e**2 + b*d**2)) - d*(a + b*x**2)**(p + 1)*(3*p + 4)/(b*e**3*(p + 1)*(2*p + 3)) + (a + b*x**2)**(p + 1)*(d + e*x)/(b*e**3*(2*p + 3)) - x*(a + b*x**2)**p*(a**2*e**4 - 2*a*b*d**2*e**2*(3*p + 4) - 2*b**2*d**4*(2*p**2 + 7*p + 6))*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*x**2/a)/(b*e**4*(1 + b*x**2/a)**p*(2*p + 3)*(a*e**2 + b*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_414():
    f = x**3*(a + b*x**2)**p/(d + e*x)**2
    F = d**3*(a + b*x**2)**(p + 1)/(e**2*(d + e*x)*(a*e**2 + b*d**2)) - d**2*(a + b*x**2)**(p + 1)*(3*a*e**2 + b*d**2*(2*p + 3))*hyper((1, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(2*e**2*(p + 1)*(a*e**2 + b*d**2)**2) - d*x*(a + b*x**2)**p*(2*a*e**2 + b*d**2*(2*p + 3))*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*x**2/a)/(e**3*(1 + b*x**2/a)**p*(a*e**2 + b*d**2)) + d*x*(a + b*x**2)**p*(3*a*e**2 + b*d**2*(2*p + 3))*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(e**3*(1 + b*x**2/a)**p*(a*e**2 + b*d**2)) + (a + b*x**2)**(p + 1)/(2*b*e**2*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_415():
    f = x**2*(a + b*x**2)**p/(d + e*x)**2
    F = -d**2*(a + b*x**2)**(p + 1)/(e*(d + e*x)*(a*e**2 + b*d**2)) + d*(a + b*x**2)**(p + 1)*(a*e**2 + b*d**2*(p + 1))*hyper((1, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(e*(p + 1)*(a*e**2 + b*d**2)**2) + x*(a + b*x**2)**p*(a*e**2 + 2*b*d**2*(p + 1))*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*x**2/a)/(e**2*(1 + b*x**2/a)**p*(a*e**2 + b*d**2)) - x*(a + b*x**2)**p*(2*a*e**2 + 2*b*d**2*(p + 1))*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(e**2*(1 + b*x**2/a)**p*(a*e**2 + b*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_416():
    f = x*(a + b*x**2)**p/(d + e*x)**2
    F = -b*d*x*(a + b*x**2)**p*(2*p + 1)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*x**2/a)/(e*(1 + b*x**2/a)**p*(a*e**2 + b*d**2)) + d*(a + b*x**2)**(p + 1)/((d + e*x)*(a*e**2 + b*d**2)) - (a + b*x**2)**(p + 1)*(a*e**2 + b*d**2*(2*p + 1))*hyper((1, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(2*(p + 1)*(a*e**2 + b*d**2)**2) + x*(a + b*x**2)**p*(a*e**2 + b*d**2*(2*p + 1))*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(d*e*(1 + b*x**2/a)**p*(a*e**2 + b*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_417():
    f = (a + b*x**2)**p/(x*(d + e*x)**2)
    F = b*e**2*(a + b*x**2)**(p + 1)*hyper((2, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/((p + 1)*(a*e**2 + b*d**2)**2) + e**2*(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(2*d**2*(p + 1)*(a*e**2 + b*d**2)) - e*x*(a + b*x**2)**p*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(d**3*(1 + b*x**2/a)**p) - e*x*(a + b*x**2)**p*appellf1(sympy.S.Half, 2, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(d**3*(1 + b*x**2/a)**p) - e**3*x**3*(a + b*x**2)**p*appellf1(sympy.S(3)/2, 2, -p, sympy.S(5)/2, e**2*x**2/d**2, -b*x**2/a)/(3*d**5*(1 + b*x**2/a)**p) - (a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*x**2/a)/(2*a*d**2*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_418():
    f = (a + b*x**2)**p/(x**2*(d + e*x)**2)
    F = -b*e**3*(a + b*x**2)**(p + 1)*hyper((2, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(d*(p + 1)*(a*e**2 + b*d**2)**2) - (a + b*x**2)**p*hyper((sympy.S(-1)/2, -p), (sympy.S.Half,), -b*x**2/a)/(d**2*x*(1 + b*x**2/a)**p) - e**3*(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(d**3*(p + 1)*(a*e**2 + b*d**2)) + 2*e**2*x*(a + b*x**2)**p*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(d**4*(1 + b*x**2/a)**p) + e**2*x*(a + b*x**2)**p*appellf1(sympy.S.Half, 2, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(d**4*(1 + b*x**2/a)**p) + e**4*x**3*(a + b*x**2)**p*appellf1(sympy.S(3)/2, 2, -p, sympy.S(5)/2, e**2*x**2/d**2, -b*x**2/a)/(3*d**6*(1 + b*x**2/a)**p) + e*(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*x**2/a)/(a*d**3*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_419():
    f = x**4*(a + b*x**2)**p/(d + e*x)**3
    F = -d**4*(a + b*x**2)**(p + 1)/(2*e**3*(d + e*x)**2*(a*e**2 + b*d**2)) + d**3*(a + b*x**2)**(p + 1)*(4*a*e**2 + b*d**2*(p + 3))/(e**3*(d + e*x)*(a*e**2 + b*d**2)**2) - d**2*(a + b*x**2)**(p + 1)*(6*a**2*e**4 + 3*a*b*d**2*e**2*(3*p + 4) + b**2*d**4*(2*p**2 + 7*p + 6))*hyper((1, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(2*e**3*(p + 1)*(a*e**2 + b*d**2)**3) - d*x*(a + b*x**2)**p*(3*a**2*e**4 + 2*a*b*d**2*e**2*(4*p + 5) + b**2*d**4*(2*p**2 + 7*p + 6))*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*x**2/a)/(e**4*(1 + b*x**2/a)**p*(a*e**2 + b*d**2)**2) + d*x*(a + b*x**2)**p*(6*a**2*e**4 + 3*a*b*d**2*e**2*(3*p + 4) + b**2*d**4*(2*p**2 + 7*p + 6))*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(e**4*(1 + b*x**2/a)**p*(a*e**2 + b*d**2)**2) + (a + b*x**2)**(p + 1)/(2*b*e**3*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_420():
    f = x**3*(a + b*x**2)**p/(d + e*x)**3
    F = d**3*(a + b*x**2)**(p + 1)/(2*e**2*(d + e*x)**2*(a*e**2 + b*d**2)) - d**2*(a + b*x**2)**(p + 1)*(3*a*e**2 + b*d**2*(p + 2))/(e**2*(d + e*x)*(a*e**2 + b*d**2)**2) + d*(a + b*x**2)**(p + 1)*(3*a**2*e**4 + a*b*d**2*e**2*(7*p + 6) + b**2*d**4*(2*p**2 + 5*p + 3))*hyper((1, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(2*e**2*(p + 1)*(a*e**2 + b*d**2)**3) + x*(a + b*x**2)**p*(a**2*e**4 + a*b*d**2*e**2*(6*p + 5) + b**2*d**4*(2*p**2 + 5*p + 3))*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*x**2/a)/(e**3*(1 + b*x**2/a)**p*(a*e**2 + b*d**2)**2) - x*(a + b*x**2)**p*(3*a**2*e**4 + a*b*d**2*e**2*(7*p + 6) + b**2*d**4*(2*p**2 + 5*p + 3))*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(e**3*(1 + b*x**2/a)**p*(a*e**2 + b*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_421():
    f = x**2*(a + b*x**2)**p/(d + e*x)**3
    F = -b*d*x*(a + b*x**2)**p*(2*p + 1)*(2*a*e**2 + b*d**2*(p + 1))*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*x**2/a)/(e**2*(1 + b*x**2/a)**p*(a*e**2 + b*d**2)**2) - d**2*(a + b*x**2)**(p + 1)/(2*e*(d + e*x)**2*(a*e**2 + b*d**2)) + d*(a + b*x**2)**(p + 1)*(2*a*e**2 + b*d**2*(p + 1))/(e*(d + e*x)*(a*e**2 + b*d**2)**2) - (a + b*x**2)**(p + 1)*(a**2*e**4 + a*b*d**2*e**2*(5*p + 2) + b**2*d**4*(2*p**2 + 3*p + 1))*hyper((1, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(2*e*(p + 1)*(a*e**2 + b*d**2)**3) + x*(a + b*x**2)**p*(a**2*e**4 + a*b*d**2*e**2*(5*p + 2) + b**2*d**4*(2*p**2 + 3*p + 1))*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(d*e**2*(1 + b*x**2/a)**p*(a*e**2 + b*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_422():
    f = x*(a + b*x**2)**p/(d + e*x)**3
    F = b*d*p*(a + b*x**2)**(p + 1)*(3*a*e**2 + b*d**2*(2*p + 1))*hyper((1, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(2*(p + 1)*(a*e**2 + b*d**2)**3) - b*p*x*(a + b*x**2)**p*(3*a*e**2 + b*d**2*(2*p + 1))*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(e*(1 + b*x**2/a)**p*(a*e**2 + b*d**2)**2) + b*x*(a + b*x**2)**p*(2*p + 1)*(a*e**2 + b*d**2*p)*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*x**2/a)/(e*(1 + b*x**2/a)**p*(a*e**2 + b*d**2)**2) + d*(a + b*x**2)**(p + 1)/((d + e*x)**2*(2*a*e**2 + 2*b*d**2)) - (a + b*x**2)**(p + 1)*(a*e**2 + b*d**2*p)/((d + e*x)*(a*e**2 + b*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_423():
    f = (a + b*x**2)**p/(d + e*x)**3
    F = -3*b**2*d**2*e*(a + b*x**2)**(p + 1)*hyper((3, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(2*(p + 1)*(a*e**2 + b*d**2)**3) + b*e*(a + b*x**2)**(p + 1)*(2*a*e**2 + b*d**2*(p + 1))*hyper((2, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(4*(p + 1)*(a*e**2 + b*d**2)**3) - d**2*e*(a + b*x**2)**(p + 1)/((d**2 - e**2*x**2)**2*(4*a*e**2 + 4*b*d**2)) + x*(a + b*x**2)**p*appellf1(sympy.S.Half, 3, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(d**3*(1 + b*x**2/a)**p) + e**2*x**3*(a + b*x**2)**p*appellf1(sympy.S(3)/2, 3, -p, sympy.S(5)/2, e**2*x**2/d**2, -b*x**2/a)/(d**5*(1 + b*x**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_424():
    f = (a + b*x**2)**p/(x*(d + e*x)**3)
    F = 3*b**2*d*e**2*(a + b*x**2)**(p + 1)*hyper((3, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(2*(p + 1)*(a*e**2 + b*d**2)**3) + b*e**2*(a + b*x**2)**(p + 1)*hyper((2, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(d*(p + 1)*(a*e**2 + b*d**2)**2) - b*e**2*(a + b*x**2)**(p + 1)*(2*a*e**2 + b*d**2*(p + 1))*hyper((2, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(4*d*(p + 1)*(a*e**2 + b*d**2)**3) + d*e**2*(a + b*x**2)**(p + 1)/((d**2 - e**2*x**2)**2*(4*a*e**2 + 4*b*d**2)) + e**2*(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(2*d**3*(p + 1)*(a*e**2 + b*d**2)) - e*x*(a + b*x**2)**p*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(d**4*(1 + b*x**2/a)**p) - e*x*(a + b*x**2)**p*appellf1(sympy.S.Half, 2, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(d**4*(1 + b*x**2/a)**p) - e*x*(a + b*x**2)**p*appellf1(sympy.S.Half, 3, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(d**4*(1 + b*x**2/a)**p) - e**3*x**3*(a + b*x**2)**p*appellf1(sympy.S(3)/2, 2, -p, sympy.S(5)/2, e**2*x**2/d**2, -b*x**2/a)/(3*d**6*(1 + b*x**2/a)**p) - e**3*x**3*(a + b*x**2)**p*appellf1(sympy.S(3)/2, 3, -p, sympy.S(5)/2, e**2*x**2/d**2, -b*x**2/a)/(d**6*(1 + b*x**2/a)**p) - (a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*x**2/a)/(2*a*d**3*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_425():
    f = (a + b*x**2)**p/(x**2*(d + e*x)**3)
    F = -3*b**2*e**3*(a + b*x**2)**(p + 1)*hyper((3, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(2*(p + 1)*(a*e**2 + b*d**2)**3) - 2*b*e**3*(a + b*x**2)**(p + 1)*hyper((2, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(d**2*(p + 1)*(a*e**2 + b*d**2)**2) + b*e**3*(a + b*x**2)**(p + 1)*(2*a*e**2 + b*d**2*(p + 1))*hyper((2, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(4*d**2*(p + 1)*(a*e**2 + b*d**2)**3) - e**3*(a + b*x**2)**(p + 1)/((d**2 - e**2*x**2)**2*(4*a*e**2 + 4*b*d**2)) - (a + b*x**2)**p*hyper((sympy.S(-1)/2, -p), (sympy.S.Half,), -b*x**2/a)/(d**3*x*(1 + b*x**2/a)**p) - 3*e**3*(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), e**2*(a + b*x**2)/(a*e**2 + b*d**2))/(2*d**4*(p + 1)*(a*e**2 + b*d**2)) + 3*e**2*x*(a + b*x**2)**p*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(d**5*(1 + b*x**2/a)**p) + 2*e**2*x*(a + b*x**2)**p*appellf1(sympy.S.Half, 2, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(d**5*(1 + b*x**2/a)**p) + e**2*x*(a + b*x**2)**p*appellf1(sympy.S.Half, 3, -p, sympy.S(3)/2, e**2*x**2/d**2, -b*x**2/a)/(d**5*(1 + b*x**2/a)**p) + 2*e**4*x**3*(a + b*x**2)**p*appellf1(sympy.S(3)/2, 2, -p, sympy.S(5)/2, e**2*x**2/d**2, -b*x**2/a)/(3*d**7*(1 + b*x**2/a)**p) + e**4*x**3*(a + b*x**2)**p*appellf1(sympy.S(3)/2, 3, -p, sympy.S(5)/2, e**2*x**2/d**2, -b*x**2/a)/(d**7*(1 + b*x**2/a)**p) + 3*e*(a + b*x**2)**(p + 1)*hyper((1, p + 1), (p + 2,), 1 + b*x**2/a)/(2*a*d**4*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_426():
    f = (g*x)**m*(a + c*x**2)**p*(d + e*x)**3
    F = 3*d*e**2*(g*x)**(m + 1)*(a + c*x**2)**(p + 1)/(c*g*(m + 2*p + 3)) - d*(g*x)**(m + 1)*(a + c*x**2)**p*(3*a*e**2*(m + 1) - c*d**2*(m + 2*p + 3))*hyper((-p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -c*x**2/a)/(c*g*(1 + c*x**2/a)**p*(m + 1)*(m + 2*p + 3)) + e**3*(g*x)**(m + 2)*(a + c*x**2)**(p + 1)/(c*g**2*(m + 2*p + 4)) - e*(g*x)**(m + 2)*(a + c*x**2)**p*(a*e**2*(m + 2) - 3*c*d**2*(m + 2*p + 4))*hyper((-p, m/2 + 1), (m/2 + 2,), -c*x**2/a)/(c*g**2*(1 + c*x**2/a)**p*(m + 2)*(m + 2*p + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_427():
    f = (g*x)**m*(a + c*x**2)**p*(d + e*x)**2
    F = 2*d*e*(g*x)**(m + 2)*(a + c*x**2)**p*hyper((-p, m/2 + 1), (m/2 + 2,), -c*x**2/a)/(g**2*(1 + c*x**2/a)**p*(m + 2)) + e**2*(g*x)**(m + 1)*(a + c*x**2)**(p + 1)/(c*g*(m + 2*p + 3)) - (g*x)**(m + 1)*(a + c*x**2)**p*(a*e**2*(m + 1) - c*d**2*(m + 2*p + 3))*hyper((-p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -c*x**2/a)/(c*g*(1 + c*x**2/a)**p*(m + 1)*(m + 2*p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_428():
    f = (g*x)**m*(a + c*x**2)**p*(d + e*x)
    F = d*(g*x)**(m + 1)*(a + c*x**2)**p*hyper((-p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -c*x**2/a)/(g*(1 + c*x**2/a)**p*(m + 1)) + e*(g*x)**(m + 2)*(a + c*x**2)**p*hyper((-p, m/2 + 1), (m/2 + 2,), -c*x**2/a)/(g**2*(1 + c*x**2/a)**p*(m + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_429():
    f = (g*x)**m*(a + c*x**2)**p
    F = (g*x)**(m + 1)*(a + c*x**2)**p*hyper((-p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -c*x**2/a)/(g*(1 + c*x**2/a)**p*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_430():
    f = (g*x)**m*(a + c*x**2)**p/(d + e*x)
    F = x*(g*x)**m*(a + c*x**2)**p*appellf1(m/2 + sympy.S.Half, 1, -p, m/2 + sympy.S(3)/2, e**2*x**2/d**2, -c*x**2/a)/(d*(1 + c*x**2/a)**p*(m + 1)) - e*x**2*(g*x)**m*(a + c*x**2)**p*appellf1(m/2 + 1, 1, -p, m/2 + 2, e**2*x**2/d**2, -c*x**2/a)/(d**2*(1 + c*x**2/a)**p*(m + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_431():
    f = (g*x)**m*(a + c*x**2)**p/(d + e*x)**2
    F = x*(g*x)**m*(a + c*x**2)**p*appellf1(m/2 + sympy.S.Half, 2, -p, m/2 + sympy.S(3)/2, e**2*x**2/d**2, -c*x**2/a)/(d**2*(1 + c*x**2/a)**p*(m + 1)) - 2*e*x**2*(g*x)**m*(a + c*x**2)**p*appellf1(m/2 + 1, 2, -p, m/2 + 2, e**2*x**2/d**2, -c*x**2/a)/(d**3*(1 + c*x**2/a)**p*(m + 2)) + e**2*x**3*(g*x)**m*(a + c*x**2)**p*appellf1(m/2 + sympy.S(3)/2, 2, -p, m/2 + sympy.S(5)/2, e**2*x**2/d**2, -c*x**2/a)/(d**4*(1 + c*x**2/a)**p*(m + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_432():
    f = (g*x)**m*(a + c*x**2)**p/(d + e*x)**3
    F = x*(g*x)**m*(a + c*x**2)**p*appellf1(m/2 + sympy.S.Half, 3, -p, m/2 + sympy.S(3)/2, e**2*x**2/d**2, -c*x**2/a)/(d**3*(1 + c*x**2/a)**p*(m + 1)) - 3*e*x**2*(g*x)**m*(a + c*x**2)**p*appellf1(m/2 + 1, 3, -p, m/2 + 2, e**2*x**2/d**2, -c*x**2/a)/(d**4*(1 + c*x**2/a)**p*(m + 2)) + 3*e**2*x**3*(g*x)**m*(a + c*x**2)**p*appellf1(m/2 + sympy.S(3)/2, 3, -p, m/2 + sympy.S(5)/2, e**2*x**2/d**2, -c*x**2/a)/(d**5*(1 + c*x**2/a)**p*(m + 3)) - e**3*x**4*(g*x)**m*(a + c*x**2)**p*appellf1(m/2 + 2, 3, -p, m/2 + 3, e**2*x**2/d**2, -c*x**2/a)/(d**6*(1 + c*x**2/a)**p*(m + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_433():
    f = x**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(d + e*x)
    F = x**2*(a/(24*c*d) - 7*d/(24*e**2))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)) + x**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*e) - sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))*(-15*a**3*e**6 - 17*a**2*c*d**2*e**4 - 25*a*c**2*d**4*e**2 + 105*c**3*d**6 - 2*c*d*e*x*(-5*a**2*e**4 - 6*a*c*d**2*e**2 + 35*c**2*d**4))/(192*c**3*d**3*e**4) + (-a*e**2 + c*d**2)*(5*a**3*e**6 + 9*a**2*c*d**2*e**4 + 15*a*c**2*d**4*e**2 + 35*c**3*d**6)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(128*c**(sympy.S(7)/2)*d**(sympy.S(7)/2)*e**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_434():
    f = x**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(d + e*x)
    F = x**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3*e) + (-2*c*d*e*x*(-a*e**2 + 5*c*d**2) + (-3*a*e**2 + 5*c*d**2)*(a*e**2 + 3*c*d**2))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(24*c**2*d**2*e**3) - (-a*e**2 + c*d**2)*(a**2*e**4 + 2*a*c*d**2*e**2 + 5*c**2*d**4)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(16*c**(sympy.S(5)/2)*d**(sympy.S(5)/2)*e**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_435():
    f = x*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(d + e*x)
    F = (-a/(4*c*d) - 3*d/(4*e**2))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)) + (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(2*c*d*e*(d + e*x)) + (-a*e**2 + c*d**2)*(a*e**2 + 3*c*d**2)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(8*c**(sympy.S(3)/2)*d**(sympy.S(3)/2)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_436():
    f = sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(d + e*x)
    F = sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/e - (-a*e**2 + c*d**2)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(2*sqrt(c)*sqrt(d)*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_437():
    f = sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(x*(d + e*x))
    F = -sqrt(a)*sqrt(e)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/sqrt(d) + sqrt(c)*sqrt(d)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/sqrt(e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_438():
    f = sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(x**2*(d + e*x))
    F = -sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(d*x) - (-a*e**2 + c*d**2)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(2*sqrt(a)*d**(sympy.S(3)/2)*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_439():
    f = sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(x**3*(d + e*x))
    F = -(-3*e/d**2 + c/(a*e))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*x) - sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(2*d*x**2) + (-a*e**2 + c*d**2)*(3*a*e**2 + c*d**2)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(8*a**(sympy.S(3)/2)*d**(sympy.S(5)/2)*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_440():
    f = sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(x**4*(d + e*x))
    F = -(-5*e/d**2 + c/(a*e))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(12*x**2) - sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3*d*x**3) + (-5*a*e**2 + 3*c*d**2)*(3*a*e**2 + c*d**2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(24*a**2*d**3*e**2*x) - (-a*e**2 + c*d**2)*(5*a**2*e**4 + 2*a*c*d**2*e**2 + c**2*d**4)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(16*a**(sympy.S(5)/2)*d**(sympy.S(7)/2)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_441():
    f = sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(x**5*(d + e*x))
    F = -(-7*e/d**2 + c/(a*e))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(24*x**3) - sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*d*x**4) + (-35*a**2*e**4 + 6*a*c*d**2*e**2 + 5*c**2*d**4)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(96*a**2*d**3*e**2*x**2) - sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))*(-105*a**3*e**6 + 25*a**2*c*d**2*e**4 + 17*a*c**2*d**4*e**2 + 15*c**3*d**6)/(192*a**3*d**4*e**3*x) + (-a*e**2 + c*d**2)*(35*a**3*e**6 + 15*a**2*c*d**2*e**4 + 9*a*c**2*d**4*e**2 + 5*c**3*d**6)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(128*a**(sympy.S(7)/2)*d**(sympy.S(9)/2)*e**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_442():
    f = x**3*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(d + e*x)
    F = x**2*(a/(20*c*d) - 3*d/(20*e**2))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2) + x**3*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(6*e) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)*(-35*a**3*e**6 - 33*a**2*c*d**2*e**4 - 21*a*c**2*d**4*e**2 + 105*c**3*d**6 - 6*c*d*e*x*(-7*a**2*e**4 - 6*a*c*d**2*e**2 + 21*c**2*d**4))/(960*c**3*d**3*e**4) + (a*e**2 + c*d**2 + 2*c*d*e*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))*(-7*a**4*e**8 - 8*a**3*c*d**2*e**6 - 6*a**2*c**2*d**4*e**4 + 21*c**4*d**8)/(512*c**4*d**4*e**5) - (-a*e**2 + c*d**2)**3*(7*a**3*e**6 + 15*a**2*c*d**2*e**4 + 21*a*c**2*d**4*e**2 + 21*c**3*d**6)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(1024*c**(sympy.S(9)/2)*d**(sympy.S(9)/2)*e**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_443():
    f = x**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(d + e*x)
    F = x**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(5*e) + (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)*(-15*a**2*e**4 - 12*a*c*d**2*e**2 + 35*c**2*d**4 - 6*c*d*e*x*(-3*a*e**2 + 7*c*d**2))/(240*c**2*d**2*e**3) - (-a*e**2 + c*d**2)*(a*e**2 + c*d**2 + 2*c*d*e*x)*(3*a**2*e**4 + 6*a*c*d**2*e**2 + 7*c**2*d**4)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(128*c**3*d**3*e**4) + (-a*e**2 + c*d**2)**3*(3*a**2*e**4 + 6*a*c*d**2*e**2 + 7*c**2*d**4)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(256*c**(sympy.S(7)/2)*d**(sympy.S(7)/2)*e**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_444():
    f = x*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(d + e*x)
    F = -(a/(8*c*d) + 5*d/(24*e**2))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2) + (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(4*c*d*e*(d + e*x)) + (-a*e**2 + c*d**2)*(3*a*e**2 + 5*c*d**2)*(a*e**2 + c*d**2 + 2*c*d*e*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(64*c**2*d**2*e**3) - (-a*e**2 + c*d**2)**3*(3*a*e**2 + 5*c*d**2)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(128*c**(sympy.S(5)/2)*d**(sympy.S(5)/2)*e**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_445():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(d + e*x)
    F = (a/(8*c*d) - d/(8*e**2))*(a*e**2 + c*d**2 + 2*c*d*e*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)) + (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(3*e) + (-a*e**2 + c*d**2)**3*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(16*c**(sympy.S(3)/2)*d**(sympy.S(3)/2)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_446():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(x*(d + e*x))
    F = -a**(sympy.S(3)/2)*sqrt(d)*e**(sympy.S(3)/2)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))) + (5*a*e**2 + c*d**2 + 2*c*d*e*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*e) - (-3*a**2*e**4 - 6*a*c*d**2*e**2 + c**2*d**4)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(8*sqrt(c)*sqrt(d)*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_447():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(x**2*(d + e*x))
    F = -sqrt(a)*sqrt(e)*(a*e**2 + 3*c*d**2)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(2*sqrt(d)) + sqrt(c)*sqrt(d)*(3*a*e**2 + c*d**2)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(2*sqrt(e)) - (a*e - c*d*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_448():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(x**3*(d + e*x))
    F = c**(sympy.S(3)/2)*d**(sympy.S(3)/2)*sqrt(e)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))) - (2*a*d*e + x*(a*e**2 + 5*c*d**2))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*d*x**2) - (-a**2*e**4 + 6*a*c*d**2*e**2 + 3*c**2*d**4)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(8*sqrt(a)*d**(sympy.S(3)/2)*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_449():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(x**4*(d + e*x))
    F = -(-e/d**2 + c/(a*e))*(2*a*d*e + x*(a*e**2 + c*d**2))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(8*x**2) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(3*d*x**3) + (-a*e**2 + c*d**2)**3*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(16*a**(sympy.S(3)/2)*d**(sympy.S(5)/2)*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_450():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(x**5*(d + e*x))
    F = -(-5*e/d**2 + 3*c/(a*e))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(24*x**3) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(4*d*x**4) + (-a*e**2 + c*d**2)*(5*a*e**2 + 3*c*d**2)*(2*a*d*e + x*(a*e**2 + c*d**2))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(64*a**2*d**3*e**2*x**2) - (-a*e**2 + c*d**2)**3*(5*a*e**2 + 3*c*d**2)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(128*a**(sympy.S(5)/2)*d**(sympy.S(7)/2)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_451():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(x**6*(d + e*x))
    F = -(-7*e/d**2 + 3*c/(a*e))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(40*x**4) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(5*d*x**5) + (-35*a**2*e**4 + 12*a*c*d**2*e**2 + 15*c**2*d**4)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(240*a**2*d**3*e**2*x**3) - (-a*e**2 + c*d**2)*(2*a*d*e + x*(a*e**2 + c*d**2))*(7*a**2*e**4 + 6*a*c*d**2*e**2 + 3*c**2*d**4)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(128*a**3*d**4*e**3*x**2) + (-a*e**2 + c*d**2)**3*(7*a**2*e**4 + 6*a*c*d**2*e**2 + 3*c**2*d**4)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(256*a**(sympy.S(7)/2)*d**(sympy.S(9)/2)*e**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_452():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(x**7*(d + e*x))
    F = -(-3*e/d**2 + c/(a*e))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(20*x**5) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(6*d*x**6) + (-21*a**2*e**4 + 6*a*c*d**2*e**2 + 7*c**2*d**4)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(160*a**2*d**3*e**2*x**4) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)*(-105*a**3*e**6 + 21*a**2*c*d**2*e**4 + 33*a*c**2*d**4*e**2 + 35*c**3*d**6)/(960*a**3*d**4*e**3*x**3) + (2*a*d*e + x*(a*e**2 + c*d**2))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))*(-21*a**4*e**8 + 6*a**2*c**2*d**4*e**4 + 8*a*c**3*d**6*e**2 + 7*c**4*d**8)/(512*a**4*d**5*e**4*x**2) - (-a*e**2 + c*d**2)**3*(21*a**3*e**6 + 21*a**2*c*d**2*e**4 + 15*a*c**2*d**4*e**2 + 7*c**3*d**6)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(1024*a**(sympy.S(9)/2)*d**(sympy.S(11)/2)*e**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_453():
    f = x**3*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(d + e*x)
    F = x**2*(5*a/(112*c*d) - 11*d/(112*e**2))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2) + x**3*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(8*e) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)*(-105*a**3*e**6 - 95*a**2*c*d**2*e**4 - 15*a*c**2*d**4*e**2 + 231*c**3*d**6 - 10*c*d*e*x*(-15*a**2*e**4 - 10*a*c*d**2*e**2 + 33*c**2*d**4))/(4480*c**3*d**3*e**4) + (-a*e**2 + c*d**2)*(a*e**2 + c*d**2 + 2*c*d*e*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)*(15*a**3*e**6 + 35*a**2*c*d**2*e**4 + 45*a*c**2*d**4*e**2 + 33*c**3*d**6)/(2048*c**4*d**4*e**5) - 3*(-a*e**2 + c*d**2)**3*(a*e**2 + c*d**2 + 2*c*d*e*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))*(15*a**3*e**6 + 35*a**2*c*d**2*e**4 + 45*a*c**2*d**4*e**2 + 33*c**3*d**6)/(16384*c**5*d**5*e**6) + 3*(-a*e**2 + c*d**2)**5*(15*a**3*e**6 + 35*a**2*c*d**2*e**4 + 45*a*c**2*d**4*e**2 + 33*c**3*d**6)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(32768*c**(sympy.S(11)/2)*d**(sympy.S(11)/2)*e**(sympy.S(13)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_454():
    f = x**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(d + e*x)
    F = x**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(7*e) + (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)*(-35*a**2*e**4 - 20*a*c*d**2*e**2 + 63*c**2*d**4 - 10*c*d*e*x*(-5*a*e**2 + 9*c*d**2))/(840*c**2*d**2*e**3) - (-a*e**2 + c*d**2)*(a*e**2 + c*d**2 + 2*c*d*e*x)*(5*a**2*e**4 + 10*a*c*d**2*e**2 + 9*c**2*d**4)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(384*c**3*d**3*e**4) + (-a*e**2 + c*d**2)**3*(a*e**2 + c*d**2 + 2*c*d*e*x)*(5*a**2*e**4 + 10*a*c*d**2*e**2 + 9*c**2*d**4)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(1024*c**4*d**4*e**5) - (-a*e**2 + c*d**2)**5*(5*a**2*e**4 + 10*a*c*d**2*e**2 + 9*c**2*d**4)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(2048*c**(sympy.S(9)/2)*d**(sympy.S(9)/2)*e**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_455():
    f = x*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(d + e*x)
    F = -(a/(12*c*d) + 7*d/(60*e**2))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2) + (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(6*c*d*e*(d + e*x)) + (-a*e**2 + c*d**2)*(5*a*e**2 + 7*c*d**2)*(a*e**2 + c*d**2 + 2*c*d*e*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(192*c**2*d**2*e**3) - (-a*e**2 + c*d**2)**3*(5*a*e**2 + 7*c*d**2)*(a*e**2 + c*d**2 + 2*c*d*e*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(512*c**3*d**3*e**4) + (-a*e**2 + c*d**2)**5*(5*a*e**2 + 7*c*d**2)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(1024*c**(sympy.S(7)/2)*d**(sympy.S(7)/2)*e**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_456():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(d + e*x)
    F = (a/(16*c*d) - d/(16*e**2))*(a*e**2 + c*d**2 + 2*c*d*e*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2) + (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(5*e) + 3*(-a*e**2 + c*d**2)**3*(a*e**2 + c*d**2 + 2*c*d*e*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(128*c**2*d**2*e**3) - 3*(-a*e**2 + c*d**2)**5*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(256*c**(sympy.S(5)/2)*d**(sympy.S(5)/2)*e**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_457():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(x*(d + e*x))
    F = -a**(sympy.S(5)/2)*d**(sympy.S(3)/2)*e**(sympy.S(5)/2)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))) + (11*a*e**2 + 3*c*d**2 + 6*c*d*e*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(24*e) - sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))*(-5*a**3*e**6 - 83*a**2*c*d**2*e**4 - 11*a*c**2*d**4*e**2 + 3*c**3*d**6 + 2*c*d*e*x*(-5*a*e**2 + c*d**2)*(a*e**2 + 3*c*d**2))/(64*c*d*e**2) + (-5*a**4*e**8 + 60*a**3*c*d**2*e**6 + 90*a**2*c**2*d**4*e**4 - 20*a*c**3*d**6*e**2 + 3*c**4*d**8)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(128*c**(sympy.S(3)/2)*d**(sympy.S(3)/2)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_458():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(x**2*(d + e*x))
    F = -a**(sympy.S(3)/2)*sqrt(d)*e**(sympy.S(3)/2)*(3*a*e**2 + 5*c*d**2)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/2 - (3*a*e - c*d*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(3*x) + sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))*(19*a**2*e**4 + 28*a*c*d**2*e**2 + c**2*d**4 + 2*c*d*e*x*(7*a*e**2 + c*d**2))/(8*e) - (-5*a**3*e**6 - 45*a**2*c*d**2*e**4 - 15*a*c**2*d**4*e**2 + c**3*d**6)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(16*sqrt(c)*sqrt(d)*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_459():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(x**3*(d + e*x))
    F = -3*sqrt(a)*sqrt(e)*(a**2*e**4 + 10*a*c*d**2*e**2 + 5*c**2*d**4)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(8*sqrt(d)) + 3*sqrt(c)*sqrt(d)*(5*a**2*e**4 + 10*a*c*d**2*e**2 + c**2*d**4)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(8*sqrt(e)) - (3*a*e*(a*e**2 + 3*c*d**2) - 3*c*d*x*(3*a*e**2 + c*d**2))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*x) - (a*e - c*d*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_460():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(x**4*(d + e*x))
    F = c**(sympy.S(3)/2)*d**(sympy.S(3)/2)*sqrt(e)*(5*a*e**2 + 3*c*d**2)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/2 - sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))*(-a**2*e**4 + 12*a*c*d**2*e**2 + 5*c**2*d**4 - 2*c*d*e*x*(a*e**2 + 7*c*d**2))/(8*d*x) - (4*a*d*e + x*(3*a*e**2 + 9*c*d**2))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(12*d*x**3) - (-a**3*e**6 + 15*a**2*c*d**2*e**4 + 45*a*c**2*d**4*e**2 + 5*c**3*d**6)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(16*sqrt(a)*d**(sympy.S(3)/2)*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_461():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(x**5*(d + e*x))
    F = c**(sympy.S(5)/2)*d**(sympy.S(5)/2)*e**(sympy.S(3)/2)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))) - (6*a*d*e + x*(3*a*e**2 + 11*c*d**2))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(24*d*x**4) - (2*a*d*e*(-a*e**2 + 5*c*d**2)*(3*a*e**2 + c*d**2) + x*(-3*a**3*e**6 + 11*a**2*c*d**2*e**4 + 83*a*c**2*d**4*e**2 + 5*c**3*d**6))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(64*a*d**2*e*x**2) + (-3*a**4*e**8 + 20*a**3*c*d**2*e**6 - 90*a**2*c**2*d**4*e**4 - 60*a*c**3*d**6*e**2 + 5*c**4*d**8)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(128*a**(sympy.S(3)/2)*d**(sympy.S(5)/2)*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_462():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(x**6*(d + e*x))
    F = -(-e/d**2 + c/(a*e))*(2*a*d*e + x*(a*e**2 + c*d**2))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(16*x**4) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(5*d*x**5) + 3*(-a*e**2 + c*d**2)**3*(2*a*d*e + x*(a*e**2 + c*d**2))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(128*a**2*d**3*e**2*x**2) - 3*(-a*e**2 + c*d**2)**5*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(256*a**(sympy.S(5)/2)*d**(sympy.S(7)/2)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_463():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(x**7*(d + e*x))
    F = -(-7*e/d**2 + 5*c/(a*e))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(60*x**5) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(6*d*x**6) + (-a*e**2 + c*d**2)*(7*a*e**2 + 5*c*d**2)*(2*a*d*e + x*(a*e**2 + c*d**2))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(192*a**2*d**3*e**2*x**4) - (-a*e**2 + c*d**2)**3*(7*a*e**2 + 5*c*d**2)*(2*a*d*e + x*(a*e**2 + c*d**2))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(512*a**3*d**4*e**3*x**2) + (-a*e**2 + c*d**2)**5*(7*a*e**2 + 5*c*d**2)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(1024*a**(sympy.S(7)/2)*d**(sympy.S(9)/2)*e**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_464():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(x**8*(d + e*x))
    F = -(-9*e/d**2 + 5*c/(a*e))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(84*x**6) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(7*d*x**7) + (-63*a**2*e**4 + 20*a*c*d**2*e**2 + 35*c**2*d**4)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(840*a**2*d**3*e**2*x**5) - (-a*e**2 + c*d**2)*(2*a*d*e + x*(a*e**2 + c*d**2))*(9*a**2*e**4 + 10*a*c*d**2*e**2 + 5*c**2*d**4)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(384*a**3*d**4*e**3*x**4) + (-a*e**2 + c*d**2)**3*(2*a*d*e + x*(a*e**2 + c*d**2))*(9*a**2*e**4 + 10*a*c*d**2*e**2 + 5*c**2*d**4)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(1024*a**4*d**5*e**4*x**2) - (-a*e**2 + c*d**2)**5*(9*a**2*e**4 + 10*a*c*d**2*e**2 + 5*c**2*d**4)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(2048*a**(sympy.S(9)/2)*d**(sympy.S(11)/2)*e**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_465():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(x**9*(d + e*x))
    F = -(-11*e/d**2 + 5*c/(a*e))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(112*x**7) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(8*d*x**8) + (-33*a**2*e**4 + 10*a*c*d**2*e**2 + 15*c**2*d**4)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(448*a**2*d**3*e**2*x**6) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)*(-231*a**3*e**6 + 15*a**2*c*d**2*e**4 + 95*a*c**2*d**4*e**2 + 105*c**3*d**6)/(4480*a**3*d**4*e**3*x**5) + (-a*e**2 + c*d**2)*(2*a*d*e + x*(a*e**2 + c*d**2))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)*(33*a**3*e**6 + 45*a**2*c*d**2*e**4 + 35*a*c**2*d**4*e**2 + 15*c**3*d**6)/(2048*a**4*d**5*e**4*x**4) - 3*(-a*e**2 + c*d**2)**3*(2*a*d*e + x*(a*e**2 + c*d**2))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))*(33*a**3*e**6 + 45*a**2*c*d**2*e**4 + 35*a*c**2*d**4*e**2 + 15*c**3*d**6)/(16384*a**5*d**6*e**5*x**2) + 3*(-a*e**2 + c*d**2)**5*(33*a**3*e**6 + 45*a**2*c*d**2*e**4 + 35*a*c**2*d**4*e**2 + 15*c**3*d**6)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(32768*a**(sympy.S(11)/2)*d**(sympy.S(13)/2)*e**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_466():
    f = x**2/((d + e*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    F = 2*d**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(e**2*(d + e*x)*(-a*e**2 + c*d**2)) + sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(c*d*e**2) - (a*e**2 + 3*c*d**2)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(2*c**(sympy.S(3)/2)*d**(sympy.S(3)/2)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_467():
    f = x/((d + e*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    F = -2*d*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(e*(d + e*x)*(-a*e**2 + c*d**2)) + atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(sqrt(c)*sqrt(d)*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_468():
    f = 1/((d + e*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    F = 2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/((d + e*x)*(-a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_469():
    f = 1/(x*(d + e*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    F = -2*e*(a*e + c*d*x)/(d*(-a*e**2 + c*d**2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) - atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(sqrt(a)*d**(sympy.S(3)/2)*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_470():
    f = 1/(x**2*(d + e*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    F = -2*e*(a*e + c*d*x)/(d*x*(-a*e**2 + c*d**2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) - (-3*a*e**2 + c*d**2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(a*d**2*e*x*(-a*e**2 + c*d**2)) + (3*a*e**2 + c*d**2)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(2*a**(sympy.S(3)/2)*d**(sympy.S(5)/2)*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_471():
    f = 1/(x**3*(d + e*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    F = -2*e*(a*e + c*d*x)/(d*x**2*(-a*e**2 + c*d**2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) - (-5*a*e**2 + c*d**2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(2*a*d**2*e*x**2*(-a*e**2 + c*d**2)) + (-5*a*e**2 + 3*c*d**2)*(3*a*e**2 + c*d**2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*a**2*d**3*e**2*x*(-a*e**2 + c*d**2)) - (15*a**2*e**4 + 6*a*c*d**2*e**2 + 3*c**2*d**4)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(8*a**(sympy.S(5)/2)*d**(sympy.S(7)/2)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_472():
    f = x**5/((d + e*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    F = -2*d*x**4*(a*e*(-a*e**2 + c*d**2) + c*d*x*(-a*e**2 + c*d**2))/(3*e*(-a*e**2 + c*d**2)**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)) - 2*x**2*(a*d*e*(-a*e**2 + c*d**2)*(-3*a**2*e**4 - 12*a*c*d**2*e**2 + 7*c**2*d**4) + x*(-a*e**2 + c*d**2)*(-3*a**3*e**6 - a**2*c*d**2*e**4 - 11*a*c**2*d**4*e**2 + 7*c**3*d**6))/(3*c*d*e**2*(-a*e**2 + c*d**2)**4*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) - sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))*(-45*a**4*e**8 + 30*a**3*c*d**2*e**6 + 36*a**2*c**2*d**4*e**4 - 190*a*c**3*d**6*e**2 + 105*c**4*d**8 - 2*c*d*e*x*(-15*a**3*e**6 + 9*a**2*c*d**2*e**4 - 61*a*c**2*d**4*e**2 + 35*c**3*d**6))/(12*c**3*d**3*e**4*(-a*e**2 + c*d**2)**3) + (15*a**2*e**4 + 30*a*c*d**2*e**2 + 35*c**2*d**4)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(8*c**(sympy.S(7)/2)*d**(sympy.S(7)/2)*e**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_473():
    f = x**4/((d + e*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    F = -2*d*x**3*(a*e*(-a*e**2 + c*d**2) + c*d*x*(-a*e**2 + c*d**2))/(3*e*(-a*e**2 + c*d**2)**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)) - 2*x*(a*d*e*(-a*e**2 + c*d**2)*(-3*a**2*e**4 - 10*a*c*d**2*e**2 + 5*c**2*d**4) + x*(-a*e**2 + c*d**2)*(-3*a**3*e**6 - a**2*c*d**2*e**4 - 9*a*c**2*d**4*e**2 + 5*c**3*d**6))/(3*c*d*e**2*(-a*e**2 + c*d**2)**4*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) + sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))*(-9*a**3*e**6 + 9*a**2*c*d**2*e**4 - 31*a*c**2*d**4*e**2 + 15*c**3*d**6)/(3*c**2*d**2*e**3*(-a*e**2 + c*d**2)**3) - (3*a*e**2 + 5*c*d**2)*atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(2*c**(sympy.S(5)/2)*d**(sympy.S(5)/2)*e**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_474():
    f = x**3/((d + e*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    F = -2*d*x**2*(a*e*(-a*e**2 + c*d**2) + c*d*x*(-a*e**2 + c*d**2))/(3*e*(-a*e**2 + c*d**2)**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)) - (2*a*d*e*(-3*a*e**2 + c*d**2)*(a*e**2 + 3*c*d**2) + 2*x*(-3*a**3*e**6 - a**2*c*d**2*e**4 - 7*a*c**2*d**4*e**2 + 3*c**3*d**6))/(3*c*d*e**2*(-a*e**2 + c*d**2)**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) + atanh((a*e**2 + c*d**2 + 2*c*d*e*x)/(2*sqrt(c)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(c**(sympy.S(3)/2)*d**(sympy.S(3)/2)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_475():
    f = x**2/((d + e*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    F = -8*a*e*(2*a*d*e + x*(a*e**2 + c*d**2))/(3*(-a*e**2 + c*d**2)**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) + 2*x**2/((d + e*x)*(-3*a*e**2 + 3*c*d**2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_476():
    f = x/((d + e*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    F = -2*d/(3*e*(d + e*x)*(-a*e**2 + c*d**2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) + (6*a*e**2 + 2*c*d**2)*(a*e**2 + c*d**2 + 2*c*d*e*x)/(3*e*(-a*e**2 + c*d**2)**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_477():
    f = 1/((d + e*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    F = -8*c*d*(a*e**2 + c*d**2 + 2*c*d*e*x)/(3*(-a*e**2 + c*d**2)**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) + 2/((d + e*x)*(-3*a*e**2 + 3*c*d**2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_478():
    f = 1/(x*(d + e*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    F = -2*e*(a*e + c*d*x)/(3*d*(-a*e**2 + c*d**2)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)) + (-6*a**3*e**6 + 14*a**2*c*d**2*e**4 + 2*a*c**2*d**4*e**2 + 6*c**3*d**6 + 2*c*d*e*x*(-a*e**2 + 3*c*d**2)*(3*a*e**2 + c*d**2))/(3*a*d**2*e*(-a*e**2 + c*d**2)**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) - atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(a**(sympy.S(3)/2)*d**(sympy.S(5)/2)*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_479():
    f = 1/(x**2*(d + e*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    F = -2*e*(a*e + c*d*x)/(3*d*x*(-a*e**2 + c*d**2)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)) + (-10*a**3*e**6 + 18*a**2*c*d**2*e**4 + 2*a*c**2*d**4*e**2 + 6*c**3*d**6 + 2*c*d*e*x*(-5*a**2*e**4 + 10*a*c*d**2*e**2 + 3*c**2*d**4))/(3*a*d**2*e*x*(-a*e**2 + c*d**2)**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) - sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))*(-15*a**3*e**6 + 31*a**2*c*d**2*e**4 - 9*a*c**2*d**4*e**2 + 9*c**3*d**6)/(3*a**2*d**3*e**2*x*(-a*e**2 + c*d**2)**3) + (5*a*e**2 + 3*c*d**2)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(2*a**(sympy.S(5)/2)*d**(sympy.S(7)/2)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_480():
    f = 1/(x**3*(d + e*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    F = -2*e*(a*e + c*d*x)/(3*d*x**2*(-a*e**2 + c*d**2)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)) + (-14*a**3*e**6 + 22*a**2*c*d**2*e**4 + 2*a*c**2*d**4*e**2 + 6*c**3*d**6 + 2*c*d*e*x*(-7*a**2*e**4 + 12*a*c*d**2*e**2 + 3*c**2*d**4))/(3*a*d**2*e*x**2*(-a*e**2 + c*d**2)**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) - sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))*(-35*a**3*e**6 + 61*a**2*c*d**2*e**4 - 9*a*c**2*d**4*e**2 + 15*c**3*d**6)/(6*a**2*d**3*e**2*x**2*(-a*e**2 + c*d**2)**3) + sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))*(-105*a**4*e**8 + 190*a**3*c*d**2*e**6 - 36*a**2*c**2*d**4*e**4 - 30*a*c**3*d**6*e**2 + 45*c**4*d**8)/(12*a**3*d**4*e**3*x*(-a*e**2 + c*d**2)**3) - (35*a**2*e**4 + 30*a*c*d**2*e**2 + 15*c**2*d**4)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(8*a**(sympy.S(7)/2)*d**(sympy.S(9)/2)*e**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_481():
    f = 1/(x**4*(d + e*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    F = -2*e*(a*e + c*d*x)/(3*d*x**3*(-a*e**2 + c*d**2)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)) + (-18*a**3*e**6 + 26*a**2*c*d**2*e**4 + 2*a*c**2*d**4*e**2 + 6*c**3*d**6 + 2*c*d*e*x*(-9*a**2*e**4 + 14*a*c*d**2*e**2 + 3*c**2*d**4))/(3*a*d**2*e*x**3*(-a*e**2 + c*d**2)**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) - sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))*(-21*a**3*e**6 + 33*a**2*c*d**2*e**4 - 3*a*c**2*d**4*e**2 + 7*c**3*d**6)/(3*a**2*d**3*e**2*x**3*(-a*e**2 + c*d**2)**3) + sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))*(-105*a**4*e**8 + 168*a**3*c*d**2*e**6 - 18*a**2*c**2*d**4*e**4 - 16*a*c**3*d**6*e**2 + 35*c**4*d**8)/(12*a**3*d**4*e**3*x**2*(-a*e**2 + c*d**2)**3) - sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))*(-315*a**5*e**10 + 525*a**4*c*d**2*e**8 - 78*a**3*c**2*d**4*e**6 - 54*a**2*c**3*d**6*e**4 - 55*a*c**4*d**8*e**2 + 105*c**5*d**10)/(24*a**4*d**5*e**4*x*(-a*e**2 + c*d**2)**3) + (105*a**3*e**6 + 105*a**2*c*d**2*e**4 + 75*a*c**2*d**4*e**2 + 35*c**3*d**6)*atanh((2*a*d*e + x*(a*e**2 + c*d**2))/(2*sqrt(a)*sqrt(d)*sqrt(e)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))))/(16*a**(sympy.S(9)/2)*d**(sympy.S(11)/2)*e**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_482():
    f = x**2/((d + e*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2))
    F = 2*x**2/((d + e*x)*(-5*a*e**2 + 5*c*d**2)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)) - (8*a*d*e*(-a*e**2 + c*d**2)*(3*a*e**2 + c*d**2) + 8*x*(-2*a**3*e**6 + a**2*c*d**2*e**4 + c**3*d**6))/(15*e*(-a*e**2 + c*d**2)**4*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)) + (a*e**2 + c*d**2 + 2*c*d*e*x)*(40*a**2*e**4 + 80*a*c*d**2*e**2 + 8*c**2*d**4)/(15*e*(-a*e**2 + c*d**2)**5*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_483():
    f = x**2/((d + e*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2))
    F = -128*c*d*(a*e**2 + c*d**2 + 2*c*d*e*x)*(7*a**2*e**4 + 14*a*c*d**2*e**2 + 3*c**2*d**4)/(105*(-a*e**2 + c*d**2)**7*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) + 2*x**2/((d + e*x)*(-7*a*e**2 + 7*c*d**2)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)) - (16*a*d*e*(2*a*e**2 + c*d**2) + 8*x*(3*a**2*e**4 + a*c*d**2*e**2 + 2*c**2*d**4))/(35*e*(-a*e**2 + c*d**2)**3*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)) + (a*e**2 + c*d**2 + 2*c*d*e*x)*(112*a**2*e**4 + 224*a*c*d**2*e**2 + 48*c**2*d**4)/(105*e*(-a*e**2 + c*d**2)**5*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_484():
    f = x**3*sqrt(x + 1)*sqrt(x**2 - x + 1)
    F = 2*x**4*sqrt(x + 1)*sqrt(x**2 - x + 1)/11 + 6*x*sqrt(x + 1)*sqrt(x**2 - x + 1)/55 - 4*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*(x + 1)**(sympy.S(3)/2)*sqrt(x**2 - x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(55*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_485():
    f = x**2*sqrt(x + 1)*sqrt(x**2 - x + 1)
    F = 2*(x + 1)**(sympy.S(3)/2)*(x**2 - x + 1)**(sympy.S(3)/2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_486():
    f = x*sqrt(x + 1)*sqrt(x**2 - x + 1)
    F = 2*x**2*sqrt(x + 1)*sqrt(x**2 - x + 1)/7 + 6*sqrt(x + 1)*sqrt(x**2 - x + 1)/(7*x + 7 + 7*sqrt(3)) - 3*3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(2 - sqrt(3))*(x + 1)**(sympy.S(3)/2)*sqrt(x**2 - x + 1)*elliptic_e(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(7*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*(x**3 + 1)) + 2*sqrt(2)*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*(x + 1)**(sympy.S(3)/2)*sqrt(x**2 - x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(7*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_487():
    f = sqrt(x + 1)*sqrt(x**2 - x + 1)
    F = 2*x*sqrt(x + 1)*sqrt(x**2 - x + 1)/5 + 2*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*(x + 1)**(sympy.S(3)/2)*sqrt(x**2 - x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(5*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_488():
    f = sqrt(x + 1)*sqrt(x**2 - x + 1)/x
    F = 2*sqrt(x + 1)*sqrt(x**2 - x + 1)/3 - 2*sqrt(x + 1)*sqrt(x**2 - x + 1)*atanh(sqrt(x**3 + 1))/(3*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_489():
    f = sqrt(x + 1)*sqrt(x**2 - x + 1)/x**2
    F = 3*sqrt(x + 1)*sqrt(x**2 - x + 1)/(x + 1 + sqrt(3)) - 3*3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(2 - sqrt(3))*(x + 1)**(sympy.S(3)/2)*sqrt(x**2 - x + 1)*elliptic_e(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(2*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*(x**3 + 1)) + sqrt(2)*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*(x + 1)**(sympy.S(3)/2)*sqrt(x**2 - x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(sqrt((x + 1)/(x + 1 + sqrt(3))**2)*(x**3 + 1)) - sqrt(x + 1)*sqrt(x**2 - x + 1)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_490():
    f = sqrt(x + 1)*sqrt(x**2 - x + 1)/x**3
    F = 3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*(x + 1)**(sympy.S(3)/2)*sqrt(x**2 - x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(2*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*(x**3 + 1)) - sqrt(x + 1)*sqrt(x**2 - x + 1)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_491():
    f = x**3*(x + 1)**(sympy.S(3)/2)*(x**2 - x + 1)**(sympy.S(3)/2)
    F = 2*x**4*sqrt(x + 1)*(x**3 + 1)*sqrt(x**2 - x + 1)/17 + 18*x**4*sqrt(x + 1)*sqrt(x**2 - x + 1)/187 + 54*x*sqrt(x + 1)*sqrt(x**2 - x + 1)/935 - 36*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*(x + 1)**(sympy.S(3)/2)*sqrt(x**2 - x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(935*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_492():
    f = x**2*(x + 1)**(sympy.S(3)/2)*(x**2 - x + 1)**(sympy.S(3)/2)
    F = 2*(x + 1)**(sympy.S(5)/2)*(x**2 - x + 1)**(sympy.S(5)/2)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_493():
    f = x*(x + 1)**(sympy.S(3)/2)*(x**2 - x + 1)**(sympy.S(3)/2)
    F = 2*x**2*sqrt(x + 1)*(x**3 + 1)*sqrt(x**2 - x + 1)/13 + 18*x**2*sqrt(x + 1)*sqrt(x**2 - x + 1)/91 + 54*sqrt(x + 1)*sqrt(x**2 - x + 1)/(91*x + 91 + 91*sqrt(3)) - 27*3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(2 - sqrt(3))*(x + 1)**(sympy.S(3)/2)*sqrt(x**2 - x + 1)*elliptic_e(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(91*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*(x**3 + 1)) + 18*sqrt(2)*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*(x + 1)**(sympy.S(3)/2)*sqrt(x**2 - x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(91*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_494():
    f = (x + 1)**(sympy.S(3)/2)*(x**2 - x + 1)**(sympy.S(3)/2)
    F = 2*x*sqrt(x + 1)*(x**3 + 1)*sqrt(x**2 - x + 1)/11 + 18*x*sqrt(x + 1)*sqrt(x**2 - x + 1)/55 + 18*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*(x + 1)**(sympy.S(3)/2)*sqrt(x**2 - x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(55*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_495():
    f = (x + 1)**(sympy.S(3)/2)*(x**2 - x + 1)**(sympy.S(3)/2)/x
    F = 2*sqrt(x + 1)*(x**3 + 1)*sqrt(x**2 - x + 1)/9 + 2*sqrt(x + 1)*sqrt(x**2 - x + 1)/3 - 2*sqrt(x + 1)*sqrt(x**2 - x + 1)*atanh(sqrt(x**3 + 1))/(3*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_496():
    f = (x + 1)**(sympy.S(3)/2)*(x**2 - x + 1)**(sympy.S(3)/2)/x**2
    F = 9*x**2*sqrt(x + 1)*sqrt(x**2 - x + 1)/7 + 27*sqrt(x + 1)*sqrt(x**2 - x + 1)/(7*x + 7 + 7*sqrt(3)) - 27*3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(2 - sqrt(3))*(x + 1)**(sympy.S(3)/2)*sqrt(x**2 - x + 1)*elliptic_e(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(14*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*(x**3 + 1)) + 9*sqrt(2)*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*(x + 1)**(sympy.S(3)/2)*sqrt(x**2 - x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(7*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*(x**3 + 1)) - sqrt(x + 1)*(x**3 + 1)*sqrt(x**2 - x + 1)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_497():
    f = (x + 1)**(sympy.S(3)/2)*(x**2 - x + 1)**(sympy.S(3)/2)/x**3
    F = 9*x*sqrt(x + 1)*sqrt(x**2 - x + 1)/10 + 9*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*(x + 1)**(sympy.S(3)/2)*sqrt(x**2 - x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(10*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*(x**3 + 1)) - sqrt(x + 1)*(x**3 + 1)*sqrt(x**2 - x + 1)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_498():
    f = x**3/(sqrt(x + 1)*sqrt(x**2 - x + 1))
    F = 2*x*(x**3 + 1)/(5*sqrt(x + 1)*sqrt(x**2 - x + 1)) - 4*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*sqrt(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(15*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**2 - x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_499():
    f = x**2/(sqrt(x + 1)*sqrt(x**2 - x + 1))
    F = 2*sqrt(x + 1)*sqrt(x**2 - x + 1)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_500():
    f = x/(sqrt(x + 1)*sqrt(x**2 - x + 1))
    F = (2*x**3 + 2)/(sqrt(x + 1)*(x + 1 + sqrt(3))*sqrt(x**2 - x + 1)) - 3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(2 - sqrt(3))*sqrt(x + 1)*elliptic_e(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**2 - x + 1)) + 2*sqrt(2)*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**2 - x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_501():
    f = 1/(sqrt(x + 1)*sqrt(x**2 - x + 1))
    F = 2*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*sqrt(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**2 - x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_502():
    f = 1/(x*sqrt(x + 1)*sqrt(x**2 - x + 1))
    F = -2*sqrt(x**3 + 1)*atanh(sqrt(x**3 + 1))/(3*sqrt(x + 1)*sqrt(x**2 - x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_503():
    f = 1/(x**2*sqrt(x + 1)*sqrt(x**2 - x + 1))
    F = (x**3 + 1)/(sqrt(x + 1)*(x + 1 + sqrt(3))*sqrt(x**2 - x + 1)) - 3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(2 - sqrt(3))*sqrt(x + 1)*elliptic_e(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(2*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**2 - x + 1)) + sqrt(2)*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**2 - x + 1)) - (x**3 + 1)/(x*sqrt(x + 1)*sqrt(x**2 - x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_504():
    f = 1/(x**3*sqrt(x + 1)*sqrt(x**2 - x + 1))
    F = -3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*sqrt(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(6*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**2 - x + 1)) + (-x**3 - 1)/(2*x**2*sqrt(x + 1)*sqrt(x**2 - x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_505():
    f = x**3/((x + 1)**(sympy.S(3)/2)*(x**2 - x + 1)**(sympy.S(3)/2))
    F = -2*x/(3*sqrt(x + 1)*sqrt(x**2 - x + 1)) + 4*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*sqrt(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(9*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**2 - x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_506():
    f = x**2/((x + 1)**(sympy.S(3)/2)*(x**2 - x + 1)**(sympy.S(3)/2))
    F = -2/(3*sqrt(x + 1)*sqrt(x**2 - x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_507():
    f = x/((x + 1)**(sympy.S(3)/2)*(x**2 - x + 1)**(sympy.S(3)/2))
    F = 2*x**2/(3*sqrt(x + 1)*sqrt(x**2 - x + 1)) - (2*x**3 + 2)/(3*sqrt(x + 1)*(x + 1 + sqrt(3))*sqrt(x**2 - x + 1)) + 3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(2 - sqrt(3))*sqrt(x + 1)*elliptic_e(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**2 - x + 1)) - 2*sqrt(2)*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(9*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**2 - x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_508():
    f = 1/((x + 1)**(sympy.S(3)/2)*(x**2 - x + 1)**(sympy.S(3)/2))
    F = 2*x/(3*sqrt(x + 1)*sqrt(x**2 - x + 1)) + 2*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*sqrt(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(9*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**2 - x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_509():
    f = 1/(x*(x + 1)**(sympy.S(3)/2)*(x**2 - x + 1)**(sympy.S(3)/2))
    F = -2*sqrt(x**3 + 1)*atanh(sqrt(x**3 + 1))/(3*sqrt(x + 1)*sqrt(x**2 - x + 1)) + 2/(3*sqrt(x + 1)*sqrt(x**2 - x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_510():
    f = 1/(x**2*(x + 1)**(sympy.S(3)/2)*(x**2 - x + 1)**(sympy.S(3)/2))
    F = (5*x**3 + 5)/(3*sqrt(x + 1)*(x + 1 + sqrt(3))*sqrt(x**2 - x + 1)) - 5*3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(2 - sqrt(3))*sqrt(x + 1)*elliptic_e(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(6*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**2 - x + 1)) + 5*sqrt(2)*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(9*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**2 - x + 1)) - (5*x**3 + 5)/(3*x*sqrt(x + 1)*sqrt(x**2 - x + 1)) + 2/(3*x*sqrt(x + 1)*sqrt(x**2 - x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_511():
    f = 1/(x**3*(x + 1)**(sympy.S(3)/2)*(x**2 - x + 1)**(sympy.S(3)/2))
    F = -7*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*sqrt(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(18*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**2 - x + 1)) - (7*x**3 + 7)/(6*x**2*sqrt(x + 1)*sqrt(x**2 - x + 1)) + 2/(3*x**2*sqrt(x + 1)*sqrt(x**2 - x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_512():
    f = x**3/((x + 1)**(sympy.S(5)/2)*(x**2 - x + 1)**(sympy.S(5)/2))
    F = 4*x/(27*sqrt(x + 1)*sqrt(x**2 - x + 1)) - 2*x/(9*sqrt(x + 1)*(x**3 + 1)*sqrt(x**2 - x + 1)) + 4*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*sqrt(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(81*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**2 - x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_513():
    f = x**2/((x + 1)**(sympy.S(5)/2)*(x**2 - x + 1)**(sympy.S(5)/2))
    F = -2/(9*(x + 1)**(sympy.S(3)/2)*(x**2 - x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_514():
    f = x/((x + 1)**(sympy.S(5)/2)*(x**2 - x + 1)**(sympy.S(5)/2))
    F = 10*x**2/(27*sqrt(x + 1)*sqrt(x**2 - x + 1)) + 2*x**2/(9*sqrt(x + 1)*(x**3 + 1)*sqrt(x**2 - x + 1)) - (10*x**3 + 10)/(27*sqrt(x + 1)*(x + 1 + sqrt(3))*sqrt(x**2 - x + 1)) + 5*3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(2 - sqrt(3))*sqrt(x + 1)*elliptic_e(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(27*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**2 - x + 1)) - 10*sqrt(2)*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(81*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**2 - x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_515():
    f = 1/((x + 1)**(sympy.S(5)/2)*(x**2 - x + 1)**(sympy.S(5)/2))
    F = 14*x/(27*sqrt(x + 1)*sqrt(x**2 - x + 1)) + 2*x/(9*sqrt(x + 1)*(x**3 + 1)*sqrt(x**2 - x + 1)) + 14*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*sqrt(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(81*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**2 - x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_516():
    f = 1/(x*(x + 1)**(sympy.S(5)/2)*(x**2 - x + 1)**(sympy.S(5)/2))
    F = -2*sqrt(x**3 + 1)*atanh(sqrt(x**3 + 1))/(3*sqrt(x + 1)*sqrt(x**2 - x + 1)) + 2/(3*sqrt(x + 1)*sqrt(x**2 - x + 1)) + 2/(9*sqrt(x + 1)*(x**3 + 1)*sqrt(x**2 - x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_517():
    f = 1/(x**2*(x + 1)**(sympy.S(5)/2)*(x**2 - x + 1)**(sympy.S(5)/2))
    F = (55*x**3 + 55)/(27*sqrt(x + 1)*(x + 1 + sqrt(3))*sqrt(x**2 - x + 1)) - 55*3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(2 - sqrt(3))*sqrt(x + 1)*elliptic_e(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(54*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**2 - x + 1)) + 55*sqrt(2)*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(81*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**2 - x + 1)) - (55*x**3 + 55)/(27*x*sqrt(x + 1)*sqrt(x**2 - x + 1)) + 22/(27*x*sqrt(x + 1)*sqrt(x**2 - x + 1)) + 2/(9*x*sqrt(x + 1)*(x**3 + 1)*sqrt(x**2 - x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_518():
    f = 1/(x**3*(x + 1)**(sympy.S(5)/2)*(x**2 - x + 1)**(sympy.S(5)/2))
    F = -91*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*sqrt(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(162*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**2 - x + 1)) - (91*x**3 + 91)/(54*x**2*sqrt(x + 1)*sqrt(x**2 - x + 1)) + 26/(27*x**2*sqrt(x + 1)*sqrt(x**2 - x + 1)) + 2/(9*x**2*sqrt(x + 1)*(x**3 + 1)*sqrt(x**2 - x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_519():
    f = x/((x - 1)**3*(4*x**2 + 5*x + 3)**2)
    F = 11*log(1 - x)/2304 - 11*log(4*x**2 + 5*x + 3)/4608 + 6023*sqrt(23)*atan(sqrt(23)*(8*x + 5)/23)/1218816 - 97/(4416 - 4416*x) + (44*x + 39)/(276*(1 - x)**2*(4*x**2 + 5*x + 3)) - 21/(736*(1 - x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_520():
    f = x**4*sqrt(d + e*x)/(a + b*x + c*x**2)
    F = -2*b*sqrt(d + e*x)*(-2*a*c + b**2)/c**4 + 2*(d + e*x)**(sympy.S(7)/2)/(7*c*e**3) - (d + e*x)**(sympy.S(5)/2)*(2*b*e + 4*c*d)/(5*c**2*e**3) + (d + e*x)**(sympy.S(3)/2)*(2*b**2*e**2 + 2*c**2*d**2 + 2*c*e*(-a*e + b*d))/(3*c**3*e**3) + sqrt(2)*(-a**2*c**2*e + 3*a*b**2*c*e - 2*a*b*c**2*d - b**4*e + b**3*c*d + (-5*a**2*b*c**2*e + 2*a**2*c**3*d + 5*a*b**3*c*e - 4*a*b**2*c**2*d - b**5*e + b**4*c*d)/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(c**(sympy.S(9)/2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + sqrt(2)*(-a**2*c**2*e + 3*a*b**2*c*e - 2*a*b*c**2*d - b**4*e + b**3*c*d - (-5*a**2*b*c**2*e + 2*a**2*c**3*d + 5*a*b**3*c*e - 4*a*b**2*c**2*d - b**5*e + b**4*c*d)/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(c**(sympy.S(9)/2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_521():
    f = x**2*sqrt(d + e*x)/(a + b*x + c*x**2)
    F = -2*b*sqrt(d + e*x)/c**2 + 2*(d + e*x)**(sympy.S(3)/2)/(3*c*e) + sqrt(2)*(a*c*e - b**2*e + b*c*d + (3*a*b*c*e - 2*a*c**2*d - b**3*e + b**2*c*d)/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(c**(sympy.S(5)/2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + sqrt(2)*(a*c*e - b**2*e + b*c*d - (3*a*b*c*e - 2*a*c**2*d - b**3*e + b**2*c*d)/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(c**(sympy.S(5)/2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_522():
    f = x*sqrt(d + e*x)/(a + b*x + c*x**2)
    F = 2*sqrt(d + e*x)/c - sqrt(2)*(2*a*c*e - b**2*e + b*c*d + sqrt(-4*a*c + b**2)*(-b*e + c*d))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(c**(sympy.S(3)/2)*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + sqrt(2)*(2*a*c*e - b**2*e + b*c*d - sqrt(-4*a*c + b**2)*(-b*e + c*d))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(c**(sympy.S(3)/2)*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_523():
    f = sqrt(d + e*x)/(a + b*x + c*x**2)
    F = -sqrt(2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(sqrt(c)*sqrt(-4*a*c + b**2)) + sqrt(2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(sqrt(c)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_524():
    f = sqrt(d + e*x)/(x*(a + b*x + c*x**2))
    F = -sqrt(2)*sqrt(c)*(-2*a*e + b*d - d*sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(a*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + sqrt(2)*sqrt(c)*(-2*a*e + b*d + d*sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(a*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) - 2*sqrt(d)*atanh(sqrt(d + e*x)/sqrt(d))/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_525():
    f = sqrt(d + e*x)/(x**2*(a + b*x + c*x**2))
    F = -sqrt(d + e*x)/(a*x) + e*atanh(sqrt(d + e*x)/sqrt(d))/(a*sqrt(d)) + sqrt(2)*sqrt(c)*(-a*(2*c*d - e*sqrt(-4*a*c + b**2)) + b**2*d - b*(a*e + d*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(a**2*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) - sqrt(2)*sqrt(c)*(-a*b*e - 2*a*c*d + b**2*d + sqrt(-4*a*c + b**2)*(-a*e + b*d))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(a**2*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) + (-2*a*e + 2*b*d)*atanh(sqrt(d + e*x)/sqrt(d))/(a**2*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_526():
    f = sqrt(d + e*x)/(x**3*(a + b*x + c*x**2))
    F = -sqrt(d + e*x)/(2*a*x**2) + 3*e*sqrt(d + e*x)/(4*a*d*x) - 3*e**2*atanh(sqrt(d + e*x)/sqrt(d))/(4*a*d**(sympy.S(3)/2)) + sqrt(d + e*x)*(-a*e + b*d)/(a**2*d*x) - e*(-a*e + b*d)*atanh(sqrt(d + e*x)/sqrt(d))/(a**2*d**(sympy.S(3)/2)) - sqrt(2)*sqrt(c)*(-a*b*(3*c*d - e*sqrt(-4*a*c + b**2)) + a*c*(2*a*e + d*sqrt(-4*a*c + b**2)) + b**3*d - b**2*(a*e + d*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(a**3*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + sqrt(2)*sqrt(c)*(-a*b*(3*c*d + e*sqrt(-4*a*c + b**2)) - a*c*(-2*a*e + d*sqrt(-4*a*c + b**2)) + b**3*d + b**2*(-a*e + d*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(a**3*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) - (-2*a*b*e - 2*a*c*d + 2*b**2*d)*atanh(sqrt(d + e*x)/sqrt(d))/(a**3*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_527():
    f = x**4*(d + e*x)**(sympy.S(3)/2)/(a + b*x + c*x**2)
    F = -2*b*(d + e*x)**(sympy.S(3)/2)*(-2*a*c + b**2)/(3*c**4) + 2*(d + e*x)**(sympy.S(9)/2)/(9*c*e**3) - (d + e*x)**(sympy.S(7)/2)*(2*b*e + 4*c*d)/(7*c**2*e**3) + (d + e*x)**(sympy.S(5)/2)*(2*b**2*e**2 + 2*c**2*d**2 + 2*c*e*(-a*e + b*d))/(5*c**3*e**3) - sqrt(d + e*x)*(-2*a**2*c**2*e + 6*a*b**2*c*e - 4*a*b*c**2*d - 2*b**4*e + 2*b**3*c*d)/c**5 + sqrt(2)*((a*c*e - b**2*e + b*c*d)*(3*a*b*c*e - 2*a*c**2*d - b**3*e + b**2*c*d) - (10*a**2*b*c**3*d*e - 2*a**2*c**3*(-a*e**2 + c*d**2) - 10*a*b**3*c**2*d*e + a*b**2*c**2*(-9*a*e**2 + 4*c*d**2) - b**6*e**2 + 2*b**5*c*d*e - b**4*c*(-6*a*e**2 + c*d**2))/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(c**(sympy.S(11)/2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + sqrt(2)*((a*c*e - b**2*e + b*c*d)*(3*a*b*c*e - 2*a*c**2*d - b**3*e + b**2*c*d) + (10*a**2*b*c**3*d*e - 2*a**2*c**3*(-a*e**2 + c*d**2) - 10*a*b**3*c**2*d*e + a*b**2*c**2*(-9*a*e**2 + 4*c*d**2) - b**6*e**2 + 2*b**5*c*d*e - b**4*c*(-6*a*e**2 + c*d**2))/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(c**(sympy.S(11)/2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_528():
    f = x**3*(d + e*x)**(sympy.S(3)/2)/(a + b*x + c*x**2)
    F = 2*(d + e*x)**(sympy.S(7)/2)/(7*c*e**2) - (d + e*x)**(sympy.S(5)/2)*(2*b*e + 2*c*d)/(5*c**2*e**2) + (d + e*x)**(sympy.S(3)/2)*(-2*a*c + 2*b**2)/(3*c**3) + sqrt(d + e*x)*(4*a*b*c*e - 2*a*c**2*d - 2*b**3*e + 2*b**2*c*d)/c**4 + sqrt(2)*(-4*a*b*c**2*d*e + a*c**2*(-a*e**2 + c*d**2) - b**4*e**2 + 2*b**3*c*d*e - b**2*c*(-3*a*e**2 + c*d**2) + (4*a**2*c**3*d*e - 8*a*b**2*c**2*d*e + a*b*c**2*(-5*a*e**2 + 3*c*d**2) - b**5*e**2 + 2*b**4*c*d*e - b**3*c*(-5*a*e**2 + c*d**2))/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(c**(sympy.S(9)/2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + sqrt(2)*(-4*a*b*c**2*d*e + a*c**2*(-a*e**2 + c*d**2) - b**4*e**2 + 2*b**3*c*d*e - b**2*c*(-3*a*e**2 + c*d**2) - (4*a**2*c**3*d*e - 8*a*b**2*c**2*d*e + a*b*c**2*(-5*a*e**2 + 3*c*d**2) - b**5*e**2 + 2*b**4*c*d*e - b**3*c*(-5*a*e**2 + c*d**2))/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(c**(sympy.S(9)/2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_529():
    f = x**2*(d + e*x)**(sympy.S(3)/2)/(a + b*x + c*x**2)
    F = -2*b*(d + e*x)**(sympy.S(3)/2)/(3*c**2) + 2*(d + e*x)**(sympy.S(5)/2)/(5*c*e) - sqrt(d + e*x)*(2*a*c*e - 2*b**2*e + 2*b*c*d)/c**3 + sqrt(2)*((-b*e + c*d)*(2*a*c*e - b**2*e + b*c*d) - (-6*a*b*c**2*d*e + 2*a*c**2*(-a*e**2 + c*d**2) - b**4*e**2 + 2*b**3*c*d*e - b**2*c*(-4*a*e**2 + c*d**2))/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(c**(sympy.S(7)/2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + sqrt(2)*((-b*e + c*d)*(2*a*c*e - b**2*e + b*c*d) + (-6*a*b*c**2*d*e + 2*a*c**2*(-a*e**2 + c*d**2) - b**4*e**2 + 2*b**3*c*d*e - b**2*c*(-4*a*e**2 + c*d**2))/sqrt(-4*a*c + b**2))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(c**(sympy.S(7)/2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_530():
    f = x*(d + e*x)**(sympy.S(3)/2)/(a + b*x + c*x**2)
    F = 2*(d + e*x)**(sympy.S(3)/2)/(3*c) + sqrt(d + e*x)*(-2*b*e + 2*c*d)/c**2 - sqrt(2)*(b**3*e**2 - b**2*e*(2*c*d - e*sqrt(-4*a*c + b**2)) + b*c*(c*d**2 - e*(3*a*e + 2*d*sqrt(-4*a*c + b**2))) - c*(a*e**2*sqrt(-4*a*c + b**2) - c*d*(4*a*e + d*sqrt(-4*a*c + b**2))))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(c**(sympy.S(5)/2)*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + sqrt(2)*(b**3*e**2 - b**2*e*(2*c*d + e*sqrt(-4*a*c + b**2)) + b*c*(c*d**2 + e*(-3*a*e + 2*d*sqrt(-4*a*c + b**2))) + c*(a*e**2*sqrt(-4*a*c + b**2) - c*d*(-4*a*e + d*sqrt(-4*a*c + b**2))))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(c**(sympy.S(5)/2)*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_531():
    f = (d + e*x)**(sympy.S(3)/2)/(a + b*x + c*x**2)
    F = 2*e*sqrt(d + e*x)/c + sqrt(2)*(b*e**2*(b + sqrt(-4*a*c + b**2)) + 2*c**2*d**2 - 2*c*e*(a*e + b*d + d*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(c**(sympy.S(3)/2)*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) - sqrt(2)*(b*e**2*(b - sqrt(-4*a*c + b**2)) + 2*c**2*d**2 - 2*c*e*(a*e + b*d - d*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(c**(sympy.S(3)/2)*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_532():
    f = (d + e*x)**(sympy.S(3)/2)/(x*(a + b*x + c*x**2))
    F = -2*d**(sympy.S(3)/2)*atanh(sqrt(d + e*x)/sqrt(d))/a - sqrt(2)*(a*e**2*sqrt(-4*a*c + b**2) + b*(a*e**2 + c*d**2) - c*d*(4*a*e + d*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(a*sqrt(c)*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) - sqrt(2)*(a*e**2*sqrt(-4*a*c + b**2) - b*(a*e**2 + c*d**2) - c*d*(-4*a*e + d*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(a*sqrt(c)*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_533():
    f = (d + e*x)**(sympy.S(3)/2)/(x**2*(a + b*x + c*x**2))
    F = sqrt(d)*e*atanh(sqrt(d + e*x)/sqrt(d))/a - d*sqrt(d + e*x)/(a*x) + sqrt(2)*sqrt(c)*(-2*a*(c*d**2 - e*(a*e + d*sqrt(-4*a*c + b**2))) + b**2*d**2 - b*d*(2*a*e + d*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(a**2*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) - sqrt(2)*sqrt(c)*(-2*a*(c*d**2 + e*(-a*e + d*sqrt(-4*a*c + b**2))) + b**2*d**2 + b*d*(-2*a*e + d*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(a**2*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) + 2*sqrt(d)*(-2*a*e + b*d)*atanh(sqrt(d + e*x)/sqrt(d))/a**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_534():
    f = (d + e*x)**(sympy.S(3)/2)/(x**3*(a + b*x + c*x**2))
    F = -d*sqrt(d + e*x)/(2*a*x**2) + 3*e*sqrt(d + e*x)/(4*a*x) - 3*e**2*atanh(sqrt(d + e*x)/sqrt(d))/(4*a*sqrt(d)) + sqrt(d + e*x)*(-2*a*e + b*d)/(a**2*x) - e*(-2*a*e + b*d)*atanh(sqrt(d + e*x)/sqrt(d))/(a**2*sqrt(d)) - sqrt(2)*sqrt(c)*(-a*b*(3*c*d**2 - e*(a*e + 2*d*sqrt(-4*a*c + b**2))) - a*(a*e**2*sqrt(-4*a*c + b**2) - c*d*(4*a*e + d*sqrt(-4*a*c + b**2))) + b**3*d**2 - b**2*d*(2*a*e + d*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(a**3*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))) + sqrt(2)*sqrt(c)*(-a*b*(3*c*d**2 + e*(-a*e + 2*d*sqrt(-4*a*c + b**2))) + a*(a*e**2*sqrt(-4*a*c + b**2) - c*d*(-4*a*e + d*sqrt(-4*a*c + b**2))) + b**3*d**2 + b**2*d*(-2*a*e + d*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*sqrt(c)*sqrt(d + e*x)/sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2))))/(a**3*sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))) - (-4*a*b*d*e - 2*a*(-a*e**2 + c*d**2) + 2*b**2*d**2)*atanh(sqrt(d + e*x)/sqrt(d))/(a**3*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_535():
    f = x**m*(e + f*x)**n/(a + b*x + c*x**2)
    F = -2*c*x**(m + 1)*(e + f*x)**n*appellf1(m + 1, 1, -n, m + 2, -2*c*x/(b + sqrt(-4*a*c + b**2)), -f*x/e)/((1 + f*x/e)**n*(b + sqrt(-4*a*c + b**2))*(m + 1)*sqrt(-4*a*c + b**2)) + 2*c*x**(m + 1)*(e + f*x)**n*appellf1(m + 1, 1, -n, m + 2, -2*c*x/(b - sqrt(-4*a*c + b**2)), -f*x/e)/((1 + f*x/e)**n*(b - sqrt(-4*a*c + b**2))*(m + 1)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_536():
    f = x**3*(e + f*x)**n/(a + b*x + c*x**2)
    F = (e + f*x)**(n + 1)*(a - b**2/c - b*(-3*a*c + b**2)/(c*sqrt(-4*a*c + b**2)))*hyper((1, n + 1), (n + 2,), 2*c*(e + f*x)/(2*c*e - f*(b + sqrt(-4*a*c + b**2))))/(c*(n + 1)*(2*c*e - f*(b + sqrt(-4*a*c + b**2)))) + (e + f*x)**(n + 1)*(a - b**2/c + b*(-3*a*c + b**2)/(c*sqrt(-4*a*c + b**2)))*hyper((1, n + 1), (n + 2,), 2*c*(e + f*x)/(2*c*e - f*(b - sqrt(-4*a*c + b**2))))/(c*(n + 1)*(2*c*e - f*(b - sqrt(-4*a*c + b**2)))) + (e + f*x)**(n + 2)/(c*f**2*(n + 2)) - (e + f*x)**(n + 1)*(b*f + c*e)/(c**2*f**2*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_537():
    f = x**2*(e + f*x)**n/(a + b*x + c*x**2)
    F = (b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 2*c*(e + f*x)/(2*c*e - f*(b - sqrt(-4*a*c + b**2))))/(c*(n + 1)*(2*c*e - f*(b - sqrt(-4*a*c + b**2)))) + (b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 2*c*(e + f*x)/(2*c*e - f*(b + sqrt(-4*a*c + b**2))))/(c*(n + 1)*(2*c*e - f*(b + sqrt(-4*a*c + b**2)))) + (e + f*x)**(n + 1)/(c*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_538():
    f = x*(e + f*x)**n/(a + b*x + c*x**2)
    F = -(e + f*x)**(n + 1)*(-b/sqrt(-4*a*c + b**2) + 1)*hyper((1, n + 1), (n + 2,), 2*c*(e + f*x)/(2*c*e - f*(b - sqrt(-4*a*c + b**2))))/((n + 1)*(2*c*e - f*(b - sqrt(-4*a*c + b**2)))) - (e + f*x)**(n + 1)*(b/sqrt(-4*a*c + b**2) + 1)*hyper((1, n + 1), (n + 2,), 2*c*(e + f*x)/(2*c*e - f*(b + sqrt(-4*a*c + b**2))))/((n + 1)*(2*c*e - f*(b + sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_539():
    f = (e + f*x)**n/(a + b*x + c*x**2)
    F = 2*c*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 2*c*(e + f*x)/(2*c*e - f*(b + sqrt(-4*a*c + b**2))))/((n + 1)*sqrt(-4*a*c + b**2)*(2*c*e - f*(b + sqrt(-4*a*c + b**2)))) - 2*c*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 2*c*(e + f*x)/(2*c*e - f*(b - sqrt(-4*a*c + b**2))))/((n + 1)*sqrt(-4*a*c + b**2)*(2*c*e - f*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_540():
    f = (e + f*x)**n/(x*(a + b*x + c*x**2))
    F = c*(e + f*x)**(n + 1)*(-b/sqrt(-4*a*c + b**2) + 1)*hyper((1, n + 1), (n + 2,), 2*c*(e + f*x)/(2*c*e - f*(b + sqrt(-4*a*c + b**2))))/(a*(n + 1)*(2*c*e - f*(b + sqrt(-4*a*c + b**2)))) + c*(e + f*x)**(n + 1)*(b/sqrt(-4*a*c + b**2) + 1)*hyper((1, n + 1), (n + 2,), 2*c*(e + f*x)/(2*c*e - f*(b - sqrt(-4*a*c + b**2))))/(a*(n + 1)*(2*c*e - f*(b - sqrt(-4*a*c + b**2)))) - (e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + f*x/e)/(a*e*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_541():
    f = (e + f*x)**n/(x**2*(a + b*x + c*x**2))
    F = f*(e + f*x)**(n + 1)*hyper((2, n + 1), (n + 2,), 1 + f*x/e)/(a*e**2*(n + 1)) + b*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + f*x/e)/(a**2*e*(n + 1)) - c*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 2*c*(e + f*x)/(2*c*e - f*(b + sqrt(-4*a*c + b**2))))/(a**2*(n + 1)*(2*c*e - f*(b + sqrt(-4*a*c + b**2)))) - c*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 2*c*(e + f*x)/(2*c*e - f*(b - sqrt(-4*a*c + b**2))))/(a**2*(n + 1)*(2*c*e - f*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_542():
    f = (d + e*x)**4*(f + g*x)**2/(d**2 - e**2*x**2)
    F = -8*d**3*(d*g + e*f)**2*log(d - e*x)/e**3 - d**2*x*(8*d**2*g**2 + 16*d*e*f*g + 7*e**2*f**2)/e**2 - d*x**2*(4*d**2*g**2 + 7*d*e*f*g + 2*e**2*f**2)/e - e**2*g**2*x**5/5 - e*g*x**4*(2*d*g + e*f)/2 - x**3*(d*g/3 + e*f/3)*(7*d*g + e*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_543():
    f = (d + e*x)**3*(f + g*x)**2/(d**2 - e**2*x**2)
    F = -4*d**2*(d*g + e*f)**2*log(d - e*x)/e**3 - d*x*(2*d*g + e*f)*(2*d*g + 3*e*f)/e**2 - e*g**2*x**4/4 - g*x**3*(3*d*g + 2*e*f)/3 - x**2*(4*d**2*g**2 + 6*d*e*f*g + e**2*f**2)/(2*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_544():
    f = (d + e*x)**2*(f + g*x)**2/(d**2 - e**2*x**2)
    F = -d*(f + g*x)**2/e - 2*d*g*x*(d*g + e*f)/e**2 - 2*d*(d*g + e*f)**2*log(d - e*x)/e**3 - (f + g*x)**3/(3*g)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_545():
    f = (d + e*x)*(f + g*x)**2/(d**2 - e**2*x**2)
    F = -(f + g*x)**2/(2*e) - g*x*(d*g + e*f)/e**2 - (d*g + e*f)**2*log(d - e*x)/e**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_546():
    f = (f + g*x)**2/(d**2 - e**2*x**2)
    F = -g**2*x/e**2 + (-d*g + e*f)**2*log(d + e*x)/(2*d*e**3) - (d*g + e*f)**2*log(d - e*x)/(2*d*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_547():
    f = (f + g*x)**2/((d + e*x)*(d**2 - e**2*x**2))
    F = -(-d*g + e*f)**2/(2*d*e**3*(d + e*x)) + (-d*g + e*f)*(3*d*g + e*f)*log(d + e*x)/(4*d**2*e**3) - (d*g + e*f)**2*log(d - e*x)/(4*d**2*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_548():
    f = (f + g*x)**2/((d + e*x)**2*(d**2 - e**2*x**2))
    F = -(-d*g + e*f)**2/(4*d*e**3*(d + e*x)**2) - (-d*g + e*f)*(3*d*g + e*f)/(4*d**2*e**3*(d + e*x)) + (d*g + e*f)**2*atanh(e*x/d)/(4*d**3*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_549():
    f = (f + g*x)**2/((d + e*x)**3*(d**2 - e**2*x**2))
    F = -(-d*g + e*f)**2/(6*d*e**3*(d + e*x)**3) - (-d*g + e*f)*(3*d*g + e*f)/(8*d**2*e**3*(d + e*x)**2) - (d*g + e*f)**2/(8*d**3*e**3*(d + e*x)) + (d*g + e*f)**2*atanh(e*x/d)/(8*d**4*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_550():
    f = (f + g*x)**2/((d + e*x)**4*(d**2 - e**2*x**2))
    F = -(-d*g + e*f)**2/(8*d*e**3*(d + e*x)**4) - (-d*g + e*f)*(3*d*g + e*f)/(12*d**2*e**3*(d + e*x)**3) - (d*g + e*f)**2/(16*d**3*e**3*(d + e*x)**2) - (d*g + e*f)**2/(16*d**4*e**3*(d + e*x)) + (d*g + e*f)**2*atanh(e*x/d)/(16*d**5*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_551():
    f = (d + e*x)**7*(f + g*x)**2/(d**2 - e**2*x**2)**2
    F = 32*d**5*(d*g + e*f)**2/(e**3*(d - e*x)) + 16*d**4*(d*g + e*f)*(9*d*g + 5*e*f)*log(d - e*x)/e**3 + d**3*x*(112*d**2*g**2 + 160*d*e*f*g + 49*e**2*f**2)/e**2 + d**2*x**2*(80*d**2*g**2 + 98*d*e*f*g + 23*e**2*f**2)/(2*e) + d*x**3*(49*d**2*g**2 + 46*d*e*f*g + 7*e**2*f**2)/3 + e**3*g**2*x**6/6 + e**2*g*x**5*(7*d*g + 2*e*f)/5 + e*x**4*(23*d**2*g**2 + 14*d*e*f*g + e**2*f**2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_552():
    f = (d + e*x)**6*(f + g*x)**2/(d**2 - e**2*x**2)**2
    F = 16*d**4*(d*g + e*f)**2/(e**3*(d - e*x)) + 32*d**3*(d*g + e*f)*(2*d*g + e*f)*log(d - e*x)/e**3 + d**2*x*(48*d**2*g**2 + 64*d*e*f*g + 17*e**2*f**2)/e**2 + d*x**2*(16*d**2*g**2 + 17*d*e*f*g + 3*e**2*f**2)/e + e**2*g**2*x**5/5 + e*g*x**4*(3*d*g + e*f)/2 + x**3*(17*d**2*g**2/3 + 4*d*e*f*g + e**2*f**2/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_553():
    f = (d + e*x)**5*(f + g*x)**2/(d**2 - e**2*x**2)**2
    F = 8*d**3*(d*g + e*f)**2/(e**3*(d - e*x)) + 4*d**2*(d*g + e*f)*(7*d*g + 3*e*f)*log(d - e*x)/e**3 + d*x*(20*d**2*g**2 + 24*d*e*f*g + 5*e**2*f**2)/e**2 + e*g**2*x**4/4 + g*x**3*(5*d*g + 2*e*f)/3 + x**2*(12*d**2*g**2 + 10*d*e*f*g + e**2*f**2)/(2*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_554():
    f = (d + e*x)**4*(f + g*x)**2/(d**2 - e**2*x**2)**2
    F = 4*d**2*(d*g + e*f)**2/(e**3*(d - e*x)) + 4*d*(d*g + e*f)*(3*d*g + e*f)*log(d - e*x)/e**3 + g**2*x**3/3 + g*x**2*(2*d*g + e*f)/e + x*(8*d**2*g**2 + 8*d*e*f*g + e**2*f**2)/e**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_555():
    f = (d + e*x)**3*(f + g*x)**2/(d**2 - e**2*x**2)**2
    F = 2*d*(d*g + e*f)**2/(e**3*(d - e*x)) + g**2*x**2/(2*e) + g*x*(3*d*g + 2*e*f)/e**2 + (d*g + e*f)*(5*d*g + e*f)*log(d - e*x)/e**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_556():
    f = (d + e*x)**2*(f + g*x)**2/(d**2 - e**2*x**2)**2
    F = g**2*x/e**2 + 2*g*(d*g + e*f)*log(d - e*x)/e**3 + (d*g + e*f)**2/(e**3*(d - e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_557():
    f = (d + e*x)*(f + g*x)**2/(d**2 - e**2*x**2)**2
    F = (d*g + e*f)**2/(2*d*e**3*(d - e*x)) - (-3*d*g + e*f)*(d*g + e*f)*log(d - e*x)/(4*d**2*e**3) + (-d*g + e*f)**2*log(d + e*x)/(4*d**2*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_558():
    f = (f + g*x)**2/(d**2 - e**2*x**2)**2
    F = (f + g*x)*(d**2*g + e**2*f*x)/(2*d**2*e**2*(d**2 - e**2*x**2)) + (-d*g + e*f)*(d*g + e*f)*atanh(e*x/d)/(2*d**3*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_559():
    f = (f + g*x)**2/((d + e*x)*(d**2 - e**2*x**2)**2)
    F = -(-d*g + e*f)**2/(8*d**2*e**3*(d + e*x)**2) - (-d**2*g**2 + e**2*f**2)/(4*d**3*e**3*(d + e*x)) + (d*g + e*f)**2/(8*d**3*e**3*(d - e*x)) + (-d*g + 3*e*f)*(d*g + e*f)*atanh(e*x/d)/(8*d**4*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_560():
    f = (f + g*x)**2/((d + e*x)**2*(d**2 - e**2*x**2)**2)
    F = -(-d*g + e*f)**2/(12*d**2*e**3*(d + e*x)**3) - (-d**2*g**2 + e**2*f**2)/(8*d**3*e**3*(d + e*x)**2) - (-d*g + 3*e*f)*(d*g + e*f)/(16*d**4*e**3*(d + e*x)) + (d*g + e*f)**2/(16*d**4*e**3*(d - e*x)) + f*(d*g + e*f)*atanh(e*x/d)/(4*d**5*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_561():
    f = (f + g*x)**2/((d + e*x)**3*(d**2 - e**2*x**2)**2)
    F = -(-d*g + e*f)**2/(16*d**2*e**3*(d + e*x)**4) - (-d**2*g**2 + e**2*f**2)/(12*d**3*e**3*(d + e*x)**3) - (-d*g + 3*e*f)*(d*g + e*f)/(32*d**4*e**3*(d + e*x)**2) - f*(d*g + e*f)/(8*d**5*e**2*(d + e*x)) + (d*g + e*f)**2/(32*d**5*e**3*(d - e*x)) + (d*g + e*f)*(d*g + 5*e*f)*atanh(e*x/d)/(32*d**6*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_562():
    f = (f + g*x)**2/((d + e*x)**4*(d**2 - e**2*x**2)**2)
    F = -(-d*g + e*f)**2/(20*d**2*e**3*(d + e*x)**5) - (-d**2*g**2 + e**2*f**2)/(16*d**3*e**3*(d + e*x)**4) - (-d*g + 3*e*f)*(d*g + e*f)/(48*d**4*e**3*(d + e*x)**3) - f*(d*g + e*f)/(16*d**5*e**2*(d + e*x)**2) - (d*g + e*f)*(d*g + 5*e*f)/(64*d**6*e**3*(d + e*x)) + (d*g + e*f)**2/(64*d**6*e**3*(d - e*x)) + (d*g + e*f)*(d*g + 3*e*f)*atanh(e*x/d)/(32*d**7*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_563():
    f = (d + e*x)**7*(f + g*x)**2/(d**2 - e**2*x**2)**3
    F = 8*d**4*(d*g + e*f)**2/(e**3*(d - e*x)**2) - 32*d**3*(d*g + e*f)*(2*d*g + e*f)/(e**3*(d - e*x)) - 8*d**2*(13*d**2*g**2 + 14*d*e*f*g + 3*e**2*f**2)*log(d - e*x)/e**3 - d*x*(56*d**2*g**2 + 48*d*e*f*g + 7*e**2*f**2)/e**2 - e*g**2*x**4/4 - g*x**3*(7*d*g + 2*e*f)/3 - x**2*(2*d*g + e*f)*(12*d*g + e*f)/(2*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_564():
    f = (d + e*x)**6*(f + g*x)**2/(d**2 - e**2*x**2)**3
    F = 4*d**3*(d*g + e*f)**2/(e**3*(d - e*x)**2) - 4*d**2*(d*g + e*f)*(7*d*g + 3*e*f)/(e**3*(d - e*x)) - 2*d*(19*d**2*g**2 + 18*d*e*f*g + 3*e**2*f**2)*log(d - e*x)/e**3 - g**2*x**3/3 - g*x**2*(3*d*g + e*f)/e - x*(18*d**2*g**2 + 12*d*e*f*g + e**2*f**2)/e**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_565():
    f = (d + e*x)**5*(f + g*x)**2/(d**2 - e**2*x**2)**3
    F = 2*d**2*(d*g + e*f)**2/(e**3*(d - e*x)**2) - 4*d*(d*g + e*f)*(3*d*g + e*f)/(e**3*(d - e*x)) - g**2*x**2/(2*e) - g*x*(5*d*g + 2*e*f)/e**2 - (13*d**2*g**2 + 10*d*e*f*g + e**2*f**2)*log(d - e*x)/e**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_566():
    f = (d + e*x)**4*(f + g*x)**2/(d**2 - e**2*x**2)**3
    F = d*(d*g + e*f)**2/(e**3*(d - e*x)**2) - g**2*x/e**2 - 2*g*(2*d*g + e*f)*log(d - e*x)/e**3 - (d*g + e*f)*(5*d*g + e*f)/(e**3*(d - e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_567():
    f = (d + e*x)**3*(f + g*x)**2/(d**2 - e**2*x**2)**3
    F = -g**2*log(d - e*x)/e**3 - 2*g*(d*g + e*f)/(e**3*(d - e*x)) + (d*g + e*f)**2/(2*e**3*(d - e*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_568():
    f = (d + e*x)**2*(f + g*x)**2/(d**2 - e**2*x**2)**3
    F = (d*g + e*f)**2/(4*d*e**3*(d - e*x)**2) + (-3*d*g + e*f)*(d*g + e*f)/(4*d**2*e**3*(d - e*x)) + (-d*g + e*f)**2*atanh(e*x/d)/(4*d**3*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_569():
    f = (d + e*x)*(f + g*x)**2/(d**2 - e**2*x**2)**3
    F = (d*g + e*f)**2/(8*d**2*e**3*(d - e*x)**2) - (-d*g + e*f)**2/(8*d**3*e**3*(d + e*x)) + (-d**2*g**2 + e**2*f**2)/(4*d**3*e**3*(d - e*x)) + (-d*g + e*f)*(d*g + 3*e*f)*atanh(e*x/d)/(8*d**4*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_570():
    f = (f + g*x)**2/(d**2 - e**2*x**2)**3
    F = (f + g*x)*(d**2*g + e**2*f*x)/(4*d**2*e**2*(d**2 - e**2*x**2)**2) + (2*d**2*f*g + x*(-d**2*g**2 + 3*e**2*f**2))/(8*d**4*e**2*(d**2 - e**2*x**2)) + (-d**2*g**2 + 3*e**2*f**2)*atanh(e*x/d)/(8*d**5*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_571():
    f = (f + g*x)**2/((d + e*x)*(d**2 - e**2*x**2)**3)
    F = -(-d*g + e*f)**2/(24*d**3*e**3*(d + e*x)**3) - (-d*g + e*f)*(d*g + 3*e*f)/(32*d**4*e**3*(d + e*x)**2) + (d*g + e*f)**2/(32*d**4*e**3*(d - e*x)**2) + f*(d*g + e*f)/(8*d**5*e**2*(d - e*x)) - (-d**2*g**2 + 3*e**2*f**2)/(16*d**5*e**3*(d + e*x)) + (-d**2*g**2 + 2*d*e*f*g + 5*e**2*f**2)*atanh(e*x/d)/(16*d**6*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_572():
    f = (f + g*x)**2/((d + e*x)**2*(d**2 - e**2*x**2)**3)
    F = -(-d*g + e*f)**2/(32*d**3*e**3*(d + e*x)**4) - (-d*g + e*f)*(d*g + 3*e*f)/(48*d**4*e**3*(d + e*x)**3) - (-d**2*g**2 + 3*e**2*f**2)/(32*d**5*e**3*(d + e*x)**2) + (d*g + e*f)**2/(64*d**5*e**3*(d - e*x)**2) - (-d**2*g**2 + 2*d*e*f*g + 5*e**2*f**2)/(32*d**6*e**3*(d + e*x)) + (d*g + e*f)*(d*g + 5*e*f)/(64*d**6*e**3*(d - e*x)) + (-d**2*g**2 + 10*d*e*f*g + 15*e**2*f**2)*atanh(e*x/d)/(64*d**7*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_573():
    f = (d + e*x)**3*(f + g*x)**5/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = g**5*x*sqrt(d**2 - e**2*x**2)/(2*e**5) + g**4*sqrt(d**2 - e**2*x**2)*(3*d*g + 5*e*f)/e**6 - g**3*(13*d**2*g**2 + 30*d*e*f*g + 20*e**2*f**2)*atan(e*x/sqrt(d**2 - e**2*x**2))/(2*e**6) + (d + e*x)**3*(d*g + e*f)**5/(5*d*e**6*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + (d + e*x)**2*(-23*d*g + 2*e*f)*(d*g + e*f)**4/(15*d**2*e**6*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (d + e*x)*(d*g + e*f)**3*(127*d**2*g**2 - 21*d*e*f*g + 2*e**2*f**2)/(15*d**3*e**6*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_574():
    f = (d + e*x)**3*(f + g*x)**4/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = g**4*sqrt(d**2 - e**2*x**2)/e**5 - g**3*(3*d*g + 4*e*f)*atan(e*x/sqrt(d**2 - e**2*x**2))/e**5 + (d + e*x)**3*(d*g + e*f)**4/(5*d*e**5*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + (d + e*x)**2*(-18*d*g + 2*e*f)*(d*g + e*f)**3/(15*d**2*e**5*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + 2*(d + e*x)*(d*g + e*f)**2*(36*d**2*g**2 - 8*d*e*f*g + e**2*f**2)/(15*d**3*e**5*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_575():
    f = (d + e*x)**3*(f + g*x)**3/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = -g**3*atan(e*x/sqrt(d**2 - e**2*x**2))/e**4 + (d + e*x)**3*(d*g + e*f)**3/(5*d*e**4*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + (d + e*x)**2*(-13*d*g + 2*e*f)*(d*g + e*f)**2/(15*d**2*e**4*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (d + e*x)*(d*g + e*f)*(32*d**2*g**2 - 11*d*e*f*g + 2*e**2*f**2)/(15*d**3*e**4*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_576():
    f = (d + e*x)**3*(f + g*x)**2/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = (d + e*x)**3*(d*g + e*f)**2/(5*d*e**3*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + (d + e*x)**2*(-8*d*g + 2*e*f)*(d*g + e*f)/(15*d**2*e**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + (d + e*x)*(7*d**2*g**2 - 6*d*e*f*g + 2*e**2*f**2)/(15*d**3*e**3*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_577():
    f = (d + e*x)**3*(f + g*x)/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = (d + e*x)**3*(d*g + e*f)/(5*d*e**2*(d**2 - e**2*x**2)**(sympy.S(5)/2)) + (d + e*x)*(-6*d*g + 4*e*f)/(15*d*e**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)) + x*(-3*d*g + 2*e*f)/(15*d**3*e*sqrt(d**2 - e**2*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_578():
    f = (d + e*x)**3/(d**2 - e**2*x**2)**(sympy.S(7)/2)
    F = sqrt(d**2 - e**2*x**2)/(5*d*e*(d - e*x)**3) + 2*sqrt(d**2 - e**2*x**2)/(15*d**2*e*(d - e*x)**2) + 2*sqrt(d**2 - e**2*x**2)/(15*d**3*e*(d - e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_579():
    f = (d + e*x)**3/((d**2 - e**2*x**2)**(sympy.S(7)/2)*(f + g*x))
    F = 4*d*(d + e*x)/((d**2 - e**2*x**2)**(sympy.S(5)/2)*(5*d*g + 5*e*f)) + g**3*atan((d**2*g + e**2*f*x)/(sqrt(d**2 - e**2*x**2)*sqrt(-d**2*g**2 + e**2*f**2)))/((d*g + e*f)**3*sqrt(-d**2*g**2 + e**2*f**2)) - (5*d*(-d*g + e*f) - e*x*(11*d*g + e*f))/(15*d*(d**2 - e**2*x**2)**(sympy.S(3)/2)*(d*g + e*f)**2) + (15*d**3*g**2 + e*x*(22*d**2*g**2 + 9*d*e*f*g + 2*e**2*f**2))/(15*d**3*sqrt(d**2 - e**2*x**2)*(d*g + e*f)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_580():
    f = (d + e*x)**3/((d**2 - e**2*x**2)**(sympy.S(7)/2)*(f + g*x)**2)
    F = 4*d*e*(d + e*x)/(5*(d**2 - e**2*x**2)**(sympy.S(5)/2)*(d*g + e*f)**2) + e*g**3*(-3*d*g + 4*e*f)*atan((d**2*g + e**2*f*x)/(sqrt(d**2 - e**2*x**2)*sqrt(-d**2*g**2 + e**2*f**2)))/((-d*g + e*f)*(d*g + e*f)**4*sqrt(-d**2*g**2 + e**2*f**2)) + g**4*sqrt(d**2 - e**2*x**2)/((f + g*x)*(-d*g + e*f)*(d*g + e*f)**4) - e*(5*d*(-3*d*g + e*f) - e*x*(21*d*g + e*f))/(15*d*(d**2 - e**2*x**2)**(sympy.S(3)/2)*(d*g + e*f)**3) + e*(45*d**3*g**2 + e*x*(57*d**2*g**2 + 14*d*e*f*g + 2*e**2*f**2))/(15*d**3*sqrt(d**2 - e**2*x**2)*(d*g + e*f)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_581():
    f = (d + e*x)**3/((d**2 - e**2*x**2)**(sympy.S(7)/2)*(f + g*x)**3)
    F = 4*d*e**2*(d + e*x)/(5*(d**2 - e**2*x**2)**(sympy.S(5)/2)*(d*g + e*f)**3) + e**2*g**3*(13*d**2*g**2 - 30*d*e*f*g + 20*e**2*f**2)*atan((d**2*g + e**2*f*x)/(sqrt(d**2 - e**2*x**2)*sqrt(-d**2*g**2 + e**2*f**2)))/(2*(-d*g + e*f)**2*(d*g + e*f)**5*sqrt(-d**2*g**2 + e**2*f**2)) + 3*e*g**4*sqrt(d**2 - e**2*x**2)*(-2*d*g + 3*e*f)/(2*(f + g*x)*(-d*g + e*f)**2*(d*g + e*f)**5) + g**4*sqrt(d**2 - e**2*x**2)/((f + g*x)**2*(-2*d*g + 2*e*f)*(d*g + e*f)**4) - e**2*(5*d*(-5*d*g + e*f) - e*x*(31*d*g + e*f))/(15*d*(d**2 - e**2*x**2)**(sympy.S(3)/2)*(d*g + e*f)**4) + e**2*(90*d**3*g**2 + e*x*(107*d**2*g**2 + 19*d*e*f*g + 2*e**2*f**2))/(15*d**3*sqrt(d**2 - e**2*x**2)*(d*g + e*f)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_582():
    f = (a + c*x**2)/((d + e*x)**(sympy.S(3)/2)*(f + g*x))
    F = 2*c*sqrt(d + e*x)/(e**2*g) - (2*a*g**2 + 2*c*f**2)*atan(sqrt(g)*sqrt(d + e*x)/sqrt(-d*g + e*f))/(g**(sympy.S(3)/2)*(-d*g + e*f)**(sympy.S(3)/2)) - (2*a*e**2 + 2*c*d**2)/(e**2*sqrt(d + e*x)*(-d*g + e*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_583():
    f = (a + c*x**2)*(d + e*x)**3/sqrt(f + g*x)
    F = 2*c*e**3*(f + g*x)**(sympy.S(11)/2)/(11*g**6) - 2*c*e**2*(f + g*x)**(sympy.S(9)/2)*(-3*d*g + 5*e*f)/(9*g**6) + 2*e*(f + g*x)**(sympy.S(7)/2)*(a*e**2*g**2 + c*(3*d**2*g**2 - 12*d*e*f*g + 10*e**2*f**2))/(7*g**6) - (f + g*x)**(sympy.S(5)/2)*(-2*d*g + 2*e*f)*(3*a*e**2*g**2 + c*(d**2*g**2 - 8*d*e*f*g + 10*e**2*f**2))/(5*g**6) + 2*(f + g*x)**(sympy.S(3)/2)*(-d*g + e*f)**2*(3*a*e*g**2 + c*f*(-2*d*g + 5*e*f))/(3*g**6) - 2*sqrt(f + g*x)*(a*g**2 + c*f**2)*(-d*g + e*f)**3/g**6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_584():
    f = (a + c*x**2)*(d + e*x)**2/sqrt(f + g*x)
    F = 2*c*e**2*(f + g*x)**(sympy.S(9)/2)/(9*g**5) - 4*c*e*(f + g*x)**(sympy.S(7)/2)*(-d*g + 2*e*f)/(7*g**5) + (f + g*x)**(sympy.S(5)/2)*(2*a*e**2*g**2 + 2*c*(d**2*g**2 - 6*d*e*f*g + 6*e**2*f**2))/(5*g**5) - (f + g*x)**(sympy.S(3)/2)*(-4*d*g + 4*e*f)*(a*e*g**2 + c*f*(-d*g + 2*e*f))/(3*g**5) + 2*sqrt(f + g*x)*(a*g**2 + c*f**2)*(-d*g + e*f)**2/g**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_585():
    f = (a + c*x**2)*(d + e*x)/sqrt(f + g*x)
    F = 2*c*e*(f + g*x)**(sympy.S(7)/2)/(7*g**4) - 2*c*(f + g*x)**(sympy.S(5)/2)*(-d*g + 3*e*f)/(5*g**4) + (f + g*x)**(sympy.S(3)/2)*(2*a*e*g**2 + 2*c*f*(-2*d*g + 3*e*f))/(3*g**4) - sqrt(f + g*x)*(a*g**2 + c*f**2)*(-2*d*g + 2*e*f)/g**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_586():
    f = (a + c*x**2)/sqrt(f + g*x)
    F = -4*c*f*(f + g*x)**(sympy.S(3)/2)/(3*g**3) + 2*c*(f + g*x)**(sympy.S(5)/2)/(5*g**3) + sqrt(f + g*x)*(2*a*g**2 + 2*c*f**2)/g**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_587():
    f = (a + c*x**2)/((d + e*x)*sqrt(f + g*x))
    F = 2*c*(f + g*x)**(sympy.S(3)/2)/(3*e*g**2) - 2*c*sqrt(f + g*x)*(d*g + e*f)/(e**2*g**2) - (2*a*e**2 + 2*c*d**2)*atanh(sqrt(e)*sqrt(f + g*x)/sqrt(-d*g + e*f))/(e**(sympy.S(5)/2)*sqrt(-d*g + e*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_588():
    f = (a + c*x**2)/((d + e*x)**2*sqrt(f + g*x))
    F = 2*c*sqrt(f + g*x)/(e**2*g) - (a + c*d**2/e**2)*sqrt(f + g*x)/((d + e*x)*(-d*g + e*f)) + (a*e**2*g + c*d*(-3*d*g + 4*e*f))*atanh(sqrt(e)*sqrt(f + g*x)/sqrt(-d*g + e*f))/(e**(sympy.S(5)/2)*(-d*g + e*f)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_589():
    f = (a + c*x**2)/((d + e*x)**3*sqrt(f + g*x))
    F = -(a + c*d**2/e**2)*sqrt(f + g*x)/((d + e*x)**2*(-2*d*g + 2*e*f)) + sqrt(f + g*x)*(3*a*e**2*g + c*d*(-5*d*g + 8*e*f))/(4*e**2*(d + e*x)*(-d*g + e*f)**2) - (3*a*e**2*g**2 + c*(3*d**2*g**2 - 8*d*e*f*g + 8*e**2*f**2))*atanh(sqrt(e)*sqrt(f + g*x)/sqrt(-d*g + e*f))/(4*e**(sympy.S(5)/2)*(-d*g + e*f)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_590():
    f = (a + c*x**2)*(d + e*x)**3/(f + g*x)**(sympy.S(3)/2)
    F = 2*c*e**3*(f + g*x)**(sympy.S(9)/2)/(9*g**6) - 2*c*e**2*(f + g*x)**(sympy.S(7)/2)*(-3*d*g + 5*e*f)/(7*g**6) + 2*e*(f + g*x)**(sympy.S(5)/2)*(a*e**2*g**2 + c*(3*d**2*g**2 - 12*d*e*f*g + 10*e**2*f**2))/(5*g**6) - (f + g*x)**(sympy.S(3)/2)*(-2*d*g + 2*e*f)*(3*a*e**2*g**2 + c*(d**2*g**2 - 8*d*e*f*g + 10*e**2*f**2))/(3*g**6) + 2*sqrt(f + g*x)*(-d*g + e*f)**2*(3*a*e*g**2 + c*f*(-2*d*g + 5*e*f))/g**6 + 2*(a*g**2 + c*f**2)*(-d*g + e*f)**3/(g**6*sqrt(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_591():
    f = (a + c*x**2)*(d + e*x)**2/(f + g*x)**(sympy.S(3)/2)
    F = 2*c*e**2*(f + g*x)**(sympy.S(7)/2)/(7*g**5) - 4*c*e*(f + g*x)**(sympy.S(5)/2)*(-d*g + 2*e*f)/(5*g**5) + (f + g*x)**(sympy.S(3)/2)*(2*a*e**2*g**2 + 2*c*(d**2*g**2 - 6*d*e*f*g + 6*e**2*f**2))/(3*g**5) - sqrt(f + g*x)*(-4*d*g + 4*e*f)*(a*e*g**2 + c*f*(-d*g + 2*e*f))/g**5 - 2*(a*g**2 + c*f**2)*(-d*g + e*f)**2/(g**5*sqrt(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_592():
    f = (a + c*x**2)*(d + e*x)/(f + g*x)**(sympy.S(3)/2)
    F = 2*c*e*(f + g*x)**(sympy.S(5)/2)/(5*g**4) - 2*c*(f + g*x)**(sympy.S(3)/2)*(-d*g + 3*e*f)/(3*g**4) + sqrt(f + g*x)*(2*a*e*g**2 + 2*c*f*(-2*d*g + 3*e*f))/g**4 + (a*g**2 + c*f**2)*(-2*d*g + 2*e*f)/(g**4*sqrt(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_593():
    f = (a + c*x**2)/(f + g*x)**(sympy.S(3)/2)
    F = -4*c*f*sqrt(f + g*x)/g**3 + 2*c*(f + g*x)**(sympy.S(3)/2)/(3*g**3) - (2*a*g**2 + 2*c*f**2)/(g**3*sqrt(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_594():
    f = (a + c*x**2)/((d + e*x)*(f + g*x)**(sympy.S(3)/2))
    F = 2*c*sqrt(f + g*x)/(e*g**2) + (2*a*g**2 + 2*c*f**2)/(g**2*sqrt(f + g*x)*(-d*g + e*f)) - (2*a*e**2 + 2*c*d**2)*atanh(sqrt(e)*sqrt(f + g*x)/sqrt(-d*g + e*f))/(e**(sympy.S(3)/2)*(-d*g + e*f)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_595():
    f = (a + c*x**2)/((d + e*x)**2*(f + g*x)**(sympy.S(3)/2))
    F = -(2*a*g**2 + 2*c*f**2)/(g*sqrt(f + g*x)*(-d*g + e*f)**2) - sqrt(f + g*x)*(a*e**2 + c*d**2)/(e*(d + e*x)*(-d*g + e*f)**2) + (3*a*e**2*g + c*d*(-d*g + 4*e*f))*atanh(sqrt(e)*sqrt(f + g*x)/sqrt(-d*g + e*f))/(e**(sympy.S(3)/2)*(-d*g + e*f)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_596():
    f = (a + c*x**2)/((d + e*x)**3*(f + g*x)**(sympy.S(3)/2))
    F = (2*a*g**2 + 2*c*f**2)/(sqrt(f + g*x)*(-d*g + e*f)**3) + sqrt(f + g*x)*(7*a*e**2*g + c*d*(-d*g + 8*e*f))/(4*e*(d + e*x)*(-d*g + e*f)**3) - sqrt(f + g*x)*(a*e**2 + c*d**2)/(2*e*(d + e*x)**2*(-d*g + e*f)**2) - (15*a*e**2*g**2 + c*(-d**2*g**2 + 8*d*e*f*g + 8*e**2*f**2))*atanh(sqrt(e)*sqrt(f + g*x)/sqrt(-d*g + e*f))/(4*e**(sympy.S(3)/2)*(-d*g + e*f)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_597():
    f = (a + c*x**2)/(sqrt(d + e*x)*sqrt(f + g*x))
    F = c*(d + e*x)**(sympy.S(3)/2)*sqrt(f + g*x)/(2*e**2*g) - c*sqrt(d + e*x)*sqrt(f + g*x)*(5*d*g + 3*e*f)/(4*e**2*g**2) + (8*a*e**2*g**2 + c*(3*d**2*g**2 + 2*d*e*f*g + 3*e**2*f**2))*atanh(sqrt(g)*sqrt(d + e*x)/(sqrt(e)*sqrt(f + g*x)))/(4*e**(sympy.S(5)/2)*g**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_598():
    f = (2*x**2 - 1)/(sqrt(x - 1)*sqrt(x + 1))
    F = x*sqrt(x - 1)*sqrt(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_599():
    f = (d + e*x)**(sympy.S(3)/2)*sqrt(f + g*x)/(a + c*x**2)
    F = sqrt(e)*(3*d*g + e*f)*atanh(sqrt(g)*sqrt(d + e*x)/(sqrt(e)*sqrt(f + g*x)))/(c*sqrt(g)) + e*sqrt(d + e*x)*sqrt(f + g*x)/c + (a*(a*e**2*g - c*d*(d*g + 2*e*f))/sqrt(c) + sqrt(-a)*(-a*e*(2*d*g + e*f) + c*d**2*f))*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f + g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d + e*sqrt(-a))))/(a*c*sqrt(sqrt(c)*d + e*sqrt(-a))*sqrt(sqrt(c)*f + g*sqrt(-a))) + (a*(a*e**2*g - c*d*(d*g + 2*e*f))/sqrt(c) - sqrt(-a)*(-a*e*(2*d*g + e*f) + c*d**2*f))*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f - g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d - e*sqrt(-a))))/(a*c*sqrt(sqrt(c)*d - e*sqrt(-a))*sqrt(sqrt(c)*f - g*sqrt(-a)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_600():
    f = sqrt(d + e*x)*sqrt(f + g*x)/(a + c*x**2)
    F = 2*sqrt(e)*sqrt(g)*atanh(sqrt(g)*sqrt(d + e*x)/(sqrt(e)*sqrt(f + g*x)))/c - (-a*e*g + sqrt(c)*sqrt(-a)*(d*g + e*f) + c*d*f)*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f + g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d + e*sqrt(-a))))/(c*sqrt(-a)*sqrt(sqrt(c)*d + e*sqrt(-a))*sqrt(sqrt(c)*f + g*sqrt(-a))) + (-a*e*g - sqrt(c)*sqrt(-a)*(d*g + e*f) + c*d*f)*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f - g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d - e*sqrt(-a))))/(c*sqrt(-a)*sqrt(sqrt(c)*d - e*sqrt(-a))*sqrt(sqrt(c)*f - g*sqrt(-a)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_601():
    f = sqrt(f + g*x)/((a + c*x**2)*sqrt(d + e*x))
    F = -sqrt(sqrt(c)*f + g*sqrt(-a))*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f + g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d + e*sqrt(-a))))/(sqrt(c)*sqrt(-a)*sqrt(sqrt(c)*d + e*sqrt(-a))) + sqrt(sqrt(c)*f - g*sqrt(-a))*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f - g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d - e*sqrt(-a))))/(sqrt(c)*sqrt(-a)*sqrt(sqrt(c)*d - e*sqrt(-a)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_602():
    f = sqrt(f + g*x)/((a + c*x**2)*(d + e*x)**(sympy.S(3)/2))
    F = -2*e*sqrt(f + g*x)/(sqrt(d + e*x)*(a*e**2 + c*d**2)) - (a*e*g - sqrt(c)*sqrt(-a)*(-d*g + e*f) + c*d*f)*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f + g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d + e*sqrt(-a))))/(sqrt(-a)*(a*e**2 + c*d**2)*sqrt(sqrt(c)*d + e*sqrt(-a))*sqrt(sqrt(c)*f + g*sqrt(-a))) + (a*e*g + sqrt(c)*sqrt(-a)*(-d*g + e*f) + c*d*f)*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f - g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d - e*sqrt(-a))))/(sqrt(-a)*(a*e**2 + c*d**2)*sqrt(sqrt(c)*d - e*sqrt(-a))*sqrt(sqrt(c)*f - g*sqrt(-a)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_603():
    f = sqrt(f + g*x)/((a + c*x**2)*(d + e*x)**(sympy.S(5)/2))
    F = sqrt(c)*(a*e*g + sqrt(c)*sqrt(-a)*(-d*g + e*f) + c*d*f)*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f - g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d - e*sqrt(-a))))/(sqrt(-a)*(a*e**2 + c*d**2)*(sqrt(c)*d - e*sqrt(-a))**(sympy.S(3)/2)*sqrt(sqrt(c)*f - g*sqrt(-a))) + 4*e*g*sqrt(f + g*x)/(sqrt(d + e*x)*(3*a*e**2 + 3*c*d**2)*(-d*g + e*f)) - 2*e*sqrt(f + g*x)/((d + e*x)**(sympy.S(3)/2)*(3*a*e**2 + 3*c*d**2)) + e*sqrt(f + g*x)*(a*e*g - sqrt(c)*sqrt(-a)*(-d*g + e*f) + c*d*f)/(sqrt(-a)*sqrt(d + e*x)*(a*e**2 + c*d**2)*(sqrt(c)*d + e*sqrt(-a))*(-d*g + e*f)) - e*sqrt(f + g*x)*(a*e*g + sqrt(c)*sqrt(-a)*(-d*g + e*f) + c*d*f)/(sqrt(-a)*sqrt(d + e*x)*(a*e**2 + c*d**2)*(sqrt(c)*d - e*sqrt(-a))*(-d*g + e*f)) + sqrt(c)*(a*sqrt(c)*(-d*g + e*f) + a*e*g*sqrt(-a) + c*d*f*sqrt(-a))*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f + g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d + e*sqrt(-a))))/(a*(a*e**2 + c*d**2)*(sqrt(c)*d + e*sqrt(-a))**(sympy.S(3)/2)*sqrt(sqrt(c)*f + g*sqrt(-a)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_604():
    f = (d + e*x)**(sympy.S(3)/2)/((a + c*x**2)*sqrt(f + g*x))
    F = 2*e**(sympy.S(3)/2)*atanh(sqrt(g)*sqrt(d + e*x)/(sqrt(e)*sqrt(f + g*x)))/(c*sqrt(g)) - (-a*e**2 + 2*sqrt(c)*d*e*sqrt(-a) + c*d**2)*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f + g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d + e*sqrt(-a))))/(c*sqrt(-a)*sqrt(sqrt(c)*d + e*sqrt(-a))*sqrt(sqrt(c)*f + g*sqrt(-a))) + (-a*e**2 - 2*sqrt(c)*d*e*sqrt(-a) + c*d**2)*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f - g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d - e*sqrt(-a))))/(c*sqrt(-a)*sqrt(sqrt(c)*d - e*sqrt(-a))*sqrt(sqrt(c)*f - g*sqrt(-a)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_605():
    f = sqrt(d + e*x)/((a + c*x**2)*sqrt(f + g*x))
    F = sqrt(sqrt(c)*d - e*sqrt(-a))*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f - g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d - e*sqrt(-a))))/(sqrt(c)*sqrt(-a)*sqrt(sqrt(c)*f - g*sqrt(-a))) - sqrt(sqrt(c)*d + e*sqrt(-a))*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f + g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d + e*sqrt(-a))))/(sqrt(c)*sqrt(-a)*sqrt(sqrt(c)*f + g*sqrt(-a)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_606():
    f = 1/((a + c*x**2)*sqrt(d + e*x)*sqrt(f + g*x))
    F = -atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f + g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d + e*sqrt(-a))))/(sqrt(-a)*sqrt(sqrt(c)*d + e*sqrt(-a))*sqrt(sqrt(c)*f + g*sqrt(-a))) + atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f - g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d - e*sqrt(-a))))/(sqrt(-a)*sqrt(sqrt(c)*d - e*sqrt(-a))*sqrt(sqrt(c)*f - g*sqrt(-a)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_607():
    f = 1/((a + c*x**2)*(d + e*x)**(sympy.S(3)/2)*sqrt(f + g*x))
    F = -sqrt(c)*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f + g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d + e*sqrt(-a))))/(sqrt(-a)*(sqrt(c)*d + e*sqrt(-a))**(sympy.S(3)/2)*sqrt(sqrt(c)*f + g*sqrt(-a))) + sqrt(c)*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f - g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d - e*sqrt(-a))))/(sqrt(-a)*(sqrt(c)*d - e*sqrt(-a))**(sympy.S(3)/2)*sqrt(sqrt(c)*f - g*sqrt(-a))) + e*sqrt(f + g*x)/(sqrt(-a)*sqrt(d + e*x)*(sqrt(c)*d + e*sqrt(-a))*(-d*g + e*f)) - e*sqrt(f + g*x)/(sqrt(-a)*sqrt(d + e*x)*(sqrt(c)*d - e*sqrt(-a))*(-d*g + e*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_608():
    f = (d + e*x)**(sympy.S(3)/2)/((a + c*x**2)*(f + g*x)**(sympy.S(3)/2))
    F = -2*sqrt(e)*(-d*g + e*f)*atanh(sqrt(g)*sqrt(d + e*x)/(sqrt(e)*sqrt(f + g*x)))/(sqrt(g)*(a*g**2 + c*f**2)) + sqrt(d + e*x)*(-2*d*g + 2*e*f)/(sqrt(f + g*x)*(a*g**2 + c*f**2)) - sqrt(e)*(a*e*g - sqrt(c)*sqrt(-a)*(-d*g + e*f) + c*d*f)*atanh(sqrt(g)*sqrt(d + e*x)/(sqrt(e)*sqrt(f + g*x)))/(sqrt(c)*sqrt(g)*sqrt(-a)*(a*g**2 + c*f**2)) + sqrt(e)*(a*e*g + sqrt(c)*sqrt(-a)*(-d*g + e*f) + c*d*f)*atanh(sqrt(g)*sqrt(d + e*x)/(sqrt(e)*sqrt(f + g*x)))/(sqrt(c)*sqrt(g)*sqrt(-a)*(a*g**2 + c*f**2)) + sqrt(sqrt(c)*d - e*sqrt(-a))*(a*e*g - sqrt(c)*sqrt(-a)*(-d*g + e*f) + c*d*f)*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f - g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d - e*sqrt(-a))))/(sqrt(c)*sqrt(-a)*(a*g**2 + c*f**2)*sqrt(sqrt(c)*f - g*sqrt(-a))) - sqrt(sqrt(c)*d + e*sqrt(-a))*(a*e*g + sqrt(c)*sqrt(-a)*(-d*g + e*f) + c*d*f)*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f + g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d + e*sqrt(-a))))/(sqrt(c)*sqrt(-a)*(a*g**2 + c*f**2)*sqrt(sqrt(c)*f + g*sqrt(-a)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_609():
    f = sqrt(d + e*x)/((a + c*x**2)*(f + g*x)**(sympy.S(3)/2))
    F = -2*g*sqrt(d + e*x)/(sqrt(f + g*x)*(a*g**2 + c*f**2)) - (a*e*g + sqrt(c)*sqrt(-a)*(-d*g + e*f) + c*d*f)*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f + g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d + e*sqrt(-a))))/(sqrt(-a)*(a*g**2 + c*f**2)*sqrt(sqrt(c)*d + e*sqrt(-a))*sqrt(sqrt(c)*f + g*sqrt(-a))) + (a*e*g - sqrt(c)*sqrt(-a)*(-d*g + e*f) + c*d*f)*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f - g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d - e*sqrt(-a))))/(sqrt(-a)*(a*g**2 + c*f**2)*sqrt(sqrt(c)*d - e*sqrt(-a))*sqrt(sqrt(c)*f - g*sqrt(-a)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_610():
    f = 1/((a + c*x**2)*sqrt(d + e*x)*(f + g*x)**(sympy.S(3)/2))
    F = -sqrt(c)*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f + g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d + e*sqrt(-a))))/(sqrt(-a)*sqrt(sqrt(c)*d + e*sqrt(-a))*(sqrt(c)*f + g*sqrt(-a))**(sympy.S(3)/2)) + sqrt(c)*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f - g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d - e*sqrt(-a))))/(sqrt(-a)*sqrt(sqrt(c)*d - e*sqrt(-a))*(sqrt(c)*f - g*sqrt(-a))**(sympy.S(3)/2)) - g*sqrt(d + e*x)/(sqrt(-a)*sqrt(f + g*x)*(sqrt(c)*f + g*sqrt(-a))*(-d*g + e*f)) + g*sqrt(d + e*x)/(sqrt(-a)*sqrt(f + g*x)*(sqrt(c)*f - g*sqrt(-a))*(-d*g + e*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_611():
    f = 1/((a + c*x**2)*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(3)/2))
    F = -c*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f + g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d + e*sqrt(-a))))/(sqrt(-a)*(sqrt(c)*d + e*sqrt(-a))**(sympy.S(3)/2)*(sqrt(c)*f + g*sqrt(-a))**(sympy.S(3)/2)) + c*atanh(sqrt(d + e*x)*sqrt(sqrt(c)*f - g*sqrt(-a))/(sqrt(f + g*x)*sqrt(sqrt(c)*d - e*sqrt(-a))))/(sqrt(-a)*(sqrt(c)*d - e*sqrt(-a))**(sympy.S(3)/2)*(sqrt(c)*f - g*sqrt(-a))**(sympy.S(3)/2)) + e/(sqrt(-a)*sqrt(d + e*x)*sqrt(f + g*x)*(sqrt(c)*d + e*sqrt(-a))*(-d*g + e*f)) - e/(sqrt(-a)*sqrt(d + e*x)*sqrt(f + g*x)*(sqrt(c)*d - e*sqrt(-a))*(-d*g + e*f)) + g*sqrt(d + e*x)*(sqrt(c)*(d*g + e*f) + 2*e*g*sqrt(-a))/(sqrt(-a)*sqrt(f + g*x)*(sqrt(c)*d + e*sqrt(-a))*(sqrt(c)*f + g*sqrt(-a))*(-d*g + e*f)**2) + g*sqrt(d + e*x)*(-sqrt(c)*(d*g + e*f) + 2*e*g*sqrt(-a))/(sqrt(-a)*sqrt(f + g*x)*(sqrt(c)*d - e*sqrt(-a))*(sqrt(c)*f - g*sqrt(-a))*(-d*g + e*f)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_612():
    f = sqrt(x)/(sqrt(x + 1)*(x**2 + 1))
    F = -(1 - I)**(sympy.S(3)/2)*atanh(sqrt(x)*sqrt(1 - I)/sqrt(x + 1))/2 - (1 + I)**(sympy.S(3)/2)*atanh(sqrt(x)*sqrt(1 + I)/sqrt(x + 1))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_613():
    f = sqrt(1 - x**2)*(f + g*x)**2/(1 - x)**4
    F = -g**2*asin(x) + 2*g**2*(x + 1)/sqrt(1 - x**2) + (f - 9*g)*(f + g)*(x + 1)**3/(15*(1 - x**2)**(sympy.S(3)/2)) + (f + g)**2*(x + 1)**4/(5*(1 - x**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_614():
    f = (-a**2*x**2 + 1)**(sympy.S(3)/2)/((c + d*x)*(-a*x + 1)**2)
    F = -sqrt(-a**2*x**2 + 1)/d - (a*c - 2*d)*asin(a*x)/d**2 + (a*c - d)**2*atan((a**2*c*x + d)/(sqrt(a**2*c**2 - d**2)*sqrt(-a**2*x**2 + 1)))/(d**2*sqrt(a**2*c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_615():
    f = (a*x + 1)**2/((c + d*x)*sqrt(-a**2*x**2 + 1))
    F = -sqrt(-a**2*x**2 + 1)/d - (a*c - 2*d)*asin(a*x)/d**2 + (a*c - d)**2*atan((a**2*c*x + d)/(sqrt(a**2*c**2 - d**2)*sqrt(-a**2*x**2 + 1)))/(d**2*sqrt(a**2*c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_616():
    f = sqrt(a + c*x**2)*(d + e*x)**3*sqrt(f + g*x)
    F = 2*e**2*sqrt(a + c*x**2)*(f + g*x)**(sympy.S(7)/2)*(-3*d*g + e*f)/(99*g**4) + 2*sqrt(a + c*x**2)*(d + e*x)**4*sqrt(f + g*x)/(11*e) + 2*e*sqrt(a + c*x**2)*(f + g*x)**(sympy.S(5)/2)*(18*a*e**2*g**2 - c*(81*d**2*g**2 - 96*d*e*f*g + 29*e**2*f**2))/(693*c*g**4) - sqrt(a + c*x**2)*(f + g*x)**(sympy.S(3)/2)*(4*a*e**2*g**2*(-231*d*g + 74*e*f) - 2*c*(-567*d**3*g**3 + 1107*d**2*e*f*g**2 - 843*d*e**2*f**2*g + 233*e**3*f**3))/(3465*c*g**4) + sqrt(a + c*x**2)*sqrt(f + g*x)*(-300*a**2*e**4*g**4 + 12*a*c*e**2*g**2*(165*d**2*g**2 - 33*d*e*f*g + 2*e**2*f**2) - 2*c**2*(315*d**4*g**4 - 798*d**3*e*f*g**3 + 1098*d**2*e**2*f**2*g**2 - 732*d*e**3*f**3*g + 187*e**4*f**4))/(3465*c**2*e*g**4) + 4*sqrt(-a)*sqrt(1 + c*x**2/a)*sqrt(f + g*x)*(3*a**2*e**2*g**4*(231*d*g + 26*e*f) - 9*a*c*g**2*(77*d**3*g**3 + 88*d**2*e*f*g**2 - 33*d*e**2*f**2*g + 6*e**3*f**3) - c**2*f**2*(-231*d**3*g**3 + 396*d**2*e*f*g**2 - 264*d*e**2*f**2*g + 64*e**3*f**3))*elliptic_e(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(3465*c**(sympy.S(3)/2)*g**5*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(a + c*x**2)) - 4*sqrt(-a)*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(1 + c*x**2/a)*(a*g**2 + c*f**2)*(75*a**2*e**3*g**4 - 3*a*c*e*g**2*(165*d**2*g**2 - 33*d*e*f*g + 2*e**2*f**2) - c**2*f*(-231*d**3*g**3 + 396*d**2*e*f*g**2 - 264*d*e**2*f**2*g + 64*e**3*f**3))*elliptic_f(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(3465*c**(sympy.S(5)/2)*g**5*sqrt(a + c*x**2)*sqrt(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_617():
    f = sqrt(a + c*x**2)*(d + e*x)**2*sqrt(f + g*x)
    F = 2*e*sqrt(a + c*x**2)*(f + g*x)**(sympy.S(5)/2)*(-3*d*g + e*f)/(63*g**3) + 2*sqrt(a + c*x**2)*(d + e*x)**3*sqrt(f + g*x)/(9*e) + sqrt(a + c*x**2)*(f + g*x)**(sympy.S(3)/2)*(28*a*e**2*g**2 - 4*c*(21*d**2*g**2 - 24*d*e*f*g + 8*e**2*f**2))/(315*c*g**3) + sqrt(a + c*x**2)*sqrt(f + g*x)*(-12*a*e**2*g**2*(-10*d*g + e*f) + 2*c*(-35*d**3*g**3 + 63*d**2*e*f*g**2 - 57*d*e**2*f**2*g + 19*e**3*f**3))/(315*c*e*g**3) - 4*sqrt(-a)*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(1 + c*x**2/a)*(a*g**2 + c*f**2)*(3*a*e*g**2*(-10*d*g + e*f) + c*f*(21*d**2*g**2 - 24*d*e*f*g + 8*e**2*f**2))*elliptic_f(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(315*c**(sympy.S(3)/2)*g**4*sqrt(a + c*x**2)*sqrt(f + g*x)) + 4*sqrt(-a)*sqrt(1 + c*x**2/a)*sqrt(f + g*x)*(21*a**2*e**2*g**4 + 3*a*c*g**2*(-21*d**2*g**2 - 16*d*e*f*g + 3*e**2*f**2) + c**2*f**2*(21*d**2*g**2 - 24*d*e*f*g + 8*e**2*f**2))*elliptic_e(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(315*c**(sympy.S(3)/2)*g**4*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(a + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_618():
    f = sqrt(a + c*x**2)*(d + e*x)*sqrt(f + g*x)
    F = 2*e*(a + c*x**2)**(sympy.S(3)/2)*sqrt(f + g*x)/(7*c) - 2*sqrt(a + c*x**2)*sqrt(f + g*x)*(5*a*e*g**2 + c*f*(-7*d*g + 4*e*f) - 3*c*g*x*(7*d*g + e*f))/(105*c*g**2) - 4*sqrt(-a)*sqrt(1 + c*x**2/a)*sqrt(f + g*x)*(a*g**2*(21*d*g + 8*e*f) + c*f**2*(-7*d*g + 4*e*f))*elliptic_e(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(105*sqrt(c)*g**3*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(a + c*x**2)) + 4*sqrt(-a)*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(1 + c*x**2/a)*(a*g**2 + c*f**2)*(5*a*e*g**2 + c*f*(-7*d*g + 4*e*f))*elliptic_f(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(105*c**(sympy.S(3)/2)*g**3*sqrt(a + c*x**2)*sqrt(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_619():
    f = sqrt(a + c*x**2)*sqrt(f + g*x)
    F = -4*f*sqrt(a + c*x**2)*sqrt(f + g*x)/(15*g) + 2*sqrt(a + c*x**2)*(f + g*x)**(sympy.S(3)/2)/(5*g) - 4*f*sqrt(-a)*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(1 + c*x**2/a)*(a*g**2 + c*f**2)*elliptic_f(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(15*sqrt(c)*g**2*sqrt(a + c*x**2)*sqrt(f + g*x)) + 4*sqrt(-a)*sqrt(1 + c*x**2/a)*sqrt(f + g*x)*(-3*a*g**2 + c*f**2)*elliptic_e(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(15*sqrt(c)*g**2*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(a + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_620():
    f = sqrt(a + c*x**2)*sqrt(f + g*x)/(d + e*x)
    F = ((Integer(2) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(3) * Symbol('e')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('e'))**(Integer(2)) * Symbol('g') * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g')))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('e'))**(Integer(2)) * Symbol('g') * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * ((Integer(2) * Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Integer(-1) * (Integer(3) * Symbol('c') * Symbol('d') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g'))))))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((Integer(3) * sympy.sqrt(Symbol('c')) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('e')) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g')) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1))))) * (((Symbol('e'))**(Integer(3)) * (((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_621():
    f = sqrt(a + c*x**2)*sqrt(f + g*x)/(d + e*x)**2
    F = (Integer(-1) * ((sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))) * ((Symbol('e') * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * (((Symbol('e'))**(Integer(2)) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((Integer(3) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f') * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * (((Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * ((Integer(2) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g')))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * (((Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g'))))))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('e')) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g')) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1))))) * (((Symbol('e'))**(Integer(3)) * (((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_622():
    f = sqrt(a + c*x**2)*sqrt(f + g*x)/(d + e*x)**3
    F = (Integer(-1) * ((sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g'))))))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(4) * Symbol('e') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * ((Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g'))))))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('e'))**(Integer(2)) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('g') * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f') * ((Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g'))))))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('e'))**(Integer(2)) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('d') * Symbol('g') * ((Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g'))))))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('e'))**(Integer(3)) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('c') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g')))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('e')) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g')) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1))))) * (((Symbol('e'))**(Integer(3)) * (((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + (((((Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g'))))))))**(Integer(2)) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('e')) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g')) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('e'))**(Integer(3)) * (((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_623():
    f = sqrt(a + c*x**2)*(d + e*x)**3/sqrt(f + g*x)
    F = -4*e**2*sqrt(a + c*x**2)*(f + g*x)**(sympy.S(5)/2)*(-3*d*g + 4*e*f)/(63*g**4) + 2*sqrt(a + c*x**2)*(d + e*x)**3*sqrt(f + g*x)/(9*g) + 4*e*sqrt(a + c*x**2)*(f + g*x)**(sympy.S(3)/2)*(7*a*e**2*g**2 + c*(42*d**2*g**2 - 111*d*e*f*g + 64*e**2*f**2))/(315*c*g**4) + sqrt(a + c*x**2)*sqrt(f + g*x)*(-36*a*e**2*g**2*(-5*d*g + 2*e*f) - 4*c*(-35*d**3*g**3 + 168*d**2*e*f*g**2 - 204*d*e**2*f**2*g + 76*e**3*f**3))/(315*c*g**4) - 4*sqrt(-a)*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(1 + c*x**2/a)*(a*g**2 + c*f**2)*(9*a*e**2*g**2*(-5*d*g + 2*e*f) - c*(-105*d**3*g**3 + 252*d**2*e*f*g**2 - 216*d*e**2*f**2*g + 64*e**3*f**3))*elliptic_f(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(315*c**(sympy.S(3)/2)*g**5*sqrt(a + c*x**2)*sqrt(f + g*x)) + 4*sqrt(-a)*sqrt(1 + c*x**2/a)*sqrt(f + g*x)*(21*a**2*e**3*g**4 - 3*a*c*e*g**2*(63*d**2*g**2 - 39*d*e*f*g + 10*e**2*f**2) - c**2*f*(-105*d**3*g**3 + 252*d**2*e*f*g**2 - 216*d*e**2*f**2*g + 64*e**3*f**3))*elliptic_e(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(315*c**(sympy.S(3)/2)*g**5*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(a + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_624():
    f = sqrt(a + c*x**2)*(d + e*x)**2/sqrt(f + g*x)
    F = -4*e*sqrt(a + c*x**2)*(f + g*x)**(sympy.S(3)/2)*(-2*d*g + 3*e*f)/(35*g**3) + 2*sqrt(a + c*x**2)*(d + e*x)**2*sqrt(f + g*x)/(7*g) + sqrt(a + c*x**2)*sqrt(f + g*x)*(20*a*e**2*g**2 + 4*c*(10*d**2*g**2 - 34*d*e*f*g + 21*e**2*f**2))/(105*c*g**3) + 4*sqrt(-a)*sqrt(1 + c*x**2/a)*sqrt(f + g*x)*(a*e*g**2*(-42*d*g + 13*e*f) + c*f*(35*d**2*g**2 - 56*d*e*f*g + 24*e**2*f**2))*elliptic_e(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(105*sqrt(c)*g**4*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(a + c*x**2)) + 4*sqrt(-a)*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(1 + c*x**2/a)*(a*g**2 + c*f**2)*(5*a*e**2*g**2 - c*(35*d**2*g**2 - 56*d*e*f*g + 24*e**2*f**2))*elliptic_f(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(105*c**(sympy.S(3)/2)*g**4*sqrt(a + c*x**2)*sqrt(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_625():
    f = sqrt(a + c*x**2)*(d + e*x)/sqrt(f + g*x)
    F = -2*sqrt(a + c*x**2)*sqrt(f + g*x)*(-5*d*g + 4*e*f - 3*e*g*x)/(15*g**2) + 4*sqrt(-a)*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(1 + c*x**2/a)*(a*g**2 + c*f**2)*(-5*d*g + 4*e*f)*elliptic_f(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(15*sqrt(c)*g**3*sqrt(a + c*x**2)*sqrt(f + g*x)) - 4*sqrt(-a)*sqrt(1 + c*x**2/a)*sqrt(f + g*x)*(3*a*e*g**2 + c*f*(-5*d*g + 4*e*f))*elliptic_e(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(15*sqrt(c)*g**3*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(a + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_626():
    f = sqrt(a + c*x**2)/sqrt(f + g*x)
    F = 4*sqrt(c)*f*sqrt(-a)*sqrt(1 + c*x**2/a)*sqrt(f + g*x)*elliptic_e(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(3*g**2*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(a + c*x**2)) + 2*sqrt(a + c*x**2)*sqrt(f + g*x)/(3*g) - 4*sqrt(-a)*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(1 + c*x**2/a)*(a*g**2 + c*f**2)*elliptic_f(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(3*sqrt(c)*g**2*sqrt(a + c*x**2)*sqrt(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_627():
    f = sqrt(a + c*x**2)/((d + e*x)*sqrt(f + g*x))
    F = ((Integer(-2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1))))) * ((Symbol('e') * Symbol('g') * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * ((Symbol('e') * Symbol('f')) + (Symbol('d') * Symbol('g'))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1))))) * (((Symbol('e'))**(Integer(2)) * Symbol('g') * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('e')) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g')) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1))))) * (((Symbol('e'))**(Integer(2)) * (((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_628():
    f = sqrt(a + c*x**2)/((d + e*x)**2*sqrt(f + g*x))
    F = (Integer(-1) * ((sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))) * ((((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((Symbol('e') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f') * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((Symbol('e') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * ((Integer(2) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * (((Symbol('e'))**(Integer(2)) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((((Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('e')) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g')) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1))))) * (((Symbol('e'))**(Integer(2)) * (((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_629():
    f = sqrt(a + c*x**2)/((d + e*x)**3*sqrt(f + g*x))
    F = (Integer(-1) * ((sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(2) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + ((((Integer(3) * Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Symbol('d') * Symbol('g'))))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(4) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(2)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1))) + ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * ((Integer(3) * Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Symbol('d') * Symbol('g'))))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((Integer(4) * Symbol('e') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(2)) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('g') * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**(Integer(2)) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f') * ((Integer(3) * Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Symbol('d') * Symbol('g'))))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((Integer(4) * Symbol('e') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(2)) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('d') * Symbol('g') * ((Integer(3) * Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Symbol('d') * Symbol('g'))))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('e'))**(Integer(2)) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(2)) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * ((Symbol('e') * Symbol('f')) + (Symbol('d') * Symbol('g'))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('e')) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g')) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1))))) * (((Symbol('e'))**(Integer(2)) * (((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g'))))))) * ((Integer(3) * Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Symbol('d') * Symbol('g'))))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('e')) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g')) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('e'))**(Integer(2)) * (((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(2)) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_630():
    f = (d + e*x)**3*sqrt(f + g*x)/sqrt(a + c*x**2)
    F = 2*e**2*sqrt(a + c*x**2)*(f + g*x)**(sympy.S(3)/2)*(11*d*g + e*f)/(35*c*g**2) + 2*e*sqrt(a + c*x**2)*(d + e*x)**2*sqrt(f + g*x)/(7*c) - 2*e*sqrt(a + c*x**2)*sqrt(f + g*x)*(25*a*e**2*g**2 + c*(-90*d**2*g**2 + 12*d*e*f*g + 7*e**2*f**2))/(105*c**2*g**2) + 2*sqrt(-a)*sqrt(1 + c*x**2/a)*sqrt(f + g*x)*(a*e**2*g**2*(189*d*g + 19*e*f) - c*(105*d**3*g**3 + 105*d**2*e*f*g**2 - 42*d*e**2*f**2*g + 8*e**3*f**3))*elliptic_e(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(105*c**(sympy.S(3)/2)*g**3*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(a + c*x**2)) - 2*e*sqrt(-a)*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(1 + c*x**2/a)*(a*g**2 + c*f**2)*(25*a*e**2*g**2 - c*(105*d**2*g**2 - 42*d*e*f*g + 8*e**2*f**2))*elliptic_f(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(105*c**(sympy.S(5)/2)*g**3*sqrt(a + c*x**2)*sqrt(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_631():
    f = (d + e*x)**2*sqrt(f + g*x)/sqrt(a + c*x**2)
    F = 2*e*sqrt(a + c*x**2)*(d + e*x)*sqrt(f + g*x)/(5*c) + 2*e*sqrt(a + c*x**2)*sqrt(f + g*x)*(7*d*g + e*f)/(15*c*g) - 4*e*sqrt(-a)*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(1 + c*x**2/a)*(a*g**2 + c*f**2)*(-5*d*g + e*f)*elliptic_f(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(15*c**(sympy.S(3)/2)*g**2*sqrt(a + c*x**2)*sqrt(f + g*x)) + 2*sqrt(-a)*sqrt(1 + c*x**2/a)*sqrt(f + g*x)*(9*a*e**2*g**2 + c*(-15*d**2*g**2 - 10*d*e*f*g + 2*e**2*f**2))*elliptic_e(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(15*c**(sympy.S(3)/2)*g**2*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(a + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_632():
    f = (d + e*x)*sqrt(f + g*x)/sqrt(a + c*x**2)
    F = 2*e*sqrt(a + c*x**2)*sqrt(f + g*x)/(3*c) - 2*sqrt(-a)*sqrt(1 + c*x**2/a)*sqrt(f + g*x)*(3*d*g + e*f)*elliptic_e(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(3*sqrt(c)*g*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(a + c*x**2)) + 2*e*sqrt(-a)*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(1 + c*x**2/a)*(a*g**2 + c*f**2)*elliptic_f(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(3*c**(sympy.S(3)/2)*g*sqrt(a + c*x**2)*sqrt(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_633():
    f = sqrt(f + g*x)/sqrt(a + c*x**2)
    F = -2*sqrt(-a)*sqrt(1 + c*x**2/a)*sqrt(f + g*x)*elliptic_e(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(sqrt(c)*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(a + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_634():
    f = sqrt(f + g*x)/(sqrt(a + c*x**2)*(d + e*x))
    F = ((Integer(-2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g') * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('c')) * Symbol('e') * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('e')) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g')) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1))))) * ((Symbol('e') * (((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_635():
    f = sqrt(f + g*x)/(sqrt(a + c*x**2)*(d + e*x)**2)
    F = (Integer(-1) * ((Symbol('e') * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))) * ((((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1))))) * ((((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f') * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1))))) * ((((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('d') * Symbol('g') * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1))))) * ((Symbol('e') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('e')) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g')) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1))))) * ((Symbol('e') * (((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_636():
    f = sqrt(f + g*x)/(sqrt(a + c*x**2)*(d + e*x)**3)
    F = ((Integer(-1) * (Symbol('e') * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2))))))) * ((Integer(2) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * ((Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Symbol('c') * Symbol('d') * ((Integer(6) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(5) * Symbol('d') * Symbol('g')))))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(4) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * ((Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Symbol('c') * Symbol('d') * ((Integer(6) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(5) * Symbol('d') * Symbol('g')))))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1))))) * ((Integer(4) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('g') * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1))))) * ((Integer(2) * Symbol('e') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f') * ((Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Symbol('c') * Symbol('d') * ((Integer(6) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(5) * Symbol('d') * Symbol('g')))))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1))))) * ((Integer(4) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('d') * Symbol('g') * ((Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Symbol('c') * Symbol('d') * ((Integer(6) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(5) * Symbol('d') * Symbol('g')))))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1))))) * ((Integer(4) * Symbol('e') * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((Symbol('c') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g')))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('e')) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g')) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1))))) * ((Symbol('e') * (((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + ((((Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Symbol('c') * Symbol('d') * ((Integer(6) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(5) * Symbol('d') * Symbol('g')))))) * ((Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g'))))))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('e')) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g')) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1))))) * ((Integer(4) * Symbol('e') * (((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_637():
    f = (f + g*x)**(sympy.S(3)/2)/(sqrt(a + c*x**2)*(d + e*x))
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g') * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((sympy.sqrt(Symbol('c')) * Symbol('e') * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((sympy.sqrt(Symbol('c')) * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(2)) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('e')) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g')) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1))))) * (((Symbol('e'))**(Integer(2)) * (((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_638():
    f = (d + e*x)**3/(sqrt(a + c*x**2)*sqrt(f + g*x))
    F = 2*e**2*sqrt(a + c*x**2)*(d + e*x)*sqrt(f + g*x)/(5*c*g) - 8*e**2*sqrt(a + c*x**2)*sqrt(f + g*x)*(-3*d*g + e*f)/(15*c*g**2) + 2*e*sqrt(-a)*sqrt(1 + c*x**2/a)*sqrt(f + g*x)*(9*a*e**2*g**2 - c*(45*d**2*g**2 - 30*d*e*f*g + 8*e**2*f**2))*elliptic_e(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(15*c**(sympy.S(3)/2)*g**3*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(a + c*x**2)) - 2*sqrt(-a)*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(1 + c*x**2/a)*(a*e**2*g**2*(-15*d*g + 7*e*f) - c*(-15*d**3*g**3 + 45*d**2*e*f*g**2 - 30*d*e**2*f**2*g + 8*e**3*f**3))*elliptic_f(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(15*c**(sympy.S(3)/2)*g**3*sqrt(a + c*x**2)*sqrt(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_639():
    f = (d + e*x)**2/(sqrt(a + c*x**2)*sqrt(f + g*x))
    F = 2*e**2*sqrt(a + c*x**2)*sqrt(f + g*x)/(3*c*g) + 4*e*sqrt(-a)*sqrt(1 + c*x**2/a)*sqrt(f + g*x)*(-3*d*g + e*f)*elliptic_e(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(3*sqrt(c)*g**2*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(a + c*x**2)) - 2*sqrt(-a)*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(1 + c*x**2/a)*(2*c*e*f*(-3*d*g + e*f) + g**2*(-a*e**2 + 3*c*d**2))*elliptic_f(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(3*c**(sympy.S(3)/2)*g**2*sqrt(a + c*x**2)*sqrt(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_640():
    f = (d + e*x)/(sqrt(a + c*x**2)*sqrt(f + g*x))
    F = -2*e*sqrt(-a)*sqrt(1 + c*x**2/a)*sqrt(f + g*x)*elliptic_e(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(sqrt(c)*g*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(a + c*x**2)) + 2*sqrt(-a)*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(1 + c*x**2/a)*(-d*g + e*f)*elliptic_f(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(sqrt(c)*g*sqrt(a + c*x**2)*sqrt(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_641():
    f = 1/(sqrt(a + c*x**2)*sqrt(f + g*x))
    F = -2*sqrt(-a)*sqrt(sqrt(c)*(f + g*x)/(sqrt(c)*f + g*sqrt(-a)))*sqrt(1 + c*x**2/a)*elliptic_f(asin(sqrt(2)*sqrt(-sqrt(c)*x/sqrt(-a) + 1)/2), -2*a*g/(-a*g + sqrt(c)*f*sqrt(-a)))/(sqrt(c)*sqrt(a + c*x**2)*sqrt(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_642():
    f = 1/(sqrt(a + c*x**2)*(d + e*x)*sqrt(f + g*x))
    F = (Integer(-2) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('e')) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g')) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1))))) * (((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_643():
    f = 1/(sqrt(a + c*x**2)*(d + e*x)**2*sqrt(f + g*x))
    F = (Integer(-1) * (((Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))) * ((((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('e') * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1))))) * ((((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('e') * Symbol('f') * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1))))) * ((((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('d') * Symbol('g') * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1))))) * ((((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((((Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g'))))))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('e')) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g')) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1))))) * (((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_644():
    f = 1/(sqrt(a + c*x**2)*(d + e*x)**3*sqrt(f + g*x))
    F = (Integer(-1) * (((Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(2) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (Symbol('e'))**(Integer(2)) * ((Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g'))))))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(4) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(2)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1))) + ((Integer(3) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('e') * ((Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g'))))))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((Integer(4) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(2)) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + ((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('g') * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((Integer(2) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('e') * Symbol('f') * ((Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g'))))))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((Integer(4) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(2)) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((Integer(3) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('d') * Symbol('g') * ((Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g'))))))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((Integer(4) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(2)) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + ((Symbol('c') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g')))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('e')) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g')) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1))))) * (((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (((Symbol('a') * (Symbol('e'))**(Integer(2)) * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g'))))))))**(Integer(2)) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('e')) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g')) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1))))) * ((Integer(4) * (((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(2)) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_645():
    f = 1/(sqrt(a + c*x**2)*(d + e*x)*(f + g*x)**(sympy.S(3)/2))
    F = ((Integer(2) * (Symbol('g'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))) * ((((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * ((Symbol('c') * (Symbol('f'))**(Integer(2))) + (Symbol('a') * (Symbol('g'))**(Integer(2)))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x)))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('g') * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * ((Symbol('c') * (Symbol('f'))**(Integer(2))) + (Symbol('a') * (Symbol('g'))**(Integer(2)))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('e') * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('e')) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g')) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1))))) * (((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_646():
    f = 1/(sqrt(a + c*x**2)*(d + e*x)*(f + g*x)**(sympy.S(5)/2))
    F = ((Integer(2) * (Symbol('g'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(3) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * ((Symbol('c') * (Symbol('f'))**(Integer(2))) + (Symbol('a') * (Symbol('g'))**(Integer(2)))) * ((Symbol('f') + (Symbol('g') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(8) * Symbol('c') * Symbol('f') * (Symbol('g'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(3) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * (((Symbol('c') * (Symbol('f'))**(Integer(2))) + (Symbol('a') * (Symbol('g'))**(Integer(2)))))**(Integer(2)) * sympy.sqrt((Symbol('f') + (Symbol('g') * x)))))**(Integer(-1))) + ((Integer(2) * Symbol('e') * (Symbol('g'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))) * (((((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(2)) * ((Symbol('c') * (Symbol('f'))**(Integer(2))) + (Symbol('a') * (Symbol('g'))**(Integer(2)))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x)))))**(Integer(-1))) + ((Integer(8) * sympy.sqrt((Integer(-1) * Symbol('a'))) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f') * Symbol('g') * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((Integer(3) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * (((Symbol('c') * (Symbol('f'))**(Integer(2))) + (Symbol('a') * (Symbol('g'))**(Integer(2)))))**(Integer(2)) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('e') * Symbol('g') * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * (((((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(2)) * ((Symbol('c') * (Symbol('f'))**(Integer(2))) + (Symbol('a') * (Symbol('g'))**(Integer(2)))) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('g') * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('c')) * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g')))))**(Integer(-1)))))) * ((Integer(3) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * ((Symbol('c') * (Symbol('f'))**(Integer(2))) + (Symbol('a') * (Symbol('g'))**(Integer(2)))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('e'))**(Integer(2)) * sympy.sqrt(((sympy.sqrt(Symbol('c')) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('e')) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g')) * (((sympy.sqrt(Symbol('c')) * Symbol('f')) + (sympy.sqrt((Integer(-1) * Symbol('a'))) * Symbol('g'))))**(Integer(-1))))) * (((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('a'))))**(Integer(-1))) + Symbol('e')) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(2)) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_647():
    f = 1/((d + e*x)*sqrt(f + g*x)*sqrt(c*x**2 + 1))
    F = Integer(-1) * ((Integer(2) * sympy.sqrt(((sympy.sqrt((Integer(-1) * Symbol('c'))) * (Symbol('f') + (Symbol('g') * x))) * (((sympy.sqrt((Integer(-1) * Symbol('c'))) * Symbol('f')) + Symbol('g')))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('e')) * (((sympy.sqrt((Integer(-1) * Symbol('c'))) * Symbol('d')) + Symbol('e')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('c'))) * x)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('g')) * (((sympy.sqrt((Integer(-1) * Symbol('c'))) * Symbol('f')) + Symbol('g')))**(Integer(-1))))) * ((((sympy.sqrt((Integer(-1) * Symbol('c'))) * Symbol('d')) + Symbol('e')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_648():
    f = 1/(sqrt(a + c*x**2)*sqrt(d + e*x)*sqrt(f + g*x))
    F = -sqrt((1 - (f + g*x)*(2*a*e*g + 2*c*d*f)/((d + e*x)*(a*g**2 + c*f**2)) + (f + g*x)**2*(a*e**2 + c*d**2)/((d + e*x)**2*(a*g**2 + c*f**2)))/(1 + (f + g*x)*sqrt(a*e**2 + c*d**2)/((d + e*x)*sqrt(a*g**2 + c*f**2)))**2)*sqrt((a + c*x**2)*(-d*g + e*f)**2/((d + e*x)**2*(a*g**2 + c*f**2)))*(1 + (f + g*x)*sqrt(a*e**2 + c*d**2)/((d + e*x)*sqrt(a*g**2 + c*f**2)))*(d + e*x)*(a*g**2 + c*f**2)**(sympy.S(1)/4)*elliptic_f(2*atan(sqrt(f + g*x)*(a*e**2 + c*d**2)**(sympy.S(1)/4)/(sqrt(d + e*x)*(a*g**2 + c*f**2)**(sympy.S(1)/4))), sympy.S.Half + (a*e*g + c*d*f)/(2*sqrt(a*e**2 + c*d**2)*sqrt(a*g**2 + c*f**2)))/(sqrt(a + c*x**2)*(a*e**2 + c*d**2)**(sympy.S(1)/4)*(-d*g + e*f)*sqrt(1 - (f + g*x)*(2*a*e*g + 2*c*d*f)/((d + e*x)*(a*g**2 + c*f**2)) + (f + g*x)**2*(a*e**2 + c*d**2)/((d + e*x)**2*(a*g**2 + c*f**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_649():
    f = 1/(sqrt(x - 1)*sqrt(x + 1)*sqrt(2*x**2 - 1))
    F = sqrt(1 - 2*x**2)*sqrt(1 - x**2)*elliptic_f(asin(x), 2)/(sqrt(x - 1)*sqrt(x + 1)*sqrt(2*x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_650():
    f = sqrt(d + e*x)*(f + g*x)**3/sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))
    F = 2*(f + g*x)**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(7*c*d*sqrt(d + e*x)) + (f + g*x)**2*(-12*a*e*g + 12*c*d*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(35*c**2*d**2*sqrt(d + e*x)) + 16*g*sqrt(d + e*x)*(-a*e*g + c*d*f)**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(35*c**3*d**3*e) - 16*(-a*e*g + c*d*f)**2*(2*a*e**2*g - c*d*(-d*g + 3*e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(35*c**4*d**4*e*sqrt(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_651():
    f = sqrt(d + e*x)*(f + g*x)**2/sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))
    F = 2*(f + g*x)**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(5*c*d*sqrt(d + e*x)) + 8*g*sqrt(d + e*x)*(-a*e*g + c*d*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(15*c**2*d**2*e) - (-8*a*e*g + 8*c*d*f)*(2*a*e**2*g - c*d*(-d*g + 3*e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(15*c**3*d**3*e*sqrt(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_652():
    f = sqrt(d + e*x)*(f + g*x)/sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))
    F = 2*g*sqrt(d + e*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3*c*d*e) - (4*a*e**2*g - 2*c*d*(-d*g + 3*e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3*c**2*d**2*e*sqrt(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_653():
    f = sqrt(d + e*x)/sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))
    F = 2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(c*d*sqrt(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_654():
    f = sqrt(d + e*x)/((f + g*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    F = 2*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(sqrt(g)*sqrt(-a*e*g + c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_655():
    f = sqrt(d + e*x)/((f + g*x)**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    F = c*d*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(sqrt(g)*(-a*e*g + c*d*f)**(sympy.S(3)/2)) + sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*(f + g*x)*(-a*e*g + c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_656():
    f = sqrt(d + e*x)/((f + g*x)**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    F = 3*c**2*d**2*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(4*sqrt(g)*(-a*e*g + c*d*f)**(sympy.S(5)/2)) + 3*c*d*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*sqrt(d + e*x)*(f + g*x)*(-a*e*g + c*d*f)**2) + sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*(f + g*x)**2*(-2*a*e*g + 2*c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_657():
    f = sqrt(d + e*x)/((f + g*x)**4*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    F = 5*c**3*d**3*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(8*sqrt(g)*(-a*e*g + c*d*f)**(sympy.S(7)/2)) + 5*c**2*d**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(8*sqrt(d + e*x)*(f + g*x)*(-a*e*g + c*d*f)**3) + 5*c*d*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(12*sqrt(d + e*x)*(f + g*x)**2*(-a*e*g + c*d*f)**2) + sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*(f + g*x)**3*(-3*a*e*g + 3*c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_658():
    f = (d + e*x)**(sympy.S(3)/2)*(f + g*x)**3/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)
    F = -2*sqrt(d + e*x)*(f + g*x)**3/(c*d*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) + 12*g*(f + g*x)**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(5*c**2*d**2*sqrt(d + e*x)) + 16*g**2*sqrt(d + e*x)*(-a*e*g + c*d*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(5*c**3*d**3*e) - 16*g*(-a*e*g + c*d*f)*(2*a*e**2*g - c*d*(-d*g + 3*e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(5*c**4*d**4*e*sqrt(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_659():
    f = (d + e*x)**(sympy.S(3)/2)*(f + g*x)**2/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)
    F = -2*sqrt(d + e*x)*(f + g*x)**2/(c*d*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) + 8*g**2*sqrt(d + e*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3*c**2*d**2*e) - 8*g*(2*a*e**2*g - c*d*(-d*g + 3*e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3*c**3*d**3*e*sqrt(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_660():
    f = (d + e*x)**(sympy.S(3)/2)*(f + g*x)/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)
    F = -(d + e*x)**(sympy.S(3)/2)*(-2*a*e*g + 2*c*d*f)/(c*d*(-a*e**2 + c*d**2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) - (4*a*e**2*g - 2*c*d*(d*g + e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(c**2*d**2*sqrt(d + e*x)*(-a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_661():
    f = (d + e*x)**(sympy.S(3)/2)/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)
    F = -2*sqrt(d + e*x)/(c*d*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_662():
    f = (d + e*x)**(sympy.S(3)/2)/((f + g*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    F = -2*sqrt(g)*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(-a*e*g + c*d*f)**(sympy.S(3)/2) - 2*sqrt(d + e*x)/((-a*e*g + c*d*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_663():
    f = (d + e*x)**(sympy.S(3)/2)/((f + g*x)**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    F = -3*c*d*sqrt(g)*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(-a*e*g + c*d*f)**(sympy.S(5)/2) - 3*g*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*(f + g*x)*(-a*e*g + c*d*f)**2) - 2*sqrt(d + e*x)/((f + g*x)*(-a*e*g + c*d*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_664():
    f = (d + e*x)**(sympy.S(3)/2)/((f + g*x)**3*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    F = -15*c**2*d**2*sqrt(g)*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(4*(-a*e*g + c*d*f)**(sympy.S(7)/2)) - 15*c*d*g*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*sqrt(d + e*x)*(f + g*x)*(-a*e*g + c*d*f)**3) - 5*g*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(2*sqrt(d + e*x)*(f + g*x)**2*(-a*e*g + c*d*f)**2) - 2*sqrt(d + e*x)/((f + g*x)**2*(-a*e*g + c*d*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_665():
    f = (d + e*x)**(sympy.S(5)/2)*(f + g*x)**3/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)
    F = -2*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**3/(3*c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)) - 4*g*sqrt(d + e*x)*(f + g*x)**2/(c**2*d**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) + 16*g**3*sqrt(d + e*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3*c**3*d**3*e) - 16*g**2*(2*a*e**2*g - c*d*(-d*g + 3*e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3*c**4*d**4*e*sqrt(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_666():
    f = (d + e*x)**(sympy.S(5)/2)*(f + g*x)**2/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)
    F = -2*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**2/(3*c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)) - 8*g*(d + e*x)**(sympy.S(3)/2)*(-a*e*g + c*d*f)/(3*c**2*d**2*(-a*e**2 + c*d**2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) - 8*g*(2*a*e**2*g - c*d*(d*g + e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3*c**3*d**3*sqrt(d + e*x)*(-a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_667():
    f = (d + e*x)**(sympy.S(5)/2)*(f + g*x)/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)
    F = -(d + e*x)**(sympy.S(5)/2)*(-2*a*e*g + 2*c*d*f)/(3*c*d*(-a*e**2 + c*d**2)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)) + sqrt(d + e*x)*(4*a*e**2*g + 2*c*d*(-3*d*g + e*f))/(3*c**2*d**2*(-a*e**2 + c*d**2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_668():
    f = (d + e*x)**(sympy.S(5)/2)/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)
    F = -2*(d + e*x)**(sympy.S(3)/2)/(3*c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_669():
    f = (d + e*x)**(sympy.S(5)/2)/((f + g*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2))
    F = 2*g**(sympy.S(3)/2)*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(-a*e*g + c*d*f)**(sympy.S(5)/2) + 2*g*sqrt(d + e*x)/((-a*e*g + c*d*f)**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) - 2*(d + e*x)**(sympy.S(3)/2)/((-3*a*e*g + 3*c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_670():
    f = (d + e*x)**(sympy.S(5)/2)/((f + g*x)**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2))
    F = 5*c*d*g**(sympy.S(3)/2)*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(-a*e*g + c*d*f)**(sympy.S(7)/2) + 5*g**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*(f + g*x)*(-a*e*g + c*d*f)**3) + 10*g*sqrt(d + e*x)/(3*(f + g*x)*(-a*e*g + c*d*f)**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) - 2*(d + e*x)**(sympy.S(3)/2)/((f + g*x)*(-3*a*e*g + 3*c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_671():
    f = (d + e*x)**(sympy.S(5)/2)/((f + g*x)**3*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2))
    F = 35*c**2*d**2*g**(sympy.S(3)/2)*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(4*(-a*e*g + c*d*f)**(sympy.S(9)/2)) + 35*c*d*g**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*sqrt(d + e*x)*(f + g*x)*(-a*e*g + c*d*f)**4) + 35*g**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(6*sqrt(d + e*x)*(f + g*x)**2*(-a*e*g + c*d*f)**3) + 14*g*sqrt(d + e*x)/(3*(f + g*x)**2*(-a*e*g + c*d*f)**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) - 2*(d + e*x)**(sympy.S(3)/2)/((f + g*x)**2*(-3*a*e*g + 3*c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_672():
    f = (f + g*x)**4*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/sqrt(d + e*x)
    F = 2*(f + g*x)**4*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(11*c*d*(d + e*x)**(sympy.S(3)/2)) + (f + g*x)**3*(-16*a*e*g + 16*c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(99*c**2*d**2*(d + e*x)**(sympy.S(3)/2)) + 32*(f + g*x)**2*(-a*e*g + c*d*f)**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(231*c**3*d**3*(d + e*x)**(sympy.S(3)/2)) + 128*g*(-a*e*g + c*d*f)**3*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(1155*c**4*d**4*e*sqrt(d + e*x)) - 128*(-a*e*g + c*d*f)**3*(2*a*e**2*g - c*d*(-3*d*g + 5*e*f))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(3465*c**5*d**5*e*(d + e*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_673():
    f = (f + g*x)**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/sqrt(d + e*x)
    F = 2*(f + g*x)**3*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(9*c*d*(d + e*x)**(sympy.S(3)/2)) + (f + g*x)**2*(-4*a*e*g + 4*c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(21*c**2*d**2*(d + e*x)**(sympy.S(3)/2)) + 16*g*(-a*e*g + c*d*f)**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(105*c**3*d**3*e*sqrt(d + e*x)) - 16*(-a*e*g + c*d*f)**2*(2*a*e**2*g - c*d*(-3*d*g + 5*e*f))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(315*c**4*d**4*e*(d + e*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_674():
    f = (f + g*x)**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/sqrt(d + e*x)
    F = 2*(f + g*x)**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(7*c*d*(d + e*x)**(sympy.S(3)/2)) + 8*g*(-a*e*g + c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(35*c**2*d**2*e*sqrt(d + e*x)) - (-8*a*e*g + 8*c*d*f)*(2*a*e**2*g - c*d*(-3*d*g + 5*e*f))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(105*c**3*d**3*e*(d + e*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_675():
    f = (f + g*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/sqrt(d + e*x)
    F = 2*g*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(5*c*d*e*sqrt(d + e*x)) - (4*a*e**2*g - 2*c*d*(-3*d*g + 5*e*f))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(15*c**2*d**2*e*(d + e*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_676():
    f = sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/sqrt(d + e*x)
    F = 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(3*c*d*(d + e*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_677():
    f = sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*(f + g*x))
    F = 2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(g*sqrt(d + e*x)) - 2*sqrt(-a*e*g + c*d*f)*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/g**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_678():
    f = sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*(f + g*x)**2)
    F = c*d*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(g**(sympy.S(3)/2)*sqrt(-a*e*g + c*d*f)) - sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(g*sqrt(d + e*x)*(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_679():
    f = sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*(f + g*x)**3)
    F = c**2*d**2*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(4*g**(sympy.S(3)/2)*(-a*e*g + c*d*f)**(sympy.S(3)/2)) + c*d*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*g*sqrt(d + e*x)*(f + g*x)*(-a*e*g + c*d*f)) - sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(2*g*sqrt(d + e*x)*(f + g*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_680():
    f = sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*(f + g*x)**4)
    F = c**3*d**3*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(8*g**(sympy.S(3)/2)*(-a*e*g + c*d*f)**(sympy.S(5)/2)) + c**2*d**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(8*g*sqrt(d + e*x)*(f + g*x)*(-a*e*g + c*d*f)**2) + c*d*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(12*g*sqrt(d + e*x)*(f + g*x)**2*(-a*e*g + c*d*f)) - sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3*g*sqrt(d + e*x)*(f + g*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_681():
    f = sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*(f + g*x)**5)
    F = 5*c**4*d**4*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(64*g**(sympy.S(3)/2)*(-a*e*g + c*d*f)**(sympy.S(7)/2)) + 5*c**3*d**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(64*g*sqrt(d + e*x)*(f + g*x)*(-a*e*g + c*d*f)**3) + 5*c**2*d**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(96*g*sqrt(d + e*x)*(f + g*x)**2*(-a*e*g + c*d*f)**2) + c*d*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(24*g*sqrt(d + e*x)*(f + g*x)**3*(-a*e*g + c*d*f)) - sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*g*sqrt(d + e*x)*(f + g*x)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_682():
    f = (f + g*x)**4*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(d + e*x)**(sympy.S(3)/2)
    F = 2*(f + g*x)**4*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(13*c*d*(d + e*x)**(sympy.S(5)/2)) + (f + g*x)**3*(-16*a*e*g + 16*c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(143*c**2*d**2*(d + e*x)**(sympy.S(5)/2)) + 32*(f + g*x)**2*(-a*e*g + c*d*f)**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(429*c**3*d**3*(d + e*x)**(sympy.S(5)/2)) + 128*g*(-a*e*g + c*d*f)**3*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(3003*c**4*d**4*e*(d + e*x)**(sympy.S(3)/2)) - 128*(-a*e*g + c*d*f)**3*(2*a*e**2*g - c*d*(-5*d*g + 7*e*f))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(15015*c**5*d**5*e*(d + e*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_683():
    f = (f + g*x)**3*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(d + e*x)**(sympy.S(3)/2)
    F = 2*(f + g*x)**3*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(11*c*d*(d + e*x)**(sympy.S(5)/2)) + (f + g*x)**2*(-4*a*e*g + 4*c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(33*c**2*d**2*(d + e*x)**(sympy.S(5)/2)) + 16*g*(-a*e*g + c*d*f)**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(231*c**3*d**3*e*(d + e*x)**(sympy.S(3)/2)) - 16*(-a*e*g + c*d*f)**2*(2*a*e**2*g - c*d*(-5*d*g + 7*e*f))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(1155*c**4*d**4*e*(d + e*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_684():
    f = (f + g*x)**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(d + e*x)**(sympy.S(3)/2)
    F = 2*(f + g*x)**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(9*c*d*(d + e*x)**(sympy.S(5)/2)) + 8*g*(-a*e*g + c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(63*c**2*d**2*e*(d + e*x)**(sympy.S(3)/2)) - (-8*a*e*g + 8*c*d*f)*(2*a*e**2*g - c*d*(-5*d*g + 7*e*f))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(315*c**3*d**3*e*(d + e*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_685():
    f = (f + g*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(d + e*x)**(sympy.S(3)/2)
    F = 2*g*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(7*c*d*e*(d + e*x)**(sympy.S(3)/2)) - (4*a*e**2*g - 2*c*d*(-5*d*g + 7*e*f))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(35*c**2*d**2*e*(d + e*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_686():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(d + e*x)**(sympy.S(3)/2)
    F = 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(5*c*d*(d + e*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_687():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/((d + e*x)**(sympy.S(3)/2)*(f + g*x))
    F = 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(3*g*(d + e*x)**(sympy.S(3)/2)) - (-2*a*e*g + 2*c*d*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(g**2*sqrt(d + e*x)) + 2*(-a*e*g + c*d*f)**(sympy.S(3)/2)*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/g**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_688():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/((d + e*x)**(sympy.S(3)/2)*(f + g*x)**2)
    F = 3*c*d*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(g**2*sqrt(d + e*x)) - 3*c*d*sqrt(-a*e*g + c*d*f)*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/g**(sympy.S(5)/2) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(g*(d + e*x)**(sympy.S(3)/2)*(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_689():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/((d + e*x)**(sympy.S(3)/2)*(f + g*x)**3)
    F = 3*c**2*d**2*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(4*g**(sympy.S(5)/2)*sqrt(-a*e*g + c*d*f)) - 3*c*d*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*g**2*sqrt(d + e*x)*(f + g*x)) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(2*g*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_690():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/((d + e*x)**(sympy.S(3)/2)*(f + g*x)**4)
    F = c**3*d**3*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(8*g**(sympy.S(5)/2)*(-a*e*g + c*d*f)**(sympy.S(3)/2)) + c**2*d**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(8*g**2*sqrt(d + e*x)*(f + g*x)*(-a*e*g + c*d*f)) - c*d*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*g**2*sqrt(d + e*x)*(f + g*x)**2) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(3*g*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_691():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/((d + e*x)**(sympy.S(3)/2)*(f + g*x)**5)
    F = 3*c**4*d**4*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(64*g**(sympy.S(5)/2)*(-a*e*g + c*d*f)**(sympy.S(5)/2)) + 3*c**3*d**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(64*g**2*sqrt(d + e*x)*(f + g*x)*(-a*e*g + c*d*f)**2) + c**2*d**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(32*g**2*sqrt(d + e*x)*(f + g*x)**2*(-a*e*g + c*d*f)) - c*d*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(8*g**2*sqrt(d + e*x)*(f + g*x)**3) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(4*g*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_692():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/((d + e*x)**(sympy.S(3)/2)*(f + g*x)**6)
    F = 3*c**5*d**5*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(128*g**(sympy.S(5)/2)*(-a*e*g + c*d*f)**(sympy.S(7)/2)) + 3*c**4*d**4*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(128*g**2*sqrt(d + e*x)*(f + g*x)*(-a*e*g + c*d*f)**3) + c**3*d**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(64*g**2*sqrt(d + e*x)*(f + g*x)**2*(-a*e*g + c*d*f)**2) + c**2*d**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(80*g**2*sqrt(d + e*x)*(f + g*x)**3*(-a*e*g + c*d*f)) - 3*c*d*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(40*g**2*sqrt(d + e*x)*(f + g*x)**4) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(5*g*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_693():
    f = (f + g*x)**4*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(d + e*x)**(sympy.S(5)/2)
    F = 2*(f + g*x)**4*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(15*c*d*(d + e*x)**(sympy.S(7)/2)) + (f + g*x)**3*(-16*a*e*g + 16*c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(195*c**2*d**2*(d + e*x)**(sympy.S(7)/2)) + 32*(f + g*x)**2*(-a*e*g + c*d*f)**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(715*c**3*d**3*(d + e*x)**(sympy.S(7)/2)) + 128*g*(-a*e*g + c*d*f)**3*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(6435*c**4*d**4*e*(d + e*x)**(sympy.S(5)/2)) - 128*(-a*e*g + c*d*f)**3*(2*a*e**2*g - c*d*(-7*d*g + 9*e*f))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(45045*c**5*d**5*e*(d + e*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_694():
    f = (f + g*x)**3*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(d + e*x)**(sympy.S(5)/2)
    F = 2*(f + g*x)**3*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(13*c*d*(d + e*x)**(sympy.S(7)/2)) + (f + g*x)**2*(-12*a*e*g + 12*c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(143*c**2*d**2*(d + e*x)**(sympy.S(7)/2)) + 16*g*(-a*e*g + c*d*f)**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(429*c**3*d**3*e*(d + e*x)**(sympy.S(5)/2)) - 16*(-a*e*g + c*d*f)**2*(2*a*e**2*g - c*d*(-7*d*g + 9*e*f))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(3003*c**4*d**4*e*(d + e*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_695():
    f = (f + g*x)**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(d + e*x)**(sympy.S(5)/2)
    F = 2*(f + g*x)**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(11*c*d*(d + e*x)**(sympy.S(7)/2)) + 8*g*(-a*e*g + c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(99*c**2*d**2*e*(d + e*x)**(sympy.S(5)/2)) - (-8*a*e*g + 8*c*d*f)*(2*a*e**2*g - c*d*(-7*d*g + 9*e*f))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(693*c**3*d**3*e*(d + e*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_696():
    f = (f + g*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(d + e*x)**(sympy.S(5)/2)
    F = 2*g*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(9*c*d*e*(d + e*x)**(sympy.S(5)/2)) - (4*a*e**2*g - 2*c*d*(-7*d*g + 9*e*f))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(63*c**2*d**2*e*(d + e*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_697():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(d + e*x)**(sympy.S(5)/2)
    F = 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(7*c*d*(d + e*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_698():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/((d + e*x)**(sympy.S(5)/2)*(f + g*x))
    F = 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(5*g*(d + e*x)**(sympy.S(5)/2)) - (-2*a*e*g + 2*c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(3*g**2*(d + e*x)**(sympy.S(3)/2)) + 2*(-a*e*g + c*d*f)**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(g**3*sqrt(d + e*x)) - 2*(-a*e*g + c*d*f)**(sympy.S(5)/2)*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/g**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_699():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/((d + e*x)**(sympy.S(5)/2)*(f + g*x)**2)
    F = 5*c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(3*g**2*(d + e*x)**(sympy.S(3)/2)) - 5*c*d*(-a*e*g + c*d*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(g**3*sqrt(d + e*x)) + 5*c*d*(-a*e*g + c*d*f)**(sympy.S(3)/2)*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/g**(sympy.S(7)/2) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(g*(d + e*x)**(sympy.S(5)/2)*(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_700():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/((d + e*x)**(sympy.S(5)/2)*(f + g*x)**3)
    F = 15*c**2*d**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*g**3*sqrt(d + e*x)) - 15*c**2*d**2*sqrt(-a*e*g + c*d*f)*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(4*g**(sympy.S(7)/2)) - 5*c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(4*g**2*(d + e*x)**(sympy.S(3)/2)*(f + g*x)) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(2*g*(d + e*x)**(sympy.S(5)/2)*(f + g*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_701():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/((d + e*x)**(sympy.S(5)/2)*(f + g*x)**4)
    F = 5*c**3*d**3*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(8*g**(sympy.S(7)/2)*sqrt(-a*e*g + c*d*f)) - 5*c**2*d**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(8*g**3*sqrt(d + e*x)*(f + g*x)) - 5*c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(12*g**2*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**2) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(3*g*(d + e*x)**(sympy.S(5)/2)*(f + g*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_702():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/((d + e*x)**(sympy.S(5)/2)*(f + g*x)**5)
    F = 5*c**4*d**4*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(64*g**(sympy.S(7)/2)*(-a*e*g + c*d*f)**(sympy.S(3)/2)) + 5*c**3*d**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(64*g**3*sqrt(d + e*x)*(f + g*x)*(-a*e*g + c*d*f)) - 5*c**2*d**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(32*g**3*sqrt(d + e*x)*(f + g*x)**2) - 5*c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(24*g**2*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**3) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(4*g*(d + e*x)**(sympy.S(5)/2)*(f + g*x)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_703():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/((d + e*x)**(sympy.S(5)/2)*(f + g*x)**6)
    F = 3*c**5*d**5*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(128*g**(sympy.S(7)/2)*(-a*e*g + c*d*f)**(sympy.S(5)/2)) + 3*c**4*d**4*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(128*g**3*sqrt(d + e*x)*(f + g*x)*(-a*e*g + c*d*f)**2) + c**3*d**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(64*g**3*sqrt(d + e*x)*(f + g*x)**2*(-a*e*g + c*d*f)) - c**2*d**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(16*g**3*sqrt(d + e*x)*(f + g*x)**3) - c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(8*g**2*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**4) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(5*g*(d + e*x)**(sympy.S(5)/2)*(f + g*x)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_704():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/((d + e*x)**(sympy.S(5)/2)*(f + g*x)**7)
    F = 5*c**6*d**6*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(512*g**(sympy.S(7)/2)*(-a*e*g + c*d*f)**(sympy.S(7)/2)) + 5*c**5*d**5*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(512*g**3*sqrt(d + e*x)*(f + g*x)*(-a*e*g + c*d*f)**3) + 5*c**4*d**4*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(768*g**3*sqrt(d + e*x)*(f + g*x)**2*(-a*e*g + c*d*f)**2) + c**3*d**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(192*g**3*sqrt(d + e*x)*(f + g*x)**3*(-a*e*g + c*d*f)) - c**2*d**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(32*g**3*sqrt(d + e*x)*(f + g*x)**4) - c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(12*g**2*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**5) - (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(6*g*(d + e*x)**(sympy.S(5)/2)*(f + g*x)**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_705():
    f = sqrt(d + e*x)*(f + g*x)**(sympy.S(5)/2)/sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))
    F = (f + g*x)**(sympy.S(5)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3*c*d*sqrt(d + e*x)) + (f + g*x)**(sympy.S(3)/2)*(-5*a*e*g + 5*c*d*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(12*c**2*d**2*sqrt(d + e*x)) + 5*sqrt(f + g*x)*(-a*e*g + c*d*f)**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(8*c**3*d**3*sqrt(d + e*x)) + 5*sqrt(d + e*x)*sqrt(a*e + c*d*x)*(-a*e*g + c*d*f)**3*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(8*c**(sympy.S(7)/2)*d**(sympy.S(7)/2)*sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_706():
    f = sqrt(d + e*x)*(f + g*x)**(sympy.S(3)/2)/sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))
    F = (f + g*x)**(sympy.S(3)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(2*c*d*sqrt(d + e*x)) + sqrt(f + g*x)*(-3*a*e*g + 3*c*d*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*c**2*d**2*sqrt(d + e*x)) + 3*sqrt(d + e*x)*sqrt(a*e + c*d*x)*(-a*e*g + c*d*f)**2*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(4*c**(sympy.S(5)/2)*d**(sympy.S(5)/2)*sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_707():
    f = sqrt(d + e*x)*sqrt(f + g*x)/sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))
    F = sqrt(f + g*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(c*d*sqrt(d + e*x)) + sqrt(d + e*x)*sqrt(a*e + c*d*x)*(-a*e*g + c*d*f)*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(c**(sympy.S(3)/2)*d**(sympy.S(3)/2)*sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_708():
    f = sqrt(d + e*x)/(sqrt(f + g*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    F = 2*sqrt(d + e*x)*sqrt(a*e + c*d*x)*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(sqrt(c)*sqrt(d)*sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_709():
    f = sqrt(d + e*x)/((f + g*x)**(sympy.S(3)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    F = 2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(f + g*x)*(-a*e*g + c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_710():
    f = sqrt(d + e*x)/((f + g*x)**(sympy.S(5)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    F = 4*c*d*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3*sqrt(d + e*x)*sqrt(f + g*x)*(-a*e*g + c*d*f)**2) + 2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*(f + g*x)**(sympy.S(3)/2)*(-3*a*e*g + 3*c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_711():
    f = sqrt(d + e*x)/((f + g*x)**(sympy.S(7)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    F = 16*c**2*d**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(15*sqrt(d + e*x)*sqrt(f + g*x)*(-a*e*g + c*d*f)**3) + 8*c*d*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(15*sqrt(d + e*x)*(f + g*x)**(sympy.S(3)/2)*(-a*e*g + c*d*f)**2) + 2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*(f + g*x)**(sympy.S(5)/2)*(-5*a*e*g + 5*c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_712():
    f = sqrt(d + e*x)/((f + g*x)**(sympy.S(9)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    F = 32*c**3*d**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(35*sqrt(d + e*x)*sqrt(f + g*x)*(-a*e*g + c*d*f)**4) + 16*c**2*d**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(35*sqrt(d + e*x)*(f + g*x)**(sympy.S(3)/2)*(-a*e*g + c*d*f)**3) + 12*c*d*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(35*sqrt(d + e*x)*(f + g*x)**(sympy.S(5)/2)*(-a*e*g + c*d*f)**2) + 2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*(f + g*x)**(sympy.S(7)/2)*(-7*a*e*g + 7*c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_713():
    f = (d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(5)/2)/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)
    F = -2*sqrt(d + e*x)*(f + g*x)**(sympy.S(5)/2)/(c*d*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) + 5*g*(f + g*x)**(sympy.S(3)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(2*c**2*d**2*sqrt(d + e*x)) + 15*g*sqrt(f + g*x)*(-a*e*g + c*d*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*c**3*d**3*sqrt(d + e*x)) + 15*sqrt(g)*sqrt(d + e*x)*sqrt(a*e + c*d*x)*(-a*e*g + c*d*f)**2*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(4*c**(sympy.S(7)/2)*d**(sympy.S(7)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_714():
    f = (d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(3)/2)/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)
    F = -2*sqrt(d + e*x)*(f + g*x)**(sympy.S(3)/2)/(c*d*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) + 3*g*sqrt(f + g*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(c**2*d**2*sqrt(d + e*x)) + 3*sqrt(g)*sqrt(d + e*x)*sqrt(a*e + c*d*x)*(-a*e*g + c*d*f)*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(c**(sympy.S(5)/2)*d**(sympy.S(5)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_715():
    f = (d + e*x)**(sympy.S(3)/2)*sqrt(f + g*x)/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)
    F = -2*sqrt(d + e*x)*sqrt(f + g*x)/(c*d*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) + 2*sqrt(g)*sqrt(d + e*x)*sqrt(a*e + c*d*x)*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(c**(sympy.S(3)/2)*d**(sympy.S(3)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_716():
    f = (d + e*x)**(sympy.S(3)/2)/(sqrt(f + g*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    F = -2*sqrt(d + e*x)*sqrt(f + g*x)/((-a*e*g + c*d*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_717():
    f = (d + e*x)**(sympy.S(3)/2)/((f + g*x)**(sympy.S(3)/2)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    F = -4*g*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(f + g*x)*(-a*e*g + c*d*f)**2) - 2*sqrt(d + e*x)/(sqrt(f + g*x)*(-a*e*g + c*d*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_718():
    f = (d + e*x)**(sympy.S(3)/2)/((f + g*x)**(sympy.S(5)/2)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    F = -16*c*d*g*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3*sqrt(d + e*x)*sqrt(f + g*x)*(-a*e*g + c*d*f)**3) - 8*g*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3*sqrt(d + e*x)*(f + g*x)**(sympy.S(3)/2)*(-a*e*g + c*d*f)**2) - 2*sqrt(d + e*x)/((f + g*x)**(sympy.S(3)/2)*(-a*e*g + c*d*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_719():
    f = (d + e*x)**(sympy.S(3)/2)/((f + g*x)**(sympy.S(7)/2)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    F = -32*c**2*d**2*g*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(5*sqrt(d + e*x)*sqrt(f + g*x)*(-a*e*g + c*d*f)**4) - 16*c*d*g*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(5*sqrt(d + e*x)*(f + g*x)**(sympy.S(3)/2)*(-a*e*g + c*d*f)**3) - 12*g*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(5*sqrt(d + e*x)*(f + g*x)**(sympy.S(5)/2)*(-a*e*g + c*d*f)**2) - 2*sqrt(d + e*x)/((f + g*x)**(sympy.S(5)/2)*(-a*e*g + c*d*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_720():
    f = (d + e*x)**(sympy.S(5)/2)*(f + g*x)**(sympy.S(5)/2)/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)
    F = -2*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(5)/2)/(3*c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)) - 10*g*sqrt(d + e*x)*(f + g*x)**(sympy.S(3)/2)/(3*c**2*d**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) + 5*g**2*sqrt(f + g*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(c**3*d**3*sqrt(d + e*x)) + 5*g**(sympy.S(3)/2)*sqrt(d + e*x)*sqrt(a*e + c*d*x)*(-a*e*g + c*d*f)*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(c**(sympy.S(7)/2)*d**(sympy.S(7)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_721():
    f = (d + e*x)**(sympy.S(5)/2)*(f + g*x)**(sympy.S(3)/2)/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)
    F = -2*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(3)/2)/(3*c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)) - 2*g*sqrt(d + e*x)*sqrt(f + g*x)/(c**2*d**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) + 2*g**(sympy.S(3)/2)*sqrt(d + e*x)*sqrt(a*e + c*d*x)*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(c**(sympy.S(5)/2)*d**(sympy.S(5)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_722():
    f = (d + e*x)**(sympy.S(5)/2)*sqrt(f + g*x)/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)
    F = -2*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(3)/2)/((-3*a*e*g + 3*c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_723():
    f = (d + e*x)**(sympy.S(5)/2)/(sqrt(f + g*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2))
    F = 4*g*sqrt(d + e*x)*sqrt(f + g*x)/(3*(-a*e*g + c*d*f)**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) - 2*(d + e*x)**(sympy.S(3)/2)*sqrt(f + g*x)/((-3*a*e*g + 3*c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_724():
    f = (d + e*x)**(sympy.S(5)/2)/((f + g*x)**(sympy.S(3)/2)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2))
    F = 16*g**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3*sqrt(d + e*x)*sqrt(f + g*x)*(-a*e*g + c*d*f)**3) + 8*g*sqrt(d + e*x)/(3*sqrt(f + g*x)*(-a*e*g + c*d*f)**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) - 2*(d + e*x)**(sympy.S(3)/2)/(sqrt(f + g*x)*(-3*a*e*g + 3*c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_725():
    f = (d + e*x)**(sympy.S(5)/2)/((f + g*x)**(sympy.S(5)/2)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2))
    F = 32*c*d*g**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3*sqrt(d + e*x)*sqrt(f + g*x)*(-a*e*g + c*d*f)**4) + 16*g**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3*sqrt(d + e*x)*(f + g*x)**(sympy.S(3)/2)*(-a*e*g + c*d*f)**3) + 4*g*sqrt(d + e*x)/((f + g*x)**(sympy.S(3)/2)*(-a*e*g + c*d*f)**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) - 2*(d + e*x)**(sympy.S(3)/2)/((f + g*x)**(sympy.S(3)/2)*(-3*a*e*g + 3*c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_726():
    f = (f + g*x)**(sympy.S(5)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/sqrt(d + e*x)
    F = (f + g*x)**(sympy.S(5)/2)*(a*e/(c*d) - f/g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(24*sqrt(d + e*x)) + (f + g*x)**(sympy.S(7)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*g*sqrt(d + e*x)) - 5*(f + g*x)**(sympy.S(3)/2)*(-a*e*g + c*d*f)**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(96*c**2*d**2*g*sqrt(d + e*x)) - 5*sqrt(f + g*x)*(-a*e*g + c*d*f)**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(64*c**3*d**3*g*sqrt(d + e*x)) - 5*sqrt(d + e*x)*sqrt(a*e + c*d*x)*(-a*e*g + c*d*f)**4*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(64*c**(sympy.S(7)/2)*d**(sympy.S(7)/2)*g**(sympy.S(3)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_727():
    f = (f + g*x)**(sympy.S(3)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/sqrt(d + e*x)
    F = (f + g*x)**(sympy.S(3)/2)*(a*e/(c*d) - f/g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(12*sqrt(d + e*x)) + (f + g*x)**(sympy.S(5)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3*g*sqrt(d + e*x)) - sqrt(f + g*x)*(-a*e*g + c*d*f)**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(8*c**2*d**2*g*sqrt(d + e*x)) - sqrt(d + e*x)*sqrt(a*e + c*d*x)*(-a*e*g + c*d*f)**3*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(8*c**(sympy.S(5)/2)*d**(sympy.S(5)/2)*g**(sympy.S(3)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_728():
    f = sqrt(f + g*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/sqrt(d + e*x)
    F = sqrt(f + g*x)*(a*e/(c*d) - f/g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*sqrt(d + e*x)) + (f + g*x)**(sympy.S(3)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(2*g*sqrt(d + e*x)) - sqrt(d + e*x)*sqrt(a*e + c*d*x)*(-a*e*g + c*d*f)**2*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(4*c**(sympy.S(3)/2)*d**(sympy.S(3)/2)*g**(sympy.S(3)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_729():
    f = sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(f + g*x))
    F = sqrt(f + g*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(g*sqrt(d + e*x)) - sqrt(d + e*x)*sqrt(a*e + c*d*x)*(-a*e*g + c*d*f)*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(sqrt(c)*sqrt(d)*g**(sympy.S(3)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_730():
    f = sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*(f + g*x)**(sympy.S(3)/2))
    F = 2*sqrt(c)*sqrt(d)*sqrt(d + e*x)*sqrt(a*e + c*d*x)*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(g**(sympy.S(3)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) - 2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(g*sqrt(d + e*x)*sqrt(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_731():
    f = sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*(f + g*x)**(sympy.S(5)/2))
    F = 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/((d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(3)/2)*(-3*a*e*g + 3*c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_732():
    f = sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*(f + g*x)**(sympy.S(7)/2))
    F = 4*c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(15*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(3)/2)*(-a*e*g + c*d*f)**2) + 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/((d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(5)/2)*(-5*a*e*g + 5*c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_733():
    f = sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*(f + g*x)**(sympy.S(9)/2))
    F = 16*c**2*d**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(105*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(3)/2)*(-a*e*g + c*d*f)**3) + 8*c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(35*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(5)/2)*(-a*e*g + c*d*f)**2) + 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/((d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(7)/2)*(-7*a*e*g + 7*c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_734():
    f = sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*(f + g*x)**(sympy.S(11)/2))
    F = 32*c**3*d**3*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(315*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(3)/2)*(-a*e*g + c*d*f)**4) + 16*c**2*d**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(105*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(5)/2)*(-a*e*g + c*d*f)**3) + 4*c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(21*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(7)/2)*(-a*e*g + c*d*f)**2) + 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/((d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(9)/2)*(-9*a*e*g + 9*c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_735():
    f = (f + g*x)**(sympy.S(3)/2)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(d + e*x)**(sympy.S(3)/2)
    F = (f + g*x)**(sympy.S(5)/2)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(4*g*(d + e*x)**(sympy.S(3)/2)) - (f + g*x)**(sympy.S(5)/2)*(-a*e*g + c*d*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(8*g**2*sqrt(d + e*x)) + (f + g*x)**(sympy.S(3)/2)*(-a*e*g + c*d*f)**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(32*c*d*g**2*sqrt(d + e*x)) + 3*sqrt(f + g*x)*(-a*e*g + c*d*f)**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(64*c**2*d**2*g**2*sqrt(d + e*x)) + 3*sqrt(d + e*x)*sqrt(a*e + c*d*x)*(-a*e*g + c*d*f)**4*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(64*c**(sympy.S(5)/2)*d**(sympy.S(5)/2)*g**(sympy.S(5)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_736():
    f = sqrt(f + g*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(d + e*x)**(sympy.S(3)/2)
    F = (f + g*x)**(sympy.S(3)/2)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(3*g*(d + e*x)**(sympy.S(3)/2)) - (f + g*x)**(sympy.S(3)/2)*(-a*e*g + c*d*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*g**2*sqrt(d + e*x)) + sqrt(f + g*x)*(-a*e*g + c*d*f)**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(8*c*d*g**2*sqrt(d + e*x)) + sqrt(d + e*x)*sqrt(a*e + c*d*x)*(-a*e*g + c*d*f)**3*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(8*c**(sympy.S(3)/2)*d**(sympy.S(3)/2)*g**(sympy.S(5)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_737():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/((d + e*x)**(sympy.S(3)/2)*sqrt(f + g*x))
    F = sqrt(f + g*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(2*g*(d + e*x)**(sympy.S(3)/2)) - sqrt(f + g*x)*(-3*a*e*g + 3*c*d*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*g**2*sqrt(d + e*x)) + 3*sqrt(d + e*x)*sqrt(a*e + c*d*x)*(-a*e*g + c*d*f)**2*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(4*sqrt(c)*sqrt(d)*g**(sympy.S(5)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_738():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/((d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(3)/2))
    F = -3*sqrt(c)*sqrt(d)*sqrt(d + e*x)*sqrt(a*e + c*d*x)*(-a*e*g + c*d*f)*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(g**(sympy.S(5)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) + 3*c*d*sqrt(f + g*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(g**2*sqrt(d + e*x)) - 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(g*(d + e*x)**(sympy.S(3)/2)*sqrt(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_739():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/((d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(5)/2))
    F = 2*c**(sympy.S(3)/2)*d**(sympy.S(3)/2)*sqrt(d + e*x)*sqrt(a*e + c*d*x)*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(g**(sympy.S(5)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) - 2*c*d*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(g**2*sqrt(d + e*x)*sqrt(f + g*x)) - 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(3*g*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_740():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/((d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(7)/2))
    F = 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/((d + e*x)**(sympy.S(5)/2)*(f + g*x)**(sympy.S(5)/2)*(-5*a*e*g + 5*c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_741():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/((d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(9)/2))
    F = 4*c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(35*(d + e*x)**(sympy.S(5)/2)*(f + g*x)**(sympy.S(5)/2)*(-a*e*g + c*d*f)**2) + 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/((d + e*x)**(sympy.S(5)/2)*(f + g*x)**(sympy.S(7)/2)*(-7*a*e*g + 7*c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_742():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/((d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(11)/2))
    F = 16*c**2*d**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(315*(d + e*x)**(sympy.S(5)/2)*(f + g*x)**(sympy.S(5)/2)*(-a*e*g + c*d*f)**3) + 8*c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(63*(d + e*x)**(sympy.S(5)/2)*(f + g*x)**(sympy.S(7)/2)*(-a*e*g + c*d*f)**2) + 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/((d + e*x)**(sympy.S(5)/2)*(f + g*x)**(sympy.S(9)/2)*(-9*a*e*g + 9*c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_743():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/((d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(13)/2))
    F = 32*c**3*d**3*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(1155*(d + e*x)**(sympy.S(5)/2)*(f + g*x)**(sympy.S(5)/2)*(-a*e*g + c*d*f)**4) + 16*c**2*d**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(231*(d + e*x)**(sympy.S(5)/2)*(f + g*x)**(sympy.S(7)/2)*(-a*e*g + c*d*f)**3) + 4*c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(33*(d + e*x)**(sympy.S(5)/2)*(f + g*x)**(sympy.S(9)/2)*(-a*e*g + c*d*f)**2) + 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/((d + e*x)**(sympy.S(5)/2)*(f + g*x)**(sympy.S(11)/2)*(-11*a*e*g + 11*c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_744():
    f = (f + g*x)**(sympy.S(3)/2)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(d + e*x)**(sympy.S(5)/2)
    F = (f + g*x)**(sympy.S(5)/2)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(5*g*(d + e*x)**(sympy.S(5)/2)) - (f + g*x)**(sympy.S(5)/2)*(-a*e*g + c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(8*g**2*(d + e*x)**(sympy.S(3)/2)) + (f + g*x)**(sympy.S(5)/2)*(-a*e*g + c*d*f)**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(16*g**3*sqrt(d + e*x)) - (f + g*x)**(sympy.S(3)/2)*(-a*e*g + c*d*f)**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(64*c*d*g**3*sqrt(d + e*x)) - 3*sqrt(f + g*x)*(-a*e*g + c*d*f)**4*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(128*c**2*d**2*g**3*sqrt(d + e*x)) - 3*sqrt(d + e*x)*sqrt(a*e + c*d*x)*(-a*e*g + c*d*f)**5*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(128*c**(sympy.S(5)/2)*d**(sympy.S(5)/2)*g**(sympy.S(7)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_745():
    f = sqrt(f + g*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(d + e*x)**(sympy.S(5)/2)
    F = (f + g*x)**(sympy.S(3)/2)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(4*g*(d + e*x)**(sympy.S(5)/2)) - (f + g*x)**(sympy.S(3)/2)*(-5*a*e*g + 5*c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(24*g**2*(d + e*x)**(sympy.S(3)/2)) + 5*(f + g*x)**(sympy.S(3)/2)*(-a*e*g + c*d*f)**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(32*g**3*sqrt(d + e*x)) - 5*sqrt(f + g*x)*(-a*e*g + c*d*f)**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(64*c*d*g**3*sqrt(d + e*x)) - 5*sqrt(d + e*x)*sqrt(a*e + c*d*x)*(-a*e*g + c*d*f)**4*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(64*c**(sympy.S(3)/2)*d**(sympy.S(3)/2)*g**(sympy.S(7)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_746():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/((d + e*x)**(sympy.S(5)/2)*sqrt(f + g*x))
    F = sqrt(f + g*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(3*g*(d + e*x)**(sympy.S(5)/2)) - sqrt(f + g*x)*(-5*a*e*g + 5*c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(12*g**2*(d + e*x)**(sympy.S(3)/2)) + 5*sqrt(f + g*x)*(-a*e*g + c*d*f)**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(8*g**3*sqrt(d + e*x)) - 5*sqrt(d + e*x)*sqrt(a*e + c*d*x)*(-a*e*g + c*d*f)**3*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(8*sqrt(c)*sqrt(d)*g**(sympy.S(7)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_747():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/((d + e*x)**(sympy.S(5)/2)*(f + g*x)**(sympy.S(3)/2))
    F = 15*sqrt(c)*sqrt(d)*sqrt(d + e*x)*sqrt(a*e + c*d*x)*(-a*e*g + c*d*f)**2*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(4*g**(sympy.S(7)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) + 5*c*d*sqrt(f + g*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(2*g**2*(d + e*x)**(sympy.S(3)/2)) - 15*c*d*sqrt(f + g*x)*(-a*e*g + c*d*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*g**3*sqrt(d + e*x)) - 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(g*(d + e*x)**(sympy.S(5)/2)*sqrt(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_748():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/((d + e*x)**(sympy.S(5)/2)*(f + g*x)**(sympy.S(5)/2))
    F = -5*c**(sympy.S(3)/2)*d**(sympy.S(3)/2)*sqrt(d + e*x)*sqrt(a*e + c*d*x)*(-a*e*g + c*d*f)*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(g**(sympy.S(7)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) + 5*c**2*d**2*sqrt(f + g*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(g**3*sqrt(d + e*x)) - 10*c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(3*g**2*(d + e*x)**(sympy.S(3)/2)*sqrt(f + g*x)) - 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(3*g*(d + e*x)**(sympy.S(5)/2)*(f + g*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_749():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/((d + e*x)**(sympy.S(5)/2)*(f + g*x)**(sympy.S(7)/2))
    F = 2*c**(sympy.S(5)/2)*d**(sympy.S(5)/2)*sqrt(d + e*x)*sqrt(a*e + c*d*x)*atanh(sqrt(g)*sqrt(a*e + c*d*x)/(sqrt(c)*sqrt(d)*sqrt(f + g*x)))/(g**(sympy.S(7)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))) - 2*c**2*d**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(g**3*sqrt(d + e*x)*sqrt(f + g*x)) - 2*c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(3)/2)/(3*g**2*(d + e*x)**(sympy.S(3)/2)*(f + g*x)**(sympy.S(3)/2)) - 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/(5*g*(d + e*x)**(sympy.S(5)/2)*(f + g*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_750():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/((d + e*x)**(sympy.S(5)/2)*(f + g*x)**(sympy.S(9)/2))
    F = 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/((d + e*x)**(sympy.S(7)/2)*(f + g*x)**(sympy.S(7)/2)*(-7*a*e*g + 7*c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_751():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/((d + e*x)**(sympy.S(5)/2)*(f + g*x)**(sympy.S(11)/2))
    F = 4*c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(63*(d + e*x)**(sympy.S(7)/2)*(f + g*x)**(sympy.S(7)/2)*(-a*e*g + c*d*f)**2) + 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/((d + e*x)**(sympy.S(7)/2)*(f + g*x)**(sympy.S(9)/2)*(-9*a*e*g + 9*c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_752():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/((d + e*x)**(sympy.S(5)/2)*(f + g*x)**(sympy.S(13)/2))
    F = 16*c**2*d**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(693*(d + e*x)**(sympy.S(7)/2)*(f + g*x)**(sympy.S(7)/2)*(-a*e*g + c*d*f)**3) + 8*c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(99*(d + e*x)**(sympy.S(7)/2)*(f + g*x)**(sympy.S(9)/2)*(-a*e*g + c*d*f)**2) + 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/((d + e*x)**(sympy.S(7)/2)*(f + g*x)**(sympy.S(11)/2)*(-11*a*e*g + 11*c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_753():
    f = (a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(5)/2)/((d + e*x)**(sympy.S(5)/2)*(f + g*x)**(sympy.S(15)/2))
    F = 32*c**3*d**3*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(3003*(d + e*x)**(sympy.S(7)/2)*(f + g*x)**(sympy.S(7)/2)*(-a*e*g + c*d*f)**4) + 16*c**2*d**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(429*(d + e*x)**(sympy.S(7)/2)*(f + g*x)**(sympy.S(9)/2)*(-a*e*g + c*d*f)**3) + 12*c*d*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/(143*(d + e*x)**(sympy.S(7)/2)*(f + g*x)**(sympy.S(11)/2)*(-a*e*g + c*d*f)**2) + 2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(sympy.S(7)/2)/((d + e*x)**(sympy.S(7)/2)*(f + g*x)**(sympy.S(13)/2)*(-13*a*e*g + 13*c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_754():
    f = (d + e*x)**m*(f + g*x)**3/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m
    F = (d + e*x)**(m - 1)*(f + g*x)**3*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(1 - m)/(c*d*(4 - m)) + (d + e*x)**(m - 1)*(f + g*x)**2*(-3*a*e*g + 3*c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(1 - m)/(c**2*d**2*(3 - m)*(4 - m)) + 6*g*(d + e*x)**m*(-a*e*g + c*d*f)**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(1 - m)/(c**3*d**3*e*(2 - m)*(3 - m)*(4 - m)) - 6*(d + e*x)**(m - 1)*(-a*e*g + c*d*f)**2*(a*e**2*g + c*d*(d*g*(1 - m) - e*f*(2 - m)))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(1 - m)/(c**4*d**4*e*(1 - m)*(2 - m)*(3 - m)*(4 - m))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_755():
    f = (d + e*x)**m*(f + g*x)**2/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m
    F = (d + e*x)**(m - 1)*(f + g*x)**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(1 - m)/(c*d*(3 - m)) + 2*g*(d + e*x)**m*(-a*e*g + c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(1 - m)/(c**2*d**2*e*(2 - m)*(3 - m)) - (d + e*x)**(m - 1)*(-2*a*e*g + 2*c*d*f)*(a*e**2*g + c*d*(d*g*(1 - m) - e*f*(2 - m)))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(1 - m)/(c**3*d**3*e*(1 - m)*(2 - m)*(3 - m))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_756():
    f = (d + e*x)**m*(f + g*x)/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m
    F = g*(d + e*x)**m*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(1 - m)/(c*d*e*(2 - m)) - (d + e*x)**(m - 1)*(a*e**2*g + c*d*(d*g*(1 - m) - e*f*(2 - m)))*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(1 - m)/(c**2*d**2*e*(1 - m)*(2 - m))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_757():
    f = (d + e*x)**m/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m
    F = (d + e*x)**(m - 1)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(1 - m)/(c*d*(1 - m))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_758():
    f = (d + e*x)**m/((f + g*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m)
    F = (d + e*x)**m*(a*e + c*d*x)*hyper((1, 1 - m), (2 - m,), -g*(a*e + c*d*x)/(-a*e*g + c*d*f))/((1 - m)*(-a*e*g + c*d*f)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_759():
    f = (d + e*x)**m/((f + g*x)**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m)
    F = c*d*(d + e*x)**m*(a*e + c*d*x)*hyper((2, 1 - m), (2 - m,), -g*(a*e + c*d*x)/(-a*e*g + c*d*f))/((1 - m)*(-a*e*g + c*d*f)**2*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_760():
    f = (d + e*x)**m/((f + g*x)**3*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m)
    F = c**2*d**2*(d + e*x)**m*(a*e + c*d*x)*hyper((3, 1 - m), (2 - m,), -g*(a*e + c*d*x)/(-a*e*g + c*d*f))/((1 - m)*(-a*e*g + c*d*f)**3*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_761():
    f = (d + e*x)**m*(f + g*x)**(sympy.S(3)/2)/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m
    F = 2*(-g*(a*e + c*d*x)/(-a*e*g + c*d*f))**m*(d + e*x)**m*(f + g*x)**(sympy.S(5)/2)*hyper((sympy.S(5)/2, m), (sympy.S(7)/2,), c*d*(f + g*x)/(-a*e*g + c*d*f))/(5*g*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_762():
    f = (d + e*x)**m*sqrt(f + g*x)/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m
    F = 2*(-g*(a*e + c*d*x)/(-a*e*g + c*d*f))**m*(d + e*x)**m*(f + g*x)**(sympy.S(3)/2)*hyper((sympy.S(3)/2, m), (sympy.S(5)/2,), c*d*(f + g*x)/(-a*e*g + c*d*f))/(3*g*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_763():
    f = (d + e*x)**m/(sqrt(f + g*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m)
    F = 2*(-g*(a*e + c*d*x)/(-a*e*g + c*d*f))**m*(d + e*x)**m*sqrt(f + g*x)*hyper((sympy.S.Half, m), (sympy.S(3)/2,), c*d*(f + g*x)/(-a*e*g + c*d*f))/(g*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_764():
    f = (d + e*x)**m/((f + g*x)**(sympy.S(3)/2)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m)
    F = -2*(-g*(a*e + c*d*x)/(-a*e*g + c*d*f))**m*(d + e*x)**m*hyper((sympy.S(-1)/2, m), (sympy.S.Half,), c*d*(f + g*x)/(-a*e*g + c*d*f))/(g*sqrt(f + g*x)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_765():
    f = (d + e*x)**m/((f + g*x)**(sympy.S(5)/2)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m)
    F = -2*(-g*(a*e + c*d*x)/(-a*e*g + c*d*f))**m*(d + e*x)**m*hyper((sympy.S(-3)/2, m), (sympy.S(-1)/2,), c*d*(f + g*x)/(-a*e*g + c*d*f))/(3*g*(f + g*x)**(sympy.S(3)/2)*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_766():
    f = (d + e*x)**m*(a*e + c*d*x)**n/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m
    F = (d + e*x)**(m - 1)*(a*e + c*d*x)**n*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**(1 - m)/(c*d*(-m + n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_767():
    f = (d + e*x)**m*(c*d**2*e*g - c*d*e**2*g*x - e*g*(a*e**2 + c*d**2))**(m - 1)/(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m
    F = -(d + e*x)**m*(-a*e**3*g - c*d*e**2*g*x)**m*log(a*e + c*d*x)/(c*d*e**2*g*(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_768():
    f = (d + e*x)**(sympy.S(3)/2)*(f + g*x)**4/sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))
    F = 2*e*(f + g*x)**5*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(11*c*d*g*sqrt(d + e*x)) - (f + g*x)**4*(20*a*e**2*g + 2*c*d*(-11*d*g + e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(99*c**2*d**2*g*sqrt(d + e*x)) - (f + g*x)**3*(-16*a*e*g + 16*c*d*f)*(10*a*e**2*g + c*d*(-11*d*g + e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(693*c**3*d**3*g*sqrt(d + e*x)) - 32*(f + g*x)**2*(-a*e*g + c*d*f)**2*(10*a*e**2*g + c*d*(-11*d*g + e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(1155*c**4*d**4*g*sqrt(d + e*x)) - 128*sqrt(d + e*x)*(-a*e*g + c*d*f)**3*(10*a*e**2*g + c*d*(-11*d*g + e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3465*c**5*d**5*e) + 128*(-a*e*g + c*d*f)**3*(2*a*e**2*g - c*d*(-d*g + 3*e*f))*(10*a*e**2*g + c*d*(-11*d*g + e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3465*c**6*d**6*e*g*sqrt(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_769():
    f = (d + e*x)**(sympy.S(3)/2)*(f + g*x)**3/sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))
    F = 2*e*(f + g*x)**4*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(9*c*d*g*sqrt(d + e*x)) - (f + g*x)**3*(16*a*e**2*g + 2*c*d*(-9*d*g + e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(63*c**2*d**2*g*sqrt(d + e*x)) - (f + g*x)**2*(-4*a*e*g + 4*c*d*f)*(8*a*e**2*g + c*d*(-9*d*g + e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(105*c**3*d**3*g*sqrt(d + e*x)) - 16*sqrt(d + e*x)*(-a*e*g + c*d*f)**2*(8*a*e**2*g + c*d*(-9*d*g + e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(315*c**4*d**4*e) + 16*(-a*e*g + c*d*f)**2*(2*a*e**2*g - c*d*(-d*g + 3*e*f))*(8*a*e**2*g + c*d*(-9*d*g + e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(315*c**5*d**5*e*g*sqrt(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_770():
    f = (d + e*x)**(sympy.S(3)/2)*(f + g*x)**2/sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))
    F = 2*e*(f + g*x)**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(7*c*d*g*sqrt(d + e*x)) - (f + g*x)**2*(12*a*e**2*g + 2*c*d*(-7*d*g + e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(35*c**2*d**2*g*sqrt(d + e*x)) - sqrt(d + e*x)*(-8*a*e*g + 8*c*d*f)*(6*a*e**2*g + c*d*(-7*d*g + e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(105*c**3*d**3*e) + (-8*a*e*g + 8*c*d*f)*(2*a*e**2*g - c*d*(-d*g + 3*e*f))*(6*a*e**2*g + c*d*(-7*d*g + e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(105*c**4*d**4*e*g*sqrt(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_771():
    f = (d + e*x)**(sympy.S(3)/2)*(f + g*x)/sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))
    F = 2*g*(d + e*x)**(sympy.S(3)/2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(5*c*d*e) - sqrt(d + e*x)*(8*a*e**2*g - 2*c*d*(-d*g + 5*e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(15*c**2*d**2*e) - (-4*a*e**2 + 4*c*d**2)*(4*a*e**2*g - c*d*(-d*g + 5*e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(15*c**3*d**3*e*sqrt(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_772():
    f = (d + e*x)**(sympy.S(3)/2)/sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))
    F = 2*sqrt(d + e*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3*c*d) + (-4*a*e**2 + 4*c*d**2)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3*c**2*d**2*sqrt(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_773():
    f = (d + e*x)**(sympy.S(3)/2)/((f + g*x)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    F = -(-2*d*g + 2*e*f)*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(g**(sympy.S(3)/2)*sqrt(-a*e*g + c*d*f)) + 2*e*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(c*d*g*sqrt(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_774():
    f = (d + e*x)**(sympy.S(3)/2)/((f + g*x)**2*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    F = -(-d*g + e*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(g*sqrt(d + e*x)*(f + g*x)*(-a*e*g + c*d*f)) - (2*a*e**2*g - c*d*(d*g + e*f))*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(g**(sympy.S(3)/2)*(-a*e*g + c*d*f)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_775():
    f = (d + e*x)**(sympy.S(3)/2)/((f + g*x)**3*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    F = -c*d*(4*a*e**2*g - c*d*(3*d*g + e*f))*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(4*g**(sympy.S(3)/2)*(-a*e*g + c*d*f)**(sympy.S(5)/2)) - (4*a*e**2*g - c*d*(3*d*g + e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(4*g*sqrt(d + e*x)*(f + g*x)*(-a*e*g + c*d*f)**2) - (-d*g + e*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(2*g*sqrt(d + e*x)*(f + g*x)**2*(-a*e*g + c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_776():
    f = (d + e*x)**(sympy.S(3)/2)/((f + g*x)**4*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2)))
    F = -c**2*d**2*(6*a*e**2*g - c*d*(5*d*g + e*f))*atan(sqrt(g)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(sqrt(d + e*x)*sqrt(-a*e*g + c*d*f)))/(8*g**(sympy.S(3)/2)*(-a*e*g + c*d*f)**(sympy.S(7)/2)) - c*d*(6*a*e**2*g - c*d*(5*d*g + e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(8*g*sqrt(d + e*x)*(f + g*x)*(-a*e*g + c*d*f)**3) - (6*a*e**2*g - c*d*(5*d*g + e*f))*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(12*g*sqrt(d + e*x)*(f + g*x)**2*(-a*e*g + c*d*f)**2) - (-d*g + e*f)*sqrt(a*d*e + c*d*e*x**2 + x*(a*e**2 + c*d**2))/(3*g*sqrt(d + e*x)*(f + g*x)**3*(-a*e*g + c*d*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_777():
    f = (a + b*x + c*x**2)**3/(sqrt(-d*x + 1)*sqrt(d*x + 1))
    F = -3*b*c**2*x**4*sqrt(-d**2*x**2 + 1)/(5*d**2) - b*x**2*sqrt(-d**2*x**2 + 1)*(30*a*c*d**2 + 5*b**2*d**2 + 12*c**2)/(15*d**4) - b*sqrt(-d**2*x**2 + 1)*(45*a**2*d**4 + 60*a*c*d**2 + 10*b**2*d**2 + 24*c**2)/(15*d**6) - c**3*x**5*sqrt(-d**2*x**2 + 1)/(6*d**2) - c*x**3*sqrt(-d**2*x**2 + 1)*(18*a*c*d**2 + 18*b**2*d**2 + 5*c**2)/(24*d**4) - x*sqrt(-d**2*x**2 + 1)*(24*a**2*c*d**4 + 24*a*b**2*d**4 + 18*a*c**2*d**2 + 18*b**2*c*d**2 + 5*c**3)/(16*d**6) + (16*a**3*d**6 + 24*a**2*c*d**4 + 24*a*b**2*d**4 + 18*a*c**2*d**2 + 18*b**2*c*d**2 + 5*c**3)*asin(d*x)/(16*d**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_778():
    f = (a + b*x + c*x**2)**2/(sqrt(-d*x + 1)*sqrt(d*x + 1))
    F = -2*b*c*x**2*sqrt(-d**2*x**2 + 1)/(3*d**2) - 2*b*(3*a*d**2 + 2*c)*sqrt(-d**2*x**2 + 1)/(3*d**4) - c**2*x**3*sqrt(-d**2*x**2 + 1)/(4*d**2) - x*(4*b**2 + c*(8*a + 3*c/d**2))*sqrt(-d**2*x**2 + 1)/(8*d**2) + (8*a**2*d**4 + 8*a*c*d**2 + 4*b**2*d**2 + 3*c**2)*asin(d*x)/(8*d**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_779():
    f = (a + b*x + c*x**2)/(sqrt(-d*x + 1)*sqrt(d*x + 1))
    F = -b*sqrt(-d**2*x**2 + 1)/d**2 - c*x*sqrt(-d**2*x**2 + 1)/(2*d**2) + (2*a*d**2 + c)*asin(d*x)/(2*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_780():
    f = 1/(sqrt(-d*x + 1)*sqrt(d*x + 1)*(a + b*x + c*x**2))
    F = sqrt(2)*c*atanh(sqrt(2)*(2*c + d**2*x*(b + sqrt(-4*a*c + b**2)))/(2*sqrt(-d**2*x**2 + 1)*sqrt(2*a*c*d**2 - b*d**2*(b + sqrt(-4*a*c + b**2)) + 2*c**2)))/(sqrt(-4*a*c + b**2)*sqrt(2*a*c*d**2 - b*d**2*(b + sqrt(-4*a*c + b**2)) + 2*c**2)) - sqrt(2)*c*atanh(sqrt(2)*(2*c + d**2*x*(b - sqrt(-4*a*c + b**2)))/(2*sqrt(-d**2*x**2 + 1)*sqrt(2*a*c*d**2 - b*d**2*(b - sqrt(-4*a*c + b**2)) + 2*c**2)))/(sqrt(-4*a*c + b**2)*sqrt(2*a*c*d**2 - b*d**2*(b - sqrt(-4*a*c + b**2)) + 2*c**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_781():
    f = 1/(sqrt(-d*x + 1)*sqrt(d*x + 1)*(a + b*x + c*x**2)**2)
    F = sqrt(2)*c*(-2*a*b**2*d**4 + 12*a*c**2*d**2 - b*d**2*(b + sqrt(-4*a*c + b**2))*(-a*d**2 + c) + 4*c**3 - 4*c*d**2*(-2*a**2*d**2 + b**2))*atanh(sqrt(2)*(2*c + d**2*x*(b + sqrt(-4*a*c + b**2)))/(2*sqrt(-d**2*x**2 + 1)*sqrt(2*a*c*d**2 - b*d**2*(b + sqrt(-4*a*c + b**2)) + 2*c**2)))/(2*(-4*a*c + b**2)**(sympy.S(3)/2)*(b**2*d**2 - (a*d**2 + c)**2)*sqrt(2*a*c*d**2 - b*d**2*(b + sqrt(-4*a*c + b**2)) + 2*c**2)) - sqrt(2)*c*(-a*b*d**4*(b + sqrt(-4*a*c + b**2)) + 12*a*c**2*d**2 + 4*c**3 - c*d**2*(-8*a**2*d**2 + 5*b**2 - b*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*(2*c + d**2*x*(b - sqrt(-4*a*c + b**2)))/(2*sqrt(-d**2*x**2 + 1)*sqrt(2*a*c*d**2 - b*d**2*(b - sqrt(-4*a*c + b**2)) + 2*c**2)))/(2*(-4*a*c + b**2)**(sympy.S(3)/2)*(b**2*d**2 - (a*d**2 + c)**2)*sqrt(2*a*c*d**2 - b*d**2*(b - sqrt(-4*a*c + b**2)) + 2*c**2)) - (b*(b**2*d**2 - c*(3*a*d**2 + c)) - c*x*(2*a*c*d**2 - b**2*d**2 + 2*c**2))*sqrt(-d**2*x**2 + 1)/((-4*a*c + b**2)*(b**2*d**2 - (a*d**2 + c)**2)*(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_782():
    f = (a + b*x + c*x**2)**3/((-d*x + 1)**(sympy.S(3)/2)*(d*x + 1)**(sympy.S(3)/2))
    F = b*c**2*x**2*sqrt(-d**2*x**2 + 1)/d**4 + b*sqrt(-d**2*x**2 + 1)*(6*a*c*d**2 + b**2*d**2 + 5*c**2)/d**6 + c**3*x**3*sqrt(-d**2*x**2 + 1)/(4*d**4) + c*x*sqrt(-d**2*x**2 + 1)*(12*a*c*d**2 + 12*b**2*d**2 + 7*c**2)/(8*d**6) + (b*d**4*(3*a**2 + 6*a*c/d**2 + b**2/d**2 + 3*c**2/d**4) + x*(a*d**2 + c)*(a**2*d**4 + 2*a*c*d**2 + 3*b**2*d**2 + c**2))/(d**6*sqrt(-d**2*x**2 + 1)) - (24*a**2*c*d**4 + 24*a*b**2*d**4 + 36*a*c**2*d**2 + 36*b**2*c*d**2 + 15*c**3)*asin(d*x)/(8*d**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_783():
    f = (a + b*x + c*x**2)**2/((-d*x + 1)**(sympy.S(3)/2)*(d*x + 1)**(sympy.S(3)/2))
    F = 2*b*c*sqrt(-d**2*x**2 + 1)/d**4 + c**2*x*sqrt(-d**2*x**2 + 1)/(2*d**4) - (2*b**2 + c*(4*a + 3*c/d**2))*asin(d*x)/(2*d**3) + (2*b*d**2*(a + c/d**2) + x*(a**2*d**4 + 2*a*c*d**2 + b**2*d**2 + c**2))/(d**4*sqrt(-d**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_784():
    f = (a + b*x + c*x**2)/((-d*x + 1)**(sympy.S(3)/2)*(d*x + 1)**(sympy.S(3)/2))
    F = -c*asin(d*x)/d**3 + (b + x*(a*d**2 + c))/(d**2*sqrt(-d**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_785():
    f = 1/((-d*x + 1)**(sympy.S(3)/2)*(d*x + 1)**(sympy.S(3)/2)*(a + b*x + c*x**2))
    F = -sqrt(2)*c*(2*a*c*d**2 - b*d**2*(b - sqrt(-4*a*c + b**2)) + 2*c**2)*atanh(sqrt(2)*(2*c + d**2*x*(b + sqrt(-4*a*c + b**2)))/(2*sqrt(-d**2*x**2 + 1)*sqrt(2*a*c*d**2 - b*d**2*(b + sqrt(-4*a*c + b**2)) + 2*c**2)))/(2*sqrt(-4*a*c + b**2)*(b**2*d**2 - (a*d**2 + c)**2)*sqrt(2*a*c*d**2 - b*d**2*(b + sqrt(-4*a*c + b**2)) + 2*c**2)) + sqrt(2)*c*(2*a*c*d**2 - b*d**2*(b + sqrt(-4*a*c + b**2)) + 2*c**2)*atanh(sqrt(2)*(2*c + d**2*x*(b - sqrt(-4*a*c + b**2)))/(2*sqrt(-d**2*x**2 + 1)*sqrt(2*a*c*d**2 - b*d**2*(b - sqrt(-4*a*c + b**2)) + 2*c**2)))/(2*sqrt(-4*a*c + b**2)*(b**2*d**2 - (a*d**2 + c)**2)*sqrt(2*a*c*d**2 - b*d**2*(b - sqrt(-4*a*c + b**2)) + 2*c**2)) + d**2*(b - x*(a*d**2 + c))/((b**2*d**2 - (a*d**2 + c)**2)*sqrt(-d**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_786():
    f = 1/((-d*x + 1)**(sympy.S(3)/2)*(d*x + 1)**(sympy.S(3)/2)*(a + b*x + c*x**2)**2)
    F = sqrt(2)*c*(3*a*b**3*d**6*(b + sqrt(-4*a*c + b**2)) + 24*a*c**4*d**2 - 2*a*c**2*d**4*(-8*a**2*d**2 + 7*b**2 + 5*b*sqrt(-4*a*c + b**2)) + b*c*d**4*(-17*a**2*b*d**2 - 11*a**2*d**2*sqrt(-4*a*c + b**2) + 2*b**3 + 2*b**2*sqrt(-4*a*c + b**2)) + 4*c**5 - c**3*d**2*(-36*a**2*d**2 + 9*b**2 - b*sqrt(-4*a*c + b**2)))*atanh(sqrt(2)*(2*c + d**2*x*(b - sqrt(-4*a*c + b**2)))/(2*sqrt(-d**2*x**2 + 1)*sqrt(2*a*c*d**2 - b*d**2*(b - sqrt(-4*a*c + b**2)) + 2*c**2)))/(2*(-4*a*c + b**2)**(sympy.S(3)/2)*sqrt(2*a*c*d**2 - b*d**2*(b - sqrt(-4*a*c + b**2)) + 2*c**2)*(a**2*d**4 + 2*a*c*d**2 - b**2*d**2 + c**2)**2) + sqrt(2)*c*(-6*a*b**4*d**8 - 24*a*c**4*d**4 - 4*b**2*c*d**6*(-7*a**2*d**2 + b**2) + b*d**4*(b + sqrt(-4*a*c + b**2))*(-11*a**2*c*d**4 + 3*a*b**2*d**4 - 10*a*c**2*d**2 + 2*b**2*c*d**2 + c**3) - 4*c**5*d**2 + 2*c**3*(-18*a**2*d**6 + 4*b**2*d**4) + 8*c**2*(-2*a**3*d**8 + 3*a*b**2*d**6))*atanh(sqrt(2)*(2*c + d**2*x*(b + sqrt(-4*a*c + b**2)))/(2*sqrt(-d**2*x**2 + 1)*sqrt(2*a*c*d**2 - b*d**2*(b + sqrt(-4*a*c + b**2)) + 2*c**2)))/(2*d**2*(-4*a*c + b**2)**(sympy.S(3)/2)*sqrt(2*a*c*d**2 - b*d**2*(b + sqrt(-4*a*c + b**2)) + 2*c**2)*(a**2*d**4 + 2*a*c*d**2 - b**2*d**2 + c**2)**2) - d**2*(b*(-11*a**2*c*d**4 + 3*a*b**2*d**4 - 10*a*c**2*d**2 + 2*b**2*c*d**2 + c**3) - x*(b**2*d**4*(a**2*d**2 + 2*b**2) + 2*c**4 - c**2*d**2*(6*a**2*d**2 + b**2) - c*(4*a**3*d**6 + 6*a*b**2*d**4)))/((-4*a*c + b**2)*sqrt(-d**2*x**2 + 1)*(a*d**2 - b*d + c)**2*(a*d**2 + b*d + c)**2) - (b*(b**2*d**2 - c*(3*a*d**2 + c)) - c*x*(2*a*c*d**2 - b**2*d**2 + 2*c**2))/((-4*a*c + b**2)*(b**2*d**2 - (a*d**2 + c)**2)*sqrt(-d**2*x**2 + 1)*(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_787():
    f = (a + c*x**2)**p*(-e*x + 1)**m*(e*x + 1)**m
    F = x*(a + c*x**2)**p*appellf1(sympy.S.Half, -m, -p, sympy.S(3)/2, e**2*x**2, -c*x**2/a)/(1 + c*x**2/a)**p
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_788():
    f = (a + c*x**2)**p*(d - e*x)**m*(d + e*x)**m
    F = x*(a + c*x**2)**p*(d - e*x)**m*(d + e*x)**m*appellf1(sympy.S.Half, -m, -p, sympy.S(3)/2, e**2*x**2/d**2, -c*x**2/a)/((1 + c*x**2/a)**p*(1 - e**2*x**2/d**2)**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_789():
    f = (a + c*x**2)**p*(d + e*x)**m*(d*f - e*f*x)**m
    F = x*(a + c*x**2)**p*(d + e*x)**m*(d*f - e*f*x)**m*appellf1(sympy.S.Half, -m, -p, sympy.S(3)/2, e**2*x**2/d**2, -c*x**2/a)/((1 + c*x**2/a)**p*(1 - e**2*x**2/d**2)**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_790():
    f = (d + e*x)**3*(f + g*x)**n*(a + 2*c*d*x + c*e*x**2)
    F = c*e**4*(f + g*x)**(n + 6)/(g**6*(n + 6)) - 5*c*e**3*(f + g*x)**(n + 5)*(-d*g + e*f)/(g**6*(n + 5)) + e**2*(f + g*x)**(n + 4)*(a*e*g**2 + c*(9*d**2*g**2 - 20*d*e*f*g + 10*e**2*f**2))/(g**6*(n + 4)) - e*(f + g*x)**(n + 3)*(-d*g + e*f)*(3*a*e*g**2 + c*(7*d**2*g**2 - 20*d*e*f*g + 10*e**2*f**2))/(g**6*(n + 3)) - (f + g*x)**(n + 1)*(a*g**2 + c*f*(-2*d*g + e*f))*(-d*g + e*f)**3/(g**6*(n + 1)) + (f + g*x)**(n + 2)*(-d*g + e*f)**2*(3*a*e*g**2 + c*(2*d**2*g**2 - 10*d*e*f*g + 5*e**2*f**2))/(g**6*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_791():
    f = (d + e*x)**2*(f + g*x)**n*(a + 2*c*d*x + c*e*x**2)
    F = c*e**3*(f + g*x)**(n + 5)/(g**5*(n + 5)) - 4*c*e**2*(f + g*x)**(n + 4)*(-d*g + e*f)/(g**5*(n + 4)) + e*(f + g*x)**(n + 3)*(a*e*g**2 + c*(5*d**2*g**2 - 12*d*e*f*g + 6*e**2*f**2))/(g**5*(n + 3)) + (f + g*x)**(n + 1)*(a*g**2 + c*f*(-2*d*g + e*f))*(-d*g + e*f)**2/(g**5*(n + 1)) - (f + g*x)**(n + 2)*(-2*d*g + 2*e*f)*(a*e*g**2 + c*(d**2*g**2 - 4*d*e*f*g + 2*e**2*f**2))/(g**5*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_792():
    f = (d + e*x)*(f + g*x)**n*(a + 2*c*d*x + c*e*x**2)
    F = c*e**2*(f + g*x)**(n + 4)/(g**4*(n + 4)) - 3*c*e*(f + g*x)**(n + 3)*(-d*g + e*f)/(g**4*(n + 3)) - (f + g*x)**(n + 1)*(a*g**2 + c*f*(-2*d*g + e*f))*(-d*g + e*f)/(g**4*(n + 1)) + (f + g*x)**(n + 2)*(a*e*g**2 + c*(2*d**2*g**2 - 6*d*e*f*g + 3*e**2*f**2))/(g**4*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_793():
    f = (f + g*x)**n*(a + 2*c*d*x + c*e*x**2)
    F = c*e*(f + g*x)**(n + 3)/(g**3*(n + 3)) - 2*c*(f + g*x)**(n + 2)*(-d*g + e*f)/(g**3*(n + 2)) + (f + g*x)**(n + 1)*(a*g**2 + c*f*(-2*d*g + e*f))/(g**3*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_794():
    f = (f + g*x)**n*(a + 2*c*d*x + c*e*x**2)/(d + e*x)
    F = c*(f + g*x)**(n + 2)/(g**2*(n + 2)) - c*(f + g*x)**(n + 1)*(-d*g + e*f)/(e*g**2*(n + 1)) + (f + g*x)**(n + 1)*(-a*e + c*d**2)*hyper((1, n + 1), (n + 2,), e*(f + g*x)/(-d*g + e*f))/(e*(n + 1)*(-d*g + e*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_795():
    f = (f + g*x)**n*(a + 2*c*d*x + c*e*x**2)/(d + e*x)**2
    F = c*(f + g*x)**(n + 1)/(e*g*(n + 1)) - g*(f + g*x)**(n + 1)*(-a*e + c*d**2)*hyper((2, n + 1), (n + 2,), e*(f + g*x)/(-d*g + e*f))/(e*(n + 1)*(-d*g + e*f)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_796():
    f = (f + g*x)**n*(a + 2*c*d*x + c*e*x**2)/(d + e*x)**3
    F = -(a - c*d**2/e)*(f + g*x)**(n + 1)/((d + e*x)**2*(-2*d*g + 2*e*f)) - g*(1 - n)*(f + g*x)**(n + 1)*(-a*e + c*d**2)/(2*e*(d + e*x)*(-d*g + e*f)**2) + (f + g*x)**(n + 1)*(a*e*g**2*n*(1 - n) - c*(d**2*g**2*(-n**2 + n + 2) - 4*d*e*f*g + 2*e**2*f**2))*hyper((1, n + 1), (n + 2,), e*(f + g*x)/(-d*g + e*f))/(2*e*(n + 1)*(-d*g + e*f)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_797():
    f = (f + g*x)**n*(a + 2*c*d*x + c*e*x**2)/(d + e*x)**4
    F = -(a - c*d**2/e)*(f + g*x)**(n + 1)/((d + e*x)**3*(-3*d*g + 3*e*f)) - g*(2 - n)*(f + g*x)**(n + 1)*(-a*e + c*d**2)/(6*e*(d + e*x)**2*(-d*g + e*f)**2) + g*(f + g*x)**(n + 1)*(a*e*g**2*(n**2 - 3*n + 2) + c*(d**2*g**2*(-n**2 + 3*n + 4) - 12*d*e*f*g + 6*e**2*f**2))*hyper((2, n + 1), (n + 2,), e*(f + g*x)/(-d*g + e*f))/(6*e*(n + 1)*(-d*g + e*f)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_798():
    f = (d + e*x)**m*(f + g*x)**n*(a + 2*c*d*x + c*e*x**2)
    F = c*(d + e*x)**(m + 2)*(f + g*x)**(n + 1)/(e*g*(m + n + 3)) - c*(d + e*x)**(m + 1)*(f + g*x)**(n + 1)*(m + 2)*(-d*g + e*f)/(e*g**2*(m + n + 2)*(m + n + 3)) + (d + e*x)**(m + 1)*(f + g*x)**n*(c*(m + 2)*(-d*g + e*f)*(d*g*(n + 1) + e*f*(m + 1)) + g*(a*e*g*(m + n + 3) - c*d*(d*g*(n + 1) + e*f*(m + 2)))*(m + n + 2))*hyper((-n, m + 1), (m + 2,), -g*(d + e*x)/(-d*g + e*f))/(e**2*g**2*(e*(f + g*x)/(-d*g + e*f))**n*(m + 1)*(m + n + 2)*(m + n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_799():
    f = (a + b*x + c*x**2)/((d + e*x)*(f + g*x))
    F = c*x/(e*g) - (a*g**2 - b*f*g + c*f**2)*log(f + g*x)/(g**2*(-d*g + e*f)) + (a*e**2 - b*d*e + c*d**2)*log(d + e*x)/(e**2*(-d*g + e*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_800():
    f = (a + b*x + c*x**2)**2/((d + e*x)*(f + g*x))
    F = c**2*x**3/(3*e*g) - c*x**2*(-2*b*e*g + c*d*g + c*e*f)/(2*e**2*g**2) - (a*g**2 - b*f*g + c*f**2)**2*log(f + g*x)/(g**4*(-d*g + e*f)) + x*(b**2*e**2*g**2 + c**2*(d**2*g**2 + d*e*f*g + e**2*f**2) - 2*c*e*g*(-a*e*g + b*d*g + b*e*f))/(e**3*g**3) + (a*e**2 - b*d*e + c*d**2)**2*log(d + e*x)/(e**4*(-d*g + e*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_801():
    f = (a + b*x + c*x**2)**3/((d + e*x)*(f + g*x))
    F = c**3*x**5/(5*e*g) - c**2*x**4*(-3*b*e*g + c*d*g + c*e*f)/(4*e**2*g**2) + c*x**3*(3*b**2*e**2*g**2 + c**2*(d**2*g**2 + d*e*f*g + e**2*f**2) - 3*c*e*g*(-a*e*g + b*d*g + b*e*f))/(3*e**3*g**3) - (a*g**2 - b*f*g + c*f**2)**3*log(f + g*x)/(g**6*(-d*g + e*f)) + x**2*(b**3*e**3*g**3 - 3*b*c*e**2*g**2*(-2*a*e*g + b*d*g + b*e*f) - c**3*(d**3*g**3 + d**2*e*f*g**2 + d*e**2*f**2*g + e**3*f**3) - 3*c**2*e*g*(a*e*g*(d*g + e*f) - b*(d**2*g**2 + d*e*f*g + e**2*f**2)))/(2*e**4*g**4) - x*(b**2*e**3*g**3*(-3*a*e*g + b*d*g + b*e*f) - c**3*(d**4*g**4 + d**3*e*f*g**3 + d**2*e**2*f**2*g**2 + d*e**3*f**3*g + e**4*f**4) - 3*c**2*e*g*(a*e*g*(d**2*g**2 + d*e*f*g + e**2*f**2) - b*(d**3*g**3 + d**2*e*f*g**2 + d*e**2*f**2*g + e**3*f**3)) - 3*c*e**2*g**2*(a**2*e**2*g**2 - 2*a*b*e*g*(d*g + e*f) + b**2*(d**2*g**2 + d*e*f*g + e**2*f**2)))/(e**5*g**5) + (a*e**2 - b*d*e + c*d**2)**3*log(d + e*x)/(e**6*(-d*g + e*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_802():
    f = 1/((d + e*x)*(f + g*x)*(a + b*x + c*x**2))
    F = e**2*log(d + e*x)/((-d*g + e*f)*(a*e**2 - b*d*e + c*d**2)) - g**2*log(f + g*x)/((-d*g + e*f)*(a*g**2 - b*f*g + c*f**2)) - (-b*e*g + c*d*g + c*e*f)*log(a + b*x + c*x**2)/((c*f**2 - g*(-a*g + b*f))*(2*a*e**2 - 2*b*d*e + 2*c*d**2)) - (b**2*e*g + 2*c**2*d*f - c*(2*a*e*g + b*d*g + b*e*f))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(sqrt(-4*a*c + b**2)*(c*f**2 - g*(-a*g + b*f))*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_803():
    f = 1/((d + e*x)*(f + g*x)*(a + b*x + c*x**2)**2)
    F = 2*c*(b**2*e*g + 2*c**2*d*f - c*(2*a*e*g + b*d*g + b*e*f))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/((-4*a*c + b**2)**(sympy.S(3)/2)*(c*f**2 - g*(-a*g + b*f))*(a*e**2 - b*d*e + c*d**2)) + e**4*log(d + e*x)/((-d*g + e*f)*(a*e**2 - b*d*e + c*d**2)**2) - g**4*log(f + g*x)/((-d*g + e*f)*(a*g**2 - b*f*g + c*f**2)**2) - (c*(d**2*g**2 + e**2*f**2) + e*g*(2*a*e*g - b*(d*g + e*f)))*(-b*e*g + c*d*g + c*e*f)*log(a + b*x + c*x**2)/(2*(c*f**2 - g*(-a*g + b*f))**2*(a*e**2 - b*d*e + c*d**2)**2) - (2*a*c**2*(d*g + e*f) + b**3*e*g - b**2*c*(d*g + e*f) + b*c*(-3*a*e*g + c*d*f) + c*x*(b**2*e*g + 2*c**2*d*f - c*(2*a*e*g + b*d*g + b*e*f)))/((-4*a*c + b**2)*(c*f**2 - g*(-a*g + b*f))*(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2)) + (b**2*e**2*g**2*(-2*a*e*g + b*d*g + b*e*f) - 2*c**3*d*f*(d**2*g**2 + d*e*f*g + e**2*f**2) - c**2*(4*a*d*e**2*f*g**2 - b*(d**3*g**3 + 5*d**2*e*f*g**2 + 5*d*e**2*f**2*g + e**3*f**3)) + 2*c*e*g*(a**2*e**2*g**2 + a*b*e*g*(d*g + e*f) - b**2*(d*g + e*f)**2))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(sqrt(-4*a*c + b**2)*(c*f**2 - g*(-a*g + b*f))**2*(a*e**2 - b*d*e + c*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_804():
    f = (d + e*x)**3*(a + b*x + c*x**2)/sqrt(f + g*x)
    F = 2*c*e**3*(f + g*x)**(sympy.S(11)/2)/(11*g**6) - 2*e**2*(f + g*x)**(sympy.S(9)/2)*(-b*e*g - 3*c*d*g + 5*c*e*f)/(9*g**6) - 2*e*(f + g*x)**(sympy.S(7)/2)*(-c*(3*d**2*g**2 - 12*d*e*f*g + 10*e**2*f**2) + e*g*(-a*e*g - 3*b*d*g + 4*b*e*f))/(7*g**6) + (f + g*x)**(sympy.S(5)/2)*(-c*(d**2*g**2 - 8*d*e*f*g + 10*e**2*f**2) + 3*e*g*(-a*e*g - b*d*g + 2*b*e*f))*(-2*d*g + 2*e*f)/(5*g**6) + 2*(f + g*x)**(sympy.S(3)/2)*(-d*g + e*f)**2*(c*f*(-2*d*g + 5*e*f) - g*(-3*a*e*g - b*d*g + 4*b*e*f))/(3*g**6) - 2*sqrt(f + g*x)*(-d*g + e*f)**3*(a*g**2 - b*f*g + c*f**2)/g**6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_805():
    f = (d + e*x)**2*(a + b*x + c*x**2)/sqrt(f + g*x)
    F = 2*c*e**2*(f + g*x)**(sympy.S(9)/2)/(9*g**5) - 2*e*(f + g*x)**(sympy.S(7)/2)*(-b*e*g - 2*c*d*g + 4*c*e*f)/(7*g**5) - (f + g*x)**(sympy.S(5)/2)*(-2*c*(d**2*g**2 - 6*d*e*f*g + 6*e**2*f**2) + 2*e*g*(-a*e*g - 2*b*d*g + 3*b*e*f))/(5*g**5) - (f + g*x)**(sympy.S(3)/2)*(-2*d*g + 2*e*f)*(2*c*f*(-d*g + 2*e*f) - g*(-2*a*e*g - b*d*g + 3*b*e*f))/(3*g**5) + 2*sqrt(f + g*x)*(-d*g + e*f)**2*(a*g**2 - b*f*g + c*f**2)/g**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_806():
    f = (d + e*x)*(a + b*x + c*x**2)/sqrt(f + g*x)
    F = 2*c*e*(f + g*x)**(sympy.S(7)/2)/(7*g**4) - (f + g*x)**(sympy.S(5)/2)*(-2*b*e*g - 2*c*d*g + 6*c*e*f)/(5*g**4) + (f + g*x)**(sympy.S(3)/2)*(2*c*f*(-2*d*g + 3*e*f) - 2*g*(-a*e*g - b*d*g + 2*b*e*f))/(3*g**4) - sqrt(f + g*x)*(-2*d*g + 2*e*f)*(a*g**2 - b*f*g + c*f**2)/g**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_807():
    f = (a + b*x + c*x**2)/sqrt(f + g*x)
    F = 2*c*(f + g*x)**(sympy.S(5)/2)/(5*g**3) - (f + g*x)**(sympy.S(3)/2)*(-2*b*g + 4*c*f)/(3*g**3) + sqrt(f + g*x)*(2*a*g**2 - 2*b*f*g + 2*c*f**2)/g**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_808():
    f = (a + b*x + c*x**2)/((d + e*x)*sqrt(f + g*x))
    F = 2*c*(f + g*x)**(sympy.S(3)/2)/(3*e*g**2) + sqrt(f + g*x)*(2*b*e*g - 2*c*(d*g + e*f))/(e**2*g**2) - (2*a*e**2 - 2*b*d*e + 2*c*d**2)*atanh(sqrt(e)*sqrt(f + g*x)/sqrt(-d*g + e*f))/(e**(sympy.S(5)/2)*sqrt(-d*g + e*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_809():
    f = (a + b*x + c*x**2)/((d + e*x)**2*sqrt(f + g*x))
    F = 2*c*sqrt(f + g*x)/(e**2*g) - (a + d*(-b*e + c*d)/e**2)*sqrt(f + g*x)/((d + e*x)*(-d*g + e*f)) + (c*d*(-3*d*g + 4*e*f) - e*(-a*e*g - b*d*g + 2*b*e*f))*atanh(sqrt(e)*sqrt(f + g*x)/sqrt(-d*g + e*f))/(e**(sympy.S(5)/2)*(-d*g + e*f)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_810():
    f = (a + b*x + c*x**2)/((d + e*x)**3*sqrt(f + g*x))
    F = -(a + d*(-b*e + c*d)/e**2)*sqrt(f + g*x)/((d + e*x)**2*(-2*d*g + 2*e*f)) + sqrt(f + g*x)*(c*d*(-5*d*g + 8*e*f) - e*(-3*a*e*g - b*d*g + 4*b*e*f))/(4*e**2*(d + e*x)*(-d*g + e*f)**2) + (-c*(3*d**2*g**2 - 8*d*e*f*g + 8*e**2*f**2) + e*g*(-3*a*e*g - b*d*g + 4*b*e*f))*atanh(sqrt(e)*sqrt(f + g*x)/sqrt(-d*g + e*f))/(4*e**(sympy.S(5)/2)*(-d*g + e*f)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_811():
    f = (d + e*x)**3*(a + b*x + c*x**2)/(f + g*x)**(sympy.S(3)/2)
    F = 2*c*e**3*(f + g*x)**(sympy.S(9)/2)/(9*g**6) - 2*e**2*(f + g*x)**(sympy.S(7)/2)*(-b*e*g - 3*c*d*g + 5*c*e*f)/(7*g**6) - 2*e*(f + g*x)**(sympy.S(5)/2)*(-c*(3*d**2*g**2 - 12*d*e*f*g + 10*e**2*f**2) + e*g*(-a*e*g - 3*b*d*g + 4*b*e*f))/(5*g**6) + (f + g*x)**(sympy.S(3)/2)*(-c*(d**2*g**2 - 8*d*e*f*g + 10*e**2*f**2) + 3*e*g*(-a*e*g - b*d*g + 2*b*e*f))*(-2*d*g + 2*e*f)/(3*g**6) + 2*sqrt(f + g*x)*(-d*g + e*f)**2*(c*f*(-2*d*g + 5*e*f) - g*(-3*a*e*g - b*d*g + 4*b*e*f))/g**6 + 2*(-d*g + e*f)**3*(a*g**2 - b*f*g + c*f**2)/(g**6*sqrt(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_812():
    f = (d + e*x)**2*(a + b*x + c*x**2)/(f + g*x)**(sympy.S(3)/2)
    F = 2*c*e**2*(f + g*x)**(sympy.S(7)/2)/(7*g**5) - 2*e*(f + g*x)**(sympy.S(5)/2)*(-b*e*g - 2*c*d*g + 4*c*e*f)/(5*g**5) - (f + g*x)**(sympy.S(3)/2)*(-2*c*(d**2*g**2 - 6*d*e*f*g + 6*e**2*f**2) + 2*e*g*(-a*e*g - 2*b*d*g + 3*b*e*f))/(3*g**5) - sqrt(f + g*x)*(-2*d*g + 2*e*f)*(2*c*f*(-d*g + 2*e*f) - g*(-2*a*e*g - b*d*g + 3*b*e*f))/g**5 - 2*(-d*g + e*f)**2*(a*g**2 - b*f*g + c*f**2)/(g**5*sqrt(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_813():
    f = (d + e*x)*(a + b*x + c*x**2)/(f + g*x)**(sympy.S(3)/2)
    F = 2*c*e*(f + g*x)**(sympy.S(5)/2)/(5*g**4) - (f + g*x)**(sympy.S(3)/2)*(-2*b*e*g - 2*c*d*g + 6*c*e*f)/(3*g**4) + sqrt(f + g*x)*(2*c*f*(-2*d*g + 3*e*f) - 2*g*(-a*e*g - b*d*g + 2*b*e*f))/g**4 + (-2*d*g + 2*e*f)*(a*g**2 - b*f*g + c*f**2)/(g**4*sqrt(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_814():
    f = (a + b*x + c*x**2)/(f + g*x)**(sympy.S(3)/2)
    F = 2*c*(f + g*x)**(sympy.S(3)/2)/(3*g**3) - sqrt(f + g*x)*(-2*b*g + 4*c*f)/g**3 - (2*a*g**2 - 2*b*f*g + 2*c*f**2)/(g**3*sqrt(f + g*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_815():
    f = (a + b*x + c*x**2)/((d + e*x)*(f + g*x)**(sympy.S(3)/2))
    F = 2*c*sqrt(f + g*x)/(e*g**2) + (2*a*g**2 - 2*b*f*g + 2*c*f**2)/(g**2*sqrt(f + g*x)*(-d*g + e*f)) - (2*a*e**2 - 2*b*d*e + 2*c*d**2)*atanh(sqrt(e)*sqrt(f + g*x)/sqrt(-d*g + e*f))/(e**(sympy.S(3)/2)*(-d*g + e*f)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_816():
    f = (a + b*x + c*x**2)/((d + e*x)**2*(f + g*x)**(sympy.S(3)/2))
    F = -(2*a*g**2 - 2*b*f*g + 2*c*f**2)/(g*sqrt(f + g*x)*(-d*g + e*f)**2) - sqrt(f + g*x)*(a*e**2 - b*d*e + c*d**2)/(e*(d + e*x)*(-d*g + e*f)**2) + (c*d*(-d*g + 4*e*f) - e*(-3*a*e*g + b*d*g + 2*b*e*f))*atanh(sqrt(e)*sqrt(f + g*x)/sqrt(-d*g + e*f))/(e**(sympy.S(3)/2)*(-d*g + e*f)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_817():
    f = (a + b*x + c*x**2)/((d + e*x)**3*(f + g*x)**(sympy.S(3)/2))
    F = (2*a*g**2 - 2*b*f*g + 2*c*f**2)/(sqrt(f + g*x)*(-d*g + e*f)**3) + sqrt(f + g*x)*(c*d*(-d*g + 8*e*f) - e*(-7*a*e*g + 3*b*d*g + 4*b*e*f))/(4*e*(d + e*x)*(-d*g + e*f)**3) - sqrt(f + g*x)*(a*e**2 - b*d*e + c*d**2)/(2*e*(d + e*x)**2*(-d*g + e*f)**2) - (c*(-d**2*g**2 + 8*d*e*f*g + 8*e**2*f**2) + 3*e*g*(5*a*e*g - b*(d*g + 4*e*f)))*atanh(sqrt(e)*sqrt(f + g*x)/sqrt(-d*g + e*f))/(4*e**(sympy.S(3)/2)*(-d*g + e*f)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_818():
    f = sqrt(x - 1)*sqrt(x + 1)/(-x**2 + x + 1)
    F = -acosh(x) + sqrt(sympy.S(-2)/5 + 2*sqrt(5)/5)*atan(sqrt(x + 1)/(sqrt(-2 + sqrt(5))*sqrt(x - 1))) + sqrt(sympy.S(2)/5 + 2*sqrt(5)/5)*atanh(sqrt(x + 1)/(sqrt(2 + sqrt(5))*sqrt(x - 1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_819():
    f = (a + b*x + c*x**2)/(sqrt(d + e*x)*sqrt(f + g*x))
    F = c*(d + e*x)**(sympy.S(3)/2)*sqrt(f + g*x)/(2*e**2*g) - sqrt(d + e*x)*sqrt(f + g*x)*(-4*b*e*g + 5*c*d*g + 3*c*e*f)/(4*e**2*g**2) + (c*(3*d**2*g**2 + 2*d*e*f*g + 3*e**2*f**2) + 4*e*g*(2*a*e*g - b*(d*g + e*f)))*atanh(sqrt(g)*sqrt(d + e*x)/(sqrt(e)*sqrt(f + g*x)))/(4*e**(sympy.S(5)/2)*g**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_820():
    f = (d + e*x)**(sympy.S(3)/2)*(a + b*x + c*x**2)/sqrt(f + g*x)
    F = c*(d + e*x)**(sympy.S(7)/2)*sqrt(f + g*x)/(4*e**2*g) - (d + e*x)**(sympy.S(5)/2)*sqrt(f + g*x)*(-8*b*e*g + 9*c*d*g + 7*c*e*f)/(24*e**2*g**2) + (d + e*x)**(sympy.S(3)/2)*sqrt(f + g*x)*(c*(3*d**2*g**2 + 10*d*e*f*g + 35*e**2*f**2) + 8*e*g*(6*a*e*g - b*(d*g + 5*e*f)))/(96*e**2*g**3) - sqrt(d + e*x)*sqrt(f + g*x)*(c*(3*d**2*g**2 + 10*d*e*f*g + 35*e**2*f**2) + 8*e*g*(6*a*e*g - b*(d*g + 5*e*f)))*(-d*g + e*f)/(64*e**2*g**4) + (c*(3*d**2*g**2 + 10*d*e*f*g + 35*e**2*f**2) + 8*e*g*(6*a*e*g - b*(d*g + 5*e*f)))*(-d*g + e*f)**2*atanh(sqrt(g)*sqrt(d + e*x)/(sqrt(e)*sqrt(f + g*x)))/(64*e**(sympy.S(5)/2)*g**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_821():
    f = sqrt(d + e*x)*(a + b*x + c*x**2)/sqrt(f + g*x)
    F = c*(d + e*x)**(sympy.S(5)/2)*sqrt(f + g*x)/(3*e**2*g) - (d + e*x)**(sympy.S(3)/2)*sqrt(f + g*x)*(-6*b*e*g + 7*c*d*g + 5*c*e*f)/(12*e**2*g**2) + sqrt(d + e*x)*sqrt(f + g*x)*(c*(d**2*g**2 + 2*d*e*f*g + 5*e**2*f**2) + 2*e*g*(4*a*e*g - b*(d*g + 3*e*f)))/(8*e**2*g**3) - (c*(d**2*g**2 + 2*d*e*f*g + 5*e**2*f**2) + 2*e*g*(4*a*e*g - b*(d*g + 3*e*f)))*(-d*g + e*f)*atanh(sqrt(g)*sqrt(d + e*x)/(sqrt(e)*sqrt(f + g*x)))/(8*e**(sympy.S(5)/2)*g**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_822():
    f = (a + b*x + c*x**2)/(sqrt(d + e*x)*sqrt(f + g*x))
    F = c*(d + e*x)**(sympy.S(3)/2)*sqrt(f + g*x)/(2*e**2*g) - sqrt(d + e*x)*sqrt(f + g*x)*(-4*b*e*g + 5*c*d*g + 3*c*e*f)/(4*e**2*g**2) + (c*(3*d**2*g**2 + 2*d*e*f*g + 3*e**2*f**2) + 4*e*g*(2*a*e*g - b*(d*g + e*f)))*atanh(sqrt(g)*sqrt(d + e*x)/(sqrt(e)*sqrt(f + g*x)))/(4*e**(sympy.S(5)/2)*g**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_823():
    f = (a + b*x + c*x**2)/((d + e*x)**(sympy.S(3)/2)*sqrt(f + g*x))
    F = c*sqrt(d + e*x)*sqrt(f + g*x)/(e**2*g) - (2*a + 2*d*(-b*e + c*d)/e**2)*sqrt(f + g*x)/(sqrt(d + e*x)*(-d*g + e*f)) - (-2*b*e*g + 3*c*d*g + c*e*f)*atanh(sqrt(g)*sqrt(d + e*x)/(sqrt(e)*sqrt(f + g*x)))/(e**(sympy.S(5)/2)*g**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_824():
    f = (a + b*x + c*x**2)/((d + e*x)**(sympy.S(5)/2)*sqrt(f + g*x))
    F = 2*c*atanh(sqrt(g)*sqrt(d + e*x)/(sqrt(e)*sqrt(f + g*x)))/(e**(sympy.S(5)/2)*sqrt(g)) - (2*a + 2*d*(-b*e + c*d)/e**2)*sqrt(f + g*x)/((d + e*x)**(sympy.S(3)/2)*(-3*d*g + 3*e*f)) + sqrt(f + g*x)*(2*c*(-4*d**2*g + 6*d*e*f) - 2*e*(-2*a*e*g - b*d*g + 3*b*e*f))/(3*e**2*sqrt(d + e*x)*(-d*g + e*f)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_825():
    f = (a + b*x + c*x**2)/((d + e*x)**(sympy.S(7)/2)*sqrt(f + g*x))
    F = -(2*a + 2*d*(-b*e + c*d)/e**2)*sqrt(f + g*x)/((d + e*x)**(sympy.S(5)/2)*(-5*d*g + 5*e*f)) + sqrt(f + g*x)*(-2*c*(3*d**2*g**2 - 10*d*e*f*g + 15*e**2*f**2) + 4*e*g*(-4*a*e*g - b*d*g + 5*b*e*f))/(15*e**2*sqrt(d + e*x)*(-d*g + e*f)**3) + sqrt(f + g*x)*(4*c*d*(-3*d*g + 5*e*f) - 2*e*(-4*a*e*g - b*d*g + 5*b*e*f))/(15*e**2*(d + e*x)**(sympy.S(3)/2)*(-d*g + e*f)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_826():
    f = (a + b*x + c*x**2)/((d + e*x)**(sympy.S(9)/2)*sqrt(f + g*x))
    F = -(2*a + 2*d*(-b*e + c*d)/e**2)*sqrt(f + g*x)/((d + e*x)**(sympy.S(7)/2)*(-7*d*g + 7*e*f)) - 4*g*sqrt(f + g*x)*(-c*(3*d**2*g**2 - 14*d*e*f*g + 35*e**2*f**2) + 4*e*g*(-6*a*e*g - b*d*g + 7*b*e*f))/(105*e**2*sqrt(d + e*x)*(-d*g + e*f)**4) + sqrt(f + g*x)*(-2*c*(3*d**2*g**2 - 14*d*e*f*g + 35*e**2*f**2) + 8*e*g*(-6*a*e*g - b*d*g + 7*b*e*f))/(105*e**2*(d + e*x)**(sympy.S(3)/2)*(-d*g + e*f)**3) + sqrt(f + g*x)*(4*c*d*(-4*d*g + 7*e*f) - 2*e*(-6*a*e*g - b*d*g + 7*b*e*f))/(35*e**2*(d + e*x)**(sympy.S(5)/2)*(-d*g + e*f)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_827():
    f = sqrt(d + e*x)*(a + b*x + c*x**2)/(e + f*x)**(sympy.S(3)/2)
    F = c*(d + e*x)**(sympy.S(3)/2)*sqrt(e + f*x)/(2*e*f**2) + (2*a + 2*e*(-b*f + c*e)/f**2)*(d + e*x)**(sympy.S(3)/2)/(sqrt(e + f*x)*(-d*f + e**2)) + sqrt(d + e*x)*sqrt(e + f*x)*(-c*(-d**2*f**2 - 6*d*e**2*f + 15*e**4) + 4*e*f*(-2*a*e*f - b*d*f + 3*b*e**2))/(4*e*f**3*(-d*f + e**2)) - (-c*(-d**2*f**2 - 6*d*e**2*f + 15*e**4) + 4*e*f*(-2*a*e*f - b*d*f + 3*b*e**2))*atanh(sqrt(f)*sqrt(d + e*x)/(sqrt(e)*sqrt(e + f*x)))/(4*e**(sympy.S(3)/2)*f**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_828():
    f = (d + e*x)**(sympy.S(3)/2)*(15*d**2 + 20*d*e*x + 8*e**2*x**2)/sqrt(a + b*x)
    F = 2*e*(a + b*x)**(sympy.S(3)/2)*(d + e*x)**(sympy.S(5)/2)/b**2 + sqrt(a + b*x)*(d + e*x)**(sympy.S(5)/2)*(-13*a*e + 17*b*d)/(3*b**2) + sqrt(a + b*x)*(d + e*x)**(sympy.S(3)/2)*(35*a**2*e**2 - 90*a*b*d*e + 73*b**2*d**2)/(12*b**3) + sqrt(a + b*x)*sqrt(d + e*x)*(-a*e + b*d)*(35*a**2*e**2 - 90*a*b*d*e + 73*b**2*d**2)/(8*b**4) + (-a*e + b*d)**2*(35*a**2*e**2 - 90*a*b*d*e + 73*b**2*d**2)*atanh(sqrt(e)*sqrt(a + b*x)/(sqrt(b)*sqrt(d + e*x)))/(8*b**(sympy.S(9)/2)*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_829():
    f = sqrt(d + e*x)*(15*d**2 + 20*d*e*x + 8*e**2*x**2)/sqrt(a + b*x)
    F = 8*e*(a + b*x)**(sympy.S(3)/2)*(d + e*x)**(sympy.S(3)/2)/(3*b**2) + sqrt(a + b*x)*(d + e*x)**(sympy.S(3)/2)*(-6*a*e + 8*b*d)/b**2 + sqrt(a + b*x)*sqrt(d + e*x)*(5*a**2*e**2 - 13*a*b*d*e + 11*b**2*d**2)/b**3 + (-a*e + b*d)*(5*a**2*e**2 - 13*a*b*d*e + 11*b**2*d**2)*atanh(sqrt(e)*sqrt(a + b*x)/(sqrt(b)*sqrt(d + e*x)))/(b**(sympy.S(7)/2)*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_830():
    f = (15*d**2 + 20*d*e*x + 8*e**2*x**2)/(sqrt(a + b*x)*sqrt(d + e*x))
    F = 4*e*(a + b*x)**(sympy.S(3)/2)*sqrt(d + e*x)/b**2 + sqrt(a + b*x)*sqrt(d + e*x)*(-10*a*e + 14*b*d)/b**2 + (6*a**2*e**2 - 16*a*b*d*e + 16*b**2*d**2)*atanh(sqrt(e)*sqrt(a + b*x)/(sqrt(b)*sqrt(d + e*x)))/(b**(sympy.S(5)/2)*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_831():
    f = (15*d**2 + 20*d*e*x + 8*e**2*x**2)/(sqrt(a + b*x)*(d + e*x)**(sympy.S(3)/2))
    F = 6*d**2*sqrt(a + b*x)/(sqrt(d + e*x)*(-a*e + b*d)) + 8*sqrt(a + b*x)*sqrt(d + e*x)/b + (-8*a*e + 16*b*d)*atanh(sqrt(e)*sqrt(a + b*x)/(sqrt(b)*sqrt(d + e*x)))/(b**(sympy.S(3)/2)*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_832():
    f = (15*d**2 + 20*d*e*x + 8*e**2*x**2)/(sqrt(a + b*x)*(d + e*x)**(sympy.S(5)/2))
    F = 2*d**2*sqrt(a + b*x)/((d + e*x)**(sympy.S(3)/2)*(-a*e + b*d)) + 4*d*sqrt(a + b*x)*(-2*a*e + 3*b*d)/(sqrt(d + e*x)*(-a*e + b*d)**2) + 16*atanh(sqrt(e)*sqrt(a + b*x)/(sqrt(b)*sqrt(d + e*x)))/(sqrt(b)*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_833():
    f = (15*d**2 + 20*d*e*x + 8*e**2*x**2)/(sqrt(a + b*x)*(d + e*x)**(sympy.S(7)/2))
    F = 6*d**2*sqrt(a + b*x)/((d + e*x)**(sympy.S(5)/2)*(-5*a*e + 5*b*d)) + 8*d*sqrt(a + b*x)*(-5*a*e + 8*b*d)/(15*(d + e*x)**(sympy.S(3)/2)*(-a*e + b*d)**2) + sqrt(a + b*x)*(240*a**2*e**2 - 560*a*b*d*e + 368*b**2*d**2)/(15*sqrt(d + e*x)*(-a*e + b*d)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_834():
    f = (15*d**2 + 20*d*e*x + 8*e**2*x**2)/(sqrt(a + b*x)*(d + e*x)**(sympy.S(9)/2))
    F = 32*b*sqrt(a + b*x)*(35*a**2*e**2 - 84*a*b*d*e + 58*b**2*d**2)/(105*sqrt(d + e*x)*(-a*e + b*d)**4) + 6*d**2*sqrt(a + b*x)/((d + e*x)**(sympy.S(7)/2)*(-7*a*e + 7*b*d)) + 4*d*sqrt(a + b*x)*(-14*a*e + 23*b*d)/(35*(d + e*x)**(sympy.S(5)/2)*(-a*e + b*d)**2) + sqrt(a + b*x)*(560*a**2*e**2 - 1344*a*b*d*e + 928*b**2*d**2)/(105*(d + e*x)**(sympy.S(3)/2)*(-a*e + b*d)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_835():
    f = (d + e*x)**(sympy.S(3)/2)/(sqrt(f + g*x)*(a + b*x + c*x**2))
    F = 2*e**(sympy.S(3)/2)*atanh(sqrt(g)*sqrt(d + e*x)/(sqrt(e)*sqrt(f + g*x)))/(c*sqrt(g)) - (2*e*(-b*e + 2*c*d) - 2*(b**2*e**2 + 2*c**2*d**2 - 2*c*e*(a*e + b*d))/sqrt(-4*a*c + b**2))*atanh(sqrt(d + e*x)*sqrt(2*c*f - g*(b + sqrt(-4*a*c + b**2)))/(sqrt(f + g*x)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))))/(c*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))*sqrt(2*c*f - g*(b + sqrt(-4*a*c + b**2)))) - (2*e*(-b*e + 2*c*d) + 2*(b**2*e**2 + 2*c**2*d**2 - 2*c*e*(a*e + b*d))/sqrt(-4*a*c + b**2))*atanh(sqrt(d + e*x)*sqrt(2*c*f - g*(b - sqrt(-4*a*c + b**2)))/(sqrt(f + g*x)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))))/(c*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))*sqrt(2*c*f - g*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_836():
    f = sqrt(d + e*x)/(sqrt(f + g*x)*(a + b*x + c*x**2))
    F = -2*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))*atanh(sqrt(d + e*x)*sqrt(2*c*f - g*(b - sqrt(-4*a*c + b**2)))/(sqrt(f + g*x)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))))/(sqrt(-4*a*c + b**2)*sqrt(2*c*f - g*(b - sqrt(-4*a*c + b**2)))) + 2*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))*atanh(sqrt(d + e*x)*sqrt(2*c*f - g*(b + sqrt(-4*a*c + b**2)))/(sqrt(f + g*x)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))))/(sqrt(-4*a*c + b**2)*sqrt(2*c*f - g*(b + sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_837():
    f = 1/(sqrt(d + e*x)*sqrt(f + g*x)*(a + b*x + c*x**2))
    F = 4*c*atanh(sqrt(d + e*x)*sqrt(2*c*f - g*(b + sqrt(-4*a*c + b**2)))/(sqrt(f + g*x)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))))/(sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))*sqrt(2*c*f - g*(b + sqrt(-4*a*c + b**2)))) - 4*c*atanh(sqrt(d + e*x)*sqrt(2*c*f - g*(b - sqrt(-4*a*c + b**2)))/(sqrt(f + g*x)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))))/(sqrt(-4*a*c + b**2)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))*sqrt(2*c*f - g*(b - sqrt(-4*a*c + b**2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_838():
    f = 1/((d + e*x)**(sympy.S(3)/2)*sqrt(f + g*x)*(a + b*x + c*x**2))
    F = 8*c**2*atanh(sqrt(d + e*x)*sqrt(2*c*f - g*(b + sqrt(-4*a*c + b**2)))/(sqrt(f + g*x)*sqrt(2*c*d - e*(b + sqrt(-4*a*c + b**2)))))/(sqrt(-4*a*c + b**2)*(2*c*d - e*(b + sqrt(-4*a*c + b**2)))**(sympy.S(3)/2)*sqrt(2*c*f - g*(b + sqrt(-4*a*c + b**2)))) - 8*c**2*atanh(sqrt(d + e*x)*sqrt(2*c*f - g*(b - sqrt(-4*a*c + b**2)))/(sqrt(f + g*x)*sqrt(2*c*d - e*(b - sqrt(-4*a*c + b**2)))))/(sqrt(-4*a*c + b**2)*(2*c*d - e*(b - sqrt(-4*a*c + b**2)))**(sympy.S(3)/2)*sqrt(2*c*f - g*(b - sqrt(-4*a*c + b**2)))) - 4*c*e*sqrt(f + g*x)/(sqrt(d + e*x)*sqrt(-4*a*c + b**2)*(2*c*d - e*(b + sqrt(-4*a*c + b**2)))*(-d*g + e*f)) + 4*c*e*sqrt(f + g*x)/(sqrt(d + e*x)*sqrt(-4*a*c + b**2)*(2*c*d - e*(b - sqrt(-4*a*c + b**2)))*(-d*g + e*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_839():
    f = (f + g*x)**3*sqrt(a + b*x + c*x**2)/(d + e*x)
    F = (-d*g + e*f)**3*sqrt(a*e**2 - b*d*e + c*d**2)*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/e**5 + g**3*(d + e*x)*(a + b*x + c*x**2)**(sympy.S(3)/2)/(4*c*e**2) + g**2*(a + b*x + c*x**2)**(sympy.S(3)/2)*(-5*b*e*g - 14*c*d*g + 24*c*e*f)/(24*c**2*e**2) + sqrt(a + b*x + c*x**2)*(5*b**3*e**3*g**3 + 16*b*c**2*e*g*(d**2*g**2 - 3*d*e*f*g + 3*e**2*f**2) - 4*b*c*e**2*g**2*(a*e*g - 2*b*d*g + 6*b*e*f) + 64*c**3*(-d*g + e*f)**3 + 2*c*e*g*x*(5*b**2*e**2*g**2 + 16*c**2*(d**2*g**2 - 3*d*e*f*g + 3*e**2*f**2) - 4*c*e*g*(a*e*g - 2*b*d*g + 6*b*e*f)))/(64*c**3*e**4) - (4*c*e*(-b*e + 2*c*d)*(5*b**2*d*e*g**3 + 16*c**2*e**2*f**3 - 4*c*d*g**2*(a*e*g - 2*b*d*g + 6*b*e*f)) - g*(-b**2*e**2 + 8*c**2*d**2 - 4*c*e*(-a*e + b*d))*(5*b**2*e**2*g**2 + 16*c**2*(d**2*g**2 - 3*d*e*f*g + 3*e**2*f**2) - 4*c*e*g*(a*e*g - 2*b*d*g + 6*b*e*f)))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(128*c**(sympy.S(7)/2)*e**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_840():
    f = (f + g*x)**2*sqrt(a + b*x + c*x**2)/(d + e*x)
    F = (-d*g + e*f)**2*sqrt(a*e**2 - b*d*e + c*d**2)*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/e**4 + g**2*(a + b*x + c*x**2)**(sympy.S(3)/2)/(3*c*e) - sqrt(a + b*x + c*x**2)*(b**2*e**2*g**2 - 2*b*c*e*g*(-d*g + 2*e*f) - 8*c**2*(-d*g + e*f)**2 - 2*c*e*g*x*(-b*e*g - 2*c*d*g + 4*c*e*f))/(8*c**2*e**3) + (-4*c*e*(-b*e + 2*c*d)*(-b*d*g**2 + 2*c*e*f**2) + g*(-b**2*e**2 + 8*c**2*d**2 - 4*c*e*(-a*e + b*d))*(-b*e*g - 2*c*d*g + 4*c*e*f))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(5)/2)*e**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_841():
    f = (f + g*x)*sqrt(a + b*x + c*x**2)/(d + e*x)
    F = (-d*g + e*f)*sqrt(a*e**2 - b*d*e + c*d**2)*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/e**3 + sqrt(a + b*x + c*x**2)*(b*e*g - 4*c*d*g + 4*c*e*f + 2*c*e*g*x)/(4*c*e**2) - (b**2*e**2*g + 8*c**2*d*(-d*g + e*f) - 4*c*e*(a*e*g - b*d*g + b*e*f))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*c**(sympy.S(3)/2)*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_842():
    f = sqrt(a + b*x + c*x**2)/(d + e*x)
    F = sqrt(a + b*x + c*x**2)/e + sqrt(a*e**2 - b*d*e + c*d**2)*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/e**2 - (-b*e + 2*c*d)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*sqrt(c)*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_843():
    f = sqrt(a + b*x + c*x**2)/((d + e*x)*(f + g*x))
    F = sqrt(c)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(e*g) - sqrt(a*g**2 - b*f*g + c*f**2)*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/(g*(-d*g + e*f)) + sqrt(a*e**2 - b*d*e + c*d**2)*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/(e*(-d*g + e*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_844():
    f = sqrt(a + b*x + c*x**2)/((d + e*x)*(f + g*x)**2)
    F = -sqrt(c)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(g*(-d*g + e*f)) - e*sqrt(a*g**2 - b*f*g + c*f**2)*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/(g*(-d*g + e*f)**2) + sqrt(a*e**2 - b*d*e + c*d**2)*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/(-d*g + e*f)**2 + sqrt(a + b*x + c*x**2)/((f + g*x)*(-d*g + e*f)) + (-b*g + 2*c*f)*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/(2*g*(-d*g + e*f)*sqrt(a*g**2 - b*f*g + c*f**2)) + e*(-b*g + 2*c*f)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*sqrt(c)*g*(-d*g + e*f)**2) - (-b*e + 2*c*d)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*sqrt(c)*(-d*g + e*f)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_845():
    f = sqrt(a + b*x + c*x**2)/((d + e*x)*(f + g*x)**3)
    F = -sqrt(c)*e*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(g*(-d*g + e*f)**2) - e**2*sqrt(a*g**2 - b*f*g + c*f**2)*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/(g*(-d*g + e*f)**3) + e*sqrt(a*e**2 - b*d*e + c*d**2)*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/(-d*g + e*f)**3 + e*sqrt(a + b*x + c*x**2)/((f + g*x)*(-d*g + e*f)**2) + e*(-b*g + 2*c*f)*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/(2*g*(-d*g + e*f)**2*sqrt(a*g**2 - b*f*g + c*f**2)) + g*(-4*a*c + b**2)*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/((-8*d*g + 8*e*f)*(a*g**2 - b*f*g + c*f**2)**(sympy.S(3)/2)) - g*sqrt(a + b*x + c*x**2)*(-2*a*g + b*f + x*(-b*g + 2*c*f))/((f + g*x)**2*(-4*d*g + 4*e*f)*(a*g**2 - b*f*g + c*f**2)) + e**2*(-b*g + 2*c*f)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*sqrt(c)*g*(-d*g + e*f)**3) - e*(-b*e + 2*c*d)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*sqrt(c)*(-d*g + e*f)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_846():
    f = sqrt(a + b*x + c*x**2)/((d + e*x)*(f + g*x)**4)
    F = -sqrt(c)*e**2*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(g*(-d*g + e*f)**3) - e**3*sqrt(a*g**2 - b*f*g + c*f**2)*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/(g*(-d*g + e*f)**4) + e**2*sqrt(a*e**2 - b*d*e + c*d**2)*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/(-d*g + e*f)**4 + e**2*sqrt(a + b*x + c*x**2)/((f + g*x)*(-d*g + e*f)**3) + e**2*(-b*g + 2*c*f)*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/(2*g*(-d*g + e*f)**3*sqrt(a*g**2 - b*f*g + c*f**2)) + e*g*(-4*a*c + b**2)*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/(8*(-d*g + e*f)**2*(a*g**2 - b*f*g + c*f**2)**(sympy.S(3)/2)) - e*g*sqrt(a + b*x + c*x**2)*(-2*a*g + b*f + x*(-b*g + 2*c*f))/(4*(f + g*x)**2*(-d*g + e*f)**2*(a*g**2 - b*f*g + c*f**2)) + g**2*(a + b*x + c*x**2)**(sympy.S(3)/2)/((f + g*x)**3*(-3*d*g + 3*e*f)*(a*g**2 - b*f*g + c*f**2)) + g*(-4*a*c + b**2)*(-b*g + 2*c*f)*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/((-16*d*g + 16*e*f)*(a*g**2 - b*f*g + c*f**2)**(sympy.S(5)/2)) - g*(-b*g + 2*c*f)*sqrt(a + b*x + c*x**2)*(-2*a*g + b*f + x*(-b*g + 2*c*f))/((f + g*x)**2*(-8*d*g + 8*e*f)*(a*g**2 - b*f*g + c*f**2)**2) + e**3*(-b*g + 2*c*f)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*sqrt(c)*g*(-d*g + e*f)**4) - e**2*(-b*e + 2*c*d)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*sqrt(c)*(-d*g + e*f)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_847():
    f = (f + g*x)**3*(a + b*x + c*x**2)**(sympy.S(3)/2)/(d + e*x)
    F = (-d*g + e*f)**3*(a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/e**7 + g**3*(d + e*x)*(a + b*x + c*x**2)**(sympy.S(5)/2)/(6*c*e**2) + g**2*(a + b*x + c*x**2)**(sympy.S(5)/2)*(-7*b*e*g - 22*c*d*g + 36*c*e*f)/(60*c**2*e**2) + (a + b*x + c*x**2)**(sympy.S(3)/2)*(7*b**3*e**3*g**3 + 24*b*c**2*e*g*(d**2*g**2 - 3*d*e*f*g + 3*e**2*f**2) - 4*b*c*e**2*g**2*(a*e*g - 3*b*d*g + 9*b*e*f) + 64*c**3*(-d*g + e*f)**3 + 2*c*e*g*x*(7*b**2*e**2*g**2 + 24*c**2*(d**2*g**2 - 3*d*e*f*g + 3*e**2*f**2) - 4*c*e*g*(a*e*g - 3*b*d*g + 9*b*e*f)))/(192*c**3*e**4) - sqrt(a + b*x + c*x**2)*(21*b**5*e**5*g**3 - 12*b**3*c*e**4*g**2*(8*a*e*g - 3*b*d*g + 9*b*e*f) - 96*b*c**3*e**2*(3*a*e*g*(d**2*g**2 - 3*d*e*f*g + 3*e**2*f**2) + 2*b*(-d*g + e*f)**3) + 24*b*c**2*e**3*g*(2*a**2*e**2*g**2 + 6*a*b*e*g*(-d*g + 3*e*f) + 3*b**2*(d**2*g**2 - 3*d*e*f*g + 3*e**2*f**2)) - 1536*c**5*d**2*(-d*g + e*f)**3 + 384*c**4*e*(-4*a*e + 5*b*d)*(-d*g + e*f)**3 + 2*c*e*x*(8*c*e*(-b*e + 2*c*d)*(7*b**2*d*e*g**3 + 24*c**2*e**2*f**3 - 4*c*d*g**2*(a*e*g - 3*b*d*g + 9*b*e*f)) - g*(7*b**2*e**2*g**2 + 24*c**2*(d**2*g**2 - 3*d*e*f*g + 3*e**2*f**2) - 4*c*e*g*(a*e*g - 3*b*d*g + 9*b*e*f))*(12*a*c*e**2 - 3*b**2*e**2 - 8*b*c*d*e + 16*c**2*d**2)))/(1536*c**4*e**6) + (4*c*e*(-b*e + 2*c*d)*(8*c*e*(-2*a*e + b*d)*(7*b**2*d*e*g**3 + 24*c**2*e**2*f**3 - 4*c*d*g**2*(a*e*g - 3*b*d*g + 9*b*e*f)) - d*g*(-4*a*c*e - 3*b**2*e + 8*b*c*d)*(7*b**2*e**2*g**2 + 24*c**2*(d**2*g**2 - 3*d*e*f*g + 3*e**2*f**2) - 4*c*e*g*(a*e*g - 3*b*d*g + 9*b*e*f))) - (8*c*e*(-b*e + 2*c*d)*(7*b**2*d*e*g**3 + 24*c**2*e**2*f**3 - 4*c*d*g**2*(a*e*g - 3*b*d*g + 9*b*e*f)) - g*(7*b**2*e**2*g**2 + 24*c**2*(d**2*g**2 - 3*d*e*f*g + 3*e**2*f**2) - 4*c*e*g*(a*e*g - 3*b*d*g + 9*b*e*f))*(12*a*c*e**2 - 3*b**2*e**2 - 8*b*c*d*e + 16*c**2*d**2))*(-b**2*e**2 + 8*c**2*d**2 - 4*c*e*(-a*e + b*d)))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(3072*c**(sympy.S(9)/2)*e**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_848():
    f = (f + g*x)**2*(a + b*x + c*x**2)**(sympy.S(3)/2)/(d + e*x)
    F = (-d*g + e*f)**2*(a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/e**6 + g**2*(a + b*x + c*x**2)**(sympy.S(5)/2)/(5*c*e) - (a + b*x + c*x**2)**(sympy.S(3)/2)*(3*b**2*e**2*g**2 - 6*b*c*e*g*(-d*g + 2*e*f) - 16*c**2*(-d*g + e*f)**2 - 6*c*e*g*x*(-b*e*g - 2*c*d*g + 4*c*e*f))/(48*c**2*e**3) + sqrt(a + b*x + c*x**2)*(3*b**4*e**4*g**2 - 6*b**2*c*e**3*g*(2*a*e*g - b*d*g + 2*b*e*f) + 8*b*c**2*e**2*(3*a*e*g*(-d*g + 2*e*f) + 2*b*(-d*g + e*f)**2) + 128*c**4*d**2*(-d*g + e*f)**2 - 32*c**3*e*(-4*a*e + 5*b*d)*(-d*g + e*f)**2 + 2*c*e*x*(-8*c*e*(-b*e + 2*c*d)*(-b*d*g**2 + 2*c*e*f**2) + g*(-3*b**2*e**2 + 16*c**2*d**2 - 4*c*e*(-3*a*e + 2*b*d))*(-b*e*g - 2*c*d*g + 4*c*e*f)))/(128*c**3*e**5) - (3*b**5*e**5*g**2 - 6*b**3*c*e**4*g*(4*a*e*g - b*d*g + 2*b*e*f) + 16*b*c**2*e**3*(3*a**2*e**2*g**2 + 3*a*b*e*g*(-d*g + 2*e*f) + b**2*(-d*g + e*f)**2) + 256*c**5*d**3*(-d*g + e*f)**2 - 384*c**4*d*e*(-a*e + b*d)*(-d*g + e*f)**2 + 96*c**3*e**2*(-a**2*e**2*g*(-d*g + 2*e*f) - 2*a*b*e*(-d*g + e*f)**2 + b**2*d*(-d*g + e*f)**2))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(256*c**(sympy.S(7)/2)*e**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_849():
    f = (f + g*x)*(a + b*x + c*x**2)**(sympy.S(3)/2)/(d + e*x)
    F = (-d*g + e*f)*(a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/e**5 + (a + b*x + c*x**2)**(sympy.S(3)/2)*(3*b*e*g - 8*c*d*g + 8*c*e*f + 6*c*e*g*x)/(24*c*e**2) - sqrt(a + b*x + c*x**2)*(3*b**3*e**3*g - 4*b*c*e**2*(3*a*e*g - 2*b*d*g + 2*b*e*f) - 64*c**3*d**2*(-d*g + e*f) + 16*c**2*e*(-4*a*e + 5*b*d)*(-d*g + e*f) + 2*c*e*x*(3*b**2*e**2*g + 16*c**2*d*(-d*g + e*f) - 4*c*e*(3*a*e*g - 2*b*d*g + 2*b*e*f)))/(64*c**2*e**4) + (3*b**4*e**4*g - 8*b**2*c*e**3*(3*a*e*g - b*d*g + b*e*f) - 128*c**4*d**3*(-d*g + e*f) + 192*c**3*d*e*(-a*e + b*d)*(-d*g + e*f) + 48*c**2*e**2*(a**2*e**2*g + 2*a*b*e*(-d*g + e*f) - b**2*d*(-d*g + e*f)))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(128*c**(sympy.S(5)/2)*e**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_850():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)/(d + e*x)
    F = (a + b*x + c*x**2)**(sympy.S(3)/2)/(3*e) + (a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/e**4 + sqrt(a + b*x + c*x**2)*(b**2*e**2 + 8*c**2*d**2 - 2*c*e*x*(-b*e + 2*c*d) - 2*c*e*(-4*a*e + 5*b*d))/(8*c*e**3) - (-b*e + 2*c*d)*(-b**2*e**2 + 8*c**2*d**2 - 4*c*e*(-3*a*e + 2*b*d))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(3)/2)*e**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_851():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)/((d + e*x)*(f + g*x))
    F = -(a*g**2 - b*f*g + c*f**2)**(sympy.S(3)/2)*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/(g**3*(-d*g + e*f)) - sqrt(a + b*x + c*x**2)*(4*c*e*f**2 - 2*c*g*x*(-d*g + e*f) - g*(-4*a*e*g - b*d*g + 5*b*e*f))/(4*e*g**2*(-d*g + e*f)) + sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2)/(e**2*(-d*g + e*f)) + (a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/(e**3*(-d*g + e*f)) + (b*g**2*(-4*a*e*g + b*d*g + 3*b*e*f) + 8*c**2*e*f**3 - 4*c*g*(-a*g*(-d*g + 3*e*f) + 3*b*e*f**2))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*sqrt(c)*e*g**3*(-d*g + e*f)) - (-b*e + 2*c*d)*(a*e**2 - b*d*e + c*d**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*sqrt(c)*e**3*(-d*g + e*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_852():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)/((d + e*x)*(f + g*x)**2)
    F = -e*(a*g**2 - b*f*g + c*f**2)**(sympy.S(3)/2)*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/(g**3*(-d*g + e*f)**2) + (a + b*x + c*x**2)**(sympy.S(3)/2)/((f + g*x)*(-d*g + e*f)) + sqrt(a + b*x + c*x**2)*(-9*b*g + 12*c*f - 6*c*g*x)/(4*g**2*(-d*g + e*f)) + (-3*b*g + 6*c*f)*sqrt(a*g**2 - b*f*g + c*f**2)*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/(2*g**3*(-d*g + e*f)) + (a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/(e**2*(-d*g + e*f)**2) - e*sqrt(a + b*x + c*x**2)*(b**2*g**2 + 8*c**2*f**2 - 2*c*g*x*(-b*g + 2*c*f) - 2*c*g*(-4*a*g + 5*b*f))/(8*c*g**2*(-d*g + e*f)**2) + sqrt(a + b*x + c*x**2)*(b**2*e**2 + 8*c**2*d**2 - 2*c*e*x*(-b*e + 2*c*d) - 2*c*e*(-4*a*e + 5*b*d))/(8*c*e*(-d*g + e*f)**2) - (3*b**2*g**2 + 24*c**2*f**2 - 12*c*g*(-a*g + 2*b*f))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*sqrt(c)*g**3*(-d*g + e*f)) + e*(-b*g + 2*c*f)*(-b**2*g**2 + 8*c**2*f**2 - 4*c*g*(-3*a*g + 2*b*f))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(3)/2)*g**3*(-d*g + e*f)**2) - (-b*e + 2*c*d)*(-b**2*e**2 + 8*c**2*d**2 - 4*c*e*(-3*a*e + 2*b*d))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(3)/2)*e**2*(-d*g + e*f)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_853():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)/((d + e*x)*(f + g*x)**3)
    F = 3*sqrt(c)*(-b*g + 2*c*f)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*g**3*(-d*g + e*f)) - e**2*(a*g**2 - b*f*g + c*f**2)**(sympy.S(3)/2)*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/(g**3*(-d*g + e*f)**3) + e*(a + b*x + c*x**2)**(sympy.S(3)/2)/((f + g*x)*(-d*g + e*f)**2) + 3*e*sqrt(a + b*x + c*x**2)*(-3*b*g + 4*c*f - 2*c*g*x)/(4*g**2*(-d*g + e*f)**2) + 3*e*(-b*g + 2*c*f)*sqrt(a*g**2 - b*f*g + c*f**2)*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/(2*g**3*(-d*g + e*f)**2) + (a + b*x + c*x**2)**(sympy.S(3)/2)/((f + g*x)**2*(-2*d*g + 2*e*f)) - sqrt(a + b*x + c*x**2)*(-3*b*g + 12*c*f + 6*c*g*x)/(4*g**2*(f + g*x)*(-d*g + e*f)) - (3*b**2*g**2 + 24*c**2*f**2 - 12*c*g*(-a*g + 2*b*f))*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/(8*g**3*(-d*g + e*f)*sqrt(a*g**2 - b*f*g + c*f**2)) + (a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/(e*(-d*g + e*f)**3) - e**2*sqrt(a + b*x + c*x**2)*(b**2*g**2 + 8*c**2*f**2 - 2*c*g*x*(-b*g + 2*c*f) - 2*c*g*(-4*a*g + 5*b*f))/(8*c*g**2*(-d*g + e*f)**3) + sqrt(a + b*x + c*x**2)*(b**2*e**2 + 8*c**2*d**2 - 2*c*e*x*(-b*e + 2*c*d) - 2*c*e*(-4*a*e + 5*b*d))/(8*c*(-d*g + e*f)**3) - 3*e*(b**2*g**2 + 8*c**2*f**2 - 4*c*g*(-a*g + 2*b*f))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*sqrt(c)*g**3*(-d*g + e*f)**2) + e**2*(-b*g + 2*c*f)*(-b**2*g**2 + 8*c**2*f**2 - 4*c*g*(-3*a*g + 2*b*f))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(3)/2)*g**3*(-d*g + e*f)**3) - (-b*e + 2*c*d)*(-b**2*e**2 + 8*c**2*d**2 - 4*c*e*(-3*a*e + 2*b*d))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(3)/2)*e*(-d*g + e*f)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_854():
    f = (a + b*x + c*x**2)**(sympy.S(5)/2)/((d + e*x)*(f + g*x))
    F = -(a*g**2 - b*f*g + c*f**2)**(sympy.S(5)/2)*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/(g**5*(-d*g + e*f)) - (a + b*x + c*x**2)**(sympy.S(3)/2)*(8*c*e*f**2 - 6*c*g*x*(-d*g + e*f) - g*(-8*a*e*g - 3*b*d*g + 11*b*e*f))/(24*e*g**2*(-d*g + e*f)) + (a + b*x + c*x**2)**(sympy.S(3)/2)*(a*e**2 - b*d*e + c*d**2)/(3*e**2*(-d*g + e*f)) + (a*e**2 - b*d*e + c*d**2)**(sympy.S(5)/2)*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/(e**5*(-d*g + e*f)) - sqrt(a + b*x + c*x**2)*(-b**2*g**3*(-8*a*e*g + 3*b*d*g + 5*b*e*f) + 64*c**3*e*f**4 - 16*c**2*e*f**2*g*(-8*a*g + 9*b*f) + 4*c*g**2*(16*a**2*e*g**2 - 3*a*b*g*(-d*g + 13*e*f) + 22*b**2*e*f**2) - 2*c*g*x*(b*g**2*(-8*a*e*g + 3*b*d*g + 5*b*e*f) + 16*c**2*e*f**3 - 4*c*g*(-a*g*(-3*d*g + 7*e*f) + 6*b*e*f**2)))/(64*c*e*g**4*(-d*g + e*f)) + sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2)*(b**2*e**2 + 8*c**2*d**2 - 2*c*e*x*(-b*e + 2*c*d) - 2*c*e*(-4*a*e + 5*b*d))/(8*c*e**4*(-d*g + e*f)) + (-b**3*g**4*(-8*a*e*g + 3*b*d*g + 5*b*e*f) - 8*b*c*g**3*(12*a**2*e*g**2 - 3*a*b*g*(d*g + 5*e*f) + 5*b**2*e*f**2) + 128*c**4*e*f**5 - 320*c**3*e*f**3*g*(-a*g + b*f) + 48*c**2*g**2*(a**2*g**2*(-d*g + 5*e*f) - 10*a*b*e*f**2*g + 5*b**2*e*f**3))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(128*c**(sympy.S(3)/2)*e*g**5*(-d*g + e*f)) - (-b*e + 2*c*d)*(a*e**2 - b*d*e + c*d**2)*(-b**2*e**2 + 8*c**2*d**2 - 4*c*e*(-3*a*e + 2*b*d))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(3)/2)*e**5*(-d*g + e*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_855():
    f = (f + g*x)**4/((d + e*x)*sqrt(a + b*x + c*x**2))
    F = (-d*g + e*f)**4*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/(e**4*sqrt(a*e**2 - b*d*e + c*d**2)) + g**4*(d + e*x)**2*sqrt(a + b*x + c*x**2)/(3*c*e**3) + g**3*(d + e*x)*sqrt(a + b*x + c*x**2)*(-5*b*e*g - 14*c*d*g + 24*c*e*f)/(12*c**2*e**3) + g**2*sqrt(a + b*x + c*x**2)*(15*b**2*e**2*g**2 + 4*c**2*(11*d**2*g**2 - 36*d*e*f*g + 36*e**2*f**2) - 4*c*e*g*(4*a*e*g - 7*b*d*g + 18*b*e*f))/(24*c**3*e**3) - g*(5*b**3*e**3*g**3 - 6*b*c*e**2*g**2*(2*a*e*g - b*d*g + 4*b*e*f) - 16*c**3*(-d**3*g**3 + 4*d**2*e*f*g**2 - 6*d*e**2*f**2*g + 4*e**3*f**3) + 8*c**2*e*g*(a*e*g*(-d*g + 4*e*f) + b*(d**2*g**2 - 4*d*e*f*g + 6*e**2*f**2)))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(7)/2)*e**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_856():
    f = (f + g*x)**3/((d + e*x)*sqrt(a + b*x + c*x**2))
    F = (-d*g + e*f)**3*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/(e**3*sqrt(a*e**2 - b*d*e + c*d**2)) + g**3*(d + e*x)*sqrt(a + b*x + c*x**2)/(2*c*e**2) + 3*g**2*sqrt(a + b*x + c*x**2)*(-b*e*g - 2*c*d*g + 4*c*e*f)/(4*c**2*e**2) + g*(3*b**2*e**2*g**2 + 8*c**2*(d**2*g**2 - 3*d*e*f*g + 3*e**2*f**2) - 4*c*e*g*(a*e*g - b*d*g + 3*b*e*f))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*c**(sympy.S(5)/2)*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_857():
    f = (f + g*x)**2/((d + e*x)*sqrt(a + b*x + c*x**2))
    F = (-d*g + e*f)**2*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/(e**2*sqrt(a*e**2 - b*d*e + c*d**2)) + g**2*sqrt(a + b*x + c*x**2)/(c*e) + g*(-b*e*g - 2*c*d*g + 4*c*e*f)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*c**(sympy.S(3)/2)*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_858():
    f = (f + g*x)/((d + e*x)*sqrt(a + b*x + c*x**2))
    F = (-d*g + e*f)*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/(e*sqrt(a*e**2 - b*d*e + c*d**2)) + g*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(sqrt(c)*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_859():
    f = 1/((d + e*x)*sqrt(a + b*x + c*x**2))
    F = atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/sqrt(a*e**2 - b*d*e + c*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_860():
    f = 1/((d + e*x)*(f + g*x)*sqrt(a + b*x + c*x**2))
    F = e*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/((-d*g + e*f)*sqrt(a*e**2 - b*d*e + c*d**2)) - g*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/((-d*g + e*f)*sqrt(a*g**2 - b*f*g + c*f**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_861():
    f = 1/((d + e*x)*(f + g*x)**2*sqrt(a + b*x + c*x**2))
    F = e**2*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/((-d*g + e*f)**2*sqrt(a*e**2 - b*d*e + c*d**2)) - e*g*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/((-d*g + e*f)**2*sqrt(a*g**2 - b*f*g + c*f**2)) + g**2*sqrt(a + b*x + c*x**2)/((f + g*x)*(-d*g + e*f)*(a*g**2 - b*f*g + c*f**2)) - g*(-b*g + 2*c*f)*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/((-2*d*g + 2*e*f)*(a*g**2 - b*f*g + c*f**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_862():
    f = 1/((d + e*x)*(f + g*x)**3*sqrt(a + b*x + c*x**2))
    F = e**3*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/((-d*g + e*f)**3*sqrt(a*e**2 - b*d*e + c*d**2)) - e**2*g*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/((-d*g + e*f)**3*sqrt(a*g**2 - b*f*g + c*f**2)) + e*g**2*sqrt(a + b*x + c*x**2)/((f + g*x)*(-d*g + e*f)**2*(a*g**2 - b*f*g + c*f**2)) - e*g*(-b*g + 2*c*f)*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/(2*(-d*g + e*f)**2*(a*g**2 - b*f*g + c*f**2)**(sympy.S(3)/2)) + 3*g**2*(-b*g + 2*c*f)*sqrt(a + b*x + c*x**2)/((f + g*x)*(-4*d*g + 4*e*f)*(a*g**2 - b*f*g + c*f**2)**2) + g**2*sqrt(a + b*x + c*x**2)/((f + g*x)**2*(-2*d*g + 2*e*f)*(a*g**2 - b*f*g + c*f**2)) - g*(3*b**2*g**2 + 8*c**2*f**2 - 4*c*g*(a*g + 2*b*f))*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/((-8*d*g + 8*e*f)*(a*g**2 - b*f*g + c*f**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_863():
    f = (f + g*x)**4/((d + e*x)*(a + b*x + c*x**2)**(sympy.S(3)/2))
    F = (-d*g + e*f)**4*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/(e**2*(a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)) - (2*a*b**3*d*g**4 + 4*a*c*(a**2*e*g**4 - 2*a*c*f*g**2*(-2*d*g + 3*e*f) + c**2*f**3*(-4*d*g + e*f)) - 2*b**2*(a**2*e*g**4 + 4*a*c*d*f*g**3 + c**2*e*f**4) + 2*b*c*(a**2*g**3*(-3*d*g + 4*e*f) + 2*a*c*f**2*g*(3*d*g + 2*e*f) + c**2*d*f**4) + 2*x*(b**3*g**4*(-a*e + b*d) - b*c*g**3*(-3*a**2*e*g - 4*a*b*(-d*g + e*f) + 4*b**2*d*f) + 2*c**4*d*f**4 + c**3*f**2*(4*a*g*(-3*d*g + 2*e*f) - b*f*(4*d*g + e*f)) + 2*c**2*g**2*(-a**2*g*(-d*g + 4*e*f) - 3*a*b*f*(-2*d*g + e*f) + 3*b**2*d*f**2)))/(c**2*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2)) + g**4*sqrt(a + b*x + c*x**2)/(c**2*e) + g**3*(-3*b*e*g - 2*c*d*g + 8*c*e*f)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*c**(sympy.S(5)/2)*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_864():
    f = (f + g*x)**3/((d + e*x)*(a + b*x + c*x**2)**(sympy.S(3)/2))
    F = (-d*g + e*f)**3*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/(e*(a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)) + (-4*a*c*(-a*g**2*(-d*g + 3*e*f) + c*f**2*(-3*d*g + e*f)) + 2*b**2*(a*d*g**3 + c*e*f**3) - 2*b*(a**2*e*g**3 + 3*a*c*f*g*(d*g + e*f) + c**2*d*f**3) - 2*x*(-b**2*g**3*(-a*e + b*d) + 2*c**3*d*f**3 + c**2*f*(6*a*g*(-d*g + e*f) - b*f*(3*d*g + e*f)) + c*g**2*(-2*a**2*e*g + 3*a*b*d*g - 3*a*b*e*f + 3*b**2*d*f)))/(c*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2)) + g**3*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(c**(sympy.S(3)/2)*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_865():
    f = (f + g*x)**2/((d + e*x)*(a + b*x + c*x**2)**(sympy.S(3)/2))
    F = (-d*g + e*f)**2*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/(a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2) + (4*a*(a*e*g**2 - c*f*(-2*d*g + e*f)) + 2*b**2*e*f**2 - 2*b*(a*g*(d*g + 2*e*f) + c*d*f**2) - 2*x*(b*g**2*(-a*e + b*d) + 2*c**2*d*f**2 + c*(2*a*g*(-d*g + 2*e*f) - b*f*(2*d*g + e*f))))/((-4*a*c + b**2)*sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_866():
    f = (f + g*x)/((d + e*x)*(a + b*x + c*x**2)**(sympy.S(3)/2))
    F = e*(-d*g + e*f)*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/(a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2) - (2*a*b*e*g - 4*a*c*d*g + 4*a*c*e*f - 2*b**2*e*f + 2*b*c*d*f + 2*c*x*(2*a*e*g - b*(d*g + e*f) + 2*c*d*f))/((-4*a*c + b**2)*sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_867():
    f = 1/((d + e*x)*(a + b*x + c*x**2)**(sympy.S(3)/2))
    F = e**2*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/(a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2) - (4*a*c*e - 2*b**2*e + 2*b*c*d + 2*c*x*(-b*e + 2*c*d))/((-4*a*c + b**2)*sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_868():
    f = 1/((d + e*x)*(f + g*x)*(a + b*x + c*x**2)**(sympy.S(3)/2))
    F = e**3*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/((-d*g + e*f)*(a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)) - 2*e*(2*a*c*e - b**2*e + b*c*d + c*x*(-b*e + 2*c*d))/((-4*a*c + b**2)*(-d*g + e*f)*sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2)) - g**3*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/((-d*g + e*f)*(a*g**2 - b*f*g + c*f**2)**(sympy.S(3)/2)) + 2*g*(2*a*c*g - b**2*g + b*c*f + c*x*(-b*g + 2*c*f))/((-4*a*c + b**2)*(-d*g + e*f)*sqrt(a + b*x + c*x**2)*(a*g**2 - b*f*g + c*f**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_869():
    f = 1/((d + e*x)*(f + g*x)**2*(a + b*x + c*x**2)**(sympy.S(3)/2))
    F = e**4*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/((-d*g + e*f)**2*(a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)) - 2*e**2*(2*a*c*e - b**2*e + b*c*d + c*x*(-b*e + 2*c*d))/((-4*a*c + b**2)*(-d*g + e*f)**2*sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2)) - e*g**3*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/((-d*g + e*f)**2*(a*g**2 - b*f*g + c*f**2)**(sympy.S(3)/2)) + 2*e*g*(2*a*c*g - b**2*g + b*c*f + c*x*(-b*g + 2*c*f))/((-4*a*c + b**2)*(-d*g + e*f)**2*sqrt(a + b*x + c*x**2)*(a*g**2 - b*f*g + c*f**2)) - 3*g**3*(-b*g + 2*c*f)*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/((-2*d*g + 2*e*f)*(a*g**2 - b*f*g + c*f**2)**(sympy.S(5)/2)) + g**2*sqrt(a + b*x + c*x**2)*(3*b**2*g**2 + 4*c**2*f**2 - 4*c*g*(2*a*g + b*f))/((f + g*x)*(-4*a*c + b**2)*(-d*g + e*f)*(a*g**2 - b*f*g + c*f**2)**2) + 2*g*(2*a*c*g - b**2*g + b*c*f + c*x*(-b*g + 2*c*f))/((f + g*x)*(-4*a*c + b**2)*(-d*g + e*f)*sqrt(a + b*x + c*x**2)*(a*g**2 - b*f*g + c*f**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_870():
    f = 1/((d + e*x)*(f + g*x)**3*(a + b*x + c*x**2)**(sympy.S(3)/2))
    F = e**5*atanh((-2*a*e + b*d + x*(-b*e + 2*c*d))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*e**2 - b*d*e + c*d**2)))/((-d*g + e*f)**3*(a*e**2 - b*d*e + c*d**2)**(sympy.S(3)/2)) - 2*e**3*(2*a*c*e - b**2*e + b*c*d + c*x*(-b*e + 2*c*d))/((-4*a*c + b**2)*(-d*g + e*f)**3*sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2)) - e**2*g**3*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/((-d*g + e*f)**3*(a*g**2 - b*f*g + c*f**2)**(sympy.S(3)/2)) + 2*e**2*g*(2*a*c*g - b**2*g + b*c*f + c*x*(-b*g + 2*c*f))/((-4*a*c + b**2)*(-d*g + e*f)**3*sqrt(a + b*x + c*x**2)*(a*g**2 - b*f*g + c*f**2)) - 3*e*g**3*(-b*g + 2*c*f)*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/(2*(-d*g + e*f)**2*(a*g**2 - b*f*g + c*f**2)**(sympy.S(5)/2)) + e*g**2*sqrt(a + b*x + c*x**2)*(3*b**2*g**2 + 4*c**2*f**2 - 4*c*g*(2*a*g + b*f))/((f + g*x)*(-4*a*c + b**2)*(-d*g + e*f)**2*(a*g**2 - b*f*g + c*f**2)**2) + 2*e*g*(2*a*c*g - b**2*g + b*c*f + c*x*(-b*g + 2*c*f))/((f + g*x)*(-4*a*c + b**2)*(-d*g + e*f)**2*sqrt(a + b*x + c*x**2)*(a*g**2 - b*f*g + c*f**2)) - 3*g**3*(5*b**2*g**2 + 16*c**2*f**2 - 4*c*g*(a*g + 4*b*f))*atanh((-2*a*g + b*f + x*(-b*g + 2*c*f))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*g**2 - b*f*g + c*f**2)))/((-8*d*g + 8*e*f)*(a*g**2 - b*f*g + c*f**2)**(sympy.S(7)/2)) + g**2*(-b*g + 2*c*f)*sqrt(a + b*x + c*x**2)*(15*b**2*g**2 + 8*c**2*f**2 - 4*c*g*(13*a*g + 2*b*f))/((f + g*x)*(-16*a*c + 4*b**2)*(-d*g + e*f)*(a*g**2 - b*f*g + c*f**2)**3) + g**2*sqrt(a + b*x + c*x**2)*(5*b**2*g**2 + 8*c**2*f**2 - 4*c*g*(3*a*g + 2*b*f))/((f + g*x)**2*(-8*a*c + 2*b**2)*(-d*g + e*f)*(a*g**2 - b*f*g + c*f**2)**2) + 2*g*(2*a*c*g - b**2*g + b*c*f + c*x*(-b*g + 2*c*f))/((f + g*x)**2*(-4*a*c + b**2)*(-d*g + e*f)*sqrt(a + b*x + c*x**2)*(a*g**2 - b*f*g + c*f**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_871():
    f = (d + e*x)**3*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)
    F = 2*(d + e*x)**4*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)/(11*e) + 2*e**2*(f + g*x)**(sympy.S(7)/2)*sqrt(a + b*x + c*x**2)*(b*e*g - 3*c*d*g + c*e*f)/(99*c*g**4) - 2*e*(f + g*x)**(sympy.S(5)/2)*sqrt(a + b*x + c*x**2)*(8*b**2*e**2*g**2 + c**2*(81*d**2*g**2 - 96*d*e*f*g + 29*e**2*f**2) + c*e*g*(-18*a*e*g - 33*b*d*g + 19*b*e*f))/(693*c**2*g**4) + (f + g*x)**(sympy.S(3)/2)*sqrt(a + b*x + c*x**2)*(96*b**3*e**3*g**3 + 2*b*c*e**2*g**2*(-157*a*e*g - 198*b*d*g + 67*b*e*f) + 2*c**3*(-567*d**3*g**3 + 1107*d**2*e*f*g**2 - 843*d*e**2*f**2*g + 233*e**3*f**3) - 2*c**2*e*g*(2*a*e*g*(-231*d*g + 74*e*f) - 3*b*(99*d**2*g**2 - 88*d*e*f*g + 24*e**2*f**2)))/(3465*c**3*g**4) + sqrt(f + g*x)*sqrt(a + b*x + c*x**2)*(-128*b**4*e**4*g**4 - 8*b**2*c*e**3*g**3*(-69*a*e*g - 66*b*d*g + 7*b*e*f) - 2*c**4*(315*d**4*g**4 - 798*d**3*e*f*g**3 + 1098*d**2*e**2*f**2*g**2 - 732*d*e**3*f**3*g + 187*e**4*f**4) + 2*c**3*e*g*(6*a*e*g*(165*d**2*g**2 - 33*d*e*f*g + 2*e**2*f**2) + b*(231*d**3*g**3 - 99*d**2*e*f*g**2 + 8*e**3*f**3)) - 6*c**2*e**2*g**2*(50*a**2*e**2*g**2 - a*b*e*g*(-297*d*g + 29*e*f) + 3*b**2*(44*d**2*g**2 - 11*d*e*f*g + e**2*f**2)))/(3465*c**4*e*g**4) + 2*sqrt(2)*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(a*g**2 - b*f*g + c*f**2)*(64*b**4*e**3*g**4 + 4*b**2*c*e**2*g**3*(-69*a*e*g - 66*b*d*g + 7*b*e*f) - 2*c**4*f*(-231*d**3*g**3 + 396*d**2*e*f*g**2 - 264*d*e**2*f**2*g + 64*e**3*f**3) - c**3*g*(6*a*e*g*(165*d**2*g**2 - 33*d*e*f*g + 2*e**2*f**2) + b*(231*d**3*g**3 - 99*d**2*e*f*g**2 + 8*e**3*f**3)) + 3*c**2*e*g**2*(50*a**2*e**2*g**2 - a*b*e*g*(-297*d*g + 29*e*f) + 3*b**2*(44*d**2*g**2 - 11*d*e*f*g + e**2*f**2)))*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(3465*c**5*g**5*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)) + sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(f + g*x)*sqrt(-4*a*c + b**2)*(128*b**5*e**3*g**5 - 8*b**3*c*e**2*g**4*(87*a*e*g + 66*b*d*g + 7*b*e*f) + b*c**2*e*g**3*(771*a**2*e**2*g**2 + 6*a*b*e*g*(396*d*g + 43*e*f) - b**2*(-792*d**2*g**2 - 264*d*e*f*g + 37*e**2*f**2)) + 2*c**5*f**2*(-231*d**3*g**3 + 396*d**2*e*f*g**2 - 264*d*e**2*f**2*g + 64*e**3*f**3) - c**4*g*(-18*a*g*(77*d**3*g**3 + 88*d**2*e*f*g**2 - 33*d*e**2*f**2*g + 6*e**3*f**3) + b*f*(-462*d**3*g**3 + 495*d**2*e*f*g**2 - 264*d*e**2*f**2*g + 56*e**3*f**3)) - c**3*g**2*(6*a**2*e**2*g**2*(231*d*g + 26*e*f) - 9*a*b*e*g*(-319*d**2*g**2 - 110*d*e*f*g + 15*e**2*f**2) + b**2*(462*d**3*g**3 + 495*d**2*e*f*g**2 - 198*d*e**2*f**2*g + 37*e**3*f**3)))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(3465*c**5*g**5*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_872():
    f = (d + e*x)**2*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)
    F = 2*(d + e*x)**3*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)/(9*e) + 2*e*(f + g*x)**(sympy.S(5)/2)*sqrt(a + b*x + c*x**2)*(b*e*g - 3*c*d*g + c*e*f)/(63*c*g**3) - (f + g*x)**(sympy.S(3)/2)*sqrt(a + b*x + c*x**2)*(12*b**2*e**2*g**2 + 4*c**2*(21*d**2*g**2 - 24*d*e*f*g + 8*e**2*f**2) + 4*c*e*g*(-7*a*e*g - 9*b*d*g + 4*b*e*f))/(315*c**2*g**3) + sqrt(f + g*x)*sqrt(a + b*x + c*x**2)*(16*b**3*e**3*g**3 + 6*b*c*e**2*g**2*(-9*a*e*g - 8*b*d*g + b*e*f) + 2*c**3*(-35*d**3*g**3 + 63*d**2*e*f*g**2 - 57*d*e**2*f**2*g + 19*e**3*f**3) - 6*c**2*e*g**2*(2*a*e*(-10*d*g + e*f) + b*d*(-7*d*g + 2*e*f)))/(315*c**3*e*g**3) - 2*sqrt(2)*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(a*g**2 - b*f*g + c*f**2)*(8*b**3*e**2*g**3 + 3*b*c*e*g**2*(-9*a*e*g - 8*b*d*g + b*e*f) - 2*c**3*f*(21*d**2*g**2 - 24*d*e*f*g + 8*e**2*f**2) - 3*c**2*g**2*(2*a*e*(-10*d*g + e*f) + b*d*(-7*d*g + 2*e*f)))*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(315*c**4*g**4*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)) - 2*sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(f + g*x)*sqrt(-4*a*c + b**2)*(8*b**4*e**2*g**4 - 4*b**2*c*e*g**3*(9*a*e*g + 6*b*d*g + b*e*f) + c**4*f**2*(21*d**2*g**2 - 24*d*e*f*g + 8*e**2*f**2) + c**3*g*(3*a*g*(-21*d**2*g**2 - 16*d*e*f*g + 3*e**2*f**2) - b*f*(21*d**2*g**2 - 15*d*e*f*g + 4*e**2*f**2)) + 3*c**2*g**2*(7*a**2*e**2*g**2 + a*b*e*g*(29*d*g + 5*e*f) - b**2*(-7*d**2*g**2 - 5*d*e*f*g + e**2*f**2)))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(315*c**4*g**4*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_873():
    f = (d + e*x)*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)
    F = 2*e*sqrt(f + g*x)*(a + b*x + c*x**2)**(sympy.S(3)/2)/(7*c) - 2*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)*(4*b**2*e*g**2 + c**2*f*(-7*d*g + 4*e*f) - 3*c*g*x*(-4*b*e*g + 7*c*d*g + c*e*f) - c*g*(-5*a*e*g + 7*b*d*g + 2*b*e*f))/(105*c**2*g**2) + 2*sqrt(2)*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(a*g**2 - b*f*g + c*f**2)*(4*b**2*e*g**2 - 2*c**2*f*(-7*d*g + 4*e*f) + c*g*(-10*a*e*g - 7*b*d*g + b*e*f))*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(105*c**3*g**3*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)) + sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(f + g*x)*sqrt(-4*a*c + b**2)*(-5*c*g*(-b*g + 2*c*f)*(7*c*d*f - e*(a*g + 3*b*f)) + (-2*b**2*g**2 + 8*c**2*f**2 - 3*c*g*(-2*a*g + b*f))*(-4*b*e*g + 7*c*d*g + c*e*f))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(105*c**3*g**3*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_874():
    f = sqrt(f + g*x)*sqrt(a + b*x + c*x**2)
    F = 2*(f + g*x)**(sympy.S(3)/2)*sqrt(a + b*x + c*x**2)/(5*g) + sqrt(f + g*x)*(2*b*g - 4*c*f)*sqrt(a + b*x + c*x**2)/(15*c*g) + 2*sqrt(2)*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(-b*g + 2*c*f)*(a*g**2 - b*f*g + c*f**2)*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(15*c**2*g**2*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)) - 2*sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(f + g*x)*sqrt(-4*a*c + b**2)*(b**2*g**2 + c**2*f**2 - c*g*(3*a*g + b*f))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(15*c**2*g**2*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_875():
    f = (d + e*x)**3*sqrt(a + b*x + c*x**2)/sqrt(f + g*x)
    F = 2*(d + e*x)**3*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)/(9*g) - 2*e**2*(f + g*x)**(sympy.S(5)/2)*sqrt(a + b*x + c*x**2)*(-b*e*g - 6*c*d*g + 8*c*e*f)/(63*c*g**4) - 2*e*(f + g*x)**(sympy.S(3)/2)*sqrt(a + b*x + c*x**2)*(6*b**2*e**2*g**2 - 2*c**2*(42*d**2*g**2 - 111*d*e*f*g + 64*e**2*f**2) + c*e*g*(-14*a*e*g - 27*b*d*g + 17*b*e*f))/(315*c**2*g**4) + sqrt(f + g*x)*sqrt(a + b*x + c*x**2)*(16*b**3*e**3*g**3 + 6*b*c*e**2*g**2*(-9*a*e*g - 12*b*d*g + 5*b*e*f) - 2*c**3*(-70*d**3*g**3 + 336*d**2*e*f*g**2 - 408*d*e**2*f**2*g + 152*e**3*f**3) - 6*c**2*e*g*(6*a*e*g*(-5*d*g + 2*e*f) - b*(21*d**2*g**2 - 24*d*e*f*g + 8*e**2*f**2)))/(315*c**3*g**4) - 2*sqrt(2)*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(a*g**2 - b*f*g + c*f**2)*(8*b**3*e**3*g**3 + 3*b*c*e**2*g**2*(-9*a*e*g - 12*b*d*g + 5*b*e*f) + 2*c**3*(-105*d**3*g**3 + 252*d**2*e*f*g**2 - 216*d*e**2*f**2*g + 64*e**3*f**3) - 3*c**2*e*g*(6*a*e*g*(-5*d*g + 2*e*f) - b*(21*d**2*g**2 - 24*d*e*f*g + 8*e**2*f**2)))*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(315*c**4*g**5*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)) - sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(f + g*x)*sqrt(-4*a*c + b**2)*(16*b**4*e**3*g**4 + 8*b**2*c*e**2*g**3*(-9*a*e*g - 9*b*d*g + 2*b*e*f) - 2*c**4*f*(-105*d**3*g**3 + 252*d**2*e*f*g**2 - 216*d*e**2*f**2*g + 64*e**3*f**3) - c**3*g*(6*a*e*g*(63*d**2*g**2 - 39*d*e*f*g + 10*e**2*f**2) - b*(-105*d**3*g**3 + 189*d**2*e*f*g**2 - 144*d*e**2*f**2*g + 40*e**3*f**3)) + 3*c**2*e*g**2*(14*a**2*e**2*g**2 - a*b*e*g*(-87*d*g + 19*e*f) + b**2*(42*d**2*g**2 - 27*d*e*f*g + 7*e**2*f**2)))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(315*c**4*g**5*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_876():
    f = (d + e*x)**2*sqrt(a + b*x + c*x**2)/sqrt(f + g*x)
    F = 2*(d + e*x)**2*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)/(7*g) - 2*e*(f + g*x)**(sympy.S(3)/2)*sqrt(a + b*x + c*x**2)*(-b*e*g - 4*c*d*g + 6*c*e*f)/(35*c*g**3) + sqrt(f + g*x)*sqrt(a + b*x + c*x**2)*(-8*b**2*e**2*g**2 + 4*c**2*(10*d**2*g**2 - 34*d*e*f*g + 21*e**2*f**2) - 4*c*e*g*(-5*a*e*g - 7*b*d*g + 4*b*e*f))/(105*c**2*g**3) + 4*sqrt(2)*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(a*g**2 - b*f*g + c*f**2)*(2*b**2*e**2*g**2 + c**2*(35*d**2*g**2 - 56*d*e*f*g + 24*e**2*f**2) + c*e*g*(-5*a*e*g - 7*b*d*g + 4*b*e*f))*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(105*c**3*g**4*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)) + sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(f + g*x)*sqrt(-4*a*c + b**2)*(8*b**3*e**2*g**3 + b*c*e*g**2*(-29*a*e*g - 28*b*d*g + 9*b*e*f) - 2*c**3*f*(35*d**2*g**2 - 56*d*e*f*g + 24*e**2*f**2) - c**2*g*(2*a*e*g*(-42*d*g + 13*e*f) - b*(35*d**2*g**2 - 42*d*e*f*g + 16*e**2*f**2)))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(105*c**3*g**4*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_877():
    f = (d + e*x)*sqrt(a + b*x + c*x**2)/sqrt(f + g*x)
    F = -2*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)*(-b*e*g - 5*c*d*g + 4*c*e*f - 3*c*e*g*x)/(15*c*g**2) - 2*sqrt(2)*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(a*g**2 - b*f*g + c*f**2)*(b*e*g - 10*c*d*g + 8*c*e*f)*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(15*c**2*g**3*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)) - sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(f + g*x)*sqrt(-4*a*c + b**2)*(2*b**2*e*g**2 - 2*c**2*f*(-5*d*g + 4*e*f) + c*g*(-6*a*e*g - 5*b*d*g + 3*b*e*f))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(15*c**2*g**3*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_878():
    f = sqrt(a + b*x + c*x**2)/sqrt(f + g*x)
    F = 2*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)/(3*g) + 4*sqrt(2)*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(a*g**2 - b*f*g + c*f**2)*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(3*c*g**2*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)) - sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(f + g*x)*sqrt(-4*a*c + b**2)*(-b*g + 2*c*f)*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(3*c*g**2*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_879():
    f = sqrt(a + b*x + c*x**2)/((d + e*x)*sqrt(f + g*x))
    F = ((sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1))))) * ((Symbol('e') * Symbol('g') * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * ((Symbol('c') * Symbol('e') * Symbol('f')) + (Symbol('c') * Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('b') * Symbol('e') * Symbol('g')))) * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1))))) * ((Symbol('c') * (Symbol('e'))**(Integer(2)) * Symbol('g') * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((Symbol('e') * ((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * (Symbol('b') * Symbol('g'))) + (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')))) * ((Integer(2) * Symbol('c') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g'))))))**(Integer(-1))), sympy.asin(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x)))) * (sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))))**(Integer(-1)))), ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('c')) * (Symbol('e'))**(Integer(2)) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_880():
    f = (d + e*x)**3*sqrt(f + g*x)/sqrt(a + b*x + c*x**2)
    F = 2*e*(d + e*x)**2*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)/(7*c) + 2*e**2*(f + g*x)**(sympy.S(3)/2)*sqrt(a + b*x + c*x**2)*(-6*b*e*g + 11*c*d*g + c*e*f)/(35*c**2*g**2) + 2*e*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)*(24*b**2*e**2*g**2 - c**2*(-90*d**2*g**2 + 12*d*e*f*g + 7*e**2*f**2) + c*e*g*(-25*a*e*g - 84*b*d*g + 13*b*e*f))/(105*c**3*g**2) - 2*sqrt(2)*e*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(a*g**2 - b*f*g + c*f**2)*(24*b**2*e**2*g**2 + c**2*(105*d**2*g**2 - 42*d*e*f*g + 8*e**2*f**2) + c*e*g*(-25*a*e*g - 84*b*d*g + 13*b*e*f))*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(105*c**4*g**3*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)) - sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(f + g*x)*sqrt(-4*a*c + b**2)*(48*b**3*e**3*g**3 - 8*b*c*e**2*g**2*(13*a*e*g + 21*b*d*g + 2*b*e*f) - c**3*(105*d**3*g**3 + 105*d**2*e*f*g**2 - 42*d*e**2*f**2*g + 8*e**3*f**3) + c**2*e*g*(a*e*g*(189*d*g + 19*e*f) - b*(-210*d**2*g**2 - 63*d*e*f*g + 9*e**2*f**2)))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(105*c**4*g**3*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_881():
    f = (d + e*x)**2*sqrt(f + g*x)/sqrt(a + b*x + c*x**2)
    F = 2*e*(d + e*x)*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)/(5*c) + 2*e*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)*(-4*b*e*g + 7*c*d*g + c*e*f)/(15*c**2*g) + 4*sqrt(2)*e*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(a*g**2 - b*f*g + c*f**2)*(2*b*e*g - 5*c*d*g + c*e*f)*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(15*c**3*g**2*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)) + sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(f + g*x)*sqrt(-4*a*c + b**2)*(8*b**2*e**2*g**2 - c**2*(-15*d**2*g**2 - 10*d*e*f*g + 2*e**2*f**2) - c*e*g*(9*a*e*g + 20*b*d*g + 3*b*e*f))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(15*c**3*g**2*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_882():
    f = (d + e*x)*sqrt(f + g*x)/sqrt(a + b*x + c*x**2)
    F = 2*e*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)/(3*c) - 2*sqrt(2)*e*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(a*g**2 - b*f*g + c*f**2)*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(3*c**2*g*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)) + sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(f + g*x)*sqrt(-4*a*c + b**2)*(-2*b*e*g + 3*c*d*g + c*e*f)*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(3*c**2*g*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_883():
    f = sqrt(f + g*x)/sqrt(a + b*x + c*x**2)
    F = sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(f + g*x)*sqrt(-4*a*c + b**2)*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(c*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_884():
    f = sqrt(f + g*x)/((d + e*x)*sqrt(a + b*x + c*x**2))
    F = ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g') * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1))))) * ((Symbol('c') * Symbol('e') * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((Symbol('e') * ((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * (Symbol('b') * Symbol('g'))) + (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')))) * ((Integer(2) * Symbol('c') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g'))))))**(Integer(-1))), sympy.asin(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x)))) * (sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))))**(Integer(-1)))), ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('c')) * Symbol('e') * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_885():
    f = sqrt(f + g*x)/((d + e*x)**2*sqrt(a + b*x + c*x**2))
    F = (Integer(-1) * ((Symbol('e') * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1))))) * ((sympy.sqrt(Integer(2)) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('f') * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1))))) * ((((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('d') * Symbol('g') * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1))))) * ((Symbol('e') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + ((sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))) * (((Symbol('e'))**(Integer(2)) * ((Symbol('b') * Symbol('f')) + (Integer(-1) * (Symbol('a') * Symbol('g'))))) + (Integer(-1) * (Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g'))))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((Symbol('e') * ((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * (Symbol('b') * Symbol('g'))) + (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')))) * ((Integer(2) * Symbol('c') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g'))))))**(Integer(-1))), sympy.asin(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x)))) * (sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))))**(Integer(-1)))), ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * Symbol('e') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_886():
    f = sqrt(f + g*x)/((d + e*x)**3*sqrt(a + b*x + c*x**2))
    F = ((Integer(-1) * (Symbol('e') * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))))) * ((Integer(2) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * ((Symbol('c') * Symbol('d') * ((Integer(6) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(5) * Symbol('d') * Symbol('g'))))) + (Integer(-1) * (Symbol('e') * ((Integer(3) * Symbol('b') * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(2) * Symbol('b') * Symbol('d') * Symbol('g'))) + (Integer(-1) * (Symbol('a') * Symbol('e') * Symbol('g'))))))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(4) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * ((Symbol('c') * Symbol('d') * ((Integer(6) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(5) * Symbol('d') * Symbol('g'))))) + (Integer(-1) * (Symbol('e') * ((Integer(3) * Symbol('b') * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(2) * Symbol('b') * Symbol('d') * Symbol('g'))) + (Integer(-1) * (Symbol('a') * Symbol('e') * Symbol('g'))))))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Integer(2)) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g') * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1))))) * ((sympy.sqrt(Integer(2)) * Symbol('e') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('f') * ((Symbol('c') * Symbol('d') * ((Integer(6) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(5) * Symbol('d') * Symbol('g'))))) + (Integer(-1) * (Symbol('e') * ((Integer(3) * Symbol('b') * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(2) * Symbol('b') * Symbol('d') * Symbol('g'))) + (Integer(-1) * (Symbol('a') * Symbol('e') * Symbol('g'))))))) * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('d') * Symbol('g') * ((Symbol('c') * Symbol('d') * ((Integer(6) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(5) * Symbol('d') * Symbol('g'))))) + (Integer(-1) * (Symbol('e') * ((Integer(3) * Symbol('b') * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(2) * Symbol('b') * Symbol('d') * Symbol('g'))) + (Integer(-1) * (Symbol('a') * Symbol('e') * Symbol('g'))))))) * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * Symbol('e') * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + ((sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))) * ((Symbol('c') * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('c') * Symbol('d') * Symbol('g'))) + (Symbol('b') * Symbol('e') * Symbol('g'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((Symbol('e') * ((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * (Symbol('b') * Symbol('g'))) + (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')))) * ((Integer(2) * Symbol('c') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g'))))))**(Integer(-1))), sympy.asin(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x)))) * (sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))))**(Integer(-1)))), ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * Symbol('e') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))) * ((Symbol('c') * Symbol('d') * ((Integer(6) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(5) * Symbol('d') * Symbol('g'))))) + (Integer(-1) * (Symbol('e') * ((Integer(3) * Symbol('b') * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(2) * Symbol('b') * Symbol('d') * Symbol('g'))) + (Integer(-1) * (Symbol('a') * Symbol('e') * Symbol('g'))))))) * ((Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g'))))) + (Integer(-1) * (Symbol('e') * ((Symbol('b') * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(2) * Symbol('b') * Symbol('d') * Symbol('g'))) + (Symbol('a') * Symbol('e') * Symbol('g')))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((Symbol('e') * ((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * (Symbol('b') * Symbol('g'))) + (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')))) * ((Integer(2) * Symbol('c') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g'))))))**(Integer(-1))), sympy.asin(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x)))) * (sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))))**(Integer(-1)))), ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * Symbol('e') * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_887():
    f = (f + g*x)**(sympy.S(3)/2)/((d + e*x)*sqrt(a + b*x + c*x**2))
    F = ((sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g') * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))))) * ((Symbol('c') * Symbol('e') * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))))) * ((Symbol('c') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((Symbol('e') * ((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * (Symbol('b') * Symbol('g'))) + (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')))) * ((Integer(2) * Symbol('c') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g'))))))**(Integer(-1))), sympy.asin(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x)))) * (sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))))**(Integer(-1)))), ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('c')) * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_888():
    f = (f + g*x)**(sympy.S(5)/2)/((d + e*x)*sqrt(a + b*x + c*x**2))
    F = ((Integer(2) * (Symbol('g'))**(Integer(2)) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(3) * Symbol('c') * Symbol('e')))**(Integer(-1))) + ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g') * ((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * (Symbol('b') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('c'))**(Integer(2)) * Symbol('e') * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + ((sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))))) * ((Symbol('c') * (Symbol('e'))**(Integer(2)) * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g') * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(2)) * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))))) * ((Symbol('c') * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g') * ((Symbol('c') * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('f') * Symbol('g'))) + (Symbol('a') * (Symbol('g'))**(Integer(2)))) * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('c'))**(Integer(2)) * Symbol('e') * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((Symbol('e') * ((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * (Symbol('b') * Symbol('g'))) + (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')))) * ((Integer(2) * Symbol('c') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g'))))))**(Integer(-1))), sympy.asin(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x)))) * (sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))))**(Integer(-1)))), ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('c')) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_889():
    f = (d + e*x)**3/(sqrt(f + g*x)*sqrt(a + b*x + c*x**2))
    F = 2*e**2*(d + e*x)*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)/(5*c*g) - 8*e**2*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)*(b*e*g - 3*c*d*g + c*e*f)/(15*c**2*g**2) + sqrt(2)*e*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(f + g*x)*sqrt(-4*a*c + b**2)*(8*b**2*e**2*g**2 + c**2*(45*d**2*g**2 - 30*d*e*f*g + 8*e**2*f**2) + c*e*g*(-9*a*e*g - 30*b*d*g + 7*b*e*f))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(15*c**3*g**3*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2)) - 2*sqrt(2)*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(4*b*e**3*g**2*(-a*g + b*f) + c**2*(-15*d**3*g**3 + 45*d**2*e*f*g**2 - 30*d*e**2*f**2*g + 8*e**3*f**3) - c*e**2*g*(a*g*(-15*d*g + 7*e*f) - 3*b*f*(-5*d*g + e*f)))*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(15*c**3*g**3*sqrt(f + g*x)*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_890():
    f = (d + e*x)**2/(sqrt(f + g*x)*sqrt(a + b*x + c*x**2))
    F = 2*e**2*sqrt(f + g*x)*sqrt(a + b*x + c*x**2)/(3*c*g) - 2*sqrt(2)*e*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(f + g*x)*sqrt(-4*a*c + b**2)*(b*e*g - 3*c*d*g + c*e*f)*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(3*c**2*g**2*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2)) + 2*sqrt(2)*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(c*(3*d**2*g**2 - 6*d*e*f*g + 2*e**2*f**2) + e**2*g*(-a*g + b*f))*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(3*c**2*g**2*sqrt(f + g*x)*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_891():
    f = (d + e*x)/(sqrt(f + g*x)*sqrt(a + b*x + c*x**2))
    F = sqrt(2)*e*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(f + g*x)*sqrt(-4*a*c + b**2)*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(c*g*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2)) - 2*sqrt(2)*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(-d*g + e*f)*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(c*g*sqrt(f + g*x)*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_892():
    f = 1/(sqrt(f + g*x)*sqrt(a + b*x + c*x**2))
    F = 2*sqrt(2)*sqrt(c*(f + g*x)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*g*sqrt(-4*a*c + b**2)/(2*c*f - g*(b + sqrt(-4*a*c + b**2))))/(c*sqrt(f + g*x)*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_893():
    f = 1/((d + e*x)*sqrt(f + g*x)*sqrt(a + b*x + c*x**2))
    F = Integer(-1) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((Symbol('e') * ((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * (Symbol('b') * Symbol('g'))) + (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')))) * ((Integer(2) * Symbol('c') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g'))))))**(Integer(-1))), sympy.asin(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x)))) * (sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))))**(Integer(-1)))), ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('c')) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_894():
    f = 1/((d + e*x)**2*sqrt(f + g*x)*sqrt(a + b*x + c*x**2))
    F = (Integer(-1) * (((Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e') * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1))))) * ((sympy.sqrt(Integer(2)) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e') * Symbol('f') * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1))))) * ((((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('d') * Symbol('g') * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(-2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1))))) * ((((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))) * ((Symbol('c') * Symbol('d') * ((Integer(2) * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(3) * Symbol('d') * Symbol('g'))))) + (Integer(-1) * (Symbol('e') * ((Symbol('b') * Symbol('e') * Symbol('f')) + (Integer(-1) * (Integer(2) * Symbol('b') * Symbol('d') * Symbol('g'))) + (Symbol('a') * Symbol('e') * Symbol('g')))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((Symbol('e') * ((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * (Symbol('b') * Symbol('g'))) + (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')))) * ((Integer(2) * Symbol('c') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g'))))))**(Integer(-1))), sympy.asin(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x)))) * (sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))))**(Integer(-1)))), ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_895():
    f = 1/((d + e*x)*(f + g*x)**(sympy.S(3)/2)*sqrt(a + b*x + c*x**2))
    F = ((Integer(2) * (Symbol('g'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * ((Symbol('c') * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('f') * Symbol('g'))) + (Symbol('a') * (Symbol('g'))**(Integer(2)))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x)))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g') * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))))) * ((((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * ((Symbol('c') * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('f') * Symbol('g'))) + (Symbol('a') * (Symbol('g'))**(Integer(2)))) * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * Symbol('e') * sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((Symbol('e') * ((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * (Symbol('b') * Symbol('g'))) + (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')))) * ((Integer(2) * Symbol('c') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g'))))))**(Integer(-1))), sympy.asin(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x)))) * (sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))))**(Integer(-1)))), ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('c')) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_896():
    f = 1/((d + e*x)*(f + g*x)**(sympy.S(5)/2)*sqrt(a + b*x + c*x**2))
    F = ((Integer(2) * (Symbol('g'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(3) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * ((Symbol('c') * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('f') * Symbol('g'))) + (Symbol('a') * (Symbol('g'))**(Integer(2)))) * ((Symbol('f') + (Symbol('g') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(4) * (Symbol('g'))**(Integer(2)) * ((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * (Symbol('b') * Symbol('g')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(3) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * (((Symbol('c') * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('f') * Symbol('g'))) + (Symbol('a') * (Symbol('g'))**(Integer(2)))))**(Integer(2)) * sympy.sqrt((Symbol('f') + (Symbol('g') * x)))))**(Integer(-1))) + ((Integer(2) * Symbol('e') * (Symbol('g'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * (((((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(2)) * ((Symbol('c') * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('f') * Symbol('g'))) + (Symbol('a') * (Symbol('g'))**(Integer(2)))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g') * ((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * (Symbol('b') * Symbol('g')))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))))) * ((Integer(3) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * (((Symbol('c') * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('f') * Symbol('g'))) + (Symbol('a') * (Symbol('g'))**(Integer(2)))))**(Integer(2)) * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e') * Symbol('g') * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_e(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))))) * (((((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(2)) * ((Symbol('c') * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('f') * Symbol('g'))) + (Symbol('a') * (Symbol('g'))**(Integer(2)))) * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g') * sympy.sqrt(((Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))**(Integer(-1))))) * sympy.elliptic_f(sympy.asin((sympy.sqrt(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))))) * ((Integer(3) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * ((Symbol('c') * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('f') * Symbol('g'))) + (Symbol('a') * (Symbol('g'))**(Integer(2)))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * (Symbol('e'))**(Integer(2)) * sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f') + (Symbol('g') * x))) * (((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((Symbol('e') * ((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * (Symbol('b') * Symbol('g'))) + (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('g')))) * ((Integer(2) * Symbol('c') * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g'))))))**(Integer(-1))), sympy.asin(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Symbol('f') + (Symbol('g') * x)))) * (sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g'))))))**(Integer(-1)))), ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(-1) * ((Integer(2) * Symbol('c') * Symbol('f')) * (Symbol('g'))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('c')) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_897():
    f = sqrt(d + e*x)/(sqrt(f + g*x)*sqrt(a + b*x + c*x**2))
    F = (sympy.sqrt(Integer(2)) * sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g'))))) * sympy.sqrt((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(2) * Symbol('c') * x))) * sympy.sqrt(((((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x))) * ((((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) * sympy.sqrt(((((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * ((Integer(2) * Symbol('a')) + ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * x))) * ((((Symbol('b') * Symbol('f')) + (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('f')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('g')))) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) * (Symbol('d') + (Symbol('e') * x)) * sympy.Function('EllipticPi')(((Symbol('e') * ((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g'))))) * ((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))) * Symbol('g')))**(Integer(-1))), sympy.asin(((sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e'))))) * sympy.sqrt((Symbol('f') + (Symbol('g') * x)))) * ((sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g'))))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))), ((((Symbol('b') * Symbol('d')) + (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('d')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('e')))) * ((Integer(2) * Symbol('c') * Symbol('f')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g'))))) * ((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))) * ((Symbol('b') * Symbol('f')) + (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('f')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('g'))))))**(Integer(-1))))) * ((sympy.sqrt(((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e'))))) * Symbol('g') * sympy.sqrt((((Integer(2) * Symbol('a') * Symbol('c')) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1))) + (Symbol('c') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_898():
    f = 1/(sqrt(d + e*x)*sqrt(f + g*x)*sqrt(a + b*x + c*x**2))
    F = -sqrt((1 - (f + g*x)*(2*a*e*g - b*(d*g + e*f) + 2*c*d*f)/((d + e*x)*(a*g**2 - b*f*g + c*f**2)) + (f + g*x)**2*(a*e**2 - b*d*e + c*d**2)/((d + e*x)**2*(c*f**2 - g*(-a*g + b*f))))/(1 + (f + g*x)*sqrt(a*e**2 - b*d*e + c*d**2)/((d + e*x)*sqrt(c*f**2 - g*(-a*g + b*f))))**2)*sqrt((-d*g + e*f)**2*(a + b*x + c*x**2)/((d + e*x)**2*(a*g**2 - b*f*g + c*f**2)))*(1 + (f + g*x)*sqrt(a*e**2 - b*d*e + c*d**2)/((d + e*x)*sqrt(c*f**2 - g*(-a*g + b*f))))*(d + e*x)*(c*f**2 - g*(-a*g + b*f))**(sympy.S(1)/4)*elliptic_f(2*atan(sqrt(f + g*x)*(a*e**2 - b*d*e + c*d**2)**(sympy.S(1)/4)/(sqrt(d + e*x)*(a*g**2 - b*f*g + c*f**2)**(sympy.S(1)/4))), sympy.S.Half + (2*a*e*g - b*(d*g + e*f) + 2*c*d*f)/(4*sqrt(c*d**2 - e*(-a*e + b*d))*sqrt(c*f**2 - g*(-a*g + b*f))))/((-d*g + e*f)*sqrt(1 - (f + g*x)*(2*a*e*g - b*(d*g + e*f) + 2*c*d*f)/((d + e*x)*(a*g**2 - b*f*g + c*f**2)) + (f + g*x)**2*(a*e**2 - b*d*e + c*d**2)/((d + e*x)**2*(c*f**2 - g*(-a*g + b*f))))*sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_899():
    f = (d + e*x)**m*(f + g*x)**2*(a + b*x + c*x**2)
    F = c*g**2*(d + e*x)**(m + 5)/(e**5*(m + 5)) + g*(d + e*x)**(m + 4)*(b*e*g - 4*c*d*g + 2*c*e*f)/(e**5*(m + 4)) + (d + e*x)**(m + 1)*(-d*g + e*f)**2*(a*e**2 - b*d*e + c*d**2)/(e**5*(m + 1)) - (d + e*x)**(m + 2)*(-d*g + e*f)*(2*c*d*(-2*d*g + e*f) - e*(2*a*e*g - 3*b*d*g + b*e*f))/(e**5*(m + 2)) + (d + e*x)**(m + 3)*(c*(6*d**2*g**2 - 6*d*e*f*g + e**2*f**2) + e*g*(a*e*g - 3*b*d*g + 2*b*e*f))/(e**5*(m + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_900():
    f = (d + e*x)**m*(f + g*x)*(a + b*x + c*x**2)
    F = c*g*(d + e*x)**(m + 4)/(e**4*(m + 4)) + (d + e*x)**(m + 1)*(-d*g + e*f)*(a*e**2 - b*d*e + c*d**2)/(e**4*(m + 1)) - (d + e*x)**(m + 2)*(c*d*(-3*d*g + 2*e*f) - e*(a*e*g - 2*b*d*g + b*e*f))/(e**4*(m + 2)) + (d + e*x)**(m + 3)*(b*e*g - 3*c*d*g + c*e*f)/(e**4*(m + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_901():
    f = (d + e*x)**m*(a + b*x + c*x**2)/(f + g*x)
    F = c*(d + e*x)**(m + 2)/(e**2*g*(m + 2)) + (d + e*x)**(m + 1)*(a*g**2 - b*f*g + c*f**2)*hyper((1, m + 1), (m + 2,), -g*(d + e*x)/(-d*g + e*f))/(g**2*(m + 1)*(-d*g + e*f)) - (d + e*x)**(m + 1)*(-b*e*g + c*d*g + c*e*f)/(e**2*g**2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_902():
    f = (d + e*x)**m*(a + b*x + c*x**2)/(f + g*x)**2
    F = c*(d + e*x)**(m + 1)/(e*g**2*(m + 1)) + (a + f*(-b*g + c*f)/g**2)*(d + e*x)**(m + 1)/((f + g*x)*(-d*g + e*f)) + (d + e*x)**(m + 1)*(c*f*(2*d*g - e*f*(m + 2)) - g*(a*e*g*m + b*(d*g - e*f*(m + 1))))*hyper((1, m + 1), (m + 2,), -g*(d + e*x)/(-d*g + e*f))/(g**2*(m + 1)*(-d*g + e*f)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_903():
    f = (d + e*x)**m*(a + b*x + c*x**2)/(f + g*x)**3
    F = (a + f*(-b*g + c*f)/g**2)*(d + e*x)**(m + 1)/((f + g*x)**2*(-2*d*g + 2*e*f)) + (d + e*x)**(m + 1)*(c*(2*d**2*g**2 - 4*d*e*f*g*(m + 1) + e**2*f**2*(m**2 + 3*m + 2)) - e*g*m*(a*e*g*(1 - m) - b*(2*d*g - e*f*(m + 1))))*hyper((1, m + 1), (m + 2,), -g*(d + e*x)/(-d*g + e*f))/(2*g**2*(m + 1)*(-d*g + e*f)**3) + (d + e*x)**(m + 1)*(c*f*(4*d*g - e*f*(m + 3)) + g*(a*e*g*(1 - m) - b*(2*d*g - e*f*(m + 1))))/(2*g**2*(f + g*x)*(-d*g + e*f)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_904():
    f = (d + e*x)**m*(f + g*x)**2*(a + b*x + c*x**2)**2
    F = c**2*g**2*(d + e*x)**(m + 7)/(e**7*(m + 7)) + 2*c*g*(d + e*x)**(m + 6)*(b*e*g - 3*c*d*g + c*e*f)/(e**7*(m + 6)) + (d + e*x)**(m + 1)*(-d*g + e*f)**2*(a*e**2 - b*d*e + c*d**2)**2/(e**7*(m + 1)) - (d + e*x)**(m + 2)*(-d*g + e*f)*(c*d*(-3*d*g + 2*e*f) - e*(a*e*g - 2*b*d*g + b*e*f))*(2*a*e**2 - 2*b*d*e + 2*c*d**2)/(e**7*(m + 2)) + (d + e*x)**(m + 3)*(c**2*d**2*(15*d**2*g**2 - 20*d*e*f*g + 6*e**2*f**2) + 2*c*e*(a*e*(6*d**2*g**2 - 6*d*e*f*g + e**2*f**2) - b*d*(10*d**2*g**2 - 12*d*e*f*g + 3*e**2*f**2)) + e**2*(a**2*e**2*g**2 + 2*a*b*e*g*(-3*d*g + 2*e*f) + b**2*(6*d**2*g**2 - 6*d*e*f*g + e**2*f**2)))/(e**7*(m + 3)) + (d + e*x)**(m + 4)*(2*b*e**2*g*(a*e*g - 2*b*d*g + b*e*f) - 4*c**2*d*(5*d**2*g**2 - 5*d*e*f*g + e**2*f**2) + 2*c*e*(2*a*e*g*(-2*d*g + e*f) + b*(10*d**2*g**2 - 8*d*e*f*g + e**2*f**2)))/(e**7*(m + 4)) + (d + e*x)**(m + 5)*(b**2*e**2*g**2 + c**2*(15*d**2*g**2 - 10*d*e*f*g + e**2*f**2) + 2*c*e*g*(a*e*g - 5*b*d*g + 2*b*e*f))/(e**7*(m + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_905():
    f = (d + e*x)**m*(f + g*x)*(a + b*x + c*x**2)**2
    F = c**2*g*(d + e*x)**(m + 6)/(e**6*(m + 6)) + c*(d + e*x)**(m + 5)*(2*b*e*g - 5*c*d*g + c*e*f)/(e**6*(m + 5)) + (d + e*x)**(m + 1)*(-d*g + e*f)*(a*e**2 - b*d*e + c*d**2)**2/(e**6*(m + 1)) - (d + e*x)**(m + 2)*(c*d*(-5*d*g + 4*e*f) - e*(a*e*g - 3*b*d*g + 2*b*e*f))*(a*e**2 - b*d*e + c*d**2)/(e**6*(m + 2)) + (d + e*x)**(m + 3)*(b*e**2*(2*a*e*g - 3*b*d*g + b*e*f) + 2*c**2*d**2*(-5*d*g + 3*e*f) + 2*c*e*(a*e*(-3*d*g + e*f) - 3*b*d*(-2*d*g + e*f)))/(e**6*(m + 3)) + (d + e*x)**(m + 4)*(b**2*e**2*g - 2*c**2*d*(-5*d*g + 2*e*f) + 2*c*e*(a*e*g - 4*b*d*g + b*e*f))/(e**6*(m + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_906():
    f = (d + e*x)**m*(a + b*x + c*x**2)**2/(f + g*x)
    F = c**2*(d + e*x)**(m + 4)/(e**4*g*(m + 4)) - c*(d + e*x)**(m + 3)*(-2*b*e*g + 3*c*d*g + c*e*f)/(e**4*g**2*(m + 3)) + (d + e*x)**(m + 1)*(a*g**2 - b*f*g + c*f**2)**2*hyper((1, m + 1), (m + 2,), -g*(d + e*x)/(-d*g + e*f))/(g**4*(m + 1)*(-d*g + e*f)) + (d + e*x)**(m + 2)*(b**2*e**2*g**2 + c**2*(3*d**2*g**2 + 2*d*e*f*g + e**2*f**2) + 2*c*e*g*(a*e*g - b*(2*d*g + e*f)))/(e**4*g**3*(m + 2)) + (d + e*x)**(m + 1)*(c*(d**2*g**2 + e**2*f**2) + e*g*(2*a*e*g - b*(d*g + e*f)))*(b*e*g - c*(d*g + e*f))/(e**4*g**4*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_907():
    f = (d + e*x)**m*(a + b*x + c*x**2)**2/(f + g*x)**2
    F = c**2*(d + e*x)**(m + 3)/(e**3*g**2*(m + 3)) - 2*c*(d + e*x)**(m + 2)*(-b*e*g + c*d*g + c*e*f)/(e**3*g**3*(m + 2)) + (d + e*x)**(m + 1)*(c*f*(4*d*g - e*f*(m + 4)) - g*(a*e*g*m + b*(2*d*g - e*f*(m + 2))))*(a*g**2 - b*f*g + c*f**2)*hyper((1, m + 1), (m + 2,), -g*(d + e*x)/(-d*g + e*f))/(g**4*(m + 1)*(-d*g + e*f)**2) + (d + e*x)**(m + 1)*(a*g**2 - b*f*g + c*f**2)**2/(g**4*(f + g*x)*(-d*g + e*f)) + (d + e*x)**(m + 1)*(b**2*e**2*g**2 + c**2*(d**2*g**2 + 2*d*e*f*g + 3*e**2*f**2) + 2*c*e*g*(a*e*g - b*(d*g + 2*e*f)))/(e**3*g**4*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_908():
    f = (d + e*x)**m*(a + b*x + c*x**2)**2/(f + g*x)**3
    F = c**2*(d + e*x)**(m + 2)/(e**2*g**3*(m + 2)) - c*(d + e*x)**(m + 1)*(-2*b*e*g + c*d*g + 3*c*e*f)/(e**2*g**4*(m + 1)) + (d + e*x)**(m + 1)*(c**2*f**2*(12*d**2*g**2 - 8*d*e*f*g*(m + 3) + e**2*f**2*(m**2 + 7*m + 12)) + 2*c*g*(a*g*(2*d**2*g**2 - 4*d*e*f*g*(m + 1) + e**2*f**2*(m**2 + 3*m + 2)) - b*f*(6*d**2*g**2 - 6*d*e*f*g*(m + 2) + e**2*f**2*(m**2 + 5*m + 6))) - g**2*(a**2*e**2*g**2*m*(1 - m) - 2*a*b*e*g*m*(2*d*g - e*f*(m + 1)) - b**2*(2*d**2*g**2 - 4*d*e*f*g*(m + 1) + e**2*f**2*(m**2 + 3*m + 2))))*hyper((1, m + 1), (m + 2,), -g*(d + e*x)/(-d*g + e*f))/(2*g**4*(m + 1)*(-d*g + e*f)**3) + (d + e*x)**(m + 1)*(c*f*(8*d*g - e*f*(m + 7)) + g*(a*e*g*(1 - m) - b*(4*d*g - e*f*(m + 3))))*(a*g**2 - b*f*g + c*f**2)/(2*g**4*(f + g*x)*(-d*g + e*f)**2) + (d + e*x)**(m + 1)*(a*g**2 - b*f*g + c*f**2)**2/(2*g**4*(f + g*x)**2*(-d*g + e*f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_909():
    f = (3*x + 2)**4*(4*x + 1)**m/(3*x**2 - 5*x + 1)
    F = 27*(4*x + 1)**(m + 3)/(64*m + 192) + 3687*(4*x + 1)**(m + 1)/(64*m + 64) + 207*(4*x + 1)**(m + 2)/(32*m + 64) - (16497 - 4893*sqrt(13))*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(13 - 2*sqrt(13)))/((338 - 52*sqrt(13))*(m + 1)) - (16497 + 4893*sqrt(13))*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(2*sqrt(13) + 13))/((52*sqrt(13) + 338)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_910():
    f = (3*x + 2)**3*(4*x + 1)**m/(3*x**2 - 5*x + 1)
    F = 9*(4*x + 1)**(m + 2)/(16*m + 32) + 123*(4*x + 1)**(m + 1)/(16*m + 16) - (1248 - 405*sqrt(13))*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(13 - 2*sqrt(13)))/((169 - 26*sqrt(13))*(m + 1)) - (1248 + 405*sqrt(13))*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(2*sqrt(13) + 13))/((26*sqrt(13) + 169)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_911():
    f = (3*x + 2)**2*(4*x + 1)**m/(3*x**2 - 5*x + 1)
    F = 3*(4*x + 1)**(m + 1)/(4*m + 4) - (351 - 141*sqrt(13))*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(13 - 2*sqrt(13)))/((338 - 52*sqrt(13))*(m + 1)) - (351 + 141*sqrt(13))*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(2*sqrt(13) + 13))/((52*sqrt(13) + 338)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_912():
    f = (3*x + 2)*(4*x + 1)**m/(3*x**2 - 5*x + 1)
    F = -(39 - 27*sqrt(13))*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(13 - 2*sqrt(13)))/((338 - 52*sqrt(13))*(m + 1)) - (39 + 27*sqrt(13))*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(2*sqrt(13) + 13))/((52*sqrt(13) + 338)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_913():
    f = (4*x + 1)**m/(3*x**2 - 5*x + 1)
    F = 3*sqrt(13)*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(13 - 2*sqrt(13)))/(13*(13 - 2*sqrt(13))*(m + 1)) - 3*sqrt(13)*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(2*sqrt(13) + 13))/(13*(2*sqrt(13) + 13)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_914():
    f = (4*x + 1)**m/((3*x + 2)*(3*x**2 - 5*x + 1))
    F = 3*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), -12*x/5 + sympy.S(-3)/5)/(85*m + 85) + (39 + 27*sqrt(13))*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(13 - 2*sqrt(13)))/((5746 - 884*sqrt(13))*(m + 1)) + (39 - 27*sqrt(13))*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(2*sqrt(13) + 13))/((884*sqrt(13) + 5746)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_915():
    f = (4*x + 1)**m/((3*x + 2)**2*(3*x**2 - 5*x + 1))
    F = 27*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), -12*x/5 + sympy.S(-3)/5)/(1445*m + 1445) + 12*(4*x + 1)**(m + 1)*hyper((2, m + 1), (m + 2,), -12*x/5 + sympy.S(-3)/5)/(425*m + 425) + (351 + 141*sqrt(13))*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(13 - 2*sqrt(13)))/((97682 - 15028*sqrt(13))*(m + 1)) + (351 - 141*sqrt(13))*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(2*sqrt(13) + 13))/((15028*sqrt(13) + 97682)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_916():
    f = (3*x + 2)**4*(4*x + 1)**m/(3*x**2 - 5*x + 1)**2
    F = (844 - 2355*x)*(4*x + 1)**(m + 1)/(117*x**2 - 195*x + 39) + 9*(4*x + 1)**(m + 1)/(4*m + 4) - (4*x + 1)**(m + 1)*(sqrt(13)*(4474*m + 1570*sqrt(13)*m + 297) + 13689)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(2*sqrt(13) + 13))/((338*sqrt(13) + 2197)*(m + 1)) - (4*x + 1)**(m + 1)*(-sqrt(13)*(-1570*sqrt(13)*m + 4474*m + 297) + 13689)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(13 - 2*sqrt(13)))/((2197 - 338*sqrt(13))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_917():
    f = (3*x + 2)**3*(4*x + 1)**m/(3*x**2 - 5*x + 1)**2
    F = (209 - 426*x)*(4*x + 1)**(m + 1)/(117*x**2 - 195*x + 39) - (4*x + 1)**(m + 1)*(sqrt(13)*(-1168*m + 568*sqrt(13)*m + 1701) + 1521)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(13 - 2*sqrt(13)))/((4394 - 676*sqrt(13))*(m + 1)) + (4*x + 1)**(m + 1)*(-7384*m + sqrt(13)*(1701 - 1168*m) - 1521)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(2*sqrt(13) + 13))/((676*sqrt(13) + 4394)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_918():
    f = (3*x + 2)**2*(4*x + 1)**m/(3*x**2 - 5*x + 1)**2
    F = (61 - 87*x)*(4*x + 1)**(m + 1)/(117*x**2 - 195*x + 39) - sqrt(13)*(4*x + 1)**(m + 1)*(-2*m*(23 - 29*sqrt(13)) + 306)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(13 - 2*sqrt(13)))/(169*(13 - 2*sqrt(13))*(m + 1)) + sqrt(13)*(4*x + 1)**(m + 1)*(-2*m*(23 + 29*sqrt(13)) + 306)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(2*sqrt(13) + 13))/(169*(2*sqrt(13) + 13)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_919():
    f = (3*x + 2)*(4*x + 1)**m/(3*x**2 - 5*x + 1)**2
    F = (20 - 21*x)*(4*x + 1)**(m + 1)/(117*x**2 - 195*x + 39) + sqrt(13)*(4*x + 1)**(m + 1)*(m*(10 - 14*sqrt(13)) + 81)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(2*sqrt(13) + 13))/(169*(2*sqrt(13) + 13)*(m + 1)) - sqrt(13)*(4*x + 1)**(m + 1)*(m*(10 + 14*sqrt(13)) + 81)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(13 - 2*sqrt(13)))/(169*(13 - 2*sqrt(13))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_920():
    f = (4*x + 1)**m/(3*x**2 - 5*x + 1)**2
    F = (7 - 6*x)*(4*x + 1)**(m + 1)/(117*x**2 - 195*x + 39) + sqrt(13)*(4*x + 1)**(m + 1)*(2*m*(4 - 2*sqrt(13)) + 18)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(2*sqrt(13) + 13))/(169*(2*sqrt(13) + 13)*(m + 1)) - sqrt(13)*(4*x + 1)**(m + 1)*(2*m*(4 + 2*sqrt(13)) + 18)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(13 - 2*sqrt(13)))/(169*(13 - 2*sqrt(13))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_921():
    f = (4*x + 1)**m/((3*x + 2)*(3*x**2 - 5*x + 1)**2)
    F = (43 - 33*x)*(4*x + 1)**(m + 1)/(1989*x**2 - 3315*x + 663) + 9*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), -12*x/5 + sympy.S(-3)/5)/(1445*m + 1445) + sqrt(13)*(4*x + 1)**(m + 1)*(m*(62 - 22*sqrt(13)) + 81)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(2*sqrt(13) + 13))/(2873*(2*sqrt(13) + 13)*(m + 1)) - sqrt(13)*(4*x + 1)**(m + 1)*(m*(62 + 22*sqrt(13)) + 81)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(13 - 2*sqrt(13)))/(2873*(13 - 2*sqrt(13))*(m + 1)) + (117 + 81*sqrt(13))*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(13 - 2*sqrt(13)))/((97682 - 15028*sqrt(13))*(m + 1)) + (117 - 81*sqrt(13))*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(2*sqrt(13) + 13))/((15028*sqrt(13) + 97682)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_922():
    f = (4*x + 1)**m/((3*x + 2)**2*(3*x**2 - 5*x + 1)**2)
    F = (268 - 195*x)*(4*x + 1)**(m + 1)/(33813*x**2 - 56355*x + 11271) + 162*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), -12*x/5 + sympy.S(-3)/5)/(24565*m + 24565) + 36*(4*x + 1)**(m + 1)*hyper((2, m + 1), (m + 2,), -12*x/5 + sympy.S(-3)/5)/(7225*m + 7225) + sqrt(13)*(4*x + 1)**(m + 1)*(m*(422 - 130*sqrt(13)) + 423)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(2*sqrt(13) + 13))/(48841*(2*sqrt(13) + 13)*(m + 1)) - sqrt(13)*(4*x + 1)**(m + 1)*(m*(422 + 130*sqrt(13)) + 423)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(13 - 2*sqrt(13)))/(48841*(13 - 2*sqrt(13))*(m + 1)) + (1053 + 576*sqrt(13))*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(13 - 2*sqrt(13)))/((830297 - 127738*sqrt(13))*(m + 1)) + (1053 - 576*sqrt(13))*(4*x + 1)**(m + 1)*hyper((1, m + 1), (m + 2,), (12*x + 3)/(2*sqrt(13) + 13))/((127738*sqrt(13) + 830297)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_923():
    f = (d + e*x)**m*(a + b*x + c*x**2)/(e + f*x)**(sympy.S(3)/2)
    F = 2*c*(d + e*x)**(m + 1)*sqrt(e + f*x)/(e*f**2*(2*m + 3)) + (2*a + 2*e*(-b*f + c*e)/f**2)*(d + e*x)**(m + 1)/(sqrt(e + f*x)*(-d*f + e**2)) + (d + e*x)**m*sqrt(e + f*x)*(2*c*(d**2*f**2 + 4*d*e**2*f*(m + 1) - 4*e**4*(m**2 + 3*m + 2)) - 2*e*f*(2*m + 3)*(a*e*f*(2*m + 1) + b*(d*f - 2*e**2*(m + 1))))*hyper((sympy.S.Half, -m), (sympy.S(3)/2,), e*(e + f*x)/(-d*f + e**2))/(e*f**3*(-f*(d + e*x)/(-d*f + e**2))**m*(2*m + 3)*(-d*f + e**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_924():
    f = (d + e*x)**m*(f + g*x)**2*sqrt(a + b*x + c*x**2)
    F = g**2*(d + e*x)**(m + 1)*(a + b*x + c*x**2)**(sympy.S(3)/2)/(c*e*(m + 4)) - g*(d + e*x)**(m + 2)*(b*e*g*(2*m + 5) + 2*c*(3*d*g - 2*e*f*(m + 4)))*sqrt(a + b*x + c*x**2)*appellf1(m + 2, sympy.S(-1)/2, sympy.S(-1)/2, m + 3, 2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))), 2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(2*c*e**3*(m + 2)*(m + 4)*sqrt(-2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))) + 1)*sqrt(-2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))) + 1)) + (d + e*x)**(m + 1)*(c*(3*d**2*g**2 - 2*d*e*f*g*(m + 4) + e**2*f**2*(m + 4)) + e*g**2*(m + 1)*(-a*e + b*d))*sqrt(a + b*x + c*x**2)*appellf1(m + 1, sympy.S(-1)/2, sympy.S(-1)/2, m + 2, 2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))), 2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(c*e**3*(m + 1)*(m + 4)*sqrt(-2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))) + 1)*sqrt(-2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_925():
    f = (d + e*x)**m*(f + g*x)*sqrt(a + b*x + c*x**2)
    F = g*(d + e*x)**(m + 2)*sqrt(a + b*x + c*x**2)*appellf1(m + 2, sympy.S(-1)/2, sympy.S(-1)/2, m + 3, 2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))), 2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(e**2*(m + 2)*sqrt(-2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))) + 1)*sqrt(-2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))) + 1)) + (d + e*x)**(m + 1)*(-d*g + e*f)*sqrt(a + b*x + c*x**2)*appellf1(m + 1, sympy.S(-1)/2, sympy.S(-1)/2, m + 2, 2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))), 2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(e**2*(m + 1)*sqrt(-2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))) + 1)*sqrt(-2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_926():
    f = (d + e*x)**m*sqrt(a + b*x + c*x**2)
    F = (d + e*x)**(m + 1)*sqrt(a + b*x + c*x**2)*appellf1(m + 1, sympy.S(-1)/2, sympy.S(-1)/2, m + 2, 2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))), 2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(e*(m + 1)*sqrt(-2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))) + 1)*sqrt(-2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_927():
    f = (d + e*x)**m*sqrt(a + b*x + c*x**2)/(f + g*x)
    F = sympy.Function('Unintegrable')(((((Symbol('d') + (Symbol('e') * x)))**(Symbol('m')) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Symbol('f') + (Symbol('g') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_928():
    f = (d + e*x)**m*(f + g*x)**2/sqrt(a + b*x + c*x**2)
    F = g**2*(d + e*x)**(m + 1)*sqrt(a + b*x + c*x**2)/(c*e*(m + 2)) - g*(d + e*x)**(m + 2)*sqrt(-2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))) + 1)*sqrt(-2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))) + 1)*(b*e*g*(2*m + 3) + c*(2*d*g - 4*e*f*(m + 2)))*appellf1(m + 2, sympy.S.Half, sympy.S.Half, m + 3, 2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))), 2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(2*c*e**3*(m + 2)**2*sqrt(a + b*x + c*x**2)) + (d + e*x)**(m + 1)*(c*(d**2*g**2 - 2*d*e*f*g*(m + 2) + e**2*f**2*(m + 2)) + e*g**2*(m + 1)*(-a*e + b*d))*sqrt(-2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))) + 1)*sqrt(-2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))) + 1)*appellf1(m + 1, sympy.S.Half, sympy.S.Half, m + 2, 2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))), 2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(c*e**3*(m + 1)*(m + 2)*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_929():
    f = (d + e*x)**m*(f + g*x)/sqrt(a + b*x + c*x**2)
    F = g*(d + e*x)**(m + 2)*sqrt(-2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))) + 1)*sqrt(-2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))) + 1)*appellf1(m + 2, sympy.S.Half, sympy.S.Half, m + 3, 2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))), 2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(e**2*(m + 2)*sqrt(a + b*x + c*x**2)) + (d + e*x)**(m + 1)*(-d*g + e*f)*sqrt(-2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))) + 1)*sqrt(-2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))) + 1)*appellf1(m + 1, sympy.S.Half, sympy.S.Half, m + 2, 2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))), 2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(e**2*(m + 1)*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_930():
    f = (d + e*x)**m/sqrt(a + b*x + c*x**2)
    F = (d + e*x)**(m + 1)*sqrt(-2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))) + 1)*sqrt(-2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))) + 1)*appellf1(m + 1, sympy.S.Half, sympy.S.Half, m + 2, 2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))), 2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(e*(m + 1)*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_931():
    f = (d + e*x)**m/((f + g*x)*sqrt(a + b*x + c*x**2))
    F = sympy.Function('Unintegrable')((((Symbol('d') + (Symbol('e') * x)))**(Symbol('m')) * (((Symbol('f') + (Symbol('g') * x)) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_932():
    f = (d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**2)
    F = c*(d + e*x)**(m + 2)*(f + g*x)**(n + 1)/(e**2*g*(m + n + 3)) + (d + e*x)**(m + 1)*(f + g*x)**(n + 1)*(b*e*g*(m + n + 3) - c*(d*g*(m + 2*n + 4) + e*f*(m + 2)))/(e**2*g**2*(m + n + 2)*(m + n + 3)) + (d + e*x)**(m + 1)*(f + g*x)**n*(g*(a*e**2*g*(m + n + 3) - c*d*(d*g*(n + 1) + e*f*(m + 2)))*(m + n + 2) - (d*g*(n + 1) + e*f*(m + 1))*(b*e*g*(m + n + 3) - c*(d*g*(m + 2*n + 4) + e*f*(m + 2))))*hyper((-n, m + 1), (m + 2,), -g*(d + e*x)/(-d*g + e*f))/(e**3*g**2*(e*(f + g*x)/(-d*g + e*f))**n*(m + 1)*(m + n + 2)*(m + n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_933():
    f = (d + e*x)**m*(f + g*x)**2*(a + b*x + c*x**2)**p
    F = g**2*(d + e*x)**(m + 1)*(a + b*x + c*x**2)**(p + 1)/(c*e*(m + 2*p + 3)) - g*(d + e*x)**(m + 2)*(b*e*g*(m + p + 2) + 2*c*(d*g*(p + 1) - e*f*(m + 2*p + 3)))*(a + b*x + c*x**2)**p*appellf1(m + 2, -p, -p, m + 3, 2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))), 2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(c*e**3*(m + 2)*(-2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))) + 1)**p*(-2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))) + 1)**p*(m + 2*p + 3)) + (d + e*x)**(m + 1)*(c*(2*d**2*g**2*(p + 1) - 2*d*e*f*g*(m + 2*p + 3) + e**2*f**2*(m + 2*p + 3)) + e*g**2*(m + 1)*(-a*e + b*d))*(a + b*x + c*x**2)**p*appellf1(m + 1, -p, -p, m + 2, 2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))), 2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(c*e**3*(m + 1)*(-2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))) + 1)**p*(-2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))) + 1)**p*(m + 2*p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_934():
    f = (d + e*x)**m*(f + g*x)*(a + b*x + c*x**2)**p
    F = g*(d + e*x)**(m + 2)*(a + b*x + c*x**2)**p*appellf1(m + 2, -p, -p, m + 3, 2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))), 2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(e**2*(m + 2)*(-2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))) + 1)**p*(-2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))) + 1)**p) + (d + e*x)**(m + 1)*(-d*g + e*f)*(a + b*x + c*x**2)**p*appellf1(m + 1, -p, -p, m + 2, 2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))), 2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(e**2*(m + 1)*(-2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))) + 1)**p*(-2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_935():
    f = (d + e*x)**m*(a + b*x + c*x**2)**p
    F = (d + e*x)**(m + 1)*(a + b*x + c*x**2)**p*appellf1(m + 1, -p, -p, m + 2, 2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))), 2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(e*(m + 1)*(-2*c*(d + e*x)/(2*c*d - e*(b - sqrt(-4*a*c + b**2))) + 1)**p*(-2*c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_936():
    f = (d + e*x)**m*(a + b*x + c*x**2)**p/(f + g*x)
    F = sympy.Function('Unintegrable')(((((Symbol('d') + (Symbol('e') * x)))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**(Symbol('p'))) * ((Symbol('f') + (Symbol('g') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_4_d_plus_e_x_pow_m_f_plus_g_x_pow_n_a_plus_b_x_plus_c_x_pow_2_pow_p_937():
    f = 1/(x**2*sqrt(1 - 1/(c**2*x**2))*sqrt(d + e*x))
    F = Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('c') * (Symbol('d') + (Symbol('e') * x))) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.Function('EllipticPi')(Integer(2), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('e')) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1))))) * ((sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))))) * x * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))
    assert integrate(f, x) == F

