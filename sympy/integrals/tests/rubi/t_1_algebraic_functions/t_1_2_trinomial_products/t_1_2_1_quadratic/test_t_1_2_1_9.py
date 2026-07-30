"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.2 Trinomial products/1.2.1 Quadratic/1.2.1.9 P(x) (d+e x)^m (a+b x+c x^2)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, C, a, b, c, d, e, f, g, h, i, j, k, l, m, p, q = symbols('A B C a b c d e f g h i j k l m p q')

def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_1():
    f = (d + e*x)**2*sqrt(d**2 - e**2*x**2)*(A + B*x + C*x**2)
    F = -C*x**3*(d**2 - e**2*x**2)**(sympy.S(3)/2)/6 + d**4*(10*A*e**2 + 4*B*d*e + 3*C*d**2)*atan(e*x/sqrt(d**2 - e**2*x**2))/(16*e**3) + d**2*x*sqrt(d**2 - e**2*x**2)*(10*A*e**2 + 4*B*d*e + 3*C*d**2)/(16*e**2) - d*(d**2 - e**2*x**2)**(sympy.S(3)/2)*(4*C*d**2 + e*(10*A*e + 7*B*d))/(15*e**3) - x**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)*(B*e + 2*C*d)/(5*e) - x*(d**2 - e**2*x**2)**(sympy.S(3)/2)*(3*C*d**2 + 2*e*(A*e + 2*B*d))/(8*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_2():
    f = (d + e*x)*sqrt(d**2 - e**2*x**2)*(A + B*x + C*x**2)
    F = -C*x**2*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(5*e) + d**3*(C*d**2 + e*(4*A*e + B*d))*atan(e*x/sqrt(d**2 - e**2*x**2))/(8*e**3) + d*x*sqrt(d**2 - e**2*x**2)*(C*d**2 + e*(4*A*e + B*d))/(8*e**2) - x*(d**2 - e**2*x**2)**(sympy.S(3)/2)*(B*e + C*d)/(4*e**2) - (d**2 - e**2*x**2)**(sympy.S(3)/2)*(2*C*d**2 + 5*e*(A*e + B*d))/(15*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_3():
    f = sqrt(d**2 - e**2*x**2)*(A + B*x + C*x**2)
    F = -B*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(3*e**2) - C*x*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(4*e**2) + d**2*(4*A*e**2 + C*d**2)*atan(e*x/sqrt(d**2 - e**2*x**2))/(8*e**3) + x*(A/2 + C*d**2/(8*e**2))*sqrt(d**2 - e**2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_4():
    f = sqrt(d**2 - e**2*x**2)*(A + B*x + C*x**2)/(d + e*x)
    F = -C*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(3*e**3) + d*(C*d**2 - e*(-2*A*e + B*d))*atan(e*x/sqrt(d**2 - e**2*x**2))/(2*e**3) + sqrt(d**2 - e**2*x**2)*(C*d**2 - e*(-2*A*e + B*d))/(2*e**3) + (d**2 - e**2*x**2)**(sympy.S(3)/2)*(-B*e + C*d)/(2*e**3*(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_5():
    f = sqrt(d**2 - e**2*x**2)*(A + B*x + C*x**2)/(d + e*x)**2
    F = -C*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(2*e**3*(d + e*x)) - (5*C*d**2 - 2*e*(-A*e + 2*B*d))*atan(e*x/sqrt(d**2 - e**2*x**2))/(2*e**3) - sqrt(d**2 - e**2*x**2)*(5*C*d**2 - 2*e*(-A*e + 2*B*d))/(2*d*e**3) - (d**2 - e**2*x**2)**(sympy.S(3)/2)*(A*e**2 - B*d*e + C*d**2)/(d*e**3*(d + e*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_6():
    f = sqrt(d**2 - e**2*x**2)*(A + B*x + C*x**2)/(d + e*x)**3
    F = -C*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(e**3*(d + e*x)**2) + (-B*e + 3*C*d)*atan(e*x/sqrt(d**2 - e**2*x**2))/e**3 + sqrt(d**2 - e**2*x**2)*(-2*B*e + 6*C*d)/(e**3*(d + e*x)) - (d**2 - e**2*x**2)**(sympy.S(3)/2)*(A*e**2 - B*d*e + C*d**2)/(3*d*e**3*(d + e*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_7():
    f = sqrt(d**2 - e**2*x**2)*(A + B*x + C*x**2)/(d + e*x)**4
    F = -C*atan(e*x/sqrt(d**2 - e**2*x**2))/e**3 - 2*C*sqrt(d**2 - e**2*x**2)/(e**3*(d + e*x)) + (d**2 - e**2*x**2)**(sympy.S(3)/2)*(-B*e + 2*C*d)/(3*d*e**3*(d + e*x)**3) - (d**2 - e**2*x**2)**(sympy.S(3)/2)*(A*e**2 - B*d*e + C*d**2)/(5*d*e**3*(d + e*x)**4) - (d**2 - e**2*x**2)**(sympy.S(3)/2)*(A*e**2 - B*d*e + C*d**2)/(15*d**2*e**3*(d + e*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_8():
    f = sqrt(d**2 - e**2*x**2)*(A + B*x + C*x**2)/(d + e*x)**5
    F = C*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(e**3*(d + e*x)**4) - (d**2 - e**2*x**2)**(sympy.S(3)/2)*(A*e**2 - B*d*e + C*d**2)/(7*d*e**3*(d + e*x)**5) - (d**2 - e**2*x**2)**(sympy.S(3)/2)*(23*C*d**2 + e*(2*A*e + 5*B*d))/(35*d**2*e**3*(d + e*x)**4) - (d**2 - e**2*x**2)**(sympy.S(3)/2)*(23*C*d**2 + e*(2*A*e + 5*B*d))/(105*d**3*e**3*(d + e*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_9():
    f = sqrt(d**2 - e**2*x**2)*(A + B*x + C*x**2)/(d + e*x)**6
    F = C*(d**2 - e**2*x**2)**(sympy.S(3)/2)/(2*e**3*(d + e*x)**5) - (d**2 - e**2*x**2)**(sympy.S(3)/2)*(A*e**2 - B*d*e + C*d**2)/(9*d*e**3*(d + e*x)**6) - (d**2 - e**2*x**2)**(sympy.S(3)/2)*(11*C*d**2 + 2*e*(A*e + 2*B*d))/(42*d**2*e**3*(d + e*x)**5) - (d**2 - e**2*x**2)**(sympy.S(3)/2)*(11*C*d**2 + 2*e*(A*e + 2*B*d))/(105*d**3*e**3*(d + e*x)**4) - (d**2 - e**2*x**2)**(sympy.S(3)/2)*(11*C*d**2 + 2*e*(A*e + 2*B*d))/(315*d**4*e**3*(d + e*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_10():
    f = (d + e*x)**3*(A + B*x + C*x**2)/sqrt(d**2 - e**2*x**2)
    F = -C*e*x**4*sqrt(d**2 - e**2*x**2)/5 + d**3*(20*A*e**2 + 15*B*d*e + 13*C*d**2)*atan(e*x/sqrt(d**2 - e**2*x**2))/(8*e**3) - d**2*sqrt(d**2 - e**2*x**2)*(55*A*e**2 + 45*B*d*e + 38*C*d**2)/(15*e**3) - d*x*sqrt(d**2 - e**2*x**2)*(12*A*e**2 + 15*B*d*e + 13*C*d**2)/(8*e**2) - x**3*sqrt(d**2 - e**2*x**2)*(B*e/4 + 3*C*d/4) - x**2*sqrt(d**2 - e**2*x**2)*(19*C*d**2 + 5*e*(A*e + 3*B*d))/(15*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_11():
    f = (d + e*x)**2*(A + B*x + C*x**2)/sqrt(d**2 - e**2*x**2)
    F = -C*x**3*sqrt(d**2 - e**2*x**2)/4 + d**2*(12*A*e**2 + 8*B*d*e + 7*C*d**2)*atan(e*x/sqrt(d**2 - e**2*x**2))/(8*e**3) - d*sqrt(d**2 - e**2*x**2)*(4*C*d**2 + e*(6*A*e + 5*B*d))/(3*e**3) - x**2*sqrt(d**2 - e**2*x**2)*(B*e + 2*C*d)/(3*e) - x*sqrt(d**2 - e**2*x**2)*(7*C*d**2 + 4*e*(A*e + 2*B*d))/(8*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_12():
    f = (d + e*x)*(A + B*x + C*x**2)/sqrt(d**2 - e**2*x**2)
    F = -C*x**2*sqrt(d**2 - e**2*x**2)/(3*e) + d*(C*d**2 + e*(2*A*e + B*d))*atan(e*x/sqrt(d**2 - e**2*x**2))/(2*e**3) - x*sqrt(d**2 - e**2*x**2)*(B*e + C*d)/(2*e**2) - sqrt(d**2 - e**2*x**2)*(2*C*d**2 + 3*e*(A*e + B*d))/(3*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_13():
    f = (A + B*x + C*x**2)/sqrt(d**2 - e**2*x**2)
    F = -B*sqrt(d**2 - e**2*x**2)/e**2 - C*x*sqrt(d**2 - e**2*x**2)/(2*e**2) + (2*A*e**2 + C*d**2)*atan(e*x/sqrt(d**2 - e**2*x**2))/(2*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_14():
    f = (A + B*x + C*x**2)/((d + e*x)*sqrt(d**2 - e**2*x**2))
    F = -C*sqrt(d**2 - e**2*x**2)/e**3 - (-B*e + C*d)*atan(e*x/sqrt(d**2 - e**2*x**2))/e**3 - sqrt(d**2 - e**2*x**2)*(A*e**2 - B*d*e + C*d**2)/(d*e**3*(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_15():
    f = (A + B*x + C*x**2)/((d + e*x)**2*sqrt(d**2 - e**2*x**2))
    F = C*atan(e*x/sqrt(d**2 - e**2*x**2))/e**3 + sqrt(d**2 - e**2*x**2)*(-B*e + 2*C*d)/(d*e**3*(d + e*x)) - sqrt(d**2 - e**2*x**2)*(A*e**2 - B*d*e + C*d**2)/(3*d*e**3*(d + e*x)**2) - sqrt(d**2 - e**2*x**2)*(A*e**2 - B*d*e + C*d**2)/(3*d**2*e**3*(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_16():
    f = (A + B*x + C*x**2)/((d + e*x)**3*sqrt(d**2 - e**2*x**2))
    F = C*sqrt(d**2 - e**2*x**2)/(e**3*(d + e*x)**2) - sqrt(d**2 - e**2*x**2)*(A*e**2 - B*d*e + C*d**2)/(5*d*e**3*(d + e*x)**3) - sqrt(d**2 - e**2*x**2)*(7*C*d**2 + e*(2*A*e + 3*B*d))/(15*d**2*e**3*(d + e*x)**2) - sqrt(d**2 - e**2*x**2)*(7*C*d**2 + e*(2*A*e + 3*B*d))/(15*d**3*e**3*(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_17():
    f = (A + B*x + C*x**2)/((d + e*x)**4*sqrt(d**2 - e**2*x**2))
    F = C*sqrt(d**2 - e**2*x**2)/(2*e**3*(d + e*x)**3) - sqrt(d**2 - e**2*x**2)*(A*e**2 - B*d*e + C*d**2)/(7*d*e**3*(d + e*x)**4) - sqrt(d**2 - e**2*x**2)*(6*A*e**2 + 8*B*d*e + 13*C*d**2)/(70*d**2*e**3*(d + e*x)**3) - sqrt(d**2 - e**2*x**2)*(6*A*e**2 + 8*B*d*e + 13*C*d**2)/(105*d**3*e**3*(d + e*x)**2) - sqrt(d**2 - e**2*x**2)*(6*A*e**2 + 8*B*d*e + 13*C*d**2)/(105*d**4*e**3*(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_18():
    f = (a + c*x**2)*(d + e*x)**3*(A + B*x + C*x**2)
    F = C*c*(d + e*x)**8/(8*e**5) - c*(d + e*x)**7*(-B*e + 4*C*d)/(7*e**5) + (d + e*x)**6*(C*a*e**2 + c*(6*C*d**2 - e*(-A*e + 3*B*d)))/(6*e**5) - (d + e*x)**5*(a*e**2*(-B*e + 2*C*d) + c*d*(4*C*d**2 - e*(-2*A*e + 3*B*d)))/(5*e**5) + (d + e*x)**4*(a*e**2 + c*d**2)*(A*e**2 - B*d*e + C*d**2)/(4*e**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_19():
    f = (a + c*x**2)*(d + e*x)**2*(A + B*x + C*x**2)
    F = C*c*(d + e*x)**7/(7*e**5) - c*(d + e*x)**6*(-B*e + 4*C*d)/(6*e**5) + (d + e*x)**5*(C*a*e**2 + c*(6*C*d**2 - e*(-A*e + 3*B*d)))/(5*e**5) - (d + e*x)**4*(a*e**2*(-B*e + 2*C*d) + c*d*(4*C*d**2 - e*(-2*A*e + 3*B*d)))/(4*e**5) + (d + e*x)**3*(a*e**2 + c*d**2)*(A*e**2 - B*d*e + C*d**2)/(3*e**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_20():
    f = (a + c*x**2)*(d + e*x)*(A + B*x + C*x**2)
    F = A*a*d*x + C*c*e*x**6/6 + a*x**2*(A*e + B*d)/2 + c*x**5*(B*e + C*d)/5 + x**4*(A*c*e/4 + B*c*d/4 + C*a*e/4) + x**3*(A*c*d/3 + B*a*e/3 + C*a*d/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_21():
    f = (a + c*x**2)*(A + B*x + C*x**2)
    F = A*a*x + B*a*x**2/2 + B*c*x**4/4 + C*c*x**5/5 + x**3*(A*c/3 + C*a/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_22():
    f = (a + c*x**2)*(A + B*x + C*x**2)/(d + e*x)
    F = C*c*x**4/(4*e) - c*x**3*(-B*e + C*d)/(3*e**2) + x**2*(C*a*e**2 + c*(C*d**2 - e*(-A*e + B*d)))/(2*e**3) - x*(a*e**2*(-B*e + C*d) + c*d*(C*d**2 - e*(-A*e + B*d)))/e**4 + (a*e**2 + c*d**2)*(A*e**2 - B*d*e + C*d**2)*log(d + e*x)/e**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_23():
    f = (a + c*x**2)*(A + B*x + C*x**2)/(d + e*x)**2
    F = C*c*x**3/(3*e**2) - c*x**2*(-B*e + 2*C*d)/(2*e**3) + x*(C*a*e**2 + c*(3*C*d**2 - e*(-A*e + 2*B*d)))/e**4 - (a*e**2*(-B*e + 2*C*d) + c*d*(4*C*d**2 - e*(-2*A*e + 3*B*d)))*log(d + e*x)/e**5 - (a*e**2 + c*d**2)*(A*e**2 - B*d*e + C*d**2)/(e**5*(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_24():
    f = (a + c*x**2)*(A + B*x + C*x**2)/(d + e*x)**3
    F = C*c*x**2/(2*e**3) - c*x*(-B*e + 3*C*d)/e**4 + (C*a*e**2 + c*(6*C*d**2 - e*(-A*e + 3*B*d)))*log(d + e*x)/e**5 + (a*e**2*(-B*e + 2*C*d) + c*d*(4*C*d**2 - e*(-2*A*e + 3*B*d)))/(e**5*(d + e*x)) - (a*e**2 + c*d**2)*(A*e**2 - B*d*e + C*d**2)/(2*e**5*(d + e*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_25():
    f = (a + c*x**2)**2*(d + e*x)**3*(A + B*x + C*x**2)
    F = A*a**2*d**3*x + C*c**2*e**3*x**10/10 + a**2*e*x**4*(3*C*d**2 + e*(A*e + 3*B*d))/4 + a*d*x**3*(A*(3*a*e**2 + 2*c*d**2) + a*d*(3*B*e + C*d))/3 + a*e*x**6*(C*a*e**2 + 2*c*(3*C*d**2 + e*(A*e + 3*B*d)))/6 + c**2*e**2*x**9*(B*e + 3*C*d)/9 + c*e*x**8*(2*C*a*e**2 + c*(3*C*d**2 + e*(A*e + 3*B*d)))/8 + c*x**7*(2*a*e**2*(B*e + 3*C*d) + c*d*(C*d**2 + 3*e*(A*e + B*d)))/7 + x**5*(A*c*d*(6*a*e**2 + c*d**2)/5 + a*(a*e**2*(B*e + 3*C*d) + 2*c*d**2*(3*B*e + C*d))/5) + d**2*(a + c*x**2)**3*(3*A*e + B*d)/(6*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_26():
    f = (a + c*x**2)**2*(d + e*x)**2*(A + B*x + C*x**2)
    F = A*a**2*d**2*x + C*c**2*e**2*x**9/9 + a**2*e*x**4*(B*e + 2*C*d)/4 + a*c*e*x**6*(B*e + 2*C*d)/3 + a*x**3*(A*(a*e**2 + 2*c*d**2) + a*d*(2*B*e + C*d))/3 + c**2*e*x**8*(B*e + 2*C*d)/8 + c*x**7*(2*C*a*e**2 + c*(C*d**2 + e*(A*e + 2*B*d)))/7 + x**5*(A*c*(2*a*e**2 + c*d**2)/5 + a*(C*a*e**2 + 2*c*d*(2*B*e + C*d))/5) + d*(a + c*x**2)**3*(2*A*e + B*d)/(6*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_27():
    f = (a + c*x**2)**2*(d + e*x)*(A + B*x + C*x**2)
    F = A*a**2*d*x + C*a**2*e*x**4/4 + C*a*c*e*x**6/3 + C*c**2*e*x**8/8 + a*x**3*(2*A*c*d + B*a*e + C*a*d)/3 + c**2*x**7*(B*e + C*d)/7 + c*x**5*(A*c*d + 2*a*(B*e + C*d))/5 + (a + c*x**2)**3*(A*e + B*d)/(6*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_28():
    f = (a + c*x**2)**2*(A + B*x + C*x**2)
    F = A*a**2*x + B*(a + c*x**2)**3/(6*c) + C*c**2*x**7/7 + a*x**3*(2*A*c + C*a)/3 + c*x**5*(A*c + 2*C*a)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_29():
    f = (a + c*x**2)**2*(A + B*x + C*x**2)/(d + e*x)
    F = C*c**2*x**6/(6*e) - c**2*x**5*(-B*e + C*d)/(5*e**2) + c*x**4*(2*C*a*e**2 + c*(C*d**2 - e*(-A*e + B*d)))/(4*e**3) - c*x**3*(2*a*e**2*(-B*e + C*d) + c*d*(C*d**2 - e*(-A*e + B*d)))/(3*e**4) + x**2*(C*a**2*e**4 + 2*a*c*e**2*(C*d**2 - e*(-A*e + B*d)) + c**2*d**2*(C*d**2 - e*(-A*e + B*d)))/(2*e**5) - x*(a**2*e**4*(-B*e + C*d) + 2*a*c*d*e**2*(C*d**2 - e*(-A*e + B*d)) + c**2*d**3*(C*d**2 - e*(-A*e + B*d)))/e**6 + (a*e**2 + c*d**2)**2*(A*e**2 - B*d*e + C*d**2)*log(d + e*x)/e**7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_30():
    f = (a + c*x**2)**2*(A + B*x + C*x**2)/(d + e*x)**2
    F = C*c**2*x**5/(5*e**2) - c**2*x**4*(-B*e + 2*C*d)/(4*e**3) + c*x**3*(2*C*a*e**2 + c*(3*C*d**2 - e*(-A*e + 2*B*d)))/(3*e**4) - c*x**2*(2*a*e**2*(-B*e + 2*C*d) + c*d*(4*C*d**2 - e*(-2*A*e + 3*B*d)))/(2*e**5) + x*(C*a**2*e**4 + 2*a*c*e**2*(3*C*d**2 - e*(-A*e + 2*B*d)) + c**2*d**2*(5*C*d**2 - e*(-3*A*e + 4*B*d)))/e**6 - (a*e**2 + c*d**2)*(a*e**2*(-B*e + 2*C*d) + c*d*(6*C*d**2 - e*(-4*A*e + 5*B*d)))*log(d + e*x)/e**7 - (a*e**2 + c*d**2)**2*(A*e**2 - B*d*e + C*d**2)/(e**7*(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_31():
    f = (a + c*x**2)**2*(A + B*x + C*x**2)/(d + e*x)**3
    F = C*c**2*x**4/(4*e**3) - c**2*x**3*(-B*e + 3*C*d)/(3*e**4) + c*x**2*(2*C*a*e**2 + c*(6*C*d**2 - e*(-A*e + 3*B*d)))/(2*e**5) - c*x*(2*a*e**2*(-B*e + 3*C*d) + c*d*(10*C*d**2 - 3*e*(-A*e + 2*B*d)))/e**6 + (C*a**2*e**4 + 2*a*c*e**2*(6*C*d**2 - e*(-A*e + 3*B*d)) + c**2*d**2*(15*C*d**2 - 2*e*(-3*A*e + 5*B*d)))*log(d + e*x)/e**7 + (a*e**2 + c*d**2)*(a*e**2*(-B*e + 2*C*d) + c*d*(6*C*d**2 - e*(-4*A*e + 5*B*d)))/(e**7*(d + e*x)) - (a*e**2 + c*d**2)**2*(A*e**2 - B*d*e + C*d**2)/(2*e**7*(d + e*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_32():
    f = (a + c*x**2)**3*(d + e*x)**3*(A + B*x + C*x**2)
    F = A*a**3*d**3*x + C*c**3*e**3*x**12/12 + a**3*e*x**4*(3*C*d**2 + e*(A*e + 3*B*d))/4 + a**2*d*x**3*(3*A*(a*e**2 + c*d**2) + a*d*(3*B*e + C*d))/3 + a**2*e*x**6*(C*a*e**2 + 3*c*(3*C*d**2 + e*(A*e + 3*B*d)))/6 + 3*a*c*e*x**8*(C*a*e**2 + c*(3*C*d**2 + e*(A*e + 3*B*d)))/8 + a*x**5*(3*A*c*d*(3*a*e**2 + c*d**2) + a*(a*e**2*(B*e + 3*C*d) + 3*c*d**2*(3*B*e + C*d)))/5 + c**3*e**2*x**11*(B*e + 3*C*d)/11 + c**2*e*x**10*(3*C*a*e**2 + c*(3*C*d**2 + e*(A*e + 3*B*d)))/10 + c**2*x**9*(3*a*e**2*(B*e + 3*C*d) + c*d*(C*d**2 + 3*e*(A*e + B*d)))/9 + c*x**7*(A*c*d*(9*a*e**2 + c*d**2) + 3*a*(a*e**2*(B*e + 3*C*d) + c*d**2*(3*B*e + C*d)))/7 + d**2*(a + c*x**2)**4*(3*A*e + B*d)/(8*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_33():
    f = (a + c*x**2)**3*(d + e*x)**2*(A + B*x + C*x**2)
    F = A*a**3*d**2*x + C*c**3*e**2*x**11/11 + a**3*e*x**4*(B*e + 2*C*d)/4 + a**2*c*e*x**6*(B*e + 2*C*d)/2 + a**2*x**3*(A*(a*e**2 + 3*c*d**2) + a*d*(2*B*e + C*d))/3 + 3*a*c**2*e*x**8*(B*e + 2*C*d)/8 + a*x**5*(3*A*c*(a*e**2 + c*d**2) + a*(C*a*e**2 + 3*c*d*(2*B*e + C*d)))/5 + c**3*e*x**10*(B*e + 2*C*d)/10 + c**2*x**9*(3*C*a*e**2 + c*(C*d**2 + e*(A*e + 2*B*d)))/9 + c*x**7*(A*c*(3*a*e**2 + c*d**2) + 3*a*(C*a*e**2 + c*d*(2*B*e + C*d)))/7 + d*(a + c*x**2)**4*(2*A*e + B*d)/(8*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_34():
    f = (a + c*x**2)**3*(d + e*x)*(A + B*x + C*x**2)
    F = A*a**3*d*x + C*a**3*e*x**4/4 + C*a**2*c*e*x**6/2 + 3*C*a*c**2*e*x**8/8 + C*c**3*e*x**10/10 + a**2*x**3*(3*A*c*d + B*a*e + C*a*d)/3 + 3*a*c*x**5*(A*c*d + B*a*e + C*a*d)/5 + c**3*x**9*(B*e + C*d)/9 + c**2*x**7*(A*c*d + 3*a*(B*e + C*d))/7 + (a + c*x**2)**4*(A*e + B*d)/(8*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_35():
    f = (a + c*x**2)**3*(A + B*x + C*x**2)
    F = A*a**3*x + B*(a + c*x**2)**4/(8*c) + C*c**3*x**9/9 + a**2*x**3*(3*A*c + C*a)/3 + 3*a*c*x**5*(A*c + C*a)/5 + c**2*x**7*(A*c + 3*C*a)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_36():
    f = (a + c*x**2)**3*(A + B*x + C*x**2)/(d + e*x)
    F = C*c**3*(d + e*x)**8/(8*e**9) - c**3*(d + e*x)**7*(-B*e + 8*C*d)/(7*e**9) + c**2*(d + e*x)**6*(3*C*a*e**2 + c*(28*C*d**2 - e*(-A*e + 7*B*d)))/(6*e**9) - c**2*(d + e*x)**5*(3*a*e**2*(-B*e + 6*C*d) + c*d*(56*C*d**2 - 3*e*(-2*A*e + 7*B*d)))/(5*e**9) + c*(d + e*x)**4*(3*C*a**2*e**4 + 3*a*c*e**2*(15*C*d**2 - e*(-A*e + 5*B*d)) + 5*c**2*d**2*(14*C*d**2 - e*(-3*A*e + 7*B*d)))/(4*e**9) - c*(d + e*x)**3*(3*a**2*e**4*(-B*e + 4*C*d) + 6*a*c*d*e**2*(10*C*d**2 - e*(-2*A*e + 5*B*d)) + c**2*d**3*(56*C*d**2 - 5*e*(-4*A*e + 7*B*d)))/(3*e**9) - x*(a*e**2 + c*d**2)**2*(a*e**2*(-B*e + 2*C*d) + c*d*(8*C*d**2 - e*(-6*A*e + 7*B*d)))/e**8 + (d + e*x)**2*(a*e**2 + c*d**2)*(C*a**2*e**4 + a*c*e**2*(17*C*d**2 - 3*e*(-A*e + 3*B*d)) + c**2*d**2*(28*C*d**2 - 3*e*(-5*A*e + 7*B*d)))/(2*e**9) + (a*e**2 + c*d**2)**3*(A*e**2 - B*d*e + C*d**2)*log(d + e*x)/e**9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_37():
    f = (a + c*x**2)**3*(A + B*x + C*x**2)/(d + e*x)**2
    F = C*c**3*x**7/(7*e**2) - c**3*x**6*(-B*e + 2*C*d)/(6*e**3) + c**2*x**5*(3*C*a*e**2 + c*(3*C*d**2 - e*(-A*e + 2*B*d)))/(5*e**4) - c**2*x**4*(3*a*e**2*(-B*e + 2*C*d) + c*d*(4*C*d**2 - e*(-2*A*e + 3*B*d)))/(4*e**5) + c*x**3*(3*C*a**2*e**4 + 3*a*c*e**2*(3*C*d**2 - e*(-A*e + 2*B*d)) + c**2*d**2*(5*C*d**2 - e*(-3*A*e + 4*B*d)))/(3*e**6) - c*x**2*(3*a**2*e**4*(-B*e + 2*C*d) + 3*a*c*d*e**2*(4*C*d**2 - e*(-2*A*e + 3*B*d)) + c**2*d**3*(6*C*d**2 - e*(-4*A*e + 5*B*d)))/(2*e**7) + x*(C*a**3*e**6 + 3*a**2*c*e**4*(3*C*d**2 - e*(-A*e + 2*B*d)) + 3*a*c**2*d**2*e**2*(5*C*d**2 - e*(-3*A*e + 4*B*d)) + c**3*d**4*(7*C*d**2 - e*(-5*A*e + 6*B*d)))/e**8 - (a*e**2 + c*d**2)**2*(a*e**2*(-B*e + 2*C*d) + c*d*(8*C*d**2 - e*(-6*A*e + 7*B*d)))*log(d + e*x)/e**9 - (a*e**2 + c*d**2)**3*(A*e**2 - B*d*e + C*d**2)/(e**9*(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_38():
    f = (a + c*x**2)**3*(A + B*x + C*x**2)/(d + e*x)**3
    F = C*c**3*x**6/(6*e**3) - c**3*x**5*(-B*e + 3*C*d)/(5*e**4) + c**2*x**4*(3*C*a*e**2 + c*(6*C*d**2 - e*(-A*e + 3*B*d)))/(4*e**5) - c**2*x**3*(3*a*e**2*(-B*e + 3*C*d) + c*d*(10*C*d**2 - 3*e*(-A*e + 2*B*d)))/(3*e**6) + c*x**2*(3*C*a**2*e**4 + 3*a*c*e**2*(6*C*d**2 - e*(-A*e + 3*B*d)) + c**2*d**2*(15*C*d**2 - 2*e*(-3*A*e + 5*B*d)))/(2*e**7) - c*x*(3*a**2*e**4*(-B*e + 3*C*d) + 3*a*c*d*e**2*(10*C*d**2 - 3*e*(-A*e + 2*B*d)) + c**2*d**3*(21*C*d**2 - 5*e*(-2*A*e + 3*B*d)))/e**8 + (a*e**2 + c*d**2)*(C*a**2*e**4 + a*c*e**2*(17*C*d**2 - 3*e*(-A*e + 3*B*d)) + c**2*d**2*(28*C*d**2 - 3*e*(-5*A*e + 7*B*d)))*log(d + e*x)/e**9 + (a*e**2 + c*d**2)**2*(a*e**2*(-B*e + 2*C*d) + c*d*(8*C*d**2 - e*(-6*A*e + 7*B*d)))/(e**9*(d + e*x)) - (a*e**2 + c*d**2)**3*(A*e**2 - B*d*e + C*d**2)/(2*e**9*(d + e*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_39():
    f = (a + b*x**2)*(-a*d + 4*b*c*x + 3*b*d*x**2)/(c + d*x)**2
    F = (a + b*x**2)**2/(c + d*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_40():
    f = (a + b*x**2)*(-a*d + b*x*(4*c + 3*d*x))/(c + d*x)**2
    F = (a + b*x**2)**2/(c + d*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_41():
    f = (a + b*x**2)**2*(-a*d + 6*b*c*x + 5*b*d*x**2)/(c + d*x)**2
    F = (a + b*x**2)**3/(c + d*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_42():
    f = (a + b*x**2)**2*(-a*d + b*x*(6*c + 5*d*x))/(c + d*x)**2
    F = (a + b*x**2)**3/(c + d*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_43():
    f = (d + e*x)**3*(A + B*x + C*x**2)/(a + c*x**2)
    F = C*e**3*x**4/(4*c) + e**2*x**3*(B*e + 3*C*d)/(3*c) - e*x**2*(C*a*e**2 - c*(3*C*d**2 + e*(A*e + 3*B*d)))/(2*c**2) - x*(a*e**2*(B*e + 3*C*d) - c*d*(C*d**2 + 3*e*(A*e + B*d)))/c**2 + (B*c*d*(-3*a*e**2 + c*d**2) + e*(A*c - C*a)*(-a*e**2 + 3*c*d**2))*log(a + c*x**2)/(2*c**3) + (A*c*d*(-3*a*e**2 + c*d**2) + a*(a*e**2*(B*e + 3*C*d) - c*d**2*(3*B*e + C*d)))*atan(sqrt(c)*x/sqrt(a))/(sqrt(a)*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_44():
    f = (d + e*x)**2*(A + B*x + C*x**2)/(a + c*x**2)
    F = C*e**2*x**3/(3*c) + e*x**2*(B*e + 2*C*d)/(2*c) - x*(C*a*e**2 - c*(C*d**2 + e*(A*e + 2*B*d)))/c**2 + (2*A*c*d*e - B*a*e**2 + B*c*d**2 - 2*C*a*d*e)*log(a + c*x**2)/(2*c**2) + (A*c*(-a*e**2 + c*d**2) + a*(C*a*e**2 - c*d*(2*B*e + C*d)))*atan(sqrt(c)*x/sqrt(a))/(sqrt(a)*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_45():
    f = (d + e*x)*(A + B*x + C*x**2)/(a + c*x**2)
    F = C*e*x**2/(2*c) + x*(B*e + C*d)/c + (A*c*e + B*c*d - C*a*e)*log(a + c*x**2)/(2*c**2) + (A*c*d - a*(B*e + C*d))*atan(sqrt(c)*x/sqrt(a))/(sqrt(a)*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_46():
    f = (A + B*x + C*x**2)/(a + c*x**2)
    F = B*log(a + c*x**2)/(2*c) + C*x/c + (A*c - C*a)*atan(sqrt(c)*x/sqrt(a))/(sqrt(a)*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_47():
    f = (A + B*x + C*x**2)/((a + c*x**2)*(d + e*x))
    F = (A*e**2 - B*d*e + C*d**2)*log(d + e*x)/(e*(a*e**2 + c*d**2)) + (-A*c*e + B*c*d + C*a*e)*log(a + c*x**2)/(2*c*(a*e**2 + c*d**2)) + (A*c*d + B*a*e - C*a*d)*atan(sqrt(c)*x/sqrt(a))/(sqrt(a)*sqrt(c)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_48():
    f = (A + B*x + C*x**2)/((a + c*x**2)*(d + e*x)**2)
    F = (-2*A*c*d*e - B*a*e**2 + B*c*d**2 + 2*C*a*d*e)*log(a + c*x**2)/(2*(a*e**2 + c*d**2)**2) - (-2*A*c*d*e - B*a*e**2 + B*c*d**2 + 2*C*a*d*e)*log(d + e*x)/(a*e**2 + c*d**2)**2 - (A*e**2 - B*d*e + C*d**2)/(e*(d + e*x)*(a*e**2 + c*d**2)) + (A*c*(-a*e**2 + c*d**2) + a*(C*a*e**2 - c*d*(-2*B*e + C*d)))*atan(sqrt(c)*x/sqrt(a))/(sqrt(a)*sqrt(c)*(a*e**2 + c*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_49():
    f = (A + B*x + C*x**2)/((a + c*x**2)*(d + e*x)**3)
    F = (B*c*d*(-3*a*e**2 + c*d**2) - e*(A*c - C*a)*(-a*e**2 + 3*c*d**2))*log(a + c*x**2)/(2*(a*e**2 + c*d**2)**3) - (B*c*d*(-3*a*e**2 + c*d**2) - e*(A*c - C*a)*(-a*e**2 + 3*c*d**2))*log(d + e*x)/(a*e**2 + c*d**2)**3 + (-2*A*c*d*e - B*a*e**2 + B*c*d**2 + 2*C*a*d*e)/((d + e*x)*(a*e**2 + c*d**2)**2) - (A*e**2 - B*d*e + C*d**2)/(2*e*(d + e*x)**2*(a*e**2 + c*d**2)) + sqrt(c)*(A*c*d*(-3*a*e**2 + c*d**2) - a*(-a*e**2*(-B*e + 3*C*d) + c*d**2*(-3*B*e + C*d)))*atan(sqrt(c)*x/sqrt(a))/(sqrt(a)*(a*e**2 + c*d**2)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_50():
    f = (d + e*x)**3*(A + B*x + C*x**2)/(a + c*x**2)**2
    F = -e*(2*C*a*e**2 - c*(3*C*d**2 + e*(A*e + 3*B*d)))*log(a + c*x**2)/(2*c**3) - (d + e*x)**3*(B*a - x*(A*c - C*a))/(2*a*c*(a + c*x**2)) - e**3*x**2*(A*c - 2*C*a)/(2*a*c**2) - 3*e**2*x*(A*c*d - a*(B*e + 3*C*d))/(2*a*c**2) + (A*c*d*(3*a*e**2 + c*d**2) - a*(3*a*e**2*(B*e + 3*C*d) - c*d**2*(3*B*e + C*d)))*atan(sqrt(c)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_51():
    f = (d + e*x)**2*(A + B*x + C*x**2)/(a + c*x**2)**2
    F = e*(B*e + 2*C*d)*log(a + c*x**2)/(2*c**2) - (d + e*x)**2*(B*a - x*(A*c - C*a))/(2*a*c*(a + c*x**2)) - e**2*x*(A*c - 3*C*a)/(2*a*c**2) + (a*e**2*(A*c - 3*C*a) + c*d*(A*c*d + 2*B*a*e + C*a*d))*atan(sqrt(c)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_52():
    f = (d + e*x)*(A + B*x + C*x**2)/(a + c*x**2)**2
    F = C*e*log(a + c*x**2)/(2*c**2) - (d + e*x)*(B*a - x*(A*c - C*a))/(2*a*c*(a + c*x**2)) + (A*c*d + B*a*e + C*a*d)*atan(sqrt(c)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_53():
    f = (A + B*x + C*x**2)/(a + c*x**2)**2
    F = -(B*a - x*(A*c - C*a))/(2*a*c*(a + c*x**2)) + (A*c + C*a)*atan(sqrt(c)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_54():
    f = (A + B*x + C*x**2)/((a + c*x**2)**2*(d + e*x))
    F = -e*(A*e**2 - B*d*e + C*d**2)*log(a + c*x**2)/(2*(a*e**2 + c*d**2)**2) + e*(A*e**2 - B*d*e + C*d**2)*log(d + e*x)/(a*e**2 + c*d**2)**2 - (a*(-A*c*e + B*c*d + C*a*e) - c*x*(A*c*d + B*a*e - C*a*d))/(2*a*c*(a + c*x**2)*(a*e**2 + c*d**2)) + (A*c*d*(3*a*e**2 + c*d**2) + a*(-B*e + C*d)*(-a*e**2 + c*d**2))*atan(sqrt(c)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*sqrt(c)*(a*e**2 + c*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_55():
    f = (A + B*x + C*x**2)/((a + c*x**2)**2*(d + e*x)**2)
    F = e*(a*e**2*(-B*e + 2*C*d) - c*d*(2*C*d**2 - e*(-4*A*e + 3*B*d)))*log(a + c*x**2)/(2*(a*e**2 + c*d**2)**3) - e*(a*e**2*(-B*e + 2*C*d) - c*d*(2*C*d**2 - e*(-4*A*e + 3*B*d)))*log(d + e*x)/(a*e**2 + c*d**2)**3 - e*(A*e**2 - B*d*e + C*d**2)/((d + e*x)*(a*e**2 + c*d**2)**2) - (a*(-2*A*c*d*e - B*a*e**2 + B*c*d**2 + 2*C*a*d*e) - x*(A*c*(-a*e**2 + c*d**2) + a*(C*a*e**2 - c*d*(-2*B*e + C*d))))/(2*a*(a + c*x**2)*(a*e**2 + c*d**2)**2) + (A*c*(-3*a**2*e**4 + 6*a*c*d**2*e**2 + c**2*d**4) + a*(C*a**2*e**4 - 6*a*c*d*e**2*(-B*e + C*d) + c**2*d**3*(-2*B*e + C*d)))*atan(sqrt(c)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*sqrt(c)*(a*e**2 + c*d**2)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_56():
    f = (A + B*x + C*x**2)/((a + c*x**2)**2*(d + e*x)**3)
    F = -e*(C*a**2*e**4 - 2*a*c*e**2*(4*C*d**2 - e*(-A*e + 3*B*d)) + c**2*d**2*(3*C*d**2 - 2*e*(-5*A*e + 3*B*d)))*log(a + c*x**2)/(2*(a*e**2 + c*d**2)**4) + e*(C*a**2*e**4 - 2*a*c*e**2*(4*C*d**2 - e*(-A*e + 3*B*d)) + c**2*d**2*(3*C*d**2 - 2*e*(-5*A*e + 3*B*d)))*log(d + e*x)/(a*e**2 + c*d**2)**4 + e*(a*e**2*(-B*e + 2*C*d) - c*d*(2*C*d**2 - e*(-4*A*e + 3*B*d)))/((d + e*x)*(a*e**2 + c*d**2)**3) - e*(A*e**2 - B*d*e + C*d**2)/(2*(d + e*x)**2*(a*e**2 + c*d**2)**2) - (a*(B*c*d*(-3*a*e**2 + c*d**2) - e*(A*c - C*a)*(-a*e**2 + 3*c*d**2)) - c*x*(A*c*d*(-3*a*e**2 + c*d**2) - a*(-a*e**2*(-B*e + 3*C*d) + c*d**2*(-3*B*e + C*d))))/(2*a*(a + c*x**2)*(a*e**2 + c*d**2)**3) + sqrt(c)*(A*c*d*(-15*a**2*e**4 + 10*a*c*d**2*e**2 + c**2*d**4) - a*(-3*a**2*e**4*(-B*e + 3*C*d) + 2*a*c*d**2*e**2*(-9*B*e + 7*C*d) - c**2*d**4*(-3*B*e + C*d)))*atan(sqrt(c)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*(a*e**2 + c*d**2)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_57():
    f = (d + e*x)**3*(A + B*x + C*x**2)/(a + c*x**2)**3
    F = C*e**3*log(a + c*x**2)/(2*c**3) - (d + e*x)**3*(B*a - x*(A*c - C*a))/(4*a*c*(a + c*x**2)**2) - (d + e*x)*(a*e*(3*A*c*d + 3*B*a*e + 5*C*a*d) - x*(3*A*c**2*d**2 - a*(4*C*a*e**2 - c*d*(3*B*e + C*d))))/(8*a**2*c**2*(a + c*x**2)) + (3*a*e**2*(A*c*d + B*a*e + 3*C*a*d) + c*d**2*(3*A*c*d + 3*B*a*e + C*a*d))*atan(sqrt(c)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_58():
    f = (d + e*x)*(A + B*x + C*x**2)/(a + c*x**2)**3
    F = -(d + e*x)*(B*a - x*(A*c - C*a))/(4*a*c*(a + c*x**2)**2) - (2*a*e*(A*c + C*a) - c*x*(3*A*c*d + B*a*e + C*a*d))/(8*a**2*c**2*(a + c*x**2)) + (3*A*c*d + B*a*e + C*a*d)*atan(sqrt(c)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_59():
    f = (A + B*x + C*x**2)/(a + c*x**2)**3
    F = -(B*a - x*(A*c - C*a))/(4*a*c*(a + c*x**2)**2) + x*(3*A*c + C*a)/(8*a**2*c*(a + c*x**2)) + (3*A*c + C*a)*atan(sqrt(c)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_60():
    f = (A + B*x + C*x**2)/((a + c*x**2)**3*(d + e*x))
    F = -e**3*(A*e**2 - B*d*e + C*d**2)*log(a + c*x**2)/(2*(a*e**2 + c*d**2)**3) + e**3*(A*e**2 - B*d*e + C*d**2)*log(d + e*x)/(a*e**2 + c*d**2)**3 - (a*(-A*c*e + B*c*d + C*a*e) - c*x*(A*c*d + B*a*e - C*a*d))/(4*a*c*(a + c*x**2)**2*(a*e**2 + c*d**2)) + (4*a**2*e*(A*e**2 - B*d*e + C*d**2) + x*(A*c*d*(7*a*e**2 + 3*c*d**2) + a*(-B*e + C*d)*(-3*a*e**2 + c*d**2)))/(8*a**2*(a + c*x**2)*(a*e**2 + c*d**2)**2) + (A*c*d*(15*a**2*e**4 + 10*a*c*d**2*e**2 + 3*c**2*d**4) + a*(-B*e + C*d)*(-3*a**2*e**4 + 6*a*c*d**2*e**2 + c**2*d**4))*atan(sqrt(c)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*sqrt(c)*(a*e**2 + c*d**2)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_61():
    f = (A + B*x + C*x**2)/((a + c*x**2)**3*(d + e*x)**2)
    F = e**3*(a*e**2*(-B*e + 2*C*d) - c*d*(4*C*d**2 - e*(-6*A*e + 5*B*d)))*log(a + c*x**2)/(2*(a*e**2 + c*d**2)**4) - e**3*(a*e**2*(-B*e + 2*C*d) - c*d*(4*C*d**2 - e*(-6*A*e + 5*B*d)))*log(d + e*x)/(a*e**2 + c*d**2)**4 - e**3*(A*e**2 - B*d*e + C*d**2)/((d + e*x)*(a*e**2 + c*d**2)**3) - (a*(-2*A*c*d*e - B*a*e**2 + B*c*d**2 + 2*C*a*d*e) - x*(A*c*(-a*e**2 + c*d**2) + a*(C*a*e**2 - c*d*(-2*B*e + C*d))))/(4*a*(a + c*x**2)**2*(a*e**2 + c*d**2)**2) - (4*a**2*e*(a*e**2*(-B*e + 2*C*d) - c*d*(2*C*d**2 - e*(-4*A*e + 3*B*d))) - x*(A*c*(-7*a**2*e**4 + 12*a*c*d**2*e**2 + 3*c**2*d**4) + a*(3*C*a**2*e**4 - 2*a*c*d*e**2*(-7*B*e + 6*C*d) + c**2*d**3*(-2*B*e + C*d))))/(8*a**2*(a + c*x**2)*(a*e**2 + c*d**2)**3) + (3*A*c*(-5*a**3*e**6 + 15*a**2*c*d**2*e**4 + 5*a*c**2*d**4*e**2 + c**3*d**6) + a*(3*C*a**3*e**6 - 3*a**2*c*d*e**4*(-10*B*e + 11*C*d) + a*c**2*d**3*e**2*(-20*B*e + 13*C*d) + c**3*d**5*(-2*B*e + C*d)))*atan(sqrt(c)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*sqrt(c)*(a*e**2 + c*d**2)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_62():
    f = (A + B*x + C*x**2)/((a + c*x**2)**3*(d + e*x)**3)
    F = -e**3*(C*a**2*e**4 - a*c*e**2*(3*A*e**2 - 9*B*d*e + 13*C*d**2) + c**2*d**2*(10*C*d**2 - 3*e*(-7*A*e + 5*B*d)))*log(a + c*x**2)/(2*(a*e**2 + c*d**2)**5) + e**3*(C*a**2*e**4 - a*c*e**2*(3*A*e**2 - 9*B*d*e + 13*C*d**2) + c**2*d**2*(10*C*d**2 - 3*e*(-7*A*e + 5*B*d)))*log(d + e*x)/(a*e**2 + c*d**2)**5 + e**3*(a*e**2*(-B*e + 2*C*d) - c*d*(4*C*d**2 - e*(-6*A*e + 5*B*d)))/((d + e*x)*(a*e**2 + c*d**2)**4) - e**3*(A*e**2 - B*d*e + C*d**2)/(2*(d + e*x)**2*(a*e**2 + c*d**2)**3) - (a*(B*c*d*(-3*a*e**2 + c*d**2) - e*(A*c - C*a)*(-a*e**2 + 3*c*d**2)) - c*x*(A*c*d*(-3*a*e**2 + c*d**2) - a*(-a*e**2*(-B*e + 3*C*d) + c*d**2*(-3*B*e + C*d))))/(4*a*(a + c*x**2)**2*(a*e**2 + c*d**2)**3) + (4*a**2*e*(C*a**2*e**4 - 2*a*c*e**2*(4*C*d**2 - e*(-A*e + 3*B*d)) + c**2*d**2*(3*C*d**2 - 2*e*(-5*A*e + 3*B*d))) + c*x*(3*A*c*d*(-11*a**2*e**4 + 6*a*c*d**2*e**2 + c**2*d**4) - a*(-7*a**2*e**4*(-B*e + 3*C*d) + 2*a*c*d**2*e**2*(-19*B*e + 13*C*d) - c**2*d**4*(-3*B*e + C*d))))/(8*a**2*(a + c*x**2)*(a*e**2 + c*d**2)**4) + sqrt(c)*(3*A*c*d*(-35*a**3*e**6 + 35*a**2*c*d**2*e**4 + 7*a*c**2*d**4*e**2 + c**3*d**6) + a*(15*a**3*e**6*(-B*e + 3*C*d) - 5*a**2*c*d**2*e**4*(-27*B*e + 25*C*d) + a*c**2*d**4*e**2*(-45*B*e + 23*C*d) + c**3*d**6*(-3*B*e + C*d)))*atan(sqrt(c)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*(a*e**2 + c*d**2)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_63():
    f = (d + e*x)**4*(A + B*x + C*x**2)/(a + c*x**2)**4
    F = -(d + e*x)**4*(B*a - x*(A*c - C*a))/(6*a*c*(a + c*x**2)**3) - (d + e*x)**3*(a*e*(A*c + 5*C*a) - c*x*(5*A*c*d + 4*B*a*e + C*a*d))/(24*a**2*c**2*(a + c*x**2)**2) - (d + e*x)*(a*e - c*d*x)*(a*e**2*(A*c + 5*C*a) + c*d*(5*A*c*d + 4*B*a*e + C*a*d))/(16*a**3*c**3*(a + c*x**2)) + (a*e**2 + c*d**2)*(a*e**2*(A*c + 5*C*a) + c*d*(5*A*c*d + 4*B*a*e + C*a*d))*atan(sqrt(c)*x/sqrt(a))/(16*a**(sympy.S(7)/2)*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_64():
    f = (d + e*x)**2*(A + B*x + C*x**2)/(a + c*x**2)**4
    F = -(d + e*x)**2*(B*a - x*(A*c - C*a))/(6*a*c*(a + c*x**2)**3) - (2*a*e*(4*A*c*d + B*a*e + 2*C*a*d) + x*(3*a*e**2*(A*c + C*a) - c*d*(5*A*c*d + 2*B*a*e + C*a*d)))/(24*a**2*c**2*(a + c*x**2)**2) + x*(a*e**2*(A*c + C*a) + c*d*(5*A*c*d + 2*B*a*e + C*a*d))/(16*a**3*c**2*(a + c*x**2)) + (a*e**2*(A*c + C*a) + c*d*(5*A*c*d + 2*B*a*e + C*a*d))*atan(sqrt(c)*x/sqrt(a))/(16*a**(sympy.S(7)/2)*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_65():
    f = (d + e*x)*(A + B*x + C*x**2)/(a + c*x**2)**4
    F = -(d + e*x)*(B*a - x*(A*c - C*a))/(6*a*c*(a + c*x**2)**3) - (2*a*e*(2*A*c + C*a) - c*x*(5*A*c*d + B*a*e + C*a*d))/(24*a**2*c**2*(a + c*x**2)**2) + x*(5*A*c*d + B*a*e + C*a*d)/(16*a**3*c*(a + c*x**2)) + (5*A*c*d + B*a*e + C*a*d)*atan(sqrt(c)*x/sqrt(a))/(16*a**(sympy.S(7)/2)*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_66():
    f = (A + B*x + C*x**2)/(a + c*x**2)**4
    F = -(B*a - x*(A*c - C*a))/(6*a*c*(a + c*x**2)**3) + x*(5*A*c + C*a)/(24*a**2*c*(a + c*x**2)**2) + x*(5*A*c + C*a)/(16*a**3*c*(a + c*x**2)) + (5*A*c + C*a)*atan(sqrt(c)*x/sqrt(a))/(16*a**(sympy.S(7)/2)*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_67():
    f = x**3*(x**2 + x + 1)/(x**2 + 1)**2
    F = -x**3/(2*x**2 + 2) + x**2/2 + 3*x/2 - log(x**2 + 1)/2 - 3*atan(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_68():
    f = x**2*(x**2 + x + 1)/(x**2 + 1)**2
    F = -x**2/(2*x**2 + 2) + x + log(x**2 + 1)/2 - atan(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_69():
    f = x*(x**2 + x + 1)/(x**2 + 1)**2
    F = -x/(2*x**2 + 2) + log(x**2 + 1)/2 + atan(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_70():
    f = (x**2 + x + 1)/(x**2 + 1)**2
    F = atan(x) - 1/(2*x**2 + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_71():
    f = (x**2 + x + 1)/(x*(x**2 + 1)**2)
    F = x/(2*x**2 + 2) + log(x) - log(x**2 + 1)/2 + atan(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_72():
    f = (x**2 + x + 1)/(x**2*(x**2 + 1)**2)
    F = log(x) - log(x**2 + 1)/2 - atan(x) + 1/(2*x**2 + 2) - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_73():
    f = (x**2 + x + 1)/(x**3*(x**2 + 1)**2)
    F = -x/(2*x**2 + 2) - log(x) + log(x**2 + 1)/2 - 3*atan(x)/2 - 1/x - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_74():
    f = (3*x**2 + 12*x + 2)/(x**2 + 4)**2
    F = -(5*x + 24)/(4*x**2 + 16) + 7*atan(x/2)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_75():
    f = sqrt(a + c*x**2)*(g + h*x)**3*(d + e*x + f*x**2)
    F = a*(a**2*h**2*(e*h + 3*f*g) - 2*a*c*g*(f*g**2 + 3*h*(d*h + e*g)) + 8*c**2*d*g**3)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(16*c**(sympy.S(5)/2)) + f*(a + c*x**2)**(sympy.S(3)/2)*(g + h*x)**4/(7*c*h) - (a + c*x**2)**(sympy.S(3)/2)*(g + h*x)**3*(-7*e*h + 3*f*g)/(42*c*h) + x*sqrt(a + c*x**2)*(a**2*h**2*(e*h + 3*f*g) - 2*a*c*g*(f*g**2 + 3*h*(d*h + e*g)) + 8*c**2*d*g**3)/(16*c**2) - (a + c*x**2)**(sympy.S(3)/2)*(g + h*x)**2*(8*a*f*h**2 + c*(3*f*g**2 - 7*h*(2*d*h + e*g)))/(70*c**2*h) + (a + c*x**2)**(sympy.S(3)/2)*(64*a**2*f*h**4 - 16*a*c*h**2*(15*f*g**2 + 7*h*(d*h + 3*e*g)) - 8*c**2*g**2*(3*f*g**2 - 7*h*(12*d*h + e*g)) - 3*c*h*x*(a*h**2*(35*e*h + 41*f*g) + 2*c*g*(3*f*g**2 - 7*h*(7*d*h + e*g))))/(840*c**3*h)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_76():
    f = sqrt(a + c*x**2)*(g + h*x)**2*(d + e*x + f*x**2)
    F = a*(a**2*f*h**2 - 2*a*c*(f*g**2 + h*(d*h + 2*e*g)) + 8*c**2*d*g**2)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(16*c**(sympy.S(5)/2)) + f*(a + c*x**2)**(sympy.S(3)/2)*(g + h*x)**3/(6*c*h) - (a + c*x**2)**(sympy.S(3)/2)*(g + h*x)**2*(-2*e*h + f*g)/(10*c*h) + x*sqrt(a + c*x**2)*(a**2*f*h**2 - 2*a*c*(f*g**2 + h*(d*h + 2*e*g)) + 8*c**2*d*g**2)/(16*c**2) - (a + c*x**2)**(sympy.S(3)/2)*(16*a*h**2*(e*h + 2*f*g) + 8*c*g*(f*g**2 - 2*h*(5*d*h + e*g)) - 3*h*x*(-2*c*g*(-2*e*h + f*g) + h**2*(-5*a*f + 10*c*d)))/(120*c**2*h)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_77():
    f = sqrt(a + c*x**2)*(g + h*x)*(d + e*x + f*x**2)
    F = a*(-a*e*h - a*f*g + 4*c*d*g)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(8*c**(sympy.S(3)/2)) + f*(a + c*x**2)**(sympy.S(3)/2)*(g + h*x)**2/(5*c*h) + x*sqrt(a + c*x**2)*(-a*(e*h + f*g) + 4*c*d*g)/(8*c) - (a + c*x**2)**(sympy.S(3)/2)*(8*a*f*h**2 + 3*c*h*x*(-5*e*h + 3*f*g) + 4*c*(3*f*g**2 - 5*h*(d*h + e*g)))/(60*c**2*h)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_78():
    f = sqrt(a + c*x**2)*(d + e*x + f*x**2)
    F = a*(-a*f + 4*c*d)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(8*c**(sympy.S(3)/2)) + e*(a + c*x**2)**(sympy.S(3)/2)/(3*c) + f*x*(a + c*x**2)**(sympy.S(3)/2)/(4*c) + x*sqrt(a + c*x**2)*(-a*f + 4*c*d)/(8*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_79():
    f = sqrt(a + c*x**2)*(d + e*x + f*x**2)/(g + h*x)
    F = sqrt(a + c*x**2)*(2*d*h**2 - 2*e*g*h + 2*f*g**2 - h*x*(-e*h + f*g))/(2*h**3) - sqrt(a*h**2 + c*g**2)*(d*h**2 - e*g*h + f*g**2)*atanh((a*h - c*g*x)/(sqrt(a + c*x**2)*sqrt(a*h**2 + c*g**2)))/h**4 + f*(a + c*x**2)**(sympy.S(3)/2)/(3*c*h) - (2*c*d*g*h**2 + (a*h**2 + 2*c*g**2)*(-e*h + f*g))*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(2*sqrt(c)*h**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_80():
    f = sqrt(a + c*x**2)*(d + e*x + f*x**2)/(g + h*x)**2
    F = -(a + c*x**2)**(sympy.S(3)/2)*(d*h**2 - e*g*h + f*g**2)/(h*(g + h*x)*(a*h**2 + c*g**2)) - sqrt(a + c*x**2)*(2*a*h**2*(-e*h + 2*f*g) + 2*c*g*(3*f*g**2 - h*(-d*h + 2*e*g)) - h*x*(a*f*h**2 + c*(3*f*g**2 - 2*h*(-d*h + e*g))))/(2*h**3*(a*h**2 + c*g**2)) + (a*h**2*(-e*h + 2*f*g) + c*g*(3*f*g**2 - h*(-d*h + 2*e*g)))*atanh((a*h - c*g*x)/(sqrt(a + c*x**2)*sqrt(a*h**2 + c*g**2)))/(h**4*sqrt(a*h**2 + c*g**2)) + (a*f*h**2 + 2*c*(3*f*g**2 - h*(-d*h + 2*e*g)))*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(2*sqrt(c)*h**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_81():
    f = sqrt(a + c*x**2)*(d + e*x + f*x**2)/(g + h*x)**3
    F = -sqrt(c)*(-e*h + 3*f*g)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/h**4 - (a + c*x**2)**(sympy.S(3)/2)*(d*h**2 - e*g*h + f*g**2)/(2*h*(g + h*x)**2*(a*h**2 + c*g**2)) + sqrt(a + c*x**2)*(h*x*(2*a*f*h**2 + c*(3*f*g**2 - h*(-d*h + e*g))) + (a*h**2 + c*g**2)*(-2*e*h + 6*f*g))/(2*h**3*(g + h*x)*(a*h**2 + c*g**2)) - (2*a**2*f*h**4 + a*c*h**2*(9*f*g**2 - h*(-d*h + 3*e*g)) + 2*c**2*g**3*(-e*h + 3*f*g))*atanh((a*h - c*g*x)/(sqrt(a + c*x**2)*sqrt(a*h**2 + c*g**2)))/(2*h**4*(a*h**2 + c*g**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_82():
    f = sqrt(a + c*x**2)*(d + e*x + f*x**2)/(g + h*x)**4
    F = sqrt(c)*f*atanh(sqrt(c)*x/sqrt(a + c*x**2))/h**4 + c*(a**2*h**4*(-e*h + 4*f*g) + a*c*g*h**2*(-d*h**2 + 5*f*g**2) + 2*c**2*f*g**5)*atanh((a*h - c*g*x)/(sqrt(a + c*x**2)*sqrt(a*h**2 + c*g**2)))/(2*h**4*(a*h**2 + c*g**2)**(sympy.S(5)/2)) - (a + c*x**2)**(sympy.S(3)/2)*(d*h**2 - e*g*h + f*g**2)/(3*h*(g + h*x)**3*(a*h**2 + c*g**2)) - sqrt(a + c*x**2)*(a**2*e*h**5 + a*c*g*h**2*(d*h**2 + 3*f*g**2) + 2*c**2*f*g**5 + h*x*(2*a**2*f*h**4 + a*c*g*h**2*(-e*h + 6*f*g) + c**2*(-d*g**2*h**2 + 3*f*g**4)))/(2*h**3*(g + h*x)**2*(a*h**2 + c*g**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_83():
    f = sqrt(a + c*x**2)*(d + e*x + f*x**2)/(g + h*x)**5
    F = -a*c*(4*a**2*f*h**2 - a*c*(f*g**2 - h*(-d*h + 5*e*g)) + 4*c**2*d*g**2)*atanh((a*h - c*g*x)/(sqrt(a + c*x**2)*sqrt(a*h**2 + c*g**2)))/(8*(a*h**2 + c*g**2)**(sympy.S(7)/2)) - sqrt(a + c*x**2)*(a*h - c*g*x)*(4*a**2*f*h**2 - a*c*(f*g**2 - h*(-d*h + 5*e*g)) + 4*c**2*d*g**2)/(8*(g + h*x)**2*(a*h**2 + c*g**2)**3) + (a + c*x**2)**(sympy.S(3)/2)*(4*a*h**2*(-e*h + 2*f*g) + c*g*(3*f*g**2 + h*(-5*d*h + e*g)))/(12*h*(g + h*x)**3*(a*h**2 + c*g**2)**2) - (a + c*x**2)**(sympy.S(3)/2)*(d*h**2 - e*g*h + f*g**2)/(4*h*(g + h*x)**4*(a*h**2 + c*g**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_84():
    f = sqrt(a + c*x**2)*(d + e*x + f*x**2)/(g + h*x)**6
    F = -a*c**2*(a**2*h**2*(-e*h + 6*f*g) - a*c*g*(f*g**2 - 3*h*(-d*h + 2*e*g)) + 4*c**2*d*g**3)*atanh((a*h - c*g*x)/(sqrt(a + c*x**2)*sqrt(a*h**2 + c*g**2)))/(8*(a*h**2 + c*g**2)**(sympy.S(9)/2)) - c*sqrt(a + c*x**2)*(a*h - c*g*x)*(a**2*h**2*(-e*h + 6*f*g) - a*c*g*(f*g**2 - 3*h*(-d*h + 2*e*g)) + 4*c**2*d*g**3)/(8*(g + h*x)**2*(a*h**2 + c*g**2)**4) - (a + c*x**2)**(sympy.S(3)/2)*(20*a**2*f*h**4 - a*c*h**2*(18*f*g**2 - h*(-8*d*h + 33*e*g)) - c**2*g**2*(3*f*g**2 + h*(-27*d*h + 2*e*g)))/(60*h*(g + h*x)**3*(a*h**2 + c*g**2)**3) + (a + c*x**2)**(sympy.S(3)/2)*(5*a*h**2*(-e*h + 2*f*g) + c*g*(3*f*g**2 + h*(-7*d*h + 2*e*g)))/(20*h*(g + h*x)**4*(a*h**2 + c*g**2)**2) - (a + c*x**2)**(sympy.S(3)/2)*(d*h**2 - e*g*h + f*g**2)/(5*h*(g + h*x)**5*(a*h**2 + c*g**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_85():
    f = (a + c*x**2)**(sympy.S(3)/2)*(g + h*x)**3*(d + e*x + f*x**2)
    F = a**2*(3*a**2*h**2*(e*h + 3*f*g) - 8*a*c*g*(f*g**2 + 3*h*(d*h + e*g)) + 48*c**2*d*g**3)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(128*c**(sympy.S(5)/2)) + a*x*sqrt(a + c*x**2)*(3*a**2*h**2*(e*h + 3*f*g) - 8*a*c*g*(f*g**2 + 3*h*(d*h + e*g)) + 48*c**2*d*g**3)/(128*c**2) + f*(a + c*x**2)**(sympy.S(5)/2)*(g + h*x)**4/(9*c*h) - (a + c*x**2)**(sympy.S(5)/2)*(g + h*x)**3*(-9*e*h + 5*f*g)/(72*c*h) + x*(a + c*x**2)**(sympy.S(3)/2)*(3*a**2*h**2*(e*h + 3*f*g) - 8*a*c*g*(f*g**2 + 3*h*(d*h + e*g)) + 48*c**2*d*g**3)/(192*c**2) + (a + c*x**2)**(sympy.S(5)/2)*(g + h*x)**2*(-3*c*g*(-9*e*h + 5*f*g) + h**2*(-32*a*f + 72*c*d))/(504*c**2*h) + (a + c*x**2)**(sympy.S(5)/2)*(128*a**2*f*h**4 - 32*a*c*h**2*(17*f*g**2 + 9*h*(d*h + 3*e*g)) - 12*c**2*g**2*(5*f*g**2 - 3*h*(64*d*h + 3*e*g)) - 5*c*h*x*(a*h**2*(63*e*h + 61*f*g) + 2*c*g*(5*f*g**2 - 9*h*(12*d*h + e*g))))/(5040*c**3*h)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_86():
    f = (a + c*x**2)**(sympy.S(3)/2)*(g + h*x)**2*(d + e*x + f*x**2)
    F = a**2*(3*a**2*f*h**2 - 8*a*c*(f*g**2 + h*(d*h + 2*e*g)) + 48*c**2*d*g**2)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(128*c**(sympy.S(5)/2)) + a*x*sqrt(a + c*x**2)*(3*a**2*f*h**2 - 8*a*c*(f*g**2 + h*(d*h + 2*e*g)) + 48*c**2*d*g**2)/(128*c**2) + f*(a + c*x**2)**(sympy.S(5)/2)*(g + h*x)**3/(8*c*h) - (a + c*x**2)**(sympy.S(5)/2)*(g + h*x)**2*(-8*e*h + 5*f*g)/(56*c*h) + x*(a + c*x**2)**(sympy.S(3)/2)*(3*a**2*f*h**2 - 8*a*c*(f*g**2 + h*(d*h + 2*e*g)) + 48*c**2*d*g**2)/(192*c**2) - (a + c*x**2)**(sympy.S(5)/2)*(96*a*h**2*(e*h + 2*f*g) + 12*c*g*(5*f*g**2 - 8*h*(7*d*h + e*g)) - 5*h*x*(-2*c*g*(-8*e*h + 5*f*g) + h**2*(-21*a*f + 56*c*d)))/(1680*c**2*h)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_87():
    f = (a + c*x**2)**(sympy.S(3)/2)*(g + h*x)*(d + e*x + f*x**2)
    F = a**2*(-a*e*h - a*f*g + 6*c*d*g)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(16*c**(sympy.S(3)/2)) + a*x*sqrt(a + c*x**2)*(-a*e*h - a*f*g + 6*c*d*g)/(16*c) + f*(a + c*x**2)**(sympy.S(5)/2)*(g + h*x)**2/(7*c*h) + x*(a + c*x**2)**(sympy.S(3)/2)*(-a*(e*h + f*g) + 6*c*d*g)/(24*c) - (a + c*x**2)**(sympy.S(5)/2)*(12*a*f*h**2 + 5*c*h*x*(-7*e*h + 5*f*g) + 6*c*(5*f*g**2 - 7*h*(d*h + e*g)))/(210*c**2*h)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_88():
    f = (a + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)
    F = a**2*(-a*f + 6*c*d)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(16*c**(sympy.S(3)/2)) + a*x*sqrt(a + c*x**2)*(-a*f + 6*c*d)/(16*c) + e*(a + c*x**2)**(sympy.S(5)/2)/(5*c) + f*x*(a + c*x**2)**(sympy.S(5)/2)/(6*c) + x*(a + c*x**2)**(sympy.S(3)/2)*(-a*f + 6*c*d)/(24*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_89():
    f = (a + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)/(g + h*x)
    F = (a + c*x**2)**(sympy.S(3)/2)*(4*d*h**2 - 4*e*g*h + 4*f*g**2 - 3*h*x*(-e*h + f*g))/(12*h**3) + sqrt(a + c*x**2)*(-h*x*(4*c*d*g*h**2 + (3*a*h**2 + 4*c*g**2)*(-e*h + f*g)) + (8*a*h**2 + 8*c*g**2)*(d*h**2 - e*g*h + f*g**2))/(8*h**5) - (a*h**2 + c*g**2)**(sympy.S(3)/2)*(d*h**2 - e*g*h + f*g**2)*atanh((a*h - c*g*x)/(sqrt(a + c*x**2)*sqrt(a*h**2 + c*g**2)))/h**6 + f*(a + c*x**2)**(sympy.S(5)/2)/(5*c*h) - (3*a**2*h**4*(-e*h + f*g) + 12*a*c*g*h**2*(f*g**2 - h*(-d*h + e*g)) + 8*c**2*g**3*(f*g**2 - h*(-d*h + e*g)))*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(8*sqrt(c)*h**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_90():
    f = (a + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)/(g + h*x)**2
    F = -(a + c*x**2)**(sympy.S(5)/2)*(d*h**2 - e*g*h + f*g**2)/(h*(g + h*x)*(a*h**2 + c*g**2)) - (a + c*x**2)**(sympy.S(3)/2)*(4*a*h**2*(-e*h + 2*f*g) + 4*c*g*(5*f*g**2 - h*(-3*d*h + 4*e*g)) - 3*h*x*(a*f*h**2 + c*(5*f*g**2 - 4*h*(-d*h + e*g))))/(12*h**3*(a*h**2 + c*g**2)) - sqrt(a + c*x**2)*(8*a*h**2*(-e*h + 2*f*g) + 8*c*g*(5*f*g**2 - h*(-3*d*h + 4*e*g)) - h*x*(3*a*f*h**2 + 12*c*d*h**2 - 16*c*e*g*h + 20*c*f*g**2))/(8*h**5) + sqrt(a*h**2 + c*g**2)*(a*h**2*(-e*h + 2*f*g) + c*g*(5*f*g**2 - h*(-3*d*h + 4*e*g)))*atanh((a*h - c*g*x)/(sqrt(a + c*x**2)*sqrt(a*h**2 + c*g**2)))/h**6 + (3*a**2*f*h**4 + 12*a*c*h**2*(3*f*g**2 - h*(-d*h + 2*e*g)) + 8*c**2*g**2*(5*f*g**2 - h*(-3*d*h + 4*e*g)))*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(8*sqrt(c)*h**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_91():
    f = (a + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)/(g + h*x)**3
    F = -sqrt(c)*(3*a*h**2*(-e*h + 3*f*g) + 2*c*g*(10*f*g**2 - 3*h*(-d*h + 2*e*g)))*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(2*h**6) - (a + c*x**2)**(sympy.S(5)/2)*(d*h**2 - e*g*h + f*g**2)/(2*h*(g + h*x)**2*(a*h**2 + c*g**2)) - (a + c*x**2)**(sympy.S(3)/2)*(-2*a*h*(-3*e*h + 7*f*g) + 2*c*g*(-3*d*h + 6*e*g - 10*f*g**2/h) - x*(2*a*f*h**2 + c*(5*f*g**2 - 3*h*(-d*h + e*g))))/(6*h**2*(g + h*x)*(a*h**2 + c*g**2)) + sqrt(a + c*x**2)*(2*a**2*f*h**4 + a*c*h**2*(19*f*g**2 - 3*h*(-d*h + 3*e*g)) + 2*c**2*g**2*(10*f*g**2 - 3*h*(-d*h + 2*e*g)) - c*h*x*(a*h**2*(-3*e*h + 7*f*g) + c*g*(10*f*g**2 - 3*h*(-d*h + 2*e*g))))/(2*h**5*(a*h**2 + c*g**2)) - (2*a**2*f*h**4 + a*c*h**2*(19*f*g**2 - 3*h*(-d*h + 3*e*g)) + 2*c**2*g**2*(10*f*g**2 - 3*h*(-d*h + 2*e*g)))*atanh((a*h - c*g*x)/(sqrt(a + c*x**2)*sqrt(a*h**2 + c*g**2)))/(2*h**6*sqrt(a*h**2 + c*g**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_92():
    f = (a + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)/(g + h*x)**4
    F = sqrt(c)*(3*a*f*h**2 + 2*c*(10*f*g**2 - h*(-d*h + 4*e*g)))*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(2*h**6) + c*(3*a**2*h**4*(-e*h + 4*f*g) + 3*a*c*g*h**2*(11*f*g**2 - h*(-d*h + 4*e*g)) + 2*c**2*g**3*(10*f*g**2 - h*(-d*h + 4*e*g)))*atanh((a*h - c*g*x)/(sqrt(a + c*x**2)*sqrt(a*h**2 + c*g**2)))/(2*h**6*(a*h**2 + c*g**2)**(sympy.S(3)/2)) - (a + c*x**2)**(sympy.S(5)/2)*(d*h**2 - e*g*h + f*g**2)/(3*h*(g + h*x)**3*(a*h**2 + c*g**2)) - (a + c*x**2)**(sympy.S(3)/2)*(-3*a*h*(-e*h + 3*f*g) + c*g*(-d*h + 4*e*g - 10*f*g**2/h) - x*(3*a*f*h**2 + c*(5*f*g**2 - 2*h*(-d*h + e*g))))/(6*h**2*(g + h*x)**2*(a*h**2 + c*g**2)) - sqrt(a + c*x**2)*(c*h*x*(3*a*h**2*(-e*h + 3*f*g) + c*g*(10*f*g**2 - h*(-d*h + 4*e*g))) + (a*h**2 + c*g**2)*(3*a*f*h**2 + 2*c*(10*f*g**2 - h*(-d*h + 4*e*g))))/(2*h**5*(g + h*x)*(a*h**2 + c*g**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_93():
    f = (a + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)/(g + h*x)**5
    F = -c**(sympy.S(3)/2)*(-e*h + 5*f*g)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/h**6 + c*sqrt(a + c*x**2)*(h*x*(12*a**2*f*h**4 + a*c*h**2*(35*f*g**2 - h*(-3*d*h + 7*e*g)) + 4*c**2*g**3*(-e*h + 5*f*g)) + (a*h**2 + c*g**2)**2*(-8*e*h + 40*f*g))/(8*h**5*(g + h*x)*(a*h**2 + c*g**2)**2) - c*(12*a**3*f*h**6 + 3*a**2*c*h**4*(25*f*g**2 - h*(-d*h + 5*e*g)) + 20*a*c**2*g**3*h**2*(-e*h + 5*f*g) + 8*c**3*g**5*(-e*h + 5*f*g))*atanh((a*h - c*g*x)/(sqrt(a + c*x**2)*sqrt(a*h**2 + c*g**2)))/(8*h**6*(a*h**2 + c*g**2)**(sympy.S(5)/2)) - (a + c*x**2)**(sympy.S(5)/2)*(d*h**2 - e*g*h + f*g**2)/(4*h*(g + h*x)**4*(a*h**2 + c*g**2)) + (a + c*x**2)**(sympy.S(3)/2)*(4*a**2*h**4*(-2*e*h + f*g) - a*c*g*h**2*(25*f*g**2 - h*(-9*d*h + 5*e*g)) - 4*c**2*g**4*(-e*h + 5*f*g) - 3*h*x*(4*a**2*f*h**4 + a*c*h**2*(17*f*g**2 - h*(-d*h + 5*e*g)) + 2*c**2*g**2*(5*f*g**2 - h*(d*h + e*g))))/(24*h**3*(g + h*x)**3*(a*h**2 + c*g**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_94():
    f = (a + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)/(g + h*x)**6
    F = c**(sympy.S(3)/2)*f*atanh(sqrt(c)*x/sqrt(a + c*x**2))/h**6 + c**2*(3*a**3*h**6*(-e*h + 6*f*g) + a**2*c*g*h**4*(-3*d*h**2 + 35*f*g**2) + 28*a*c**2*f*g**5*h**2 + 8*c**3*f*g**7)*atanh((a*h - c*g*x)/(sqrt(a + c*x**2)*sqrt(a*h**2 + c*g**2)))/(8*h**6*(a*h**2 + c*g**2)**(sympy.S(7)/2)) - c*sqrt(a + c*x**2)*(-a**3*h**6*(-3*e*h + 2*f*g) + a**2*c*g*h**4*(3*d*h**2 + 13*f*g**2) + 20*a*c**2*f*g**5*h**2 + 8*c**3*f*g**7 + h*x*(8*a**3*f*h**6 + a**2*c*g*h**4*(-3*e*h + 34*f*g) + a*c**2*g**2*h**2*(-3*d*h**2 + 35*f*g**2) + 12*c**3*f*g**6))/(8*h**5*(g + h*x)**2*(a*h**2 + c*g**2)**3) - (a + c*x**2)**(sympy.S(5)/2)*(d*h**2 - e*g*h + f*g**2)/(5*h*(g + h*x)**5*(a*h**2 + c*g**2)) - (a + c*x**2)**(sympy.S(3)/2)*(-a**2*h**4*(-3*e*h + 2*f*g) + a*c*g*h**2*(3*d*h**2 + 5*f*g**2) + 4*c**2*f*g**5 + h*x*(4*a**2*f*h**4 + a*c*g*h**2*(-3*e*h + 14*f*g) + c**2*(-3*d*g**2*h**2 + 7*f*g**4)))/(12*h**3*(g + h*x)**4*(a*h**2 + c*g**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_95():
    f = (a + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)/(g + h*x)**7
    F = -a**2*c**2*(6*a**2*f*h**2 - a*c*(f*g**2 - h*(-d*h + 7*e*g)) + 6*c**2*d*g**2)*atanh((a*h - c*g*x)/(sqrt(a + c*x**2)*sqrt(a*h**2 + c*g**2)))/(16*(a*h**2 + c*g**2)**(sympy.S(9)/2)) - a*c*sqrt(a + c*x**2)*(a*h - c*g*x)*(6*a**2*f*h**2 - a*c*(f*g**2 - h*(-d*h + 7*e*g)) + 6*c**2*d*g**2)/(16*(g + h*x)**2*(a*h**2 + c*g**2)**4) - (a + c*x**2)**(sympy.S(3)/2)*(a*h - c*g*x)*(6*a**2*f*h**2 - a*c*(f*g**2 - h*(-d*h + 7*e*g)) + 6*c**2*d*g**2)/(24*(g + h*x)**4*(a*h**2 + c*g**2)**3) + (a + c*x**2)**(sympy.S(5)/2)*(6*a*h**2*(-e*h + 2*f*g) + c*g*(5*f*g**2 + h*(-7*d*h + e*g)))/(30*h*(g + h*x)**5*(a*h**2 + c*g**2)**2) - (a + c*x**2)**(sympy.S(5)/2)*(d*h**2 - e*g*h + f*g**2)/(6*h*(g + h*x)**6*(a*h**2 + c*g**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_96():
    f = (a + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)/(g + h*x)**8
    F = -a**2*c**3*(a**2*h**2*(-e*h + 8*f*g) - a*c*g*(f*g**2 - h*(-3*d*h + 8*e*g)) + 6*c**2*d*g**3)*atanh((a*h - c*g*x)/(sqrt(a + c*x**2)*sqrt(a*h**2 + c*g**2)))/(16*(a*h**2 + c*g**2)**(sympy.S(11)/2)) - a*c**2*sqrt(a + c*x**2)*(a*h - c*g*x)*(a**2*h**2*(-e*h + 8*f*g) - a*c*g*(f*g**2 - h*(-3*d*h + 8*e*g)) + 6*c**2*d*g**3)/(16*(g + h*x)**2*(a*h**2 + c*g**2)**5) - c*(a + c*x**2)**(sympy.S(3)/2)*(a*h - c*g*x)*(a**2*h**2*(-e*h + 8*f*g) - a*c*g*(f*g**2 - h*(-3*d*h + 8*e*g)) + 6*c**2*d*g**3)/(24*(g + h*x)**4*(a*h**2 + c*g**2)**4) - (a + c*x**2)**(sympy.S(5)/2)*(42*a**2*f*h**4 - a*c*h**2*(26*f*g**2 - h*(-12*d*h + 61*e*g)) - c**2*g**2*(5*f*g**2 + h*(-51*d*h + 2*e*g)))/(210*h*(g + h*x)**5*(a*h**2 + c*g**2)**3) + (a + c*x**2)**(sympy.S(5)/2)*(7*a*h**2*(-e*h + 2*f*g) + c*g*(5*f*g**2 + h*(-9*d*h + 2*e*g)))/(42*h*(g + h*x)**6*(a*h**2 + c*g**2)**2) - (a + c*x**2)**(sympy.S(5)/2)*(d*h**2 - e*g*h + f*g**2)/(7*h*(g + h*x)**7*(a*h**2 + c*g**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_97():
    f = (a + c*x**2)**(sympy.S(5)/2)*(A + B*x + C*x**2)
    F = B*(a + c*x**2)**(sympy.S(7)/2)/(7*c) + C*x*(a + c*x**2)**(sympy.S(7)/2)/(8*c) + 5*a**3*(8*A*c - C*a)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(128*c**(sympy.S(3)/2)) + 5*a**2*x*sqrt(a + c*x**2)*(8*A*c - C*a)/(128*c) + 5*a*x*(a + c*x**2)**(sympy.S(3)/2)*(8*A*c - C*a)/(192*c) + x*(a + c*x**2)**(sympy.S(5)/2)*(8*A*c - C*a)/(48*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_98():
    f = (g + h*x)**3*(d + e*x + f*x**2)/sqrt(a + c*x**2)
    F = f*sqrt(a + c*x**2)*(g + h*x)**4/(5*c*h) - sqrt(a + c*x**2)*(g + h*x)**3*(-5*e*h + f*g)/(20*c*h) + sqrt(a + c*x**2)*(g + h*x)**2*(-3*c*g*(-5*e*h + f*g) + h**2*(-16*a*f + 20*c*d))/(60*c**2*h) + sqrt(a + c*x**2)*(64*a**2*f*h**4 - 16*a*c*h**2*(13*f*g**2 + 5*h*(d*h + 3*e*g)) - 4*c**2*g**2*(3*f*g**2 - 5*h*(16*d*h + 3*e*g)) - c*h*x*(a*h**2*(45*e*h + 71*f*g) + 2*c*g*(3*f*g**2 - 5*h*(10*d*h + 3*e*g))))/(120*c**3*h) + (3*a**2*h**2*(e*h + 3*f*g) - 4*a*c*g*(f*g**2 + 3*h*(d*h + e*g)) + 8*c**2*d*g**3)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(8*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_99():
    f = (g + h*x)**2*(d + e*x + f*x**2)/sqrt(a + c*x**2)
    F = f*sqrt(a + c*x**2)*(g + h*x)**3/(4*c*h) - sqrt(a + c*x**2)*(g + h*x)**2*(-4*e*h + f*g)/(12*c*h) - sqrt(a + c*x**2)*(16*a*h**2*(e*h + 2*f*g) + 4*c*g*(f*g**2 - 4*h*(3*d*h + e*g)) - h*x*(-2*c*g*(-4*e*h + f*g) + h**2*(-9*a*f + 12*c*d)))/(24*c**2*h) + (3*a**2*f*h**2 - 4*a*c*(f*g**2 + h*(d*h + 2*e*g)) + 8*c**2*d*g**2)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(8*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_100():
    f = (g + h*x)*(d + e*x + f*x**2)/sqrt(a + c*x**2)
    F = f*sqrt(a + c*x**2)*(g + h*x)**2/(3*c*h) - sqrt(a + c*x**2)*(4*a*f*h**2 + c*h*x*(-3*e*h + f*g) + 2*c*(f*g**2 - 3*h*(d*h + e*g)))/(6*c**2*h) + (-a*(e*h + f*g) + 2*c*d*g)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(2*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_101():
    f = (d + e*x + f*x**2)/sqrt(a + c*x**2)
    F = e*sqrt(a + c*x**2)/c + f*x*sqrt(a + c*x**2)/(2*c) + (-a*f + 2*c*d)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(2*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_102():
    f = (d + e*x + f*x**2)/(sqrt(a + c*x**2)*(g + h*x))
    F = -(d*h**2 - e*g*h + f*g**2)*atanh((a*h - c*g*x)/(sqrt(a + c*x**2)*sqrt(a*h**2 + c*g**2)))/(h**2*sqrt(a*h**2 + c*g**2)) + f*sqrt(a + c*x**2)/(c*h) - (-e*h + f*g)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(sqrt(c)*h**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_103():
    f = (d + e*x + f*x**2)/(sqrt(a + c*x**2)*(g + h*x)**2)
    F = -sqrt(a + c*x**2)*(d*h**2 - e*g*h + f*g**2)/(h*(g + h*x)*(a*h**2 + c*g**2)) + (a*h**2*(-e*h + 2*f*g) + c*(-d*g*h**2 + f*g**3))*atanh((a*h - c*g*x)/(sqrt(a + c*x**2)*sqrt(a*h**2 + c*g**2)))/(h**2*(a*h**2 + c*g**2)**(sympy.S(3)/2)) + f*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(sqrt(c)*h**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_104():
    f = (d + e*x + f*x**2)/(sqrt(a + c*x**2)*(g + h*x)**3)
    F = -(2*a**2*f*h**2 - a*c*(f*g**2 - h*(-d*h + 3*e*g)) + 2*c**2*d*g**2)*atanh((a*h - c*g*x)/(sqrt(a + c*x**2)*sqrt(a*h**2 + c*g**2)))/(2*(a*h**2 + c*g**2)**(sympy.S(5)/2)) + sqrt(a + c*x**2)*(2*a*h**2*(-e*h + 2*f*g) + c*g*(f*g**2 + h*(-3*d*h + e*g)))/(2*h*(g + h*x)*(a*h**2 + c*g**2)**2) - sqrt(a + c*x**2)*(d*h**2 - e*g*h + f*g**2)/(2*h*(g + h*x)**2*(a*h**2 + c*g**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_105():
    f = (g + h*x)**3*(d + e*x + f*x**2)/(a + c*x**2)**(sympy.S(3)/2)
    F = -(3*a*h**2*(e*h + 3*f*g) - 2*c*g*(f*g**2 + 3*h*(d*h + e*g)))*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(2*c**(sympy.S(5)/2)) - (g + h*x)**3*(a*e - x*(-a*f + c*d))/(a*c*sqrt(a + c*x**2)) - h*sqrt(a + c*x**2)*(g + h*x)**2*(-4*a*f + 3*c*d)/(3*a*c**2) - h*sqrt(a + c*x**2)*(16*a**2*f*h**2 - 4*a*c*(7*f*g**2 + 3*h*(d*h + 3*e*g)) + 12*c**2*d*g**2 + c*h*x*(-9*a*e*h - 11*a*f*g + 6*c*d*g))/(6*a*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_106():
    f = (g + h*x)**2*(d + e*x + f*x**2)/(a + c*x**2)**(sympy.S(3)/2)
    F = (2*c*g*(2*e*h + f*g) + h**2*(-3*a*f + 2*c*d))*atanh(sqrt(c)*x/sqrt(a + c*x**2))/(2*c**(sympy.S(5)/2)) - (g + h*x)**2*(a*e - x*(-a*f + c*d))/(a*c*sqrt(a + c*x**2)) - h*sqrt(a + c*x**2)*(-4*a*(e*h + 2*f*g) + 4*c*d*g + h*x*(-3*a*f + 2*c*d))/(2*a*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_107():
    f = (g + h*x)*(d + e*x + f*x**2)/(a + c*x**2)**(sympy.S(3)/2)
    F = (e*h + f*g)*atanh(sqrt(c)*x/sqrt(a + c*x**2))/c**(sympy.S(3)/2) - (g + h*x)*(a*e - x*(-a*f + c*d))/(a*c*sqrt(a + c*x**2)) - h*sqrt(a + c*x**2)*(-2*a*f + c*d)/(a*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_108():
    f = (d + e*x + f*x**2)/(a + c*x**2)**(sympy.S(3)/2)
    F = f*atanh(sqrt(c)*x/sqrt(a + c*x**2))/c**(sympy.S(3)/2) - (a*e - x*(-a*f + c*d))/(a*c*sqrt(a + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_109():
    f = (d + e*x + f*x**2)/((a + c*x**2)**(sympy.S(3)/2)*(g + h*x))
    F = -(d*h**2 - e*g*h + f*g**2)*atanh((a*h - c*g*x)/(sqrt(a + c*x**2)*sqrt(a*h**2 + c*g**2)))/(a*h**2 + c*g**2)**(sympy.S(3)/2) - (a*(a*f*h - c*d*h + c*e*g) - c*x*(a*e*h - a*f*g + c*d*g))/(a*c*sqrt(a + c*x**2)*(a*h**2 + c*g**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_110():
    f = (d + e*x + f*x**2)/((a + c*x**2)**(sympy.S(3)/2)*(g + h*x)**2)
    F = -h*sqrt(a + c*x**2)*(d*h**2 - e*g*h + f*g**2)/((g + h*x)*(a*h**2 + c*g**2)**2) + (a*h**2*(-e*h + 2*f*g) - c*g*(f*g**2 - h*(-3*d*h + 2*e*g)))*atanh((a*h - c*g*x)/(sqrt(a + c*x**2)*sqrt(a*h**2 + c*g**2)))/(a*h**2 + c*g**2)**(sympy.S(5)/2) - (a*(a*h*(-e*h + 2*f*g) + c*g*(-2*d*h + e*g)) - x*(a**2*f*h**2 - a*c*(f*g**2 - h*(-d*h + 2*e*g)) + c**2*d*g**2))/(a*sqrt(a + c*x**2)*(a*h**2 + c*g**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_111():
    f = (d + e*x + f*x**2)/((a + c*x**2)**(sympy.S(3)/2)*(g + h*x)**3)
    F = h*sqrt(a + c*x**2)*(2*a*h**2*(-e*h + 2*f*g) - c*g*(3*f*g**2 - h*(-7*d*h + 5*e*g)))/(2*(g + h*x)*(a*h**2 + c*g**2)**3) - h*sqrt(a + c*x**2)*(d*h**2 - e*g*h + f*g**2)/(2*(g + h*x)**2*(a*h**2 + c*g**2)**2) - (2*a**2*f*h**4 - a*c*h**2*(3*d*h**2 - 9*e*g*h + 11*f*g**2) + 2*c**2*g**2*(6*d*h**2 - 3*e*g*h + f*g**2))*atanh((a*h - c*g*x)/(sqrt(a + c*x**2)*sqrt(a*h**2 + c*g**2)))/(2*(a*h**2 + c*g**2)**(sympy.S(7)/2)) + (a*(a**2*f*h**3 - a*c*h*(3*f*g**2 - h*(-d*h + 3*e*g)) - c**2*g**2*(-3*d*h + e*g)) + c*x*(a**2*h**2*(-e*h + 3*f*g) - a*c*g*(f*g**2 - 3*h*(-d*h + e*g)) + c**2*d*g**3))/(a*sqrt(a + c*x**2)*(a*h**2 + c*g**2)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_112():
    f = (A + B*x + C*x**2)/(a + c*x**2)**(sympy.S(5)/2)
    F = -(B*a - x*(A*c - C*a))/(3*a*c*(a + c*x**2)**(sympy.S(3)/2)) + x*(2*A*c + C*a)/(3*a**2*c*sqrt(a + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_113():
    f = (A + B*x + C*x**2)/(a + c*x**2)**(sympy.S(7)/2)
    F = -(B*a - x*(A*c - C*a))/(5*a*c*(a + c*x**2)**(sympy.S(5)/2)) + x*(4*A*c + C*a)/(15*a**2*c*(a + c*x**2)**(sympy.S(3)/2)) + x*(8*A*c + 2*C*a)/(15*a**3*c*sqrt(a + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_114():
    f = (A + B*x + C*x**2)/(a + c*x**2)**(sympy.S(9)/2)
    F = -(B*a - x*(A*c - C*a))/(7*a*c*(a + c*x**2)**(sympy.S(7)/2)) + x*(6*A*c + C*a)/(35*a**2*c*(a + c*x**2)**(sympy.S(5)/2)) + x*(24*A*c + 4*C*a)/(105*a**3*c*(a + c*x**2)**(sympy.S(3)/2)) + x*(48*A*c + 8*C*a)/(105*a**4*c*sqrt(a + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_115():
    f = (2*x + 1)**3*(4*x**2 + 3*x + 1)/sqrt(3*x**2 + 2)
    F = 2*(2*x + 1)**4*sqrt(3*x**2 + 2)/15 + 13*(2*x + 1)**3*sqrt(3*x**2 + 2)/60 - 19*(2*x + 1)**2*sqrt(3*x**2 + 2)/540 - (691*x/270 + sympy.S(3937)/810)*sqrt(3*x**2 + 2) + 5*sqrt(3)*asinh(sqrt(6)*x/2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_116():
    f = (2*x + 1)**2*(4*x**2 + 3*x + 1)/sqrt(3*x**2 + 2)
    F = -(x/9 + sympy.S(61)/27)*sqrt(3*x**2 + 2) + (2*x + 1)**3*sqrt(3*x**2 + 2)/6 + 5*(2*x + 1)**2*sqrt(3*x**2 + 2)/18 - sqrt(3)*asinh(sqrt(6)*x/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_117():
    f = (2*x + 1)*(4*x**2 + 3*x + 1)/sqrt(3*x**2 + 2)
    F = (7*x/9 + sympy.S(7)/27)*sqrt(3*x**2 + 2) + 2*(2*x + 1)**2*sqrt(3*x**2 + 2)/9 - 7*sqrt(3)*asinh(sqrt(6)*x/2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_118():
    f = (4*x**2 + 3*x + 1)/((2*x + 1)*sqrt(3*x**2 + 2))
    F = 2*sqrt(3*x**2 + 2)/3 + sqrt(3)*asinh(sqrt(6)*x/2)/6 - sqrt(11)*atanh(sqrt(11)*(4 - 3*x)/(11*sqrt(3*x**2 + 2)))/22
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_119():
    f = (4*x**2 + 3*x + 1)/((2*x + 1)**2*sqrt(3*x**2 + 2))
    F = sqrt(3)*asinh(sqrt(6)*x/2)/3 + 4*sqrt(11)*atanh(sqrt(11)*(4 - 3*x)/(11*sqrt(3*x**2 + 2)))/121 - sqrt(3*x**2 + 2)/(22*x + 11)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_120():
    f = (4*x**2 + 3*x + 1)/((2*x + 1)**3*sqrt(3*x**2 + 2))
    F = -103*sqrt(11)*atanh(sqrt(11)*(4 - 3*x)/(11*sqrt(3*x**2 + 2)))/1331 + 13*sqrt(3*x**2 + 2)/(484*x + 242) - sqrt(3*x**2 + 2)/(22*(2*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_121():
    f = (2*x + 1)**3*(4*x**2 + 3*x + 1)/(3*x**2 + 2)**(sympy.S(3)/2)
    F = 32*x**2*sqrt(3*x**2 + 2)/27 + 4*x*sqrt(3*x**2 + 2) + (279*x + 398)/(54*sqrt(3*x**2 + 2)) + 292*sqrt(3*x**2 + 2)/81 - 38*sqrt(3)*asinh(sqrt(6)*x/2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_122():
    f = (2*x + 1)**2*(4*x**2 + 3*x + 1)/(3*x**2 + 2)**(sympy.S(3)/2)
    F = 8*x*sqrt(3*x**2 + 2)/9 + (70 - 47*x)/(18*sqrt(3*x**2 + 2)) + 28*sqrt(3*x**2 + 2)/9 + 4*sqrt(3)*asinh(sqrt(6)*x/2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_123():
    f = (2*x + 1)*(4*x**2 + 3*x + 1)/(3*x**2 + 2)**(sympy.S(3)/2)
    F = (2 - 51*x)/(18*sqrt(3*x**2 + 2)) + 8*sqrt(3*x**2 + 2)/9 + 10*sqrt(3)*asinh(sqrt(6)*x/2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_124():
    f = (4*x**2 + 3*x + 1)/((2*x + 1)*(3*x**2 + 2)**(sympy.S(3)/2))
    F = (21*x - 38)/(66*sqrt(3*x**2 + 2)) - 2*sqrt(11)*atanh(sqrt(11)*(4 - 3*x)/(11*sqrt(3*x**2 + 2)))/121
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_125():
    f = (4*x**2 + 3*x + 1)/((2*x + 1)**2*(3*x**2 + 2)**(sympy.S(3)/2))
    F = (97*x - 10)/(242*sqrt(3*x**2 + 2)) + 4*sqrt(11)*atanh(sqrt(11)*(4 - 3*x)/(11*sqrt(3*x**2 + 2)))/1331 - 4*sqrt(3*x**2 + 2)/(242*x + 121)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_126():
    f = (4*x**2 + 3*x + 1)/((2*x + 1)**3*(3*x**2 + 2)**(sympy.S(3)/2))
    F = (351*x + 358)/(2662*sqrt(3*x**2 + 2)) - 322*sqrt(11)*atanh(sqrt(11)*(4 - 3*x)/(11*sqrt(3*x**2 + 2)))/14641 + 2*sqrt(3*x**2 + 2)/(2662*x + 1331) - 2*sqrt(3*x**2 + 2)/(121*(2*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_127():
    f = (2*x + 1)**3*(4*x**2 + 3*x + 1)/(3*x**2 + 2)**(sympy.S(5)/2)
    F = (279*x + 398)/(162*(3*x**2 + 2)**(sympy.S(3)/2)) - (465*x + 152)/(54*sqrt(3*x**2 + 2)) + 32*sqrt(3*x**2 + 2)/27 + 8*sqrt(3)*asinh(sqrt(6)*x/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_128():
    f = (2*x + 1)**2*(4*x**2 + 3*x + 1)/(3*x**2 + 2)**(sympy.S(5)/2)
    F = (70 - 47*x)/(54*(3*x**2 + 2)**(sympy.S(3)/2)) - (59*x + 168)/(54*sqrt(3*x**2 + 2)) + 16*sqrt(3)*asinh(sqrt(6)*x/2)/27
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_129():
    f = (2*x + 1)*(4*x**2 + 3*x + 1)/(3*x**2 + 2)**(sympy.S(5)/2)
    F = (2 - 51*x)/(54*(3*x**2 + 2)**(sympy.S(3)/2)) - (16 - 13*x)/(18*sqrt(3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_130():
    f = (4*x**2 + 3*x + 1)/((2*x + 1)*(3*x**2 + 2)**(sympy.S(5)/2))
    F = (21*x - 38)/(198*(3*x**2 + 2)**(sympy.S(3)/2)) + (95*x + 24)/(726*sqrt(3*x**2 + 2)) - 8*sqrt(11)*atanh(sqrt(11)*(4 - 3*x)/(11*sqrt(3*x**2 + 2)))/1331
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_131():
    f = (4*x**2 + 3*x + 1)/((2*x + 1)**2*(3*x**2 + 2)**(sympy.S(5)/2))
    F = (97*x - 10)/(726*(3*x**2 + 2)**(sympy.S(3)/2)) + (887*x + 24)/(7986*sqrt(3*x**2 + 2)) - 32*sqrt(11)*atanh(sqrt(11)*(4 - 3*x)/(11*sqrt(3*x**2 + 2)))/14641 - 16*sqrt(3*x**2 + 2)/(2662*x + 1331)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_132():
    f = (4*x**2 + 3*x + 1)/((2*x + 1)**3*(3*x**2 + 2)**(sympy.S(5)/2))
    F = (351*x + 358)/(7986*(3*x**2 + 2)**(sympy.S(3)/2)) + (2133*x + 1216)/(29282*sqrt(3*x**2 + 2)) - 1216*sqrt(11)*atanh(sqrt(11)*(4 - 3*x)/(11*sqrt(3*x**2 + 2)))/161051 - 8*sqrt(3*x**2 + 2)/(2662*x + 1331) - 8*sqrt(3*x**2 + 2)/(1331*(2*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_133():
    f = (a + c*x**2)**p*(g + h*x)**m*(d + e*x + f*x**2)
    F = -(a + c*x**2)**p*(g + h*x)**(m + 2)*(-e*h*(m + 2*p + 3) + 2*f*g*(p + 1))*appellf1(m + 2, -p, -p, m + 3, (g + h*x)/(g - h*sqrt(-a)/sqrt(c)), (g + h*x)/(g + h*sqrt(-a)/sqrt(c)))/(h**3*(m + 2)*(-(g + h*x)/(g - h*sqrt(-a)/sqrt(c)) + 1)**p*(-(g + h*x)/(g + h*sqrt(-a)/sqrt(c)) + 1)**p*(m + 2*p + 3)) + f*(a + c*x**2)**(p + 1)*(g + h*x)**(m + 1)/(c*h*(m + 2*p + 3)) - (a + c*x**2)**p*(g + h*x)**(m + 1)*(a*f*h**2*(m + 1) - c*(2*f*g**2*(p + 1) - h*(-d*h + e*g)*(m + 2*p + 3)))*appellf1(m + 1, -p, -p, m + 2, (g + h*x)/(g - h*sqrt(-a)/sqrt(c)), (g + h*x)/(g + h*sqrt(-a)/sqrt(c)))/(c*h**3*(m + 1)*(-(g + h*x)/(g - h*sqrt(-a)/sqrt(c)) + 1)**p*(-(g + h*x)/(g + h*sqrt(-a)/sqrt(c)) + 1)**p*(m + 2*p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_134():
    f = sqrt(a + c*x**2)*(g + h*x)**m*(d + e*x + f*x**2)
    F = -sqrt(a + c*x**2)*(g + h*x)**(m + 2)*(-e*h*(m + 4) + 3*f*g)*appellf1(m + 2, sympy.S(-1)/2, sympy.S(-1)/2, m + 3, (g + h*x)/(g - h*sqrt(-a)/sqrt(c)), (g + h*x)/(g + h*sqrt(-a)/sqrt(c)))/(h**3*(m + 2)*(m + 4)*sqrt(-(g + h*x)/(g - h*sqrt(-a)/sqrt(c)) + 1)*sqrt(-(g + h*x)/(g + h*sqrt(-a)/sqrt(c)) + 1)) + f*(a + c*x**2)**(sympy.S(3)/2)*(g + h*x)**(m + 1)/(c*h*(m + 4)) - sqrt(a + c*x**2)*(g + h*x)**(m + 1)*(a*f*h**2*(m + 1) - c*(3*f*g**2 - h*(m + 4)*(-d*h + e*g)))*appellf1(m + 1, sympy.S(-1)/2, sympy.S(-1)/2, m + 2, (g + h*x)/(g - h*sqrt(-a)/sqrt(c)), (g + h*x)/(g + h*sqrt(-a)/sqrt(c)))/(c*h**3*(m + 1)*(m + 4)*sqrt(-(g + h*x)/(g - h*sqrt(-a)/sqrt(c)) + 1)*sqrt(-(g + h*x)/(g + h*sqrt(-a)/sqrt(c)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_135():
    f = (a + c*x**2)**p*(g + h*x)**(-2*p - 3)*(d + e*x + f*x**2)
    F = -f*(a + c*x**2)**p*appellf1(-2*p, -p, -p, 1 - 2*p, (g + h*x)/(g - h*sqrt(-a)/sqrt(c)), (g + h*x)/(g + h*sqrt(-a)/sqrt(c)))/(2*h**3*p*(g + h*x)**(2*p)*(-(g + h*x)/(g - h*sqrt(-a)/sqrt(c)) + 1)**p*(-(g + h*x)/(g + h*sqrt(-a)/sqrt(c)) + 1)**p) - (a + c*x**2)**(p + 1)*(g + h*x)**(-2*p - 2)*(d*h**2 - e*g*h + f*g**2)/(2*h*(p + 1)*(a*h**2 + c*g**2)) + (a + c*x**2)**p*(g + h*x)**(-2*p - 1)*(-sqrt(c)*x + sqrt(-a))*(a*h**2*(-e*h + 2*f*g) + c*(-d*g*h**2 + f*g**3))*hyper((-p, -2*p - 1), (-2*p,), 2*sqrt(c)*sqrt(-a)*(g + h*x)/((sqrt(c)*g - h*sqrt(-a))*(-sqrt(c)*x + sqrt(-a))))/(h**2*(-(sqrt(c)*g + h*sqrt(-a))*(sqrt(c)*x + sqrt(-a))/((sqrt(c)*g - h*sqrt(-a))*(-sqrt(c)*x + sqrt(-a))))**p*(2*p + 1)*(a*h**2 + c*g**2)*(sqrt(c)*g + h*sqrt(-a)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_136():
    f = (d + e*x)**m*(c*e*g*x**2 + f*(b*e - c*d) + x*(b*e*g - c*d*g + c*e*f))*(b*d*e + b*e**2*x - c*d**2 + c*e**2*x**2)**p
    F = g*(d + e*x)**(m - 1)*(b*e**2*x + c*e**2*x**2 - d*(-b*e + c*d))**(p + 2)/(c*e**2*(m + 2*p + 3)) - (c*(d + e*x)/(-b*e + 2*c*d))**(-m - p)*(d + e*x)**m*(b*e*g*(m + p + 1) + c*(d*g*(1 - m) - e*f*(m + 2*p + 3)))*(-b*e + c*d - c*e*x)**2*(b*e**2*x + c*e**2*x**2 - d*(-b*e + c*d))**p*hyper((p + 2, -m - p), (p + 3,), (-b*e + c*d - c*e*x)/(-b*e + 2*c*d))/(c**2*e**2*(p + 2)*(m + 2*p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_137():
    f = (A + C*x**2)*(a + b*x + c*x**2)**4
    F = A*a**4*x + 2*A*a**3*b*x**2 + 2*C*b*c**3*x**10/5 + C*c**4*x**11/11 + a**2*x**3*(4*A*a*c + 6*A*b**2 + C*a**2)/3 + a*b*x**4*(A*(3*a*c + b**2) + C*a**2) + b*c*x**8*(A*c**2 + C*(3*a*c + b**2))/2 + 2*b*x**6*(A*c + C*a)*(3*a*c + b**2)/3 + c**2*x**9*(A*c**2 + 4*C*a*c + 6*C*b**2)/9 + x**7*(2*A*c**2*(2*a*c + 3*b**2)/7 + C*(6*a**2*c**2 + 12*a*b**2*c + b**4)/7) + x**5*(A*(6*a**2*c**2 + 12*a*b**2*c + b**4)/5 + 2*C*a**2*(2*a*c + 3*b**2)/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_138():
    f = (A + C*x**2)*(a + b*x + c*x**2)**3
    F = A*a**3*x + 3*A*a**2*b*x**2/2 + 3*C*b*c**2*x**8/8 + C*c**3*x**9/9 + a*x**3*(3*A*(a*c + b**2) + C*a**2)/3 + b*x**6*(3*A*c**2 + C*(6*a*c + b**2))/6 + b*x**4*(A*(6*a*c + b**2) + 3*C*a**2)/4 + c*x**7*(A*c**2 + C*(3*a*c + 3*b**2))/7 + x**5*(A*c + C*a)*(3*a*c/5 + 3*b**2/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_139():
    f = (A + C*x**2)*(a + b*x + c*x**2)**2
    F = A*a**2*x + A*a*b*x**2 + C*b*c*x**6/3 + C*c**2*x**7/7 + b*x**4*(A*c + C*a)/2 + x**5*(A*c**2/5 + C*(2*a*c + b**2)/5) + x**3*(A*(2*a*c + b**2)/3 + C*a**2/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_140():
    f = (A + C*x**2)*(a + b*x + c*x**2)
    F = A*a*x + A*b*x**2/2 + C*b*x**4/4 + C*c*x**5/5 + x**3*(A*c/3 + C*a/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_141():
    f = (A + C*x**2)/(a + b*x + c*x**2)
    F = -C*b*log(a + b*x + c*x**2)/(2*c**2) + C*x/c - (2*A*c**2 + C*(-2*a*c + b**2))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_142():
    f = (A + C*x**2)/(a + b*x + c*x**2)**2
    F = (4*A*c + 4*C*a)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) - (b*c*(A + C*a/c) + x*(2*A*c**2 + C*(-2*a*c + b**2)))/(c*(-4*a*c + b**2)*(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_143():
    f = (A + C*x**2)/(a + b*x + c*x**2)**3
    F = (b + 2*c*x)*(6*A*c + 2*C*a + C*b**2/c)/(2*(-4*a*c + b**2)**2*(a + b*x + c*x**2)) - (12*A*c**2 + 2*C*(2*a*c + b**2))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2) - (b*c*(A + C*a/c) + x*(2*A*c**2 + C*(-2*a*c + b**2)))/(2*c*(-4*a*c + b**2)*(a + b*x + c*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_144():
    f = (A + C*x**2)/(a + b*x + c*x**2)**4
    F = 8*c*(5*A*c**2 + C*(a*c + b**2))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(7)/2) + (b + 2*c*x)*(5*A*c + C*(a + b**2/c))/(3*(-4*a*c + b**2)**2*(a + b*x + c*x**2)**2) - (b + 2*c*x)*(10*A*c**2 + 2*C*(a*c + b**2))/((-4*a*c + b**2)**3*(a + b*x + c*x**2)) - (b*c*(A + C*a/c) + x*(2*A*c**2 + C*(-2*a*c + b**2)))/(3*c*(-4*a*c + b**2)*(a + b*x + c*x**2)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_145():
    f = (d + e*x)**3*(f + g*x + h*x**2)/(a + b*x + c*x**2)
    F = e**3*h*x**4/(4*c) + e**2*x**3*(-b*e*h + 3*c*d*h + c*e*g)/(3*c**2) + e*x**2*(b**2*e**2*h + c**2*(3*d**2*h + 3*d*e*g + e**2*f) - c*e*(a*e*h + 3*b*d*h + b*e*g))/(2*c**3) - x*(b**3*e**3*h - b*c*e**2*(2*a*e*h + 3*b*d*h + b*e*g) - c**3*d*(d**2*h + 3*d*e*g + 3*e**2*f) + c**2*e*(a*e*(3*d*h + e*g) + b*(3*d**2*h + 3*d*e*g + e**2*f)))/c**4 + (b**4*e**3*h - b**2*c*e**2*(3*a*e*h + 3*b*d*h + b*e*g) + c**4*d**2*(d*g + 3*e*f) - c**3*(a*e*(3*d**2*h + 3*d*e*g + e**2*f) + b*d*(d**2*h + 3*d*e*g + 3*e**2*f)) + c**2*e*(a**2*e**2*h + 2*a*b*e*(3*d*h + e*g) + b**2*(3*d**2*h + 3*d*e*g + e**2*f)))*log(a + b*x + c*x**2)/(2*c**5) - (-b**5*e**3*h + b**3*c*e**2*(5*a*e*h + 3*b*d*h + b*e*g) - b*c**2*e*(5*a**2*e**2*h + 4*a*b*e*(3*d*h + e*g) + b**2*(3*d**2*h + 3*d*e*g + e**2*f)) + 2*c**5*d**3*f - c**4*d*(2*a*(d**2*h + 3*d*e*g + 3*e**2*f) + b*d*(d*g + 3*e*f)) + c**3*(2*a**2*e**2*(3*d*h + e*g) + 3*a*b*e*(3*d**2*h + 3*d*e*g + e**2*f) + b**2*d*(d**2*h + 3*d*e*g + 3*e**2*f)))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**5*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_146():
    f = (d + e*x)**2*(f + g*x + h*x**2)/(a + b*x + c*x**2)
    F = e**2*h*x**3/(3*c) + e*x**2*(-b*e*h + 2*c*d*h + c*e*g)/(2*c**2) + x*(b**2*e**2*h + c**2*(d**2*h + 2*d*e*g + e**2*f) - c*e*(a*e*h + 2*b*d*h + b*e*g))/c**3 + (-b**3*e**2*h + b*c*e*(2*a*e*h + 2*b*d*h + b*e*g) + c**3*d*(d*g + 2*e*f) - c**2*(a*e*(2*d*h + e*g) + b*(d**2*h + 2*d*e*g + e**2*f)))*log(a + b*x + c*x**2)/(2*c**4) - (b**4*e**2*h - b**2*c*e*(4*a*e*h + 2*b*d*h + b*e*g) + 2*c**4*d**2*f - c**3*(2*a*(d**2*h + 2*d*e*g + e**2*f) + b*d*(d*g + 2*e*f)) + c**2*(2*a**2*e**2*h + 3*a*b*e*(2*d*h + e*g) + b**2*(d**2*h + 2*d*e*g + e**2*f)))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**4*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_147():
    f = (d + e*x)*(f + g*x + h*x**2)/(a + b*x + c*x**2)
    F = e*h*x**2/(2*c) + x*(-b*e*h + c*d*h + c*e*g)/c**2 + (b**2*e*h + c**2*(d*g + e*f) - c*(a*e*h + b*d*h + b*e*g))*log(a + b*x + c*x**2)/(2*c**3) - (-b**3*e*h + b*c*(3*a*e*h + b*d*h + b*e*g) + 2*c**3*d*f - c**2*(2*a*d*h + 2*a*e*g + b*d*g + b*e*f))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**3*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_148():
    f = (f + g*x + h*x**2)/(a + b*x + c*x**2)
    F = h*x/c + (-b*h + c*g)*log(a + b*x + c*x**2)/(2*c**2) - (-2*a*c*h + b**2*h - b*c*g + 2*c**2*f)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_149():
    f = (f + g*x + h*x**2)/((d + e*x)*(a + b*x + c*x**2))
    F = (d**2*h - d*e*g + e**2*f)*log(d + e*x)/(e*(a*e**2 - b*d*e + c*d**2)) - (-a*e*h + b*d*h - c*d*g + c*e*f)*log(a + b*x + c*x**2)/(2*c*(a*e**2 - b*d*e + c*d**2)) - (b*h*(-a*e + b*d) + 2*c**2*d*f - c*(2*a*d*h - 2*a*e*g + b*d*g + b*e*f))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c*sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_150():
    f = (f + g*x + h*x**2)/((d + e*x)**2*(a + b*x + c*x**2))
    F = (a*e*(-2*d*h + e*g) - b*(-d**2*h + e**2*f) + c*d*(-d*g + 2*e*f))*log(d + e*x)/(a*e**2 - b*d*e + c*d**2)**2 - (a*e*(-2*d*h + e*g) - b*(-d**2*h + e**2*f) + c*d*(-d*g + 2*e*f))*log(a + b*x + c*x**2)/(2*(a*e**2 - b*d*e + c*d**2)**2) - (2*a**2*e**2*h - a*b*e*(2*d*h + e*g) + b**2*(d**2*h + e**2*f) + 2*c**2*d**2*f - c*(2*a*(d**2*h - 2*d*e*g + e**2*f) + b*d*(d*g + 2*e*f)))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2)**2) - (d**2*h - d*e*g + e**2*f)/(e*(d + e*x)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_151():
    f = (f + g*x + h*x**2)/((d + e*x)**3*(a + b*x + c*x**2))
    F = (c**2*d**2*(-d*g + 3*e*f) - c*(a*e*(3*d**2*h - 3*d*e*g + e**2*f) + b*(-d**3*h + 3*d*e**2*f)) + e**3*(a**2*h - a*b*g + b**2*f))*log(d + e*x)/(a*e**2 - b*d*e + c*d**2)**3 - (c**2*d**2*(-d*g + 3*e*f) - c*(a*e*(3*d**2*h - 3*d*e*g + e**2*f) + b*(-d**3*h + 3*d*e**2*f)) + e**3*(a**2*h - a*b*g + b**2*f))*log(a + b*x + c*x**2)/(2*(a*e**2 - b*d*e + c*d**2)**3) - (-b*e**3*(a**2*h - a*b*g + b**2*f) + 2*c**3*d**3*f - c**2*d*(2*a*(d**2*h - 3*d*e*g + 3*e**2*f) + b*d*(d*g + 3*e*f)) - c*(2*a**2*e**2*(-3*d*h + e*g) - 3*a*b*e*(-d**2*h - d*e*g + e**2*f) - b**2*(d**3*h + 3*d*e**2*f)))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2)**3) - (a*e*(-2*d*h + e*g) - b*(-d**2*h + e**2*f) + c*d*(-d*g + 2*e*f))/((d + e*x)*(a*e**2 - b*d*e + c*d**2)**2) - (d**2*h - d*e*g + e**2*f)/(2*e*(d + e*x)**2*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_152():
    f = (d + e*x)**2*(f + g*x + h*x**2)/(a + b*x + c*x**2)**2
    F = (d + e*x)**2*(c*(2*a*g - b*(a*h/c + f)) - x*(-2*a*c*h + b**2*h - b*c*g + 2*c**2*f))/(c*(-4*a*c + b**2)*(a + b*x + c*x**2)) + e**2*x*(-6*a*c*h + 2*b**2*h - b*c*g + 2*c**2*f)/(c**2*(-4*a*c + b**2)) + e*(-2*b*e*h + 2*c*d*h + c*e*g)*log(a + b*x + c*x**2)/(2*c**3) + (-6*a*c**2*e*(2*a*e*h + 2*b*d*h + b*e*g) - 2*b**4*e**2*h + b**2*c*e*(12*a*e*h + 2*b*d*h + b*e*g) + 4*c**4*d**2*f - c**3*(-4*a*(d**2*h + 2*d*e*g + e**2*f) + 2*b*d*(d*g + 2*e*f)))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**3*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_153():
    f = (d + e*x)*(f + g*x + h*x**2)/(a + b*x + c*x**2)**2
    F = (d + e*x)*(c*(2*a*g - b*(a*h/c + f)) - x*(-2*a*c*h + b**2*h - b*c*g + 2*c**2*f))/(c*(-4*a*c + b**2)*(a + b*x + c*x**2)) + e*h*log(a + b*x + c*x**2)/(2*c**2) + (-6*a*b*c*e*h + b**3*e*h + 4*c**3*d*f - 2*c**2*(-2*a*(d*h + e*g) + b*(d*g + e*f)))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**2*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_154():
    f = (f + g*x + h*x**2)/(a + b*x + c*x**2)**2
    F = (4*a*h - 2*b*g + 4*c*f)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) + (c*(2*a*g - b*(a*h/c + f)) - x*(-2*a*c*h + b**2*h - b*c*g + 2*c**2*f))/(c*(-4*a*c + b**2)*(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_155():
    f = (f + g*x + h*x**2)/((d + e*x)*(a + b*x + c*x**2)**2)
    F = e*(d**2*h - d*e*g + e**2*f)*log(d + e*x)/(a*e**2 - b*d*e + c*d**2)**2 - e*(d**2*h - d*e*g + e**2*f)*log(a + b*x + c*x**2)/(2*(a*e**2 - b*d*e + c*d**2)**2) + (-2*a*(-a*e*h - c*d*g + c*e*f) + b**2*e*f - b*(a*d*h + a*e*g + c*d*f) - x*(b*h*(-a*e + b*d) + 2*c**2*d*f - c*(2*a*d*h - 2*a*e*g + b*d*g + b*e*f)))/((-4*a*c + b**2)*(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2)) + (b*e*(-2*a**2*e**2*h + 4*a*b*d*e*h + b**2*(-d**2*h - d*e*g + e**2*f)) + 4*c**3*d**3*f - 2*c**2*d*(-2*a*(d**2*h - d*e*g + 3*e**2*f) + b*d*(d*g + 3*e*f)) + 2*c*e*(2*a**2*e*(-d*h + e*g) - a*b*(d**2*h + d*e*g + 3*e**2*f) + 2*b**2*d**2*g))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/((-4*a*c + b**2)**(sympy.S(3)/2)*(a*e**2 - b*d*e + c*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_156():
    f = (f + g*x + h*x**2)/((d + e*x)**2*(a + b*x + c*x**2)**2)
    F = -e*(-c*d*(2*d**2*h - 3*d*e*g + 4*e**2*f) + e**2*(2*a*d*h - a*e*g - b*d*g + 2*b*e*f))*log(d + e*x)/(a*e**2 - b*d*e + c*d**2)**3 + e*(-c*d*(2*d**2*h - 3*d*e*g + 4*e**2*f) + e**2*(2*a*d*h - a*e*g - b*d*g + 2*b*e*f))*log(a + b*x + c*x**2)/(2*(a*e**2 - b*d*e + c*d**2)**3) - e*(d**2*h - d*e*g + e**2*f)/((d + e*x)*(a*e**2 - b*d*e + c*d**2)**2) - (2*a*c*(a*e*(-2*d*h + e*g) + c*d*(-d*g + 2*e*f)) + b**3*e**2*f - b**2*e*(a*e*g + 2*c*d*f) + b*(a**2*e**2*h - a*c*(-d**2*h - 2*d*e*g + 3*e**2*f) + c**2*d**2*f) + c*x*(2*a**2*e**2*h - a*b*e*(2*d*h + e*g) + b**2*(d**2*h + e**2*f) + 2*c**2*d**2*f - c*(2*a*(d**2*h - 2*d*e*g + e**2*f) + b*d*(d*g + 2*e*f))))/((-4*a*c + b**2)*(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2)**2) + (-b**3*e**3*(2*a*d*h - a*e*g - b*d*g + 2*b*e*f) + 4*c**4*d**4*f - 2*c**3*d**2*(-2*a*(d**2*h - 2*d*e*g + 6*e**2*f) + b*d*(d*g + 4*e*f)) - 6*c**2*e*(2*a**2*e*(2*d**2*h - 2*d*e*g + e**2*f) + 4*a*b*d*e**2*f - b**2*d**3*g) - c*e*(-4*a**3*e**3*h + 6*a**2*b*e**3*g - 6*a*b**2*e*(2*d**2*h - d*e*g + 2*e**2*f) - b**3*d*(-2*d**2*h - 3*d*e*g + 4*e**2*f)))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/((-4*a*c + b**2)**(sympy.S(3)/2)*(a*e**2 - b*d*e + c*d**2)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_157():
    f = x**3*(x**2 + x + 1)/(x**2 - x + 1)**2
    F = x**2/2 + 3*x + (4 - 2*x)/(3*x**2 - 3*x + 3) + 2*log(x**2 - x + 1) + 10*sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_158():
    f = x**2*(x**2 + x + 1)/(x**2 - x + 1)**2
    F = x + (2 - 4*x)/(3*x**2 - 3*x + 3) + 3*log(x**2 - x + 1)/2 - 7*sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_159():
    f = x*(x**2 + x + 1)/(x**2 - x + 1)**2
    F = -(2*x + 2)/(3*x**2 - 3*x + 3) + log(x**2 - x + 1)/2 - 11*sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_160():
    f = (x**2 + x + 1)/(x**2 - x + 1)**2
    F = -(4 - 2*x)/(3*x**2 - 3*x + 3) - 10*sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_161():
    f = (x**2 + x + 1)/(x*(x**2 - x + 1)**2)
    F = -(2 - 4*x)/(3*x**2 - 3*x + 3) + log(x) - log(x**2 - x + 1)/2 - 11*sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_162():
    f = (x**2 + x + 1)/(x**2*(x**2 - x + 1)**2)
    F = (2*x + 2)/(3*x**2 - 3*x + 3) + 3*log(x) - 3*log(x**2 - x + 1)/2 - 7*sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/9 - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_163():
    f = (x**2 + x + 1)/(x**3*(x**2 - x + 1)**2)
    F = (4 - 2*x)/(3*x**2 - 3*x + 3) + 4*log(x) - 2*log(x**2 - x + 1) + 10*sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/9 - 3/x - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_164():
    f = (1 - x**2)/(x**2 + x + 1)**2
    F = x/(x**2 + x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_165():
    f = (x**2 + 1)/(x**2 + x + 1)
    F = x - log(x**2 + x + 1)/2 + sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_166():
    f = (x**2 - 1)/(x**2 - 6*x + 25)
    F = x + 3*log(x**2 - 6*x + 25) - 2*atan(x/4 + sympy.S(-3)/4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_167():
    f = (3*x**2 - 10)/(x**2 - 4*x + 4)
    F = 3*x + 12*log(2 - x) + 2/(2 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_168():
    f = (x**2 + 8)/(x**2 - 5*x + 6)
    F = x - 12*log(2 - x) + 17*log(3 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_169():
    f = (x**2 + 3*x - 4)/(x**2 - 2*x - 8)
    F = x + 4*log(4 - x) + log(x + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_170():
    f = (4*x**2 + 5*x + 7)/(4*x**2 + 4*x + 5)
    F = x + log(4*x**2 + 4*x + 5)/8 + 3*atan(x + sympy.S.Half)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_171():
    f = (x**2 - x + 2)/(x**2 + 2*x - 5)
    F = x - (sympy.S(3)/2 + 5*sqrt(6)/6)*log(x + 1 + sqrt(6)) - (sympy.S(3)/2 - 5*sqrt(6)/6)*log(x - sqrt(6) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_172():
    f = (3*x**2 + 4*x + 1)/(2*x**2 + 7*x + 4)**2
    F = -(3*x + 2)/(4*x**2 + 14*x + 8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_173():
    f = (x**2 + x + 1)/(x**2 + 2*x + 3)**2
    F = (1 - x)/(4*x**2 + 8*x + 12) + 3*sqrt(2)*atan(sqrt(2)*(x + 1)/2)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_174():
    f = (5*x**2 + 2*x - 1)/(x**2 + x + 1)**4
    F = -x/(x**2 + x + 1)**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_175():
    f = (A + C*x**2)*(a + b*x + c*x**2)**(sympy.S(5)/2)
    F = -9*C*b*(a + b*x + c*x**2)**(sympy.S(7)/2)/(112*c**2) + C*x*(a + b*x + c*x**2)**(sympy.S(7)/2)/(8*c) + (b + 2*c*x)*(a + b*x + c*x**2)**(sympy.S(5)/2)*(32*A*c**2 - 4*C*a*c + 9*C*b**2)/(384*c**3) - (b + 2*c*x)*(-20*a*c + 5*b**2)*(a + b*x + c*x**2)**(sympy.S(3)/2)*(32*A*c**2 - 4*C*a*c + 9*C*b**2)/(6144*c**4) + 5*(b + 2*c*x)*(-4*a*c + b**2)**2*sqrt(a + b*x + c*x**2)*(32*A*c**2 - 4*C*a*c + 9*C*b**2)/(16384*c**5) - 5*(-4*a*c + b**2)**3*(32*A*c**2 - 4*C*a*c + 9*C*b**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(32768*c**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_176():
    f = (A + C*x**2)*(a + b*x + c*x**2)**(sympy.S(3)/2)
    F = -7*C*b*(a + b*x + c*x**2)**(sympy.S(5)/2)/(60*c**2) + C*x*(a + b*x + c*x**2)**(sympy.S(5)/2)/(6*c) + (b + 2*c*x)*(a + b*x + c*x**2)**(sympy.S(3)/2)*(24*A*c**2 - 4*C*a*c + 7*C*b**2)/(192*c**3) - (b + 2*c*x)*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)*(24*A*c**2 - 4*C*a*c + 7*C*b**2)/(512*c**4) + (-4*a*c + b**2)**2*(24*A*c**2 - 4*C*a*c + 7*C*b**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(1024*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_177():
    f = (A + C*x**2)*sqrt(a + b*x + c*x**2)
    F = -5*C*b*(a + b*x + c*x**2)**(sympy.S(3)/2)/(24*c**2) + C*x*(a + b*x + c*x**2)**(sympy.S(3)/2)/(4*c) + (b + 2*c*x)*sqrt(a + b*x + c*x**2)*(16*A*c**2 - 4*C*a*c + 5*C*b**2)/(64*c**3) - (-4*a*c + b**2)*(16*A*c**2 - 4*C*a*c + 5*C*b**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(128*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_178():
    f = (A + C*x**2)/sqrt(a + b*x + c*x**2)
    F = -3*C*b*sqrt(a + b*x + c*x**2)/(4*c**2) + C*x*sqrt(a + b*x + c*x**2)/(2*c) + (8*A*c**2 - 4*C*a*c + 3*C*b**2)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_179():
    f = (A + C*x**2)/(a + b*x + c*x**2)**(sympy.S(3)/2)
    F = C*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/c**(sympy.S(3)/2) - (2*b*c*(A + C*a/c) + 2*x*(2*A*c**2 + C*(-2*a*c + b**2)))/(c*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_180():
    f = (A + C*x**2)/(a + b*x + c*x**2)**(sympy.S(5)/2)
    F = (b + 2*c*x)*(16*A*c + 8*C*a + 2*C*b**2/c)/(3*(-4*a*c + b**2)**2*sqrt(a + b*x + c*x**2)) - (2*b*c*(A + C*a/c) + 2*x*(2*A*c**2 + C*(-2*a*c + b**2)))/(3*c*(-4*a*c + b**2)*(a + b*x + c*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_181():
    f = (A + C*x**2)/(a + b*x + c*x**2)**(sympy.S(7)/2)
    F = (b + 2*c*x)*(32*A*c + 8*C*a + 6*C*b**2/c)/(15*(-4*a*c + b**2)**2*(a + b*x + c*x**2)**(sympy.S(3)/2)) - (b + 2*c*x)*(256*A*c**2 + 64*C*a*c + 48*C*b**2)/(15*(-4*a*c + b**2)**3*sqrt(a + b*x + c*x**2)) - (2*b*c*(A + C*a/c) + 2*x*(2*A*c**2 + C*(-2*a*c + b**2)))/(5*c*(-4*a*c + b**2)*(a + b*x + c*x**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_182():
    f = (A + C*x**2)/(a + b*x + c*x**2)**(sympy.S(9)/2)
    F = 256*c*(b + 2*c*x)*(24*A*c**2 + 4*C*a*c + 5*C*b**2)/(105*(-4*a*c + b**2)**4*sqrt(a + b*x + c*x**2)) + (b + 2*c*x)*(48*A*c + 8*C*a + 10*C*b**2/c)/(35*(-4*a*c + b**2)**2*(a + b*x + c*x**2)**(sympy.S(5)/2)) - (b + 2*c*x)*(768*A*c**2 + 128*C*a*c + 160*C*b**2)/(105*(-4*a*c + b**2)**3*(a + b*x + c*x**2)**(sympy.S(3)/2)) - (2*b*c*(A + C*a/c) + 2*x*(2*A*c**2 + C*(-2*a*c + b**2)))/(7*c*(-4*a*c + b**2)*(a + b*x + c*x**2)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_183():
    f = (g + h*x)**3*sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2)
    F = f*(g + h*x)**4*(a + b*x + c*x**2)**(sympy.S(3)/2)/(7*c*h) - (g + h*x)**3*(a + b*x + c*x**2)**(sympy.S(3)/2)*(11*b*f*h - 14*c*e*h + 6*c*f*g)/(84*c**2*h) + (g + h*x)**2*(a + b*x + c*x**2)**(sympy.S(3)/2)*(33*b**2*f*h**2 - 4*c**2*(3*f*g**2 - 7*h*(2*d*h + e*g)) - 2*c*h*(16*a*f*h + 21*b*e*h + 8*b*f*g))/(280*c**3*h) + (a + b*x + c*x**2)**(sympy.S(3)/2)*(1155*b**4*f*h**4 - 42*b**2*c*h**3*(78*a*f*h + 35*b*(e*h + 3*f*g)) - 128*c**4*g**2*(3*f*g**2 - 7*h*(12*d*h + e*g)) - 16*c**3*h*(16*a*h*(15*f*g**2 + 7*h*(d*h + 3*e*g)) + b*g*(17*f*g**2 + 21*h*(25*d*h + 19*e*g))) + 8*c**2*h**2*(128*a**2*f*h**2 + 343*a*b*h*(e*h + 3*f*g) + b**2*(537*f*g**2 + 245*h*(d*h + 3*e*g))) - 6*c*h*x*(231*b**3*f*h**3 - 6*b*c*h**2*(74*a*f*h + 49*b*e*h + 59*b*f*g) + 16*c**3*g*(3*f*g**2 - 7*h*(7*d*h + e*g)) + 8*c**2*h*(a*h*(35*e*h + 41*f*g) + b*(5*f*g**2 + 7*h*(7*d*h + 9*e*g)))))/(13440*c**5*h) + (b + 2*c*x)*sqrt(a + b*x + c*x**2)*(-33*b**5*f*h**3 + 6*b**3*c*h**2*(20*a*f*h + 7*b*(e*h + 3*f*g)) - 8*b*c**2*h*(10*a**2*f*h**2 + 14*a*b*h*(e*h + 3*f*g) + 7*b**2*(d*h**2 + 3*e*g*h + 3*f*g**2)) + 256*c**5*d*g**3 - 64*c**4*g*(a*(f*g**2 + 3*h*(d*h + e*g)) + 2*b*g*(3*d*h + e*g)) + 16*c**3*(2*a**2*h**2*(e*h + 3*f*g) + 6*a*b*h*(3*f*g**2 + h*(d*h + 3*e*g)) + 5*b**2*g*(f*g**2 + 3*h*(d*h + e*g))))/(1024*c**6) - (-4*a*c + b**2)*(-33*b**5*f*h**3 + 6*b**3*c*h**2*(20*a*f*h + 7*b*(e*h + 3*f*g)) - 8*b*c**2*h*(10*a**2*f*h**2 + 14*a*b*h*(e*h + 3*f*g) + 7*b**2*(d*h**2 + 3*e*g*h + 3*f*g**2)) + 256*c**5*d*g**3 - 64*c**4*g*(a*(f*g**2 + 3*h*(d*h + e*g)) + 2*b*g*(3*d*h + e*g)) + 16*c**3*(2*a**2*h**2*(e*h + 3*f*g) + 6*a*b*h*(3*f*g**2 + h*(d*h + 3*e*g)) + 5*b**2*g*(f*g**2 + 3*h*(d*h + e*g))))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2048*c**(sympy.S(13)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_184():
    f = (g + h*x)**2*sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2)
    F = f*(g + h*x)**3*(a + b*x + c*x**2)**(sympy.S(3)/2)/(6*c*h) - (g + h*x)**2*(a + b*x + c*x**2)**(sympy.S(3)/2)*(3*b*f*h - 4*c*e*h + 2*c*f*g)/(20*c**2*h) - (a + b*x + c*x**2)**(sympy.S(3)/2)*(105*b**3*f*h**3 - 28*b*c*h**2*(7*a*f*h + 5*b*(e*h + 2*f*g)) + 64*c**3*g*(f*g**2 - 2*h*(5*d*h + e*g)) + 8*c**2*h*(16*a*h*(e*h + 2*f*g) + b*(7*f*g**2 + 25*h*(d*h + 2*e*g))) - 6*c*h*x*(21*b**2*f*h**2 - 8*c**2*(f*g**2 - h*(5*d*h + 2*e*g)) - 4*c*h*(5*a*f*h + 7*b*e*h + 2*b*f*g)))/(960*c**4*h) + (b + 2*c*x)*sqrt(a + b*x + c*x**2)*(21*b**4*f*h**2 - 28*b**2*c*h*(2*a*f*h + b*e*h + 2*b*f*g) + 128*c**4*d*g**2 - 32*c**3*(a*(d*h**2 + 2*e*g*h + f*g**2) + 2*b*g*(2*d*h + e*g)) + 8*c**2*(2*a**2*f*h**2 + 6*a*b*h*(e*h + 2*f*g) + 5*b**2*(d*h**2 + 2*e*g*h + f*g**2)))/(512*c**5) - (-4*a*c + b**2)*(21*b**4*f*h**2 - 28*b**2*c*h*(2*a*f*h + b*e*h + 2*b*f*g) + 128*c**4*d*g**2 - 32*c**3*(a*(d*h**2 + 2*e*g*h + f*g**2) + 2*b*g*(2*d*h + e*g)) + 8*c**2*(2*a**2*f*h**2 + 6*a*b*h*(e*h + 2*f*g) + 5*b**2*(d*h**2 + 2*e*g*h + f*g**2)))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(1024*c**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_185():
    f = (g + h*x)*sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2)
    F = f*(g + h*x)**2*(a + b*x + c*x**2)**(sympy.S(3)/2)/(5*c*h) + (a + b*x + c*x**2)**(sympy.S(3)/2)*(35*b**2*f*h**2 - 16*c**2*(3*f*g**2 - 5*h*(d*h + e*g)) - 6*c*h*x*(7*b*f*h - 10*c*e*h + 6*c*f*g) - 2*c*h*(16*a*f*h + 25*b*(e*h + f*g)))/(240*c**3*h) + (b + 2*c*x)*sqrt(a + b*x + c*x**2)*(-7*b**3*f*h + 2*b*c*(6*a*f*h + 5*b*(e*h + f*g)) + 32*c**3*d*g - 8*c**2*(a*e*h + a*f*g + 2*b*d*h + 2*b*e*g))/(128*c**4) - (-4*a*c + b**2)*(-7*b**3*f*h + 2*b*c*(6*a*f*h + 5*b*(e*h + f*g)) + 32*c**3*d*g - 8*c**2*(a*e*h + a*f*g + 2*b*d*h + 2*b*e*g))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(256*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_186():
    f = sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2)
    F = f*x*(a + b*x + c*x**2)**(sympy.S(3)/2)/(4*c) + (-5*b*f + 8*c*e)*(a + b*x + c*x**2)**(sympy.S(3)/2)/(24*c**2) + (b + 2*c*x)*sqrt(a + b*x + c*x**2)*(-4*a*c*f + 5*b**2*f - 8*b*c*e + 16*c**2*d)/(64*c**3) - (-4*a*c + b**2)*(5*b**2*f + 16*c**2*d - 4*c*(a*f + 2*b*e))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(128*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_187():
    f = sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2)/(g + h*x)
    F = sqrt(a*h**2 - b*g*h + c*g**2)*(d*h**2 - e*g*h + f*g**2)*atanh((-2*a*h + b*g + x*(-b*h + 2*c*g))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*h**2 - b*g*h + c*g**2)))/h**4 + f*(a + b*x + c*x**2)**(sympy.S(3)/2)/(3*c*h) - sqrt(a + b*x + c*x**2)*(2*c*h*x*(b*f*h - 2*c*e*h + 2*c*f*g) + 4*c*h*(b*f*g - 2*c*d*h) - (-b*h + 4*c*g)*(b*f*h - 2*c*e*h + 2*c*f*g))/(8*c**2*h**3) + (4*c*h*(-b*h + 2*c*g)*(b*f*g - 2*c*d*h) - (-b**2*h**2 + 8*c**2*g**2 - 4*c*h*(-a*h + b*g))*(b*f*h - 2*c*e*h + 2*c*f*g))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(5)/2)*h**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_188():
    f = sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2)/(g + h*x)**2
    F = -(f*g**2 - h*(-d*h + e*g))*(a + b*x + c*x**2)**(sympy.S(3)/2)/(h*(g + h*x)*(a*h**2 - b*g*h + c*g**2)) - (2*c*g*(3*f*g**2 - h*(-d*h + 2*e*g)) + h*(2*a*h*(-e*h + 2*f*g) - b*(d*h**2 - 3*e*g*h + 5*f*g**2)))*atanh((-2*a*h + b*g + x*(-b*h + 2*c*g))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*h**2 - b*g*h + c*g**2)))/(2*h**4*sqrt(a*h**2 - b*g*h + c*g**2)) - sqrt(a + b*x + c*x**2)*(b*f*h**2*(-a*h + b*g) + 4*c**2*g*(3*f*g**2 - h*(-d*h + 2*e*g)) + 2*c*h**2*x*(-a*f*h + b*f*g - 2*c*d*h + 2*c*e*g - 3*c*f*g**2/h) + c*h*(4*a*h*(-e*h + 2*f*g) - b*(4*d*h**2 - 8*e*g*h + 13*f*g**2)))/(4*c*h**3*(a*h**2 - b*g*h + c*g**2)) - (b**2*f*h**2 - 8*c**2*(3*f*g**2 - h*(-d*h + 2*e*g)) + 4*c*h*(-a*f*h - b*e*h + 2*b*f*g))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*c**(sympy.S(3)/2)*h**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_189():
    f = sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2)/(g + h*x)**3
    F = -(f*g**2 - h*(-d*h + e*g))*(a + b*x + c*x**2)**(sympy.S(3)/2)/(2*h*(g + h*x)**2*(a*h**2 - b*g*h + c*g**2)) + sqrt(a + b*x + c*x**2)*(4*a*h*(-e*h + 3*f*g) - b*(-d*h**2 - 3*e*g*h + 11*f*g**2) + 4*c*g**2*(-e*h + 3*f*g)/h - 2*h*x*(-2*a*f*h + 2*b*f*g - c*d*h + c*e*g - 3*c*f*g**2/h))/(4*h**2*(g + h*x)*(a*h**2 - b*g*h + c*g**2)) + (8*c**2*g**3*(-e*h + 3*f*g) - 4*c*h*(-a*h*(d*h**2 - 3*e*g*h + 9*f*g**2) + b*g**2*(-3*e*h + 10*f*g)) + h**2*(8*a**2*f*h**2 - 4*a*b*h*(-e*h + 6*f*g) + b**2*(15*f*g**2 - h*(d*h + 3*e*g))))*atanh((-2*a*h + b*g + x*(-b*h + 2*c*g))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*h**2 - b*g*h + c*g**2)))/(8*h**4*(a*h**2 - b*g*h + c*g**2)**(sympy.S(3)/2)) - (-b*f*h - 2*c*e*h + 6*c*f*g)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*sqrt(c)*h**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_190():
    f = sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2)/(g + h*x)**4
    F = sqrt(c)*f*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/h**4 - (f*g**2 - h*(-d*h + e*g))*(a + b*x + c*x**2)**(sympy.S(3)/2)/(3*h*(g + h*x)**3*(a*h**2 - b*g*h + c*g**2)) - sqrt(a + b*x + c*x**2)*(8*c**2*f*g**5 - 2*c*g*h*(-2*a*d*h**3 - 6*a*f*g**2*h + b*d*g*h**2 + 7*b*f*g**3) + h**2*(4*a**2*e*h**3 - 2*a*b*h*(d*h**2 + 2*e*g*h + 3*f*g**2) + b**2*g*(d*h**2 + e*g*h + 5*f*g**2)) + h*x*(4*c**2*(-d*g**2*h**2 + 3*f*g**4) + 2*c*g*h*(2*a*h*(-e*h + 6*f*g) - b*(12*f*g**2 - h*(2*d*h + e*g))) + h**2*(8*a**2*f*h**2 - 2*a*b*h*(-e*h + 10*f*g) + b**2*(11*f*g**2 - h*(d*h + e*g)))))/(8*h**3*(g + h*x)**2*(a*h**2 - b*g*h + c*g**2)**2) - (-b*h**3*(8*a**2*f*h**2 - 2*a*b*h*(e*h + 6*f*g) + b**2*(d*h**2 + e*g*h + 5*f*g**2)) + 16*c**3*f*g**5 - 8*c**2*g*h*(a*d*h**3 - 5*a*f*g**2*h + 5*b*f*g**3) + 2*c*h**2*(4*a**2*h**2*(-e*h + 4*f*g) - 2*a*b*h*(-d*h**2 - e*g*h + 15*f*g**2) + b**2*(d*g*h**2 + 15*f*g**3)))*atanh((-2*a*h + b*g + x*(-b*h + 2*c*g))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*h**2 - b*g*h + c*g**2)))/(16*h**4*(a*h**2 - b*g*h + c*g**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_191():
    f = sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2)/(g + h*x)**5
    F = -(-4*a*c + b**2)*(16*a**2*f*h**2 - 8*a*b*h*(e*h + 2*f*g) + b**2*(5*d*h**2 + 3*e*g*h + 5*f*g**2) + 16*c**2*d*g**2 - 4*c*(a*(d*h**2 - 5*e*g*h + f*g**2) + 2*b*g*(2*d*h + e*g)))*atanh((-2*a*h + b*g + x*(-b*h + 2*c*g))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*h**2 - b*g*h + c*g**2)))/(128*(a*h**2 - b*g*h + c*g**2)**(sympy.S(7)/2)) + sqrt(a + b*x + c*x**2)*(-2*a*h + b*g + x*(-b*h + 2*c*g))*(16*a**2*f*h**2 - 8*a*b*h*(e*h + 2*f*g) + b**2*(5*d*h**2 + 3*e*g*h + 5*f*g**2) + 16*c**2*d*g**2 - 4*c*(a*(d*h**2 - 5*e*g*h + f*g**2) + 2*b*g*(2*d*h + e*g)))/(64*(g + h*x)**2*(a*h**2 - b*g*h + c*g**2)**3) + (2*c*g*(3*f*g**2 + h*(-5*d*h + e*g)) + h*(8*a*h*(-e*h + 2*f*g) - b*(-5*d*h**2 - 3*e*g*h + 11*f*g**2)))*(a + b*x + c*x**2)**(sympy.S(3)/2)/(24*h*(g + h*x)**3*(a*h**2 - b*g*h + c*g**2)**2) - (f*g**2 - h*(-d*h + e*g))*(a + b*x + c*x**2)**(sympy.S(3)/2)/(4*h*(g + h*x)**4*(a*h**2 - b*g*h + c*g**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_192():
    f = sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2)/(g + h*x)**6
    F = -(-4*a*c + b**2)*(-b*h*(16*a**2*f*h**2 - 2*a*b*h*(5*e*h + 6*f*g) + b**2*(7*d*h**2 + 3*e*g*h + 3*f*g**2)) + 32*c**3*d*g**3 - 8*c**2*g*(a*(3*d*h**2 - 6*e*g*h + f*g**2) + 2*b*g*(3*d*h + e*g)) + 2*c*(4*a**2*h**2*(-e*h + 6*f*g) - 6*a*b*h*(-d*h**2 + 3*e*g*h + 3*f*g**2) + b**2*g*(15*d*h**2 + 6*e*g*h + 5*f*g**2)))*atanh((-2*a*h + b*g + x*(-b*h + 2*c*g))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*h**2 - b*g*h + c*g**2)))/(256*(a*h**2 - b*g*h + c*g**2)**(sympy.S(9)/2)) + sqrt(a + b*x + c*x**2)*(-2*a*h + b*g + x*(-b*h + 2*c*g))*(-b*h*(16*a**2*f*h**2 - 2*a*b*h*(5*e*h + 6*f*g) + b**2*(7*d*h**2 + 3*e*g*h + 3*f*g**2)) + 32*c**3*d*g**3 - 8*c**2*g*(a*(3*d*h**2 - 6*e*g*h + f*g**2) + 2*b*g*(3*d*h + e*g)) + 2*c*(4*a**2*h**2*(-e*h + 6*f*g) - 6*a*b*h*(-d*h**2 + 3*e*g*h + 3*f*g**2) + b**2*g*(15*d*h**2 + 6*e*g*h + 5*f*g**2)))/(128*(g + h*x)**2*(a*h**2 - b*g*h + c*g**2)**4) + (a + b*x + c*x**2)**(sympy.S(3)/2)*(4*c**2*g**2*(3*f*g**2 + h*(-27*d*h + 2*e*g)) - 2*c*h*(-2*a*h*(8*d*h**2 - 33*e*g*h + 18*f*g**2) + b*g*(-54*d*h**2 - 21*e*g*h + 16*f*g**2)) - 5*h**2*(16*a**2*f*h**2 - 2*a*b*h*(5*e*h + 6*f*g) + b**2*(7*d*h**2 + 3*e*g*h + 3*f*g**2)))/(240*h*(g + h*x)**3*(a*h**2 - b*g*h + c*g**2)**3) + (2*c*g*(3*f*g**2 + h*(-7*d*h + 2*e*g)) + h*(10*a*h*(-e*h + 2*f*g) - b*(-7*d*h**2 - 3*e*g*h + 13*f*g**2)))*(a + b*x + c*x**2)**(sympy.S(3)/2)/(40*h*(g + h*x)**4*(a*h**2 - b*g*h + c*g**2)**2) - (f*g**2 - h*(-d*h + e*g))*(a + b*x + c*x**2)**(sympy.S(3)/2)/(5*h*(g + h*x)**5*(a*h**2 - b*g*h + c*g**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_193():
    f = (g + h*x)**3*(a + b*x + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)
    F = f*(g + h*x)**4*(a + b*x + c*x**2)**(sympy.S(5)/2)/(9*c*h) - (g + h*x)**3*(a + b*x + c*x**2)**(sympy.S(5)/2)*(13*b*f*h - 18*c*e*h + 10*c*f*g)/(144*c**2*h) + (g + h*x)**2*(a + b*x + c*x**2)**(sympy.S(5)/2)*(143*b**2*f*h**2 - 12*c**2*(5*f*g**2 - 3*h*(8*d*h + 3*e*g)) - 2*c*h*(64*a*f*h + 99*b*e*h + 24*b*f*g))/(2016*c**3*h) + (a + b*x + c*x**2)**(sympy.S(5)/2)*(3003*b**4*f*h**4 - 198*b**2*c*h**3*(38*a*f*h + 21*b*(e*h + 3*f*g)) - 192*c**4*g**2*(5*f*g**2 - 3*h*(64*d*h + 3*e*g)) - 16*c**3*h*(32*a*h*(17*f*g**2 + 9*h*(d*h + 3*e*g)) + b*g*(13*f*g**2 + 9*h*(196*d*h + 141*e*g))) + 8*c**2*h**2*(256*a**2*f*h**2 + 837*a*b*h*(e*h + 3*f*g) + b**2*(1553*f*g**2 + 756*h*(d*h + 3*e*g))) - 10*c*h*x*(429*b**3*f*h**3 - 22*b*c*h**2*(34*a*f*h + 27*b*e*h + 29*b*f*g) + 16*c**3*g*(5*f*g**2 - 9*h*(12*d*h + e*g)) + 8*c**2*h*(a*h*(63*e*h + 61*f*g) + 3*b*(f*g**2 + 6*h*(6*d*h + 7*e*g)))))/(80640*c**5*h) + (b + 2*c*x)*(a + b*x + c*x**2)**(sympy.S(3)/2)*(-143*b**5*f*h**3 + 22*b**3*c*h**2*(20*a*f*h + 9*b*(e*h + 3*f*g)) - 48*b*c**2*h*(5*a**2*f*h**2 + 9*a*b*h*(e*h + 3*f*g) + 6*b**2*(d*h**2 + 3*e*g*h + 3*f*g**2)) + 1536*c**5*d*g**3 - 256*c**4*g*(a*(f*g**2 + 3*h*(d*h + e*g)) + 3*b*g*(3*d*h + e*g)) + 32*c**3*(3*a**2*h**2*(e*h + 3*f*g) + 12*a*b*h*(3*f*g**2 + h*(d*h + 3*e*g)) + 14*b**2*g*(f*g**2 + 3*h*(d*h + e*g))))/(12288*c**6) - (b + 2*c*x)*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)*(-143*b**5*f*h**3 + 22*b**3*c*h**2*(20*a*f*h + 9*b*(e*h + 3*f*g)) - 48*b*c**2*h*(5*a**2*f*h**2 + 9*a*b*h*(e*h + 3*f*g) + 6*b**2*(d*h**2 + 3*e*g*h + 3*f*g**2)) + 1536*c**5*d*g**3 - 256*c**4*g*(a*(f*g**2 + 3*h*(d*h + e*g)) + 3*b*g*(3*d*h + e*g)) + 32*c**3*(3*a**2*h**2*(e*h + 3*f*g) + 12*a*b*h*(3*f*g**2 + h*(d*h + 3*e*g)) + 14*b**2*g*(f*g**2 + 3*h*(d*h + e*g))))/(32768*c**7) + (-4*a*c + b**2)**2*(-143*b**5*f*h**3 + 22*b**3*c*h**2*(20*a*f*h + 9*b*(e*h + 3*f*g)) - 48*b*c**2*h*(5*a**2*f*h**2 + 9*a*b*h*(e*h + 3*f*g) + 6*b**2*(d*h**2 + 3*e*g*h + 3*f*g**2)) + 1536*c**5*d*g**3 - 256*c**4*g*(a*(f*g**2 + 3*h*(d*h + e*g)) + 3*b*g*(3*d*h + e*g)) + 32*c**3*(3*a**2*h**2*(e*h + 3*f*g) + 12*a*b*h*(3*f*g**2 + h*(d*h + 3*e*g)) + 14*b**2*g*(f*g**2 + 3*h*(d*h + e*g))))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(65536*c**(sympy.S(15)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_194():
    f = (g + h*x)**2*(a + b*x + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)
    F = f*(g + h*x)**3*(a + b*x + c*x**2)**(sympy.S(5)/2)/(8*c*h) - (g + h*x)**2*(a + b*x + c*x**2)**(sympy.S(5)/2)*(11*b*f*h - 16*c*e*h + 10*c*f*g)/(112*c**2*h) - (a + b*x + c*x**2)**(sympy.S(5)/2)*(693*b**3*f*h**3 - 36*b*c*h**2*(31*a*f*h + 28*b*(e*h + 2*f*g)) + 96*c**3*g*(5*f*g**2 - 8*h*(7*d*h + e*g)) + 8*c**2*h*(96*a*h*(e*h + 2*f*g) + b*(31*f*g**2 + 196*h*(d*h + 2*e*g))) - 10*c*h*x*(99*b**2*f*h**2 - 8*c**2*(5*f*g**2 - 4*h*(7*d*h + 2*e*g)) - 12*c*h*(7*a*f*h + 2*b*(6*e*h + f*g))))/(13440*c**4*h) + (b + 2*c*x)*(a + b*x + c*x**2)**(sympy.S(3)/2)*(99*b**4*f*h**2 - 72*b**2*c*h*(3*a*f*h + 2*b*e*h + 4*b*f*g) + 768*c**4*d*g**2 - 128*c**3*(a*(d*h**2 + 2*e*g*h + f*g**2) + 3*b*g*(2*d*h + e*g)) + 16*c**2*(3*a**2*f*h**2 + 12*a*b*h*(e*h + 2*f*g) + 14*b**2*(d*h**2 + 2*e*g*h + f*g**2)))/(6144*c**5) - (b + 2*c*x)*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)*(99*b**4*f*h**2 - 72*b**2*c*h*(3*a*f*h + 2*b*e*h + 4*b*f*g) + 768*c**4*d*g**2 - 128*c**3*(a*(d*h**2 + 2*e*g*h + f*g**2) + 3*b*g*(2*d*h + e*g)) + 16*c**2*(3*a**2*f*h**2 + 12*a*b*h*(e*h + 2*f*g) + 14*b**2*(d*h**2 + 2*e*g*h + f*g**2)))/(16384*c**6) + (-4*a*c + b**2)**2*(99*b**4*f*h**2 - 72*b**2*c*h*(3*a*f*h + 2*b*e*h + 4*b*f*g) + 768*c**4*d*g**2 - 128*c**3*(a*(d*h**2 + 2*e*g*h + f*g**2) + 3*b*g*(2*d*h + e*g)) + 16*c**2*(3*a**2*f*h**2 + 12*a*b*h*(e*h + 2*f*g) + 14*b**2*(d*h**2 + 2*e*g*h + f*g**2)))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(32768*c**(sympy.S(13)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_195():
    f = (g + h*x)*(a + b*x + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)
    F = f*(g + h*x)**2*(a + b*x + c*x**2)**(sympy.S(5)/2)/(7*c*h) + (a + b*x + c*x**2)**(sympy.S(5)/2)*(63*b**2*f*h**2 - 24*c**2*(5*f*g**2 - 7*h*(d*h + e*g)) - 10*c*h*x*(9*b*f*h - 14*c*e*h + 10*c*f*g) - 2*c*h*(24*a*f*h + 49*b*(e*h + f*g)))/(840*c**3*h) + (b + 2*c*x)*(a + b*x + c*x**2)**(sympy.S(3)/2)*(-9*b**3*f*h + 2*b*c*(6*a*f*h + 7*b*(e*h + f*g)) + 48*c**3*d*g - 8*c**2*(a*e*h + a*f*g + 3*b*d*h + 3*b*e*g))/(384*c**4) - (b + 2*c*x)*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)*(-9*b**3*f*h + 2*b*c*(6*a*f*h + 7*b*(e*h + f*g)) + 48*c**3*d*g - 8*c**2*(a*e*h + a*f*g + 3*b*d*h + 3*b*e*g))/(1024*c**5) + (-4*a*c + b**2)**2*(-9*b**3*f*h + 2*b*c*(6*a*f*h + 7*b*(e*h + f*g)) + 48*c**3*d*g - 8*c**2*(a*e*h + a*f*g + 3*b*d*h + 3*b*e*g))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2048*c**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_196():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)
    F = f*x*(a + b*x + c*x**2)**(sympy.S(5)/2)/(6*c) + (-7*b*f + 12*c*e)*(a + b*x + c*x**2)**(sympy.S(5)/2)/(60*c**2) + (b + 2*c*x)*(a + b*x + c*x**2)**(sympy.S(3)/2)*(-4*a*c*f + 7*b**2*f - 12*b*c*e + 24*c**2*d)/(192*c**3) - (b + 2*c*x)*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)*(7*b**2*f + 24*c**2*d - 4*c*(a*f + 3*b*e))/(512*c**4) + (-4*a*c + b**2)**2*(7*b**2*f + 24*c**2*d - 4*c*(a*f + 3*b*e))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(1024*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_197():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)/(g + h*x)
    F = (f*g**2 - h*(-d*h + e*g))*(a*h**2 - b*g*h + c*g**2)**(sympy.S(3)/2)*atanh((-2*a*h + b*g + x*(-b*h + 2*c*g))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*h**2 - b*g*h + c*g**2)))/h**6 + f*(a + b*x + c*x**2)**(sympy.S(5)/2)/(5*c*h) - (a + b*x + c*x**2)**(sympy.S(3)/2)*(6*c*h*x*(b*f*h - 2*c*e*h + 2*c*f*g) + 8*c*h*(b*f*g - 2*c*d*h) - (-3*b*h + 8*c*g)*(b*f*h - 2*c*e*h + 2*c*f*g))/(48*c**2*h**3) + sqrt(a + b*x + c*x**2)*(3*b**4*f*h**4 + 6*b**2*c*h**3*(-2*a*f*h - b*e*h + b*f*g) - 8*b*c**2*h**2*(3*a*h*(-e*h + f*g) - 2*b*(d*h**2 - e*g*h + f*g**2)) + 128*c**4*g**2*(f*g**2 - h*(-d*h + e*g)) - 32*c**3*h*(-4*a*h + 5*b*g)*(f*g**2 - h*(-d*h + e*g)) + 2*c*h*x*(8*c*h*(-b*h + 2*c*g)*(b*f*g - 2*c*d*h) - (-3*b**2*h**2 + 16*c**2*g**2 - 4*c*h*(-3*a*h + 2*b*g))*(b*f*h - 2*c*e*h + 2*c*f*g)))/(128*c**3*h**5) - (4*c*h*(-b*h + 2*c*g)*(8*c*h*(-2*a*h + b*g)*(b*f*g - 2*c*d*h) - g*(-4*a*c*h - 3*b**2*h + 8*b*c*g)*(b*f*h - 2*c*e*h + 2*c*f*g)) - (8*c*h*(-b*h + 2*c*g)*(b*f*g - 2*c*d*h) - (-3*b**2*h**2 + 16*c**2*g**2 - 4*c*h*(-3*a*h + 2*b*g))*(b*f*h - 2*c*e*h + 2*c*f*g))*(-b**2*h**2 + 8*c**2*g**2 - 4*c*h*(-a*h + b*g)))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(256*c**(sympy.S(7)/2)*h**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_198():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)/(g + h*x)**2
    F = -(f*g**2 - h*(-d*h + e*g))*(a + b*x + c*x**2)**(sympy.S(5)/2)/(h*(g + h*x)*(a*h**2 - b*g*h + c*g**2)) - (2*c*g*(5*f*g**2 - h*(-3*d*h + 4*e*g)) + h*(2*a*h*(-e*h + 2*f*g) - b*(3*d*h**2 - 5*e*g*h + 7*f*g**2)))*sqrt(a*h**2 - b*g*h + c*g**2)*atanh((-2*a*h + b*g + x*(-b*h + 2*c*g))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*h**2 - b*g*h + c*g**2)))/(2*h**6) - (a + b*x + c*x**2)**(sympy.S(3)/2)*(3*b*f*h**2*(-a*h + b*g) + 8*c**2*g*(5*f*g**2 - h*(-3*d*h + 4*e*g)) + 6*c*h**2*x*(-a*f*h + b*f*g - 4*c*d*h + 4*c*e*g - 5*c*f*g**2/h) + c*h*(8*a*h*(-e*h + 2*f*g) - b*(43*f*g**2 - 8*h*(-3*d*h + 4*e*g))))/(24*c*h**3*(a*h**2 - b*g*h + c*g**2)) - sqrt(a + b*x + c*x**2)*(3*b**3*f*h**3 + 4*b*c*h**2*(-3*a*f*h - 2*b*e*h + 4*b*f*g) + 64*c**3*g*(5*f*g**2 - h*(-3*d*h + 4*e*g)) + 16*c**2*h*(4*a*h*(-e*h + 2*f*g) - b*(9*d*h**2 - 14*e*g*h + 19*f*g**2)) + 2*c*h*x*(3*b**2*f*h**2 - 16*c**2*(5*f*g**2 - h*(-3*d*h + 4*e*g)) + 4*c*h*(-3*a*f*h - 2*b*e*h + 4*b*f*g)))/(64*c**2*h**5) + (3*b**4*f*h**4 + 8*b**2*c*h**3*(-3*a*f*h - b*e*h + 2*b*f*g) + 128*c**4*g**2*(5*f*g**2 - h*(-3*d*h + 4*e*g)) + 192*c**3*h*(a*h*(d*h**2 - 2*e*g*h + 3*f*g**2) - b*g*(2*d*h**2 - 3*e*g*h + 4*f*g**2)) + 48*c**2*h**2*(a**2*f*h**2 - 2*a*b*h*(-e*h + 2*f*g) + b**2*(d*h**2 - 2*e*g*h + 3*f*g**2)))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(128*c**(sympy.S(5)/2)*h**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_199():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)/(g + h*x)**3
    F = -(f*g**2 - h*(-d*h + e*g))*(a + b*x + c*x**2)**(sympy.S(5)/2)/(2*h*(g + h*x)**2*(a*h**2 - b*g*h + c*g**2)) - (a + b*x + c*x**2)**(sympy.S(3)/2)*(-4*a*h*(-3*e*h + 7*f*g) + b*(31*f*g**2 - 3*h*(-d*h + 5*e*g)) + 4*c*g*(-3*d*h + 6*e*g - 10*f*g**2/h) + 2*h*x*(-2*a*f*h + 2*b*f*g - 3*c*d*h + 3*c*e*g - 5*c*f*g**2/h))/(12*h**2*(g + h*x)*(a*h**2 - b*g*h + c*g**2)) + (8*c**2*g**2*(10*f*g**2 - 3*h*(-d*h + 2*e*g)) + 4*c*h*(a*h*(3*d*h**2 - 9*e*g*h + 19*f*g**2) - b*g*(6*d*h**2 - 15*e*g*h + 28*f*g**2)) + h**2*(8*a**2*f*h**2 - 4*a*b*h*(-3*e*h + 10*f*g) + b**2*(35*f*g**2 - 3*h*(-d*h + 5*e*g))))*atanh((-2*a*h + b*g + x*(-b*h + 2*c*g))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*h**2 - b*g*h + c*g**2)))/(8*h**6*sqrt(a*h**2 - b*g*h + c*g**2)) - sqrt(a + b*x + c*x**2)*(b**2*f*h**3*(-a*h + b*g) - 8*c**3*g**2*(10*f*g**2 - 3*h*(-d*h + 2*e*g)) - 2*c**2*h*(2*a*h*(3*d*h**2 - 9*e*g*h + 19*f*g**2) - 3*b*g*(5*d*h**2 - 12*e*g*h + 22*f*g**2)) - c*h**2*(8*a**2*f*h**2 - 18*a*b*h*(-e*h + 3*f*g) + b**2*(53*f*g**2 - 6*h*(-d*h + 4*e*g))) + 2*c*h*x*(b*f*h**2*(-a*h + b*g) + 2*c**2*g*(10*f*g**2 - 3*h*(-d*h + 2*e*g)) + c*h*(2*a*h*(-3*e*h + 7*f*g) - 3*b*(d*h**2 - 3*e*g*h + 6*f*g**2))))/(8*c*h**5*(a*h**2 - b*g*h + c*g**2)) - (b**3*f*h**3 + 6*b*c*h**2*(-2*a*f*h - b*e*h + 3*b*f*g) + 16*c**3*g*(10*f*g**2 - 3*h*(-d*h + 2*e*g)) + 24*c**2*h*(a*h*(-e*h + 3*f*g) - b*(d*h**2 - 3*e*g*h + 6*f*g**2)))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(3)/2)*h**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_200():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)/(g + h*x)**4
    F = -(f*g**2 - h*(-d*h + e*g))*(a + b*x + c*x**2)**(sympy.S(5)/2)/(3*h*(g + h*x)**3*(a*h**2 - b*g*h + c*g**2)) - (a + b*x + c*x**2)**(sympy.S(3)/2)*(-6*a*h*(-e*h + 3*f*g) + b*(17*f*g**2 - h*(d*h + 5*e*g)) + 2*c*g*(-d*h + 4*e*g - 10*f*g**2/h) + 2*h*x*(-3*a*f*h + 3*b*f*g - 2*c*d*h + 2*c*e*g - 5*c*f*g**2/h))/(12*h**2*(g + h*x)**2*(a*h**2 - b*g*h + c*g**2)) - sqrt(a + b*x + c*x**2)*(8*c**2*g**2*(10*f*g**2 - h*(-d*h + 4*e*g)) - 2*c*h*(-2*a*h*(2*d*h**2 - 8*e*g*h + 23*f*g**2) + 3*b*g*(d*h**2 - 6*e*g*h + 18*f*g**2)) + h**2*(12*a**2*f*h**2 - 6*a*b*h*(-e*h + 7*f*g) + b**2*(29*f*g**2 - h*(d*h + 5*e*g))) + 2*h*x*(3*b*f*h**2*(-a*h + b*g) + 2*c**2*g*(10*f*g**2 - h*(-d*h + 4*e*g)) + c*h*(6*a*h*(-e*h + 3*f*g) - b*(d*h**2 - 7*e*g*h + 22*f*g**2))))/(8*h**5*(g + h*x)*(a*h**2 - b*g*h + c*g**2)) - (-b*h**3*(24*a**2*f*h**2 - 6*a*b*h*(-e*h + 10*f*g) + b**2*(-d*h**2 - 5*e*g*h + 35*f*g**2)) + 16*c**3*g**3*(10*f*g**2 - h*(-d*h + 4*e*g)) - 24*c**2*g*h*(-a*h*(d*h**2 - 4*e*g*h + 11*f*g**2) + b*g*(d*h**2 - 5*e*g*h + 14*f*g**2)) + 6*c*h**2*(4*a**2*h**2*(-e*h + 4*f*g) - 2*a*b*h*(d*h**2 - 7*e*g*h + 25*f*g**2) + b**2*g*(d*h**2 - 10*e*g*h + 35*f*g**2)))*atanh((-2*a*h + b*g + x*(-b*h + 2*c*g))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*h**2 - b*g*h + c*g**2)))/(16*h**6*(a*h**2 - b*g*h + c*g**2)**(sympy.S(3)/2)) + (3*b**2*f*h**2 + 8*c**2*(10*f*g**2 - h*(-d*h + 4*e*g)) - 12*c*h*(-a*f*h - b*e*h + 4*b*f*g))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*sqrt(c)*h**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_201():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)/(g + h*x)**5
    F = -sqrt(c)*(-3*b*f*h - 2*c*e*h + 10*c*f*g)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*h**6) - (f*g**2 - h*(-d*h + e*g))*(a + b*x + c*x**2)**(sympy.S(5)/2)/(4*h*(g + h*x)**4*(a*h**2 - b*g*h + c*g**2)) - (a + b*x + c*x**2)**(sympy.S(3)/2)*(16*c**2*g**4*(-e*h + 5*f*g) - 4*c*g*h*(-a*h*(9*d*h**2 - 5*e*g*h + 25*f*g**2) + b*g*(3*d*h**2 - 5*e*g*h + 31*f*g**2)) - h**2*(16*a**2*h**2*(-2*e*h + f*g) + 4*a*b*h*(3*d*h**2 + 7*e*g*h + 7*f*g**2) - b**2*g*(3*d*h**2 + 5*e*g*h + 35*f*g**2)) + 3*h*x*(8*c**2*g**2*(5*f*g**2 - h*(d*h + e*g)) - 4*c*h*(-a*h*(d*h**2 - 5*e*g*h + 17*f*g**2) + 2*b*g*(-d*h**2 - 2*e*g*h + 9*f*g**2)) + h**2*(16*a**2*f*h**2 - 8*a*b*h*(-e*h + 6*f*g) + b**2*(-3*d*h**2 - 5*e*g*h + 29*f*g**2))))/(96*h**3*(g + h*x)**3*(a*h**2 - b*g*h + c*g**2)**2) + sqrt(a + b*x + c*x**2)*(-b*h**3*(48*a**2*f*h**2 - 8*a*b*h*(e*h + 10*f*g) + b**2*(3*d*h**2 + 5*e*g*h + 35*f*g**2)) + 64*c**3*g**4*(-e*h + 5*f*g) - 16*c**2*g**2*h*(-8*a*h*(-e*h + 5*f*g) + b*g*(-7*e*h + 41*f*g)) + 4*c*h**2*(16*a**2*h**2*(-e*h + 5*f*g) - a*b*h*(-3*d*h**2 - 25*e*g*h + 173*f*g**2) + 2*b**2*g**2*(-5*e*h + 46*f*g)) + 2*c*h*x*(16*c**2*g**3*(-e*h + 5*f*g) - 4*c*h*(-a*h*(35*f*g**2 - h*(-3*d*h + 7*e*g)) + 6*b*g**2*(-e*h + 6*f*g)) + h**2*(48*a**2*f*h**2 - 8*a*b*h*(-e*h + 14*f*g) + b**2*(61*f*g**2 - h*(3*d*h + 5*e*g)))))/(64*h**5*(g + h*x)*(a*h**2 - b*g*h + c*g**2)**2) + (b**2*h**4*(48*a**2*f*h**2 - 8*a*b*h*(e*h + 10*f*g) + b**2*(3*d*h**2 + 5*e*g*h + 35*f*g**2)) + 128*c**4*g**5*(-e*h + 5*f*g) - 64*c**3*g**3*h*(-5*a*h*(-e*h + 5*f*g) + b*g*(-5*e*h + 28*f*g)) - 48*c**2*h**2*(-a**2*h**2*(d*h**2 - 5*e*g*h + 25*f*g**2) + 10*a*b*g**2*h*(-e*h + 6*f*g) - 5*b**2*g**3*(-e*h + 7*f*g)) + 8*c*h**3*(24*a**3*f*h**3 - 12*a**2*b*h**2*(-e*h + 10*f*g) + 3*a*b**2*h*(-d*h**2 - 5*e*g*h + 55*f*g**2) - 5*b**3*g**2*(-e*h + 14*f*g)))*atanh((-2*a*h + b*g + x*(-b*h + 2*c*g))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*h**2 - b*g*h + c*g**2)))/(128*h**6*(a*h**2 - b*g*h + c*g**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_202():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)/(g + h*x)**6
    F = c**(sympy.S(3)/2)*f*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/h**6 - (f*g**2 - h*(-d*h + e*g))*(a + b*x + c*x**2)**(sympy.S(5)/2)/(5*h*(g + h*x)**5*(a*h**2 - b*g*h + c*g**2)) - (a + b*x + c*x**2)**(sympy.S(3)/2)*(16*c**2*f*g**5 - 2*c*g*h*(-6*a*d*h**3 - 10*a*f*g**2*h + 3*b*d*g*h**2 + 13*b*f*g**3) - h**2*(4*a**2*h**2*(-3*e*h + 2*f*g) + 2*a*b*h*(f*g**2 + 3*h*(d*h + 2*e*g)) - b**2*g*(7*f*g**2 + 3*h*(d*h + e*g))) + h*x*(4*c**2*(-3*d*g**2*h**2 + 7*f*g**4) + 2*c*g*h*(2*a*h*(-3*e*h + 14*f*g) - b*(-6*d*h**2 - 3*e*g*h + 28*f*g**2)) + h**2*(16*a**2*f*h**2 - 2*a*b*h*(-3*e*h + 22*f*g) + b**2*(25*f*g**2 - 3*h*(d*h + e*g)))))/(48*h**3*(g + h*x)**4*(a*h**2 - b*g*h + c*g**2)**2) - sqrt(a + b*x + c*x**2)*(-b*h**4*(-2*a*h + b*g)*(16*a**2*f*h**2 - 2*a*b*h*(3*e*h + 10*f*g) + b**2*(7*f*g**2 + 3*h*(d*h + e*g))) + 128*c**4*f*g**7 - 32*c**3*f*g**5*h*(-10*a*h + 11*b*g) + 8*c**2*g*h**2*(2*a**2*h**2*(3*d*h**2 + 13*f*g**2) - a*b*g*h*(3*d*h**2 + 65*f*g**2) + 38*b**2*f*g**4) - 2*c*h**3*(8*a**3*h**3*(-3*e*h + 2*f*g) + 4*a**2*b*h**2*(3*d*h**2 + 6*e*g*h + 5*f*g**2) - 2*a*b**2*g**2*h*(3*e*h + 34*f*g) + b**3*(-3*d*g**2*h**2 + 35*f*g**4)) + h*x*(128*c*f*(c*g**2 - h*(-a*h + b*g))**3 + (-b*h + 2*c*g)*(-b*h**3*(16*a**2*f*h**2 - 2*a*b*h*(3*e*h + 10*f*g) + b**2*(7*f*g**2 + 3*h*(d*h + e*g))) + 32*c**3*f*g**5 - 8*c**2*g*h*(3*a*d*h**3 - 11*a*f*g**2*h + 10*b*f*g**3) + 2*c*h**2*(4*a**2*h**2*(-3*e*h + 10*f*g) - 6*a*b*h*(-d*h**2 - e*g*h + 11*f*g**2) + b**2*(3*d*g*h**2 + 29*f*g**3)))))/(128*h**5*(g + h*x)**2*(a*h**2 - b*g*h + c*g**2)**3) - (b**3*h**5*(16*a**2*f*h**2 - 2*a*b*h*(3*e*h + 10*f*g) + b**2*(7*f*g**2 + 3*h*(d*h + e*g))) - 2*b*c*h**4*(96*a**3*f*h**3 - 24*a**2*b*h**2*(e*h + 8*f*g) + 4*a*b**2*h*(35*f*g**2 + 3*h*(d*h + e*g)) - b**3*(-3*d*g*h**2 + 35*f*g**3)) + 256*c**5*f*g**7 - 896*c**4*f*g**5*h*(-a*h + b*g) + 32*c**3*g*h**2*(a**2*h**2*(-3*d*h**2 + 35*f*g**2) - 70*a*b*f*g**3*h + 35*b**2*f*g**4) - 16*c**2*h**3*(-6*a**3*h**3*(-e*h + 6*f*g) + 3*a**2*b*h**2*(-d*h**2 - e*g*h + 35*f*g**2) - 3*a*b**2*g*h*(d*h**2 + 35*f*g**2) + 35*b**3*f*g**4))*atanh((-2*a*h + b*g + x*(-b*h + 2*c*g))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*h**2 - b*g*h + c*g**2)))/(256*h**6*(a*h**2 - b*g*h + c*g**2)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_203():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)/(g + h*x)**7
    F = (-4*a*c + b**2)**2*(24*a**2*f*h**2 - 12*a*b*h*(e*h + 2*f*g) + b**2*(7*d*h**2 + 5*e*g*h + 7*f*g**2) + 24*c**2*d*g**2 - 4*c*(a*(d*h**2 - 7*e*g*h + f*g**2) + 3*b*g*(2*d*h + e*g)))*atanh((-2*a*h + b*g + x*(-b*h + 2*c*g))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*h**2 - b*g*h + c*g**2)))/(1024*(a*h**2 - b*g*h + c*g**2)**(sympy.S(9)/2)) - (-4*a*c + b**2)*sqrt(a + b*x + c*x**2)*(-2*a*h + b*g + x*(-b*h + 2*c*g))*(24*a**2*f*h**2 - 12*a*b*h*(e*h + 2*f*g) + b**2*(7*d*h**2 + 5*e*g*h + 7*f*g**2) + 24*c**2*d*g**2 - 4*c*(a*(d*h**2 - 7*e*g*h + f*g**2) + 3*b*g*(2*d*h + e*g)))/(512*(g + h*x)**2*(a*h**2 - b*g*h + c*g**2)**4) + (a + b*x + c*x**2)**(sympy.S(3)/2)*(-2*a*h + b*g + x*(-b*h + 2*c*g))*(24*a**2*f*h**2 - 12*a*b*h*(e*h + 2*f*g) + b**2*(7*d*h**2 + 5*e*g*h + 7*f*g**2) + 24*c**2*d*g**2 - 4*c*(a*(d*h**2 - 7*e*g*h + f*g**2) + 3*b*g*(2*d*h + e*g)))/(192*(g + h*x)**4*(a*h**2 - b*g*h + c*g**2)**3) + (2*c*g*(5*f*g**2 + h*(-7*d*h + e*g)) + h*(12*a*h*(-e*h + 2*f*g) - b*(-7*d*h**2 - 5*e*g*h + 17*f*g**2)))*(a + b*x + c*x**2)**(sympy.S(5)/2)/(60*h*(g + h*x)**5*(a*h**2 - b*g*h + c*g**2)**2) - (f*g**2 - h*(-d*h + e*g))*(a + b*x + c*x**2)**(sympy.S(5)/2)/(6*h*(g + h*x)**6*(a*h**2 - b*g*h + c*g**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_204():
    f = (a + b*x + c*x**2)**(sympy.S(3)/2)*(d + e*x + f*x**2)/(g + h*x)**8
    F = (-4*a*c + b**2)**2*(-b*h*(24*a**2*f*h**2 - 2*a*b*h*(7*e*h + 10*f*g) + b**2*(9*d*h**2 + 5*e*g*h + 5*f*g**2)) + 48*c**3*d*g**3 - 8*c**2*g*(a*(3*d*h**2 - 8*e*g*h + f*g**2) + 3*b*g*(3*d*h + e*g)) + 2*c*(4*a**2*h**2*(-e*h + 8*f*g) - 2*a*b*h*(-3*d*h**2 + 13*e*g*h + 13*f*g**2) + b**2*g*(21*d*h**2 + 10*e*g*h + 7*f*g**2)))*atanh((-2*a*h + b*g + x*(-b*h + 2*c*g))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*h**2 - b*g*h + c*g**2)))/(2048*(a*h**2 - b*g*h + c*g**2)**(sympy.S(11)/2)) - (-4*a*c + b**2)*sqrt(a + b*x + c*x**2)*(-2*a*h + b*g + x*(-b*h + 2*c*g))*(-b*h*(24*a**2*f*h**2 - 2*a*b*h*(7*e*h + 10*f*g) + b**2*(9*d*h**2 + 5*e*g*h + 5*f*g**2)) + 48*c**3*d*g**3 - 8*c**2*g*(a*(3*d*h**2 - 8*e*g*h + f*g**2) + 3*b*g*(3*d*h + e*g)) + 2*c*(4*a**2*h**2*(-e*h + 8*f*g) - 2*a*b*h*(-3*d*h**2 + 13*e*g*h + 13*f*g**2) + b**2*g*(21*d*h**2 + 10*e*g*h + 7*f*g**2)))/(1024*(g + h*x)**2*(a*h**2 - b*g*h + c*g**2)**5) + (a + b*x + c*x**2)**(sympy.S(3)/2)*(-2*a*h + b*g + x*(-b*h + 2*c*g))*(-b*h*(24*a**2*f*h**2 - 2*a*b*h*(7*e*h + 10*f*g) + b**2*(9*d*h**2 + 5*e*g*h + 5*f*g**2)) + 48*c**3*d*g**3 - 8*c**2*g*(a*(3*d*h**2 - 8*e*g*h + f*g**2) + 3*b*g*(3*d*h + e*g)) + 2*c*(4*a**2*h**2*(-e*h + 8*f*g) - 2*a*b*h*(-3*d*h**2 + 13*e*g*h + 13*f*g**2) + b**2*g*(21*d*h**2 + 10*e*g*h + 7*f*g**2)))/(384*(g + h*x)**4*(a*h**2 - b*g*h + c*g**2)**4) + (a + b*x + c*x**2)**(sympy.S(5)/2)*(4*c**2*g**2*(5*f*g**2 + h*(-51*d*h + 2*e*g)) - 2*c*h*(-2*a*h*(12*d*h**2 - 61*e*g*h + 26*f*g**2) + 3*b*g*(-34*d*h**2 - 15*e*g*h + 8*f*g**2)) - 7*h**2*(24*a**2*f*h**2 - 2*a*b*h*(7*e*h + 10*f*g) + b**2*(9*d*h**2 + 5*e*g*h + 5*f*g**2)))/(840*h*(g + h*x)**5*(a*h**2 - b*g*h + c*g**2)**3) + (2*c*g*(5*f*g**2 + h*(-9*d*h + 2*e*g)) + h*(14*a*h*(-e*h + 2*f*g) - b*(-9*d*h**2 - 5*e*g*h + 19*f*g**2)))*(a + b*x + c*x**2)**(sympy.S(5)/2)/(84*h*(g + h*x)**6*(a*h**2 - b*g*h + c*g**2)**2) - (f*g**2 - h*(-d*h + e*g))*(a + b*x + c*x**2)**(sympy.S(5)/2)/(7*h*(g + h*x)**7*(a*h**2 - b*g*h + c*g**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_205():
    f = (2*x + 1)**3*sqrt(3*x**2 - x + 2)*(4*x**2 + 3*x + 1)
    F = (5393 - 32358*x)*sqrt(3*x**2 - x + 2)/15552 + 2*(2*x + 1)**4*(3*x**2 - x + 2)**(sympy.S(3)/2)/21 + 67*(2*x + 1)**3*(3*x**2 - x + 2)**(sympy.S(3)/2)/378 + 17*(2*x + 1)**2*(3*x**2 - x + 2)**(sympy.S(3)/2)/105 - (26982*x + 75295)*(3*x**2 - x + 2)**(sympy.S(3)/2)/68040 + 124039*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/93312
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_206():
    f = (2*x + 1)**2*sqrt(3*x**2 - x + 2)*(4*x**2 + 3*x + 1)
    F = (235 - 1410*x)*sqrt(3*x**2 - x + 2)/1296 + (2*x + 1)**3*(3*x**2 - x + 2)**(sympy.S(3)/2)/9 + (2*x + 1)**2*(3*x**2 - x + 2)**(sympy.S(3)/2)/5 + (306*x + 25)*(3*x**2 - x + 2)**(sympy.S(3)/2)/810 + 5405*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/7776
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_207():
    f = (2*x + 1)*sqrt(3*x**2 - x + 2)*(4*x**2 + 3*x + 1)
    F = (19 - 114*x)*sqrt(3*x**2 - x + 2)/2592 + 2*(2*x + 1)**2*(3*x**2 - x + 2)**(sympy.S(3)/2)/15 + (738*x + 745)*(3*x**2 - x + 2)**(sympy.S(3)/2)/1620 + 437*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/15552
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_208():
    f = sqrt(3*x**2 - x + 2)*(4*x**2 + 3*x + 1)/(2*x + 1)
    F = (30*x + 13)*sqrt(3*x**2 - x + 2)/72 + 2*(3*x**2 - x + 2)**(sympy.S(3)/2)/9 - 43*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/432 - sqrt(13)*atanh(sqrt(13)*(9 - 8*x)/(26*sqrt(3*x**2 - x + 2)))/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_209():
    f = sqrt(3*x**2 - x + 2)*(4*x**2 + 3*x + 1)/(2*x + 1)**2
    F = -(67 - 96*x)*sqrt(3*x**2 - x + 2)/156 - 11*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/18 + 17*sqrt(13)*atanh(sqrt(13)*(9 - 8*x)/(26*sqrt(3*x**2 - x + 2)))/104 - (3*x**2 - x + 2)**(sympy.S(3)/2)/(26*x + 13)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_210():
    f = sqrt(3*x**2 - x + 2)*(4*x**2 + 3*x + 1)/(2*x + 1)**3
    F = (110*x + 77)*sqrt(3*x**2 - x + 2)/(208*x + 104) + 11*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/24 - 803*sqrt(13)*atanh(sqrt(13)*(9 - 8*x)/(26*sqrt(3*x**2 - x + 2)))/2704 - (3*x**2 - x + 2)**(sympy.S(3)/2)/(26*(2*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_211():
    f = (2*x + 1)**3*(3*x**2 - x + 2)**(sympy.S(3)/2)*(4*x**2 + 3*x + 1)
    F = 77*x**3*(3*x**2 - x + 2)**(sympy.S(5)/2)/81 + 913*x**2*(3*x**2 - x + 2)**(sympy.S(5)/2)/486 - (3113 - 64350*x)*(3*x**2 - x + 2)**(sympy.S(5)/2)/58320 + (54593 - 327558*x)*(3*x**2 - x + 2)**(sympy.S(3)/2)/559872 + (1255639 - 7533834*x)*sqrt(3*x**2 - x + 2)/4478976 + 2*(2*x + 1)**4*(3*x**2 - x + 2)**(sympy.S(5)/2)/27 + 28879697*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/26873856
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_212():
    f = (2*x + 1)**2*(3*x**2 - x + 2)**(sympy.S(3)/2)*(4*x**2 + 3*x + 1)
    F = (91 - 546*x)*(3*x**2 - x + 2)**(sympy.S(3)/2)/3456 + (2093 - 12558*x)*sqrt(3*x**2 - x + 2)/27648 + (2*x + 1)**3*(3*x**2 - x + 2)**(sympy.S(5)/2)/12 + 8*(2*x + 1)**2*(3*x**2 - x + 2)**(sympy.S(5)/2)/63 + (650*x + 377)*(3*x**2 - x + 2)**(sympy.S(5)/2)/2520 + 48139*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/165888
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_213():
    f = (2*x + 1)*(3*x**2 - x + 2)**(sympy.S(3)/2)*(4*x**2 + 3*x + 1)
    F = -(71 - 426*x)*(3*x**2 - x + 2)**(sympy.S(3)/2)/2592 + 2*(2*x + 1)**2*(3*x**2 - x + 2)**(sympy.S(5)/2)/21 + (102*x + 109)*(3*x**2 - x + 2)**(sympy.S(5)/2)/378 + (9798*x - 1633)*sqrt(3*x**2 - x + 2)/20736 - 37559*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/124416
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_214():
    f = (3*x**2 - x + 2)**(sympy.S(3)/2)*(4*x**2 + 3*x + 1)/(2*x + 1)
    F = (30*x + 7)*(3*x**2 - x + 2)**(sympy.S(3)/2)/144 + (402*x + 869)*sqrt(3*x**2 - x + 2)/1152 + 2*(3*x**2 - x + 2)**(sympy.S(5)/2)/15 + 2203*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/6912 - 13*sqrt(13)*atanh(sqrt(13)*(9 - 8*x)/(26*sqrt(3*x**2 - x + 2)))/32
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_215():
    f = (3*x**2 - x + 2)**(sympy.S(3)/2)*(4*x**2 + 3*x + 1)/(2*x + 1)**2
    F = -(23 - 38*x)*(3*x**2 - x + 2)**(sympy.S(3)/2)/104 - (349 - 294*x)*sqrt(3*x**2 - x + 2)/192 - 2327*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/1152 + 25*sqrt(13)*atanh(sqrt(13)*(9 - 8*x)/(26*sqrt(3*x**2 - x + 2)))/32 - (3*x**2 - x + 2)**(sympy.S(5)/2)/(26*x + 13)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_216():
    f = (3*x**2 - x + 2)**(sympy.S(3)/2)*(4*x**2 + 3*x + 1)/(2*x + 1)**3
    F = (1858 - 771*x)*sqrt(3*x**2 - x + 2)/624 + (122*x + 151)*(3*x**2 - x + 2)**(sympy.S(3)/2)/(624*x + 312) + 1519*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/576 - 1153*sqrt(13)*atanh(sqrt(13)*(9 - 8*x)/(26*sqrt(3*x**2 - x + 2)))/832 - (3*x**2 - x + 2)**(sympy.S(5)/2)/(26*(2*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_217():
    f = (2*x + 1)**3*(3*x**2 - x + 2)**(sympy.S(5)/2)*(4*x**2 + 3*x + 1)
    F = (5089 - 30534*x)*(3*x**2 - x + 2)**(sympy.S(5)/2)/155520 - (26353 - 21350*x)*(3*x**2 - x + 2)**(sympy.S(7)/2)/498960 + (117047 - 702282*x)*(3*x**2 - x + 2)**(sympy.S(3)/2)/1492992 + (2692081 - 16152486*x)*sqrt(3*x**2 - x + 2)/11943936 + 2*(2*x + 1)**4*(3*x**2 - x + 2)**(sympy.S(7)/2)/33 + 29*(2*x + 1)**3*(3*x**2 - x + 2)**(sympy.S(7)/2)/330 + 133*(2*x + 1)**2*(3*x**2 - x + 2)**(sympy.S(7)/2)/1485 + 61917863*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/71663616
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_218():
    f = (2*x + 1)**2*(3*x**2 - x + 2)**(sympy.S(5)/2)*(4*x**2 + 3*x + 1)
    F = -(293 - 1758*x)*(3*x**2 - x + 2)**(sympy.S(5)/2)/58320 - (6739 - 40434*x)*(3*x**2 - x + 2)**(sympy.S(3)/2)/559872 + (2*x + 1)**3*(3*x**2 - x + 2)**(sympy.S(7)/2)/15 + 37*(2*x + 1)**2*(3*x**2 - x + 2)**(sympy.S(7)/2)/405 + (3430*x + 2731)*(3*x**2 - x + 2)**(sympy.S(7)/2)/17010 + (929982*x - 154997)*sqrt(3*x**2 - x + 2)/4478976 - 3564931*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/26873856
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_219():
    f = (2*x + 1)*(3*x**2 - x + 2)**(sympy.S(5)/2)*(4*x**2 + 3*x + 1)
    F = -(445 - 2670*x)*(3*x**2 - x + 2)**(sympy.S(5)/2)/15552 - (51175 - 307050*x)*(3*x**2 - x + 2)**(sympy.S(3)/2)/746496 + 2*(2*x + 1)**2*(3*x**2 - x + 2)**(sympy.S(7)/2)/27 + (122*x + 137)*(3*x**2 - x + 2)**(sympy.S(7)/2)/648 + (7062150*x - 1177025)*sqrt(3*x**2 - x + 2)/5971968 - 27071575*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/35831808
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_220():
    f = (3*x**2 - x + 2)**(sympy.S(5)/2)*(4*x**2 + 3*x + 1)/(2*x + 1)
    F = (221999 - 17850*x)*sqrt(3*x**2 - x + 2)/82944 + (150*x + 29)*(3*x**2 - x + 2)**(sympy.S(5)/2)/1080 + (2154*x + 2449)*(3*x**2 - x + 2)**(sympy.S(3)/2)/10368 + 2*(3*x**2 - x + 2)**(sympy.S(7)/2)/21 + 944521*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/497664 - 169*sqrt(13)*atanh(sqrt(13)*(9 - 8*x)/(26*sqrt(3*x**2 - x + 2)))/128
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_221():
    f = (3*x**2 - x + 2)**(sympy.S(5)/2)*(4*x**2 + 3*x + 1)/(2*x + 1)**2
    F = -(407 - 660*x)*(3*x**2 - x + 2)**(sympy.S(5)/2)/2340 - (737 - 858*x)*(3*x**2 - x + 2)**(sympy.S(3)/2)/864 + (33990*x - 51997)*sqrt(3*x**2 - x + 2)/6912 - 315623*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/41472 + 429*sqrt(13)*atanh(sqrt(13)*(9 - 8*x)/(26*sqrt(3*x**2 - x + 2)))/128 - (3*x**2 - x + 2)**(sympy.S(7)/2)/(26*x + 13)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_222():
    f = (3*x**2 - x + 2)**(sympy.S(5)/2)*(4*x**2 + 3*x + 1)/(2*x + 1)**3
    F = (1227 - 838*x)*(3*x**2 - x + 2)**(sympy.S(3)/2)/832 + (21317 - 10470*x)*sqrt(3*x**2 - x + 2)/1536 + (134*x + 257)*(3*x**2 - x + 2)**(sympy.S(5)/2)/(1040*x + 520) + 118423*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/9216 - 1631*sqrt(13)*atanh(sqrt(13)*(9 - 8*x)/(26*sqrt(3*x**2 - x + 2)))/256 - (3*x**2 - x + 2)**(sympy.S(7)/2)/(26*(2*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_223():
    f = (g + h*x)**3*(d + e*x + f*x**2)/sqrt(a + b*x + c*x**2)
    F = f*(g + h*x)**4*sqrt(a + b*x + c*x**2)/(5*c*h) - (g + h*x)**3*(9*b*f*h + 2*c*(-5*e*h + f*g))*sqrt(a + b*x + c*x**2)/(40*c**2*h) + (g + h*x)**2*sqrt(a + b*x + c*x**2)*(63*b**2*f*h**2 - c**2*(12*f*g**2 - 20*h*(4*d*h + 3*e*g)) - 2*c*h*(32*a*f*h + 35*b*e*h + 24*b*f*g))/(240*c**3*h) + sqrt(a + b*x + c*x**2)*(945*b**4*f*h**4 - 210*b**2*c*h**3*(14*a*f*h + 5*b*(e*h + 3*f*g)) - 64*c**4*g**2*(3*f*g**2 - 5*h*(16*d*h + 3*e*g)) - 16*c**3*h*(16*a*h*(13*f*g**2 + 5*h*(d*h + 3*e*g)) + b*g*(39*f*g**2 + 5*h*(54*d*h + 47*e*g))) + 8*c**2*h**2*(128*a**2*f*h**2 + 275*a*b*h*(e*h + 3*f*g) + 3*b**2*(129*f*g**2 + 50*h*(d*h + 3*e*g))) - 2*c*h*x*(315*b**3*f*h**3 - 14*b*c*h**2*(46*a*f*h + 25*b*e*h + 39*b*f*g) + 16*c**3*g*(3*f*g**2 - 5*h*(10*d*h + 3*e*g)) + 8*c**2*h*(a*h*(45*e*h + 71*f*g) + b*(50*d*h**2 + 80*e*g*h + 21*f*g**2))))/(1920*c**5*h) + (-63*b**5*f*h**3 + 70*b**3*c*h**2*(4*a*f*h + b*e*h + 3*b*f*g) - 80*b*c**2*h*(3*a**2*f*h**2 + 3*a*b*h*(e*h + 3*f*g) + b**2*(d*h**2 + 3*e*g*h + 3*f*g**2)) + 256*c**5*d*g**3 - 128*c**4*g*(a*(f*g**2 + 3*h*(d*h + e*g)) + b*g*(3*d*h + e*g)) + 96*c**3*(a**2*h**2*(e*h + 3*f*g) + 2*a*b*h*(3*f*g**2 + h*(d*h + 3*e*g)) + b**2*g*(f*g**2 + 3*h*(d*h + e*g))))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(256*c**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_224():
    f = (g + h*x)**2*(d + e*x + f*x**2)/sqrt(a + b*x + c*x**2)
    F = f*(g + h*x)**3*sqrt(a + b*x + c*x**2)/(4*c*h) - (g + h*x)**2*sqrt(a + b*x + c*x**2)*(7*b*f*h - 8*c*e*h + 2*c*f*g)/(24*c**2*h) - sqrt(a + b*x + c*x**2)*(105*b**3*f*h**3 - 20*b*c*h**2*(11*a*f*h + 6*b*(e*h + 2*f*g)) + 32*c**3*g*(f*g**2 - 4*h*(3*d*h + e*g)) + 8*c**2*h*(16*a*h*(e*h + 2*f*g) + b*(11*f*g**2 + 18*h*(d*h + 2*e*g))) - 2*c*h*x*(35*b**2*f*h**2 - 8*c**2*(f*g**2 - 2*h*(3*d*h + 2*e*g)) - 4*c*h*(9*a*f*h + 10*b*e*h + 6*b*f*g)))/(192*c**4*h) + (35*b**4*f*h**2 - 40*b**2*c*h*(3*a*f*h + b*e*h + 2*b*f*g) + 128*c**4*d*g**2 - 64*c**3*(a*(d*h**2 + 2*e*g*h + f*g**2) + b*g*(2*d*h + e*g)) + 48*c**2*(a**2*f*h**2 + 2*a*b*h*(e*h + 2*f*g) + b**2*(d*h**2 + 2*e*g*h + f*g**2)))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(128*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_225():
    f = (g + h*x)*(d + e*x + f*x**2)/sqrt(a + b*x + c*x**2)
    F = f*(g + h*x)**2*sqrt(a + b*x + c*x**2)/(3*c*h) + sqrt(a + b*x + c*x**2)*(15*b**2*f*h**2 - 8*c**2*(f*g**2 - 3*h*(d*h + e*g)) - 2*c*h*x*(5*b*f*h - 6*c*e*h + 2*c*f*g) - 2*c*h*(8*a*f*h + 9*b*(e*h + f*g)))/(24*c**3*h) + (-5*b**3*f*h + 6*b*c*(2*a*f*h + b*e*h + b*f*g) + 16*c**3*d*g - 8*c**2*(a*e*h + a*f*g + b*d*h + b*e*g))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_226():
    f = (d + e*x + f*x**2)/sqrt(a + b*x + c*x**2)
    F = f*x*sqrt(a + b*x + c*x**2)/(2*c) + (-3*b*f + 4*c*e)*sqrt(a + b*x + c*x**2)/(4*c**2) + (3*b**2*f + 8*c**2*d - 4*c*(a*f + b*e))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_227():
    f = (d + e*x + f*x**2)/((g + h*x)*sqrt(a + b*x + c*x**2))
    F = (f*g**2 - h*(-d*h + e*g))*atanh((-2*a*h + b*g + x*(-b*h + 2*c*g))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*h**2 - b*g*h + c*g**2)))/(h**2*sqrt(a*h**2 - b*g*h + c*g**2)) + f*sqrt(a + b*x + c*x**2)/(c*h) - (b*f*h - 2*c*e*h + 2*c*f*g)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*c**(sympy.S(3)/2)*h**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_228():
    f = (d + e*x + f*x**2)/((g + h*x)**2*sqrt(a + b*x + c*x**2))
    F = -(f*g**2 - h*(-d*h + e*g))*sqrt(a + b*x + c*x**2)/(h*(g + h*x)*(a*h**2 - b*g*h + c*g**2)) - (2*c*(-d*g*h**2 + f*g**3) + h*(2*a*h*(-e*h + 2*f*g) - b*(-d*h**2 - e*g*h + 3*f*g**2)))*atanh((-2*a*h + b*g + x*(-b*h + 2*c*g))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*h**2 - b*g*h + c*g**2)))/(2*h**2*(a*h**2 - b*g*h + c*g**2)**(sympy.S(3)/2)) + f*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(sqrt(c)*h**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_229():
    f = (d + e*x + f*x**2)/((g + h*x)**3*sqrt(a + b*x + c*x**2))
    F = (8*a**2*f*h**2 - 4*a*b*h*(e*h + 2*f*g) + b**2*(3*d*h**2 + e*g*h + 3*f*g**2) + 8*c**2*d*g**2 - 4*c*(a*(d*h**2 - 3*e*g*h + f*g**2) + b*g*(2*d*h + e*g)))*atanh((-2*a*h + b*g + x*(-b*h + 2*c*g))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*h**2 - b*g*h + c*g**2)))/(8*(a*h**2 - b*g*h + c*g**2)**(sympy.S(5)/2)) + (2*c*g*(f*g**2 + h*(-3*d*h + e*g)) + h*(4*a*h*(-e*h + 2*f*g) - b*(-3*d*h**2 - e*g*h + 5*f*g**2)))*sqrt(a + b*x + c*x**2)/(4*h*(g + h*x)*(a*h**2 - b*g*h + c*g**2)**2) - (f*g**2 - h*(-d*h + e*g))*sqrt(a + b*x + c*x**2)/(2*h*(g + h*x)**2*(a*h**2 - b*g*h + c*g**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_230():
    f = (g + h*x)**3*(d + e*x + f*x**2)/(a + b*x + c*x**2)**(sympy.S(3)/2)
    F = (g + h*x)**3*(2*c*(2*a*e - b*(a*f/c + d)) - 2*x*(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d))/(c*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)) + h*(g + h*x)**2*sqrt(a + b*x + c*x**2)*(-16*a*c*f + 7*b**2*f - 6*b*c*e + 12*c**2*d)/(3*c**2*(-4*a*c + b**2)) + h*sqrt(a + b*x + c*x**2)*(105*b**4*f*h**2 - 10*b**2*c*h*(46*a*f*h + 9*b*(e*h + 3*f*g)) + 192*c**4*d*g**2 - 16*c**3*(4*a*(3*d*h**2 + 9*e*g*h + 7*f*g**2) + 3*b*g*(3*d*h + 2*e*g)) + 8*c**2*(32*a**2*f*h**2 + 39*a*b*h*(e*h + 3*f*g) + b**2*(20*f*g**2 + 9*h*(d*h + 3*e*g))) + 2*c*h*x*(-35*b**3*f*h + 2*b*c*(58*a*f*h + 15*b*e*h + 17*b*f*g) + 48*c**3*d*g - 8*c**2*(9*a*e*h + 11*a*f*g + 3*b*d*h + 3*b*e*g)))/(24*c**4*(-4*a*c + b**2)) - (35*b**3*f*h**3 - 30*b*c*h**2*(2*a*f*h + b*e*h + 3*b*f*g) - 16*c**3*g*(f*g**2 + 3*h*(d*h + e*g)) + 24*c**2*h*(a*h*(e*h + 3*f*g) + b*(d*h**2 + 3*e*g*h + 3*f*g**2)))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_231():
    f = (g + h*x)**2*(d + e*x + f*x**2)/(a + b*x + c*x**2)**(sympy.S(3)/2)
    F = (g + h*x)**2*(2*c*(2*a*e - b*(a*f/c + d)) - 2*x*(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d))/(c*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)) + h*sqrt(a + b*x + c*x**2)*(-15*b**3*f*h + 4*b*c*(13*a*f*h + 3*b*e*h + 6*b*f*g) + 32*c**3*d*g - 8*c**2*(4*a*e*h + 8*a*f*g + b*d*h + 2*b*e*g) + 2*c*h*x*(-12*a*c*f + 5*b**2*f - 4*b*c*e + 8*c**2*d))/(4*c**3*(-4*a*c + b**2)) + (15*b**2*f*h**2 + 8*c**2*(f*g**2 + h*(d*h + 2*e*g)) - 12*c*h*(a*f*h + b*e*h + 2*b*f*g))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_232():
    f = (g + h*x)*(d + e*x + f*x**2)/(a + b*x + c*x**2)**(sympy.S(3)/2)
    F = (g + h*x)*(2*c*(2*a*e - b*(a*f/c + d)) - 2*x*(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d))/(c*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)) + h*sqrt(a + b*x + c*x**2)*(-8*a*c*f + 3*b**2*f - 2*b*c*e + 4*c**2*d)/(c**2*(-4*a*c + b**2)) - (3*b*f*h - 2*c*(e*h + f*g))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_233():
    f = (d + e*x + f*x**2)/(a + b*x + c*x**2)**(sympy.S(3)/2)
    F = (2*c*(2*a*e - b*(a*f/c + d)) - 2*x*(-2*a*c*f + b**2*f - b*c*e + 2*c**2*d))/(c*(-4*a*c + b**2)*sqrt(a + b*x + c*x**2)) + f*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/c**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_234():
    f = (d + e*x + f*x**2)/((g + h*x)*(a + b*x + c*x**2)**(sympy.S(3)/2))
    F = (f*g**2 - h*(-d*h + e*g))*atanh((-2*a*h + b*g + x*(-b*h + 2*c*g))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*h**2 - b*g*h + c*g**2)))/(a*h**2 - b*g*h + c*g**2)**(sympy.S(3)/2) + (4*a*(a*f*h - c*d*h + c*e*g) + 2*b**2*d*h - 2*b*(a*e*h + a*f*g + c*d*g) - 2*x*(b*f*(-a*h + b*g) + 2*c**2*d*g - c*(-2*a*e*h + 2*a*f*g + b*d*h + b*e*g)))/((-4*a*c + b**2)*sqrt(a + b*x + c*x**2)*(a*h**2 - b*g*h + c*g**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_235():
    f = (d + e*x + f*x**2)/((g + h*x)**2*(a + b*x + c*x**2)**(sympy.S(3)/2))
    F = -h*(f*g**2 - h*(-d*h + e*g))*sqrt(a + b*x + c*x**2)/((g + h*x)*(a*h**2 - b*g*h + c*g**2)**2) + (2*c*g*(f*g**2 - h*(-3*d*h + 2*e*g)) - h*(2*a*h*(-e*h + 2*f*g) - b*(-3*d*h**2 + e*g*h + f*g**2)))*atanh((-2*a*h + b*g + x*(-b*h + 2*c*g))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*h**2 - b*g*h + c*g**2)))/(2*(a*h**2 - b*g*h + c*g**2)**(sympy.S(5)/2)) - (-4*a*c*(a*h*(-e*h + 2*f*g) + c*g*(-2*d*h + e*g)) + 2*b**3*d*h**2 - 2*b**2*h*(a*e*h + 2*c*d*g) + 2*b*(a**2*f*h**2 + a*c*(-3*d*h**2 + 2*e*g*h + f*g**2) + c**2*d*g**2) + 2*c*x*(2*a**2*f*h**2 - a*b*h*(e*h + 2*f*g) + b**2*(d*h**2 + f*g**2) + 2*c**2*d*g**2 - c*(2*a*(d*h**2 - 2*e*g*h + f*g**2) + b*g*(2*d*h + e*g))))/((-4*a*c + b**2)*sqrt(a + b*x + c*x**2)*(a*h**2 - b*g*h + c*g**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_236():
    f = (d + e*x + f*x**2)/((g + h*x)**3*(a + b*x + c*x**2)**(sympy.S(3)/2))
    F = -h*(2*c*g*(3*f*g**2 - h*(-7*d*h + 5*e*g)) - h*(4*a*h*(-e*h + 2*f*g) - b*(-7*d*h**2 + 3*e*g*h + f*g**2)))*sqrt(a + b*x + c*x**2)/(4*(g + h*x)*(a*h**2 - b*g*h + c*g**2)**3) - h*(f*g**2 - h*(-d*h + e*g))*sqrt(a + b*x + c*x**2)/(2*(g + h*x)**2*(a*h**2 - b*g*h + c*g**2)**2) + (8*c**2*g**2*(6*d*h**2 - 3*e*g*h + f*g**2) - 4*c*h*(a*h*(3*d*h**2 - 9*e*g*h + 11*f*g**2) - b*g*(2*f*g**2 + 3*h*(-4*d*h + e*g))) + h**2*(8*a**2*f*h**2 + 4*a*b*h*(-3*e*h + 2*f*g) - b**2*(f*g**2 + 3*h*(-5*d*h + e*g))))*atanh((-2*a*h + b*g + x*(-b*h + 2*c*g))/(2*sqrt(a + b*x + c*x**2)*sqrt(a*h**2 - b*g*h + c*g**2)))/(8*(a*h**2 - b*g*h + c*g**2)**(sympy.S(7)/2)) + (-4*a*c*(a**2*f*h**3 - a*c*h*(d*h**2 - 3*e*g*h + 3*f*g**2) - c**2*g**2*(-3*d*h + e*g)) + 2*b**4*d*h**3 - 2*b**3*h**2*(a*e*h + 3*c*d*g) + 2*b**2*h*(a**2*f*h**2 + a*c*h*(-4*d*h + 3*e*g) + 3*c**2*d*g**2) - 2*b*c*(3*a**2*h**2*(-e*h + f*g) + a*c*g*(-9*d*h**2 + 3*e*g*h + f*g**2) + c**2*d*g**3) - 2*c*x*(-b*h**3*(a**2*f - a*b*e + b**2*d) + 2*c**3*d*g**3 - c**2*g*(2*a*(3*d*h**2 - 3*e*g*h + f*g**2) + b*g*(3*d*h + e*g)) + c*(2*a**2*h**2*(-e*h + 3*f*g) - 3*a*b*h*(-d*h**2 + e*g*h + f*g**2) + b**2*(3*d*g*h**2 + f*g**3))))/((-4*a*c + b**2)*sqrt(a + b*x + c*x**2)*(a*h**2 - b*g*h + c*g**2)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_237():
    f = (2*x + 1)**3*(4*x**2 + 3*x + 1)/sqrt(3*x**2 - x + 2)
    F = 2*(2*x + 1)**4*sqrt(3*x**2 - x + 2)/15 + 19*(2*x + 1)**3*sqrt(3*x**2 - x + 2)/60 + 44*(2*x + 1)**2*sqrt(3*x**2 - x + 2)/135 - (6298*x + 24897)*sqrt(3*x**2 - x + 2)/3240 + 9211*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/3888
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_238():
    f = (2*x + 1)**2*(4*x**2 + 3*x + 1)/sqrt(3*x**2 - x + 2)
    F = (2*x + 1)**3*sqrt(3*x**2 - x + 2)/6 + 11*(2*x + 1)**2*sqrt(3*x**2 - x + 2)/27 + (286*x - 429)*sqrt(3*x**2 - x + 2)/324 + 4147*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/1944
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_239():
    f = (2*x + 1)*(4*x**2 + 3*x + 1)/sqrt(3*x**2 - x + 2)
    F = 2*(2*x + 1)**2*sqrt(3*x**2 - x + 2)/9 + (62*x + 69)*sqrt(3*x**2 - x + 2)/54 + 251*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/324
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_240():
    f = (4*x**2 + 3*x + 1)/((2*x + 1)*sqrt(3*x**2 - x + 2))
    F = 2*sqrt(3*x**2 - x + 2)/3 - 5*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/18 - sqrt(13)*atanh(sqrt(13)*(9 - 8*x)/(26*sqrt(3*x**2 - x + 2)))/26
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_241():
    f = (4*x**2 + 3*x + 1)/((2*x + 1)**2*sqrt(3*x**2 - x + 2))
    F = -sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/3 + 9*sqrt(13)*atanh(sqrt(13)*(9 - 8*x)/(26*sqrt(3*x**2 - x + 2)))/338 - sqrt(3*x**2 - x + 2)/(26*x + 13)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_242():
    f = (4*x**2 + 3*x + 1)/((2*x + 1)**3*sqrt(3*x**2 - x + 2))
    F = -581*sqrt(13)*atanh(sqrt(13)*(9 - 8*x)/(26*sqrt(3*x**2 - x + 2)))/8788 + 7*sqrt(3*x**2 - x + 2)/(338*x + 169) - sqrt(3*x**2 - x + 2)/(26*(2*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_243():
    f = (2*x + 1)**3*(4*x**2 + 3*x + 1)/(3*x**2 - x + 2)**(sympy.S(3)/2)
    F = 32*x**2*sqrt(3*x**2 - x + 2)/27 + 412*x*sqrt(3*x**2 - x + 2)/81 + (25678 - 7742*x)/(1863*sqrt(3*x**2 - x + 2)) + 746*sqrt(3*x**2 - x + 2)/81 + 353*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/243
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_244():
    f = (2*x + 1)**2*(4*x**2 + 3*x + 1)/(3*x**2 - x + 2)**(sympy.S(3)/2)
    F = 8*x*sqrt(3*x**2 - x + 2)/9 + (2498 - 4546*x)/(621*sqrt(3*x**2 - x + 2)) + 112*sqrt(3*x**2 - x + 2)/27 - 64*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/27
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_245():
    f = (2*x + 1)*(4*x**2 + 3*x + 1)/(3*x**2 - x + 2)**(sympy.S(3)/2)
    F = -(734*x + 146)/(207*sqrt(3*x**2 - x + 2)) + 8*sqrt(3*x**2 - x + 2)/9 - 14*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_246():
    f = (4*x**2 + 3*x + 1)/((2*x + 1)*(3*x**2 - x + 2)**(sympy.S(3)/2))
    F = (154*x - 202)/(299*sqrt(3*x**2 - x + 2)) - 2*sqrt(13)*atanh(sqrt(13)*(9 - 8*x)/(26*sqrt(3*x**2 - x + 2)))/169
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_247():
    f = (4*x**2 + 3*x + 1)/((2*x + 1)**2*(3*x**2 - x + 2)**(sympy.S(3)/2))
    F = -(394 - 1674*x)/(3887*sqrt(3*x**2 - x + 2)) + 2*sqrt(13)*atanh(sqrt(13)*(9 - 8*x)/(26*sqrt(3*x**2 - x + 2)))/2197 - 4*sqrt(3*x**2 - x + 2)/(338*x + 169)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_248():
    f = (4*x**2 + 3*x + 1)/((2*x + 1)**3*(3*x**2 - x + 2)**(sympy.S(3)/2))
    F = (7386*x + 4726)/(50531*sqrt(3*x**2 - x + 2)) - 487*sqrt(13)*atanh(sqrt(13)*(9 - 8*x)/(26*sqrt(3*x**2 - x + 2)))/28561 - 4*sqrt(3*x**2 - x + 2)/(4394*x + 2197) - 2*sqrt(3*x**2 - x + 2)/(169*(2*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_249():
    f = (2*x + 1)**3*(4*x**2 + 3*x + 1)/(3*x**2 - x + 2)**(sympy.S(5)/2)
    F = (25678 - 7742*x)/(5589*(3*x**2 - x + 2)**(sympy.S(3)/2)) - (1182720*x + 1002652)/(128547*sqrt(3*x**2 - x + 2)) + 32*sqrt(3*x**2 - x + 2)/27 - 296*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/81
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_250():
    f = (2*x + 1)**2*(4*x**2 + 3*x + 1)/(3*x**2 - x + 2)**(sympy.S(5)/2)
    F = (2498 - 4546*x)/(1863*(3*x**2 - x + 2)**(sympy.S(3)/2)) - (186056 - 11784*x)/(42849*sqrt(3*x**2 - x + 2)) - 16*sqrt(3)*asinh(sqrt(23)*(1 - 6*x)/23)/27
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_251():
    f = (2*x + 1)*(4*x**2 + 3*x + 1)/(3*x**2 - x + 2)**(sympy.S(5)/2)
    F = -(15556 - 17160*x)/(14283*sqrt(3*x**2 - x + 2)) - (734*x + 146)/(621*(3*x**2 - x + 2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_252():
    f = (4*x**2 + 3*x + 1)/((2*x + 1)*(3*x**2 - x + 2)**(sympy.S(5)/2))
    F = -(2764 - 54672*x)/(268203*sqrt(3*x**2 - x + 2)) + (154*x - 202)/(897*(3*x**2 - x + 2)**(sympy.S(3)/2)) - 8*sqrt(13)*atanh(sqrt(13)*(9 - 8*x)/(26*sqrt(3*x**2 - x + 2)))/2197
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_253():
    f = (4*x**2 + 3*x + 1)/((2*x + 1)**2*(3*x**2 - x + 2)**(sympy.S(5)/2))
    F = -(394 - 1674*x)/(11661*(3*x**2 - x + 2)**(sympy.S(3)/2)) - (20184 - 159192*x)/(1162213*sqrt(3*x**2 - x + 2)) - 56*sqrt(13)*atanh(sqrt(13)*(9 - 8*x)/(26*sqrt(3*x**2 - x + 2)))/28561 - 16*sqrt(3*x**2 - x + 2)/(4394*x + 2197)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_254():
    f = (4*x**2 + 3*x + 1)/((2*x + 1)**3*(3*x**2 - x + 2)**(sympy.S(5)/2))
    F = (7386*x + 4726)/(151593*(3*x**2 - x + 2)**(sympy.S(3)/2)) + (1242312*x + 309252)/(15108769*sqrt(3*x**2 - x + 2)) - 2084*sqrt(13)*atanh(sqrt(13)*(9 - 8*x)/(26*sqrt(3*x**2 - x + 2)))/371293 - 144*sqrt(3*x**2 - x + 2)/(57122*x + 28561) - 8*sqrt(3*x**2 - x + 2)/(2197*(2*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_255():
    f = (d + e*x + f*x**2)/((g + h*x)*(b*g*h + b*h**2*x - c*g**2 + c*h**2*x**2)**(sympy.S(3)/2))
    F = (2*d*h**2 - 2*e*g*h + 2*f*g**2)/(3*h**3*(g + h*x)*(-b*h + 2*c*g)*sqrt(b*h**2*x + c*h**2*x**2 - g*(-b*h + c*g))) - f/(c*h**3*sqrt(b*h**2*x + c*h**2*x**2 - g*(-b*h + c*g))) + (b + 2*c*x)*(-3*b**2*f*h**2 + 6*b*c*e*h**2 + 4*c**2*(f*g**2 - h*(2*d*h + e*g)))/(3*c*h**2*(-b*h + 2*c*g)**3*sqrt(b*h**2*x + c*h**2*x**2 - g*(-b*h + c*g)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_256():
    f = sqrt(d + e*x)*(A + B*x + C*x**2)*sqrt(a + b*x + c*x**2)
    F = 2*C*(d + e*x)**(sympy.S(3)/2)*(a + b*x + c*x**2)**(sympy.S(3)/2)/(9*c*e) - sqrt(d + e*x)*(a + b*x + c*x**2)**(sympy.S(3)/2)*(-6*B*c*e + 4*C*b*e + 4*C*c*d)/(21*c**2*e) + 2*sqrt(d + e*x)*sqrt(a + b*x + c*x**2)*(8*C*b**3*e**3 - 3*b*c*e**2*(4*B*b*e - C*a*e + C*b*d) + c**3*d*(8*C*d**2 - 3*e*(-7*A*e + 4*B*d)) + 3*c**2*e*(a*e*(-5*B*e + C*d) - b*(-7*A*e**2 - 2*B*d*e + C*d**2)) + 3*c*e*x*(8*C*b**2*e**2 - c**2*(2*C*d**2 - 3*e*(7*A*e + B*d)) - c*e*(12*B*b*e + 7*C*a*e + C*b*d)))/(315*c**3*e**3) - 2*sqrt(2)*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2)*(8*C*b**3*e**3 + 3*b*c*e**2*(-4*B*b*e - 9*C*a*e + C*b*d) - 2*c**3*d*(8*C*d**2 - 3*e*(-7*A*e + 4*B*d)) - 3*c**2*e**2*(-7*A*b*e - 10*B*a*e + B*b*d + 2*C*a*d))*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(315*c**4*e**4*sqrt(d + e*x)*sqrt(a + b*x + c*x**2)) + sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(d + e*x)*sqrt(-4*a*c + b**2)*(-5*c*e*(-b*e + 2*c*d)*(6*C*b**2*d*e + b*(2*C*a*e**2 - c*d*(9*B*e + C*d)) + c*e*(21*A*c*d - 3*B*a*e - 5*C*a*d)) + (-2*b**2*e**2 + 8*c**2*d**2 - 3*c*e*(-2*a*e + b*d))*(8*C*b**2*e**2 - c**2*(2*C*d**2 - 3*e*(7*A*e + B*d)) - c*e*(12*B*b*e + 7*C*a*e + C*b*d)))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(315*c**4*e**4*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_257():
    f = (A + B*x + C*x**2)*sqrt(a + b*x + c*x**2)/sqrt(d + e*x)
    F = 2*C*sqrt(d + e*x)*(a + b*x + c*x**2)**(sympy.S(3)/2)/(7*c*e) - 2*sqrt(d + e*x)*sqrt(a + b*x + c*x**2)*(3*c*e*x*(-7*B*c*e + 4*C*b*e + 6*C*c*d) + 5*c*e*(-7*A*c*e + C*a*e + 3*C*b*d) - (-b*e + 4*c*d)*(-7*B*c*e + 4*C*b*e + 6*C*c*d))/(105*c**2*e**3) + 2*sqrt(2)*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2)*(4*C*b**2*e**2 + c**2*(48*C*d**2 - 14*e*(-5*A*e + 4*B*d)) + c*e*(-7*B*b*e - 10*C*a*e + 8*C*b*d))*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(105*c**3*e**4*sqrt(d + e*x)*sqrt(a + b*x + c*x**2)) + sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(d + e*x)*sqrt(-4*a*c + b**2)*(5*c*e*(-b*e + 2*c*d)*(-7*A*c*e + C*a*e + 3*C*b*d) - (-2*b**2*e**2 + 8*c**2*d**2 - 3*c*e*(-2*a*e + b*d))*(-7*B*c*e + 4*C*b*e + 6*C*c*d))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(105*c**3*e**4*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_258():
    f = (A + B*x + C*x**2)*sqrt(a + b*x + c*x**2)/(d + e*x)**(sympy.S(3)/2)
    F = -(2*C*d**2 - 2*e*(-A*e + B*d))*(a + b*x + c*x**2)**(sympy.S(3)/2)/(e*sqrt(d + e*x)*(a*e**2 - b*d*e + c*d**2)) - 2*sqrt(d + e*x)*sqrt(a + b*x + c*x**2)*(C*b*e**2*(-a*e + b*d) + c**2*d*(24*C*d**2 - 5*e*(-3*A*e + 4*B*d)) + 3*c*e**2*x*(-5*A*c*e + 5*B*c*d - C*a*e + C*b*d - 6*C*c*d**2/e) + c*e*(a*e*(-5*B*e + 9*C*d) - 5*b*(3*A*e**2 - 4*B*d*e + 5*C*d**2)))/(15*c*e**3*(a*e**2 - b*d*e + c*d**2)) + 2*sqrt(2)*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(C*b*e**2*(-a*e + b*d) - 2*c**2*d*(24*C*d**2 - 5*e*(-3*A*e + 4*B*d)) - c*e*(2*a*e*(-5*B*e + 9*C*d) - b*(32*C*d**2 - 5*e*(-3*A*e + 5*B*d))))*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(15*c**2*e**4*sqrt(d + e*x)*sqrt(a + b*x + c*x**2)) - sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(d + e*x)*sqrt(-4*a*c + b**2)*(2*C*b**2*e**2 - c**2*(48*C*d**2 - 10*e*(-3*A*e + 4*B*d)) + c*e*(-5*B*b*e - 6*C*a*e + 8*C*b*d))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(15*c**2*e**4*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_259():
    f = (A + B*x + C*x**2)*sqrt(a + b*x + c*x**2)/(d + e*x)**(sympy.S(5)/2)
    F = -(2*C*d**2 - 2*e*(-A*e + B*d))*(a + b*x + c*x**2)**(sympy.S(3)/2)/(3*e*(d + e*x)**(sympy.S(3)/2)*(a*e**2 - b*d*e + c*d**2)) - sqrt(a + b*x + c*x**2)*(-2*c*d*(8*C*d**2 - e*(-A*e + 4*B*d)) + 2*e**2*x*(-A*c*e + B*c*d - C*a*e + C*b*d - 2*C*c*d**2/e) + 2*e*(-3*B*e + 7*C*d)*(-a*e + b*d))/(3*e**3*sqrt(d + e*x)*(a*e**2 - b*d*e + c*d**2)) + sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(d + e*x)*sqrt(-4*a*c + b**2)*(6*c*(b*d*(-B*e + C*d) + e*(A*c*d + B*a*e - C*a*d)) + (-b*e + 8*c*d)*(-A*c*e + B*c*d - C*a*e + C*b*d - 2*C*c*d**2/e))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(3*c*e**3*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2)) - 2*sqrt(2)*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(-2*c*(8*C*d**2 - e*(-A*e + 4*B*d)) + e*(-3*B*b*e - 2*C*a*e + 8*C*b*d))*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(3*c*e**4*sqrt(d + e*x)*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_260():
    f = (A + B*x + C*x**2)*sqrt(a + b*x + c*x**2)/(d + e*x)**(sympy.S(7)/2)
    F = -(2*C*d**2 - 2*e*(-A*e + B*d))*(a + b*x + c*x**2)**(sympy.S(3)/2)/(5*e*(d + e*x)**(sympy.S(5)/2)*(a*e**2 - b*d*e + c*d**2)) - sqrt(a + b*x + c*x**2)*(2*c**2*d**3*(24*C*d**2 - e*(A*e + 4*B*d)) - 2*c*d*e*(-a*e*(7*A*e**2 - 7*B*d*e + 37*C*d**2) + b*d*(A*e**2 - 6*B*d*e + 41*C*d**2)) + 2*e**2*(15*C*b**2*d**3 + 5*a**2*e**2*(B*e + C*d) - a*b*e*(2*A*e**2 + 3*B*d*e + 22*C*d**2)) + 2*e*x*(5*c**2*d**2*(6*C*d**2 - e*(A*e + B*d)) - c*e*(-a*e*(3*A*e**2 - 13*B*d*e + 53*C*d**2) + 5*b*d*(-A*e**2 - 2*B*d*e + 11*C*d**2)) + e**2*(15*C*a**2*e**2 - 5*a*b*e*(-B*e + 8*C*d) + b**2*(-2*A*e**2 - 3*B*d*e + 23*C*d**2))))/(15*e**3*(d + e*x)**(sympy.S(3)/2)*(a*e**2 - b*d*e + c*d**2)**2) + sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(d + e*x)*sqrt(-4*a*c + b**2)*(2*c**2*d**2*(24*C*d**2 - e*(A*e + 4*B*d)) - c*e*(-2*a*e*(3*A*e**2 - 8*B*d*e + 43*C*d**2) + b*d*(-2*A*e**2 - 13*B*d*e + 88*C*d**2)) + e**2*(30*C*a**2*e**2 - 5*a*b*e*(-B*e + 14*C*d) + b**2*(-2*A*e**2 - 3*B*d*e + 38*C*d**2)))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(15*e**4*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2)**2) - 2*sqrt(2)*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(15*C*b*e**2*(-a*e + b*d) + 2*c**2*d*(24*C*d**2 - e*(A*e + 4*B*d)) + c*e*(10*a*e*(-B*e + 5*C*d) - b*(-A*e**2 - 9*B*d*e + 64*C*d**2)))*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(15*c*e**4*sqrt(d + e*x)*sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_261():
    f = (A + B*x + C*x**2)*sqrt(a + b*x + c*x**2)/(d + e*x)**(sympy.S(9)/2)
    F = -(2*C*d**2 - 2*e*(-A*e + B*d))*(a + b*x + c*x**2)**(sympy.S(3)/2)/(7*e*(d + e*x)**(sympy.S(7)/2)*(a*e**2 - b*d*e + c*d**2)) + sqrt(a + b*x + c*x**2)*(-2*b*e**3*(35*C*a**2*e**2 - 14*a*b*e*(B*e + 3*C*d) + b**2*(8*A*e**2 + 6*B*d*e + 15*C*d**2)) + 4*c**3*d**3*(24*C*d**2 + e*(3*A*e + 4*B*d)) + 2*c**2*d*e*(2*a*e*(69*C*d**2 + e*(-29*A*e + 15*B*d)) - b*d*(128*C*d**2 + e*(9*A*e + 19*B*d))) + 2*c*e**2*(14*a**2*e**2*(-3*B*e + 11*C*d) - a*b*e*(237*C*d**2 + e*(-29*A*e + B*d)) + b**2*d*(103*C*d**2 + e*(19*A*e + 9*B*d))))/(105*e**3*sqrt(d + e*x)*(a*e**2 - b*d*e + c*d**2)**3) - sqrt(a + b*x + c*x**2)*(2*c**2*d**3*(24*C*d**2 + e*(3*A*e + 4*B*d)) - 2*c*d*e*(-a*e*(19*A*e**2 + 9*B*d*e + 33*C*d**2) + b*d*(15*A*e**2 + 6*B*d*e + 43*C*d**2)) - 2*e**2*(7*a**2*e**2*(-3*B*e + C*d) + a*b*e*(12*A*e**2 + 23*B*d*e + 12*C*d**2) - b**2*d*(8*A*e**2 + 6*B*d*e + 15*C*d**2)) + 2*e*x*(7*c**2*d**2*(6*C*d**2 + e*(-3*A*e + B*d)) + c*e*(a*e*(-5*A*e**2 - 9*B*d*e + 93*C*d**2) - b*(-21*A*d*e**2 + 91*C*d**3)) + e**2*(35*C*a**2*e**2 - 7*a*b*e*(-B*e + 12*C*d) + b**2*(-4*A*e**2 - 3*B*d*e + 45*C*d**2))))/(105*e**3*(d + e*x)**(sympy.S(5)/2)*(a*e**2 - b*d*e + c*d**2)**2) + 2*sqrt(2)*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(2*c**2*d**2*(24*C*d**2 + e*(3*A*e + 4*B*d)) + c*e*(2*a*e*(51*C*d**2 + e*(-5*A*e + 12*B*d)) - b*d*(104*C*d**2 + 3*e*(2*A*e + 5*B*d))) + e**2*(70*C*a**2*e**2 - 7*a*b*e*(B*e + 18*C*d) + b**2*(60*C*d**2 + e*(4*A*e + 3*B*d))))*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(105*e**4*sqrt(d + e*x)*sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2)**2) - sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(d + e*x)*sqrt(-4*a*c + b**2)*(-b*e**3*(35*C*a**2*e**2 - 14*a*b*e*(B*e + 3*C*d) + b**2*(8*A*e**2 + 6*B*d*e + 15*C*d**2)) + 2*c**3*d**3*(24*C*d**2 + e*(3*A*e + 4*B*d)) + c**2*d*e*(2*a*e*(69*C*d**2 + e*(-29*A*e + 15*B*d)) - b*d*(128*C*d**2 + e*(9*A*e + 19*B*d))) + c*e**2*(14*a**2*e**2*(-3*B*e + 11*C*d) - a*b*e*(237*C*d**2 + e*(-29*A*e + B*d)) + b**2*d*(103*C*d**2 + e*(19*A*e + 9*B*d))))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(105*e**4*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_262():
    f = (A + B*x + C*x**2)*sqrt(a + b*x + c*x**2)/(d + e*x)**(sympy.S(11)/2)
    F = -(2*C*d**2 - 2*e*(-A*e + B*d))*(a + b*x + c*x**2)**(sympy.S(3)/2)/(9*e*(d + e*x)**(sympy.S(9)/2)*(a*e**2 - b*d*e + c*d**2)) + sqrt(a + b*x + c*x**2)*(4*b**2*e**4*(21*C*a**2*e**2 - 6*a*b*e*(2*B*e + 3*C*d) + b**2*(8*A*e**2 + 4*B*d*e + 5*C*d**2)) + 4*c**4*d**4*(8*C*d**2 + e*(5*A*e + 4*B*d)) + 2*c**3*d**2*e*(6*a*e*(-34*A*e**2 + 8*B*d*e + 11*C*d**2) - b*d*(56*C*d**2 + 5*e*(4*A*e + 5*B*d))) - 12*c**2*e**2*(-a**2*e**2*(7*A*e**2 - 36*B*d*e + 30*C*d**2) + a*b*d*e*(-34*A*e**2 - 5*B*d*e + 30*C*d**2) - b**2*d**2*(11*A*e**2 + 3*B*d*e + 11*C*d**2)) - 2*c*e**3*(126*C*a**3*e**3 - 3*a**2*b*e**2*(29*B*e + 12*C*d) - 6*a*b**2*e*(-12*A*e**2 + 7*B*d*e + 5*C*d**2) + b**3*d*(56*A*e**2 + 25*B*d*e + 20*C*d**2)))/(315*e**3*sqrt(d + e*x)*(a*e**2 - b*d*e + c*d**2)**4) + sqrt(a + b*x + c*x**2)*(-2*b*e**3*(21*C*a**2*e**2 - 6*a*b*e*(2*B*e + 3*C*d) + b**2*(8*A*e**2 + 4*B*d*e + 5*C*d**2)) + 4*c**3*d**3*(8*C*d**2 + e*(5*A*e + 4*B*d)) + 6*c**2*d*e*(2*a*e*(-9*A*e**2 + 7*B*d*e + 9*C*d**2) - b*d*(5*A*e**2 + 7*B*d*e + 16*C*d**2)) + 6*c*e**2*(2*a**2*e**2*(-5*B*e + 17*C*d) - a*b*e*(-9*A*e**2 + 5*B*d*e + 41*C*d**2) + b**2*d*(7*A*e**2 + 3*B*d*e + 15*C*d**2)))/(315*e**3*(d + e*x)**(sympy.S(3)/2)*(a*e**2 - b*d*e + c*d**2)**3) - sqrt(a + b*x + c*x**2)*(2*c**2*d**3*(8*C*d**2 + e*(5*A*e + 4*B*d)) - 2*c*d*e*(-a*e*(13*A*e**2 + 11*B*d*e + 7*C*d**2) + 3*b*d*(5*A*e**2 + 2*B*d*e + 5*C*d**2)) - 2*e**2*(3*a**2*e**2*(-5*B*e + 3*C*d) - a*b*e*(-10*A*e**2 - 17*B*d*e + 2*C*d**2) - b**2*d*(8*A*e**2 + 4*B*d*e + 5*C*d**2)) + 2*e*x*(3*c**2*d**2*(6*C*d**2 + e*(-5*A*e + 3*B*d)) + c*e*(a*e*(-7*A*e**2 + B*d*e + 47*C*d**2) - 3*b*d*(-5*A*e**2 + 2*B*d*e + 15*C*d**2)) + e**2*(21*C*a**2*e**2 - 3*a*b*e*(-B*e + 16*C*d) + b**2*(25*C*d**2 - e*(2*A*e + B*d)))))/(105*e**3*(d + e*x)**(sympy.S(7)/2)*(a*e**2 - b*d*e + c*d**2)**2) + 2*sqrt(2)*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(-b*e**3*(21*C*a**2*e**2 - 6*a*b*e*(2*B*e + 3*C*d) + b**2*(8*A*e**2 + 4*B*d*e + 5*C*d**2)) + 2*c**3*d**3*(8*C*d**2 + e*(5*A*e + 4*B*d)) + 3*c**2*d*e*(2*a*e*(-9*A*e**2 + 7*B*d*e + 9*C*d**2) - b*d*(5*A*e**2 + 7*B*d*e + 16*C*d**2)) + 3*c*e**2*(2*a**2*e**2*(-5*B*e + 17*C*d) - a*b*e*(-9*A*e**2 + 5*B*d*e + 41*C*d**2) + b**2*d*(7*A*e**2 + 3*B*d*e + 15*C*d**2)))*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(315*e**4*sqrt(d + e*x)*sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2)**3) - sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(d + e*x)*sqrt(-4*a*c + b**2)*(2*b**2*e**4*(21*C*a**2*e**2 - 6*a*b*e*(2*B*e + 3*C*d) + b**2*(8*A*e**2 + 4*B*d*e + 5*C*d**2)) + 2*c**4*d**4*(8*C*d**2 + e*(5*A*e + 4*B*d)) + c**3*d**2*e*(6*a*e*(-34*A*e**2 + 8*B*d*e + 11*C*d**2) - b*d*(56*C*d**2 + 5*e*(4*A*e + 5*B*d))) - 6*c**2*e**2*(-a**2*e**2*(7*A*e**2 - 36*B*d*e + 30*C*d**2) + a*b*d*e*(-34*A*e**2 - 5*B*d*e + 30*C*d**2) - b**2*d**2*(11*A*e**2 + 3*B*d*e + 11*C*d**2)) - c*e**3*(126*C*a**3*e**3 - 3*a**2*b*e**2*(29*B*e + 12*C*d) - 6*a*b**2*e*(-12*A*e**2 + 7*B*d*e + 5*C*d**2) + b**3*d*(56*A*e**2 + 25*B*d*e + 20*C*d**2)))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(315*e**4*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_263():
    f = (d + e*x)**(sympy.S(3)/2)*(A + B*x + C*x**2)/sqrt(a + b*x + c*x**2)
    F = 2*C*(d + e*x)**(sympy.S(5)/2)*sqrt(a + b*x + c*x**2)/(7*c*e) - (d + e*x)**(sympy.S(3)/2)*sqrt(a + b*x + c*x**2)*(-14*B*c*e + 12*C*b*e + 4*C*c*d)/(35*c**2*e) + sqrt(d + e*x)*sqrt(a + b*x + c*x**2)*(48*C*b**2*e**2 - 2*c**2*(6*C*d**2 - 7*e*(5*A*e + 3*B*d)) - 2*c*e*(28*B*b*e + 25*C*a*e + 15*C*b*d))/(105*c**3*e) - 2*sqrt(2)*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2)*(24*C*b**2*e**2 - c**2*(6*C*d**2 - 7*e*(5*A*e + 3*B*d)) - c*e*(28*B*b*e + 25*C*a*e + 15*C*b*d))*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(105*c**4*e**2*sqrt(d + e*x)*sqrt(a + b*x + c*x**2)) - sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(d + e*x)*sqrt(-4*a*c + b**2)*(48*C*b**3*e**3 - 8*b*c*e**2*(7*B*b*e + 13*C*a*e + 9*C*b*d) + c**3*d*(6*C*d**2 - 7*e*(20*A*e + 3*B*d)) + c**2*e*(a*e*(63*B*e + 82*C*d) + b*(70*A*e**2 + 91*B*d*e + 12*C*d**2)))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(105*c**4*e**2*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_264():
    f = sqrt(d + e*x)*(A + B*x + C*x**2)/sqrt(a + b*x + c*x**2)
    F = 2*C*(d + e*x)**(sympy.S(3)/2)*sqrt(a + b*x + c*x**2)/(5*c*e) - sqrt(d + e*x)*sqrt(a + b*x + c*x**2)*(-10*B*c*e + 8*C*b*e + 4*C*c*d)/(15*c**2*e) + 2*sqrt(2)*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2)*(-5*B*c*e + 4*C*b*e + 2*C*c*d)*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(15*c**3*e**2*sqrt(d + e*x)*sqrt(a + b*x + c*x**2)) + sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(d + e*x)*sqrt(-4*a*c + b**2)*(8*C*b**2*e**2 - c**2*(2*C*d**2 - 5*e*(3*A*e + B*d)) - c*e*(10*B*b*e + 9*C*a*e + 3*C*b*d))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(15*c**3*e**2*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_265():
    f = (A + B*x + C*x**2)/(sqrt(d + e*x)*sqrt(a + b*x + c*x**2))
    F = 2*C*sqrt(d + e*x)*sqrt(a + b*x + c*x**2)/(3*c*e) + 2*sqrt(2)*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(C*e*(-a*e + b*d) + c*(2*C*d**2 - 3*e*(-A*e + B*d)))*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(3*c**2*e**2*sqrt(d + e*x)*sqrt(a + b*x + c*x**2)) - sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(d + e*x)*sqrt(-4*a*c + b**2)*(-3*B*c*e + 2*C*b*e + 2*C*c*d)*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(3*c**2*e**2*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_266():
    f = (A + B*x + C*x**2)/((d + e*x)**(sympy.S(3)/2)*sqrt(a + b*x + c*x**2))
    F = -(2*C*d**2 - 2*e*(-A*e + B*d))*sqrt(a + b*x + c*x**2)/(e*sqrt(d + e*x)*(a*e**2 - b*d*e + c*d**2)) - 2*sqrt(2)*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*(-B*e + 2*C*d)*sqrt(-4*a*c + b**2)*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(c*e**2*sqrt(d + e*x)*sqrt(a + b*x + c*x**2)) - sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(d + e*x)*sqrt(-4*a*c + b**2)*(C*e*(-a*e + b*d) - c*(2*C*d**2 - e*(-A*e + B*d)))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(c*e**2*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_267():
    f = (A + B*x + C*x**2)/((d + e*x)**(sympy.S(5)/2)*sqrt(a + b*x + c*x**2))
    F = (2*c*d*(2*C*d**2 + e*(-4*A*e + B*d)) + 2*e*(3*a*e*(-B*e + 2*C*d) - b*(-2*A*e**2 - B*d*e + 4*C*d**2)))*sqrt(a + b*x + c*x**2)/(3*e*sqrt(d + e*x)*(a*e**2 - b*d*e + c*d**2)**2) - (2*C*d**2 - 2*e*(-A*e + B*d))*sqrt(a + b*x + c*x**2)/(3*e*(d + e*x)**(sympy.S(3)/2)*(a*e**2 - b*d*e + c*d**2)) - sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(d + e*x)*sqrt(-4*a*c + b**2)*(c*d*(2*C*d**2 + e*(-4*A*e + B*d)) + e*(3*a*e*(-B*e + 2*C*d) - b*(-2*A*e**2 - B*d*e + 4*C*d**2)))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(3*e**2*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2)**2) - 2*sqrt(2)*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(3*C*e*(-a*e + b*d) - c*(2*C*d**2 + e*(-A*e + B*d)))*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(3*c*e**2*sqrt(d + e*x)*sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_268():
    f = (A + B*x + C*x**2)/((d + e*x)**(sympy.S(7)/2)*sqrt(a + b*x + c*x**2))
    F = sqrt(a + b*x + c*x**2)*(2*c**2*d**2*(2*C*d**2 + e*(-23*A*e + 3*B*d)) - 2*c*e*(-a*e*(9*A*e**2 - 29*B*d*e + 19*C*d**2) + b*d*(-23*A*e**2 - 7*B*d*e + 7*C*d**2)) - 2*e**2*(15*C*a**2*e**2 - 10*a*b*e*(B*e + C*d) + b**2*(8*A*e**2 + 2*B*d*e + 3*C*d**2)))/(15*e*sqrt(d + e*x)*(a*e**2 - b*d*e + c*d**2)**3) + (2*c*d*(2*C*d**2 + e*(-8*A*e + 3*B*d)) + 2*e*(5*a*e*(-B*e + 2*C*d) - b*(-4*A*e**2 - B*d*e + 6*C*d**2)))*sqrt(a + b*x + c*x**2)/(15*e*(d + e*x)**(sympy.S(3)/2)*(a*e**2 - b*d*e + c*d**2)**2) - (2*C*d**2 - 2*e*(-A*e + B*d))*sqrt(a + b*x + c*x**2)/(5*e*(d + e*x)**(sympy.S(5)/2)*(a*e**2 - b*d*e + c*d**2)) + 2*sqrt(2)*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(c*d*(2*C*d**2 + e*(-8*A*e + 3*B*d)) + e*(5*a*e*(-B*e + 2*C*d) - b*(-4*A*e**2 - B*d*e + 6*C*d**2)))*elliptic_f(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(15*e**2*sqrt(d + e*x)*sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2)**2) - sqrt(2)*sqrt(-c*(a + b*x + c*x**2)/(-4*a*c + b**2))*sqrt(d + e*x)*sqrt(-4*a*c + b**2)*(c**2*d**2*(2*C*d**2 + e*(-23*A*e + 3*B*d)) - c*e*(-a*e*(9*A*e**2 - 29*B*d*e + 19*C*d**2) + b*d*(-23*A*e**2 - 7*B*d*e + 7*C*d**2)) - e**2*(15*C*a**2*e**2 - 10*a*b*e*(B*e + C*d) + b**2*(8*A*e**2 + 2*B*d*e + 3*C*d**2)))*elliptic_e(asin(sqrt(2)*sqrt((b + 2*c*x + sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))/2), -2*e*sqrt(-4*a*c + b**2)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))/(15*e**2*sqrt(c*(d + e*x)/(2*c*d - e*(b + sqrt(-4*a*c + b**2))))*sqrt(a + b*x + c*x**2)*(a*e**2 - b*d*e + c*d**2)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_269():
    f = (g + h*x)**m*(a + b*x + c*x**2)**p*(d + e*x + f*x**2)
    F = f*(g + h*x)**(m + 1)*(a + b*x + c*x**2)**(p + 1)/(c*h*(m + 2*p + 3)) + (g + h*x)**(m + 1)*(c*(2*f*g**2*(p + 1) - h*(-d*h + e*g)*(m + 2*p + 3)) + f*h*(m + 1)*(-a*h + b*g))*(a + b*x + c*x**2)**p*appellf1(m + 1, -p, -p, m + 2, 2*c*(g + h*x)/(2*c*g - h*(b - sqrt(-4*a*c + b**2))), 2*c*(g + h*x)/(2*c*g - h*(b + sqrt(-4*a*c + b**2))))/(c*h**3*(m + 1)*(-2*c*(g + h*x)/(2*c*g - h*(b - sqrt(-4*a*c + b**2))) + 1)**p*(-2*c*(g + h*x)/(2*c*g - h*(b + sqrt(-4*a*c + b**2))) + 1)**p*(m + 2*p + 3)) - (g + h*x)**(m + 2)*(b*f*h*(m + p + 2) + c*(-e*h*(m + 2*p + 3) + 2*f*g*(p + 1)))*(a + b*x + c*x**2)**p*appellf1(m + 2, -p, -p, m + 3, 2*c*(g + h*x)/(2*c*g - h*(b - sqrt(-4*a*c + b**2))), 2*c*(g + h*x)/(2*c*g - h*(b + sqrt(-4*a*c + b**2))))/(c*h**3*(m + 2)*(-2*c*(g + h*x)/(2*c*g - h*(b - sqrt(-4*a*c + b**2))) + 1)**p*(-2*c*(g + h*x)/(2*c*g - h*(b + sqrt(-4*a*c + b**2))) + 1)**p*(m + 2*p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_270():
    f = (g + h*x)**m*sqrt(a + b*x + c*x**2)*(d + e*x + f*x**2)
    F = f*(g + h*x)**(m + 1)*(a + b*x + c*x**2)**(sympy.S(3)/2)/(c*h*(m + 4)) + (g + h*x)**(m + 1)*(c*(3*f*g**2 - h*(m + 4)*(-d*h + e*g)) + f*h*(m + 1)*(-a*h + b*g))*sqrt(a + b*x + c*x**2)*appellf1(m + 1, sympy.S(-1)/2, sympy.S(-1)/2, m + 2, 2*c*(g + h*x)/(2*c*g - h*(b - sqrt(-4*a*c + b**2))), 2*c*(g + h*x)/(2*c*g - h*(b + sqrt(-4*a*c + b**2))))/(c*h**3*(m + 1)*(m + 4)*sqrt(-2*c*(g + h*x)/(2*c*g - h*(b - sqrt(-4*a*c + b**2))) + 1)*sqrt(-2*c*(g + h*x)/(2*c*g - h*(b + sqrt(-4*a*c + b**2))) + 1)) - (g + h*x)**(m + 2)*(b*f*h*(2*m + 5) + c*(-2*e*h*(m + 4) + 6*f*g))*sqrt(a + b*x + c*x**2)*appellf1(m + 2, sympy.S(-1)/2, sympy.S(-1)/2, m + 3, 2*c*(g + h*x)/(2*c*g - h*(b - sqrt(-4*a*c + b**2))), 2*c*(g + h*x)/(2*c*g - h*(b + sqrt(-4*a*c + b**2))))/(2*c*h**3*(m + 2)*(m + 4)*sqrt(-2*c*(g + h*x)/(2*c*g - h*(b - sqrt(-4*a*c + b**2))) + 1)*sqrt(-2*c*(g + h*x)/(2*c*g - h*(b + sqrt(-4*a*c + b**2))) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_271():
    f = (g + h*x)**(-2*p - 3)*(a + b*x + c*x**2)**p*(d + e*x + f*x**2)
    F = -f*(a + b*x + c*x**2)**p*appellf1(-2*p, -p, -p, 1 - 2*p, 2*c*(g + h*x)/(2*c*g - h*(b - sqrt(-4*a*c + b**2))), 2*c*(g + h*x)/(2*c*g - h*(b + sqrt(-4*a*c + b**2))))/(2*h**3*p*(g + h*x)**(2*p)*(-2*c*(g + h*x)/(2*c*g - h*(b - sqrt(-4*a*c + b**2))) + 1)**p*(-2*c*(g + h*x)/(2*c*g - h*(b + sqrt(-4*a*c + b**2))) + 1)**p) - (g + h*x)**(-2*p - 2)*(f*g**2 - h*(-d*h + e*g))*(a + b*x + c*x**2)**(p + 1)/(2*h*(p + 1)*(a*h**2 - b*g*h + c*g**2)) - (g + h*x)**(-2*p - 1)*(2*c*(-d*g*h**2 + f*g**3) + h*(2*a*h*(-e*h + 2*f*g) - b*(-d*h**2 - e*g*h + 3*f*g**2)))*(a + b*x + c*x**2)**p*(b + 2*c*x - sqrt(-4*a*c + b**2))*hyper((-p, -2*p - 1), (-2*p,), -4*c*(g + h*x)*sqrt(-4*a*c + b**2)/((2*c*g - h*(b + sqrt(-4*a*c + b**2)))*(b + 2*c*x - sqrt(-4*a*c + b**2))))/(2*h**2*((2*c*g - h*(b - sqrt(-4*a*c + b**2)))*(b + 2*c*x + sqrt(-4*a*c + b**2))/((2*c*g - h*(b + sqrt(-4*a*c + b**2)))*(b + 2*c*x - sqrt(-4*a*c + b**2))))**p*(2*p + 1)*(2*c*g - h*(b - sqrt(-4*a*c + b**2)))*(a*h**2 - b*g*h + c*g**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_272():
    f = (d + f*x**2)**p*(2*b*f**2*x*(2*p + 3) + 2*c*d*f + 2*c*f**2*x**2*(2*p + 3))
    F = b*f*(d + f*x**2)**(p + 1)*(2*p + 3)/(p + 1) + 2*c*f*x*(d + f*x**2)**(p + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_273():
    f = (d + e*x + f*x**2)**p*(2*c*d*f - c*e**2*p - 2*c*e**2 + 2*c*f**2*x**2*(2*p + 3))
    F = -c*e*(p + 2)*(d + e*x + f*x**2)**(p + 1)/(p + 1) + 2*c*f*x*(d + e*x + f*x**2)**(p + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_274():
    f = (d + e*x + f*x**2)**p*(2*b*e*f*p + 3*b*e*f + 2*b*f**2*x*(2*p + 3) + 2*c*d*f - c*e**2*p - 2*c*e**2 + 2*c*f**2*x**2*(2*p + 3))
    F = 2*c*f*x*(d + e*x + f*x**2)**(p + 1) - (-b*f*(2*p + 3) + c*e*(p + 2))*(d + e*x + f*x**2)**(p + 1)/(p + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_275():
    f = (d + e*x)**3*(a + b*x + c*x**2)**5*(17*c*e**2*x**3 + d*(5*a*e + 6*b*d) + e*x**2*(11*b*e + 29*c*d) + x*(5*a*e**2 + 17*b*d*e + 12*c*d**2))
    F = (d + e*x)**5*(a + b*x + c*x**2)**6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_276():
    f = (x**3 + x**2)/(x**2 + x - 2)
    F = x**2/2 + 2*log(1 - x)/3 + 4*log(x + 2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_277():
    f = x**2*(d + e*x + f*x**2 + g*x**3)/sqrt(a + b*x + c*x**2)
    F = g*x**4*sqrt(a + b*x + c*x**2)/(5*c) + x**3*(-9*b*g + 10*c*f)*sqrt(a + b*x + c*x**2)/(40*c**2) + x**2*sqrt(a + b*x + c*x**2)*(-64*a*c*g + 63*b**2*g - 70*b*c*f + 80*c**2*e)/(240*c**3) - sqrt(a + b*x + c*x**2)*(256*a*c**2*(-4*a*g + 5*c*e) - 945*b**4*g + 1050*b**3*c*f - 60*b**2*c*(-49*a*g + 20*c*e) + 40*b*c**2*(-55*a*f + 36*c*d) - 2*c*x*(-315*b**3*g + 14*b*c*(46*a*g + 25*b*f) + 480*c**3*d - 40*c**2*(9*a*f + 10*b*e)))/(1920*c**5) + (48*a*b*c**2*(-5*a*g + 4*c*e) - 32*a*c**3*(-3*a*f + 4*c*d) - 63*b**5*g + 70*b**4*c*f - 40*b**3*c*(-7*a*g + 2*c*e) + 48*b**2*c**2*(-5*a*f + 2*c*d))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(256*c**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_278():
    f = x*(d + e*x + f*x**2 + g*x**3)/sqrt(a + b*x + c*x**2)
    F = g*x**3*sqrt(a + b*x + c*x**2)/(4*c) + x**2*(-7*b*g + 8*c*f)*sqrt(a + b*x + c*x**2)/(24*c**2) + sqrt(a + b*x + c*x**2)*(-105*b**3*g + 20*b*c*(11*a*g + 6*b*f) + 192*c**3*d - 16*c**2*(8*a*f + 9*b*e) + 2*c*x*(-36*a*c*g + 35*b**2*g - 40*b*c*f + 48*c**2*e))/(192*c**4) - (16*a*c**2*(-3*a*g + 4*c*e) - 35*b**4*g + 40*b**3*c*f - 24*b**2*c*(-5*a*g + 2*c*e) + 32*b*c**2*(-3*a*f + 2*c*d))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(128*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_279():
    f = (d + e*x + f*x**2 + g*x**3)/sqrt(a + b*x + c*x**2)
    F = g*x**2*sqrt(a + b*x + c*x**2)/(3*c) + x*(-5*b*g + 6*c*f)*sqrt(a + b*x + c*x**2)/(12*c**2) + sqrt(a + b*x + c*x**2)*(-16*a*c*g + 15*b**2*g - 18*b*c*f + 24*c**2*e)/(24*c**3) + (-5*b**3*g + 6*b*c*(2*a*g + b*f) + 16*c**3*d - 8*c**2*(a*f + b*e))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(16*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_280():
    f = (d + e*x + f*x**2 + g*x**3)/(x*sqrt(a + b*x + c*x**2))
    F = g*x*sqrt(a + b*x + c*x**2)/(2*c) + (-3*b*g + 4*c*f)*sqrt(a + b*x + c*x**2)/(4*c**2) + (3*b**2*g + 8*c**2*e - 4*c*(a*g + b*f))*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(8*c**(sympy.S(5)/2)) - d*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_281():
    f = (d + e*x + f*x**2 + g*x**3)/(x**2*sqrt(a + b*x + c*x**2))
    F = g*sqrt(a + b*x + c*x**2)/c + (-b*g + 2*c*f)*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/(2*c**(sympy.S(3)/2)) - d*sqrt(a + b*x + c*x**2)/(a*x) + (-2*a*e + b*d)*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(2*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_282():
    f = (d + e*x + f*x**2 + g*x**3)/(x**3*sqrt(a + b*x + c*x**2))
    F = g*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/sqrt(c) - d*sqrt(a + b*x + c*x**2)/(2*a*x**2) + (-4*a*e + 3*b*d)*sqrt(a + b*x + c*x**2)/(4*a**2*x) - (8*a**2*f - 4*a*b*e - 4*a*c*d + 3*b**2*d)*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(8*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_283():
    f = (d + e*x + f*x**2 + g*x**3)/(x**4*sqrt(a + b*x + c*x**2))
    F = -d*sqrt(a + b*x + c*x**2)/(3*a*x**3) + (-6*a*e + 5*b*d)*sqrt(a + b*x + c*x**2)/(12*a**2*x**2) - sqrt(a + b*x + c*x**2)*(24*a**2*f - 18*a*b*e - 16*a*c*d + 15*b**2*d)/(24*a**3*x) + (8*a**2*(-2*a*g + c*e) - 6*a*b**2*e - 4*a*b*(-2*a*f + 3*c*d) + 5*b**3*d)*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(16*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_284():
    f = (d + e*x + f*x**2 + g*x**3)/(x**5*sqrt(a + b*x + c*x**2))
    F = -d*sqrt(a + b*x + c*x**2)/(4*a*x**4) + (-8*a*e + 7*b*d)*sqrt(a + b*x + c*x**2)/(24*a**2*x**3) - sqrt(a + b*x + c*x**2)*(48*a**2*f - 40*a*b*e - 36*a*c*d + 35*b**2*d)/(96*a**3*x**2) + sqrt(a + b*x + c*x**2)*(64*a**2*(-3*a*g + 2*c*e) - 120*a*b**2*e - 4*a*b*(-36*a*f + 55*c*d) + 105*b**3*d)/(192*a**4*x) - (32*a**2*b*(-2*a*g + 3*c*e) + 16*a**2*c*(-4*a*f + 3*c*d) - 40*a*b**3*e - 24*a*b**2*(-2*a*f + 5*c*d) + 35*b**4*d)*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(128*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_285():
    f = (d + e*x + f*x**2 + g*x**3)/(x**6*sqrt(a + b*x + c*x**2))
    F = -d*sqrt(a + b*x + c*x**2)/(5*a*x**5) + (-10*a*e + 9*b*d)*sqrt(a + b*x + c*x**2)/(40*a**2*x**4) - sqrt(a + b*x + c*x**2)*(80*a**2*f - 70*a*b*e - 64*a*c*d + 63*b**2*d)/(240*a**3*x**3) + sqrt(a + b*x + c*x**2)*(120*a**2*(-4*a*g + 3*c*e) - 350*a*b**2*e - 4*a*b*(-100*a*f + 161*c*d) + 315*b**3*d)/(960*a**4*x**2) - sqrt(a + b*x + c*x**2)*(40*a**2*b*(-36*a*g + 55*c*e) + 256*a**2*c*(-5*a*f + 4*c*d) - 1050*a*b**3*e - 60*a*b**2*(-20*a*f + 49*c*d) + 945*b**4*d)/(1920*a**5*x) + (-32*a**3*c*(-4*a*g + 3*c*e) + 48*a**2*b**2*(-2*a*g + 5*c*e) + 48*a**2*b*c*(-4*a*f + 5*c*d) - 70*a*b**4*e - 40*a*b**3*(-2*a*f + 7*c*d) + 63*b**5*d)*atanh((2*a + b*x)/(2*sqrt(a)*sqrt(a + b*x + c*x**2)))/(256*a**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_286():
    f = (d + e*x)**3*(5*x**2 + 2*x + 3)*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)
    F = 2*(d + e*x)**10/e**7 - (d + e*x)**9*(120*d + 17*e)/(9*e**7) + (d + e*x)**8*(300*d**2 + 85*d*e + 17*e**2)/(8*e**7) - (d + e*x)**7*(400*d**3 + 170*d**2*e + 68*d*e**2 + 4*e**3)/(7*e**7) + (d + e*x)**6*(300*d**4 + 170*d**3*e + 102*d**2*e**2 + 12*d*e**3 + 21*e**4)/(6*e**7) - (d + e*x)**5*(120*d**5 + 85*d**4*e + 68*d**3*e**2 + 12*d**2*e**3 + 42*d*e**4 - 7*e**5)/(5*e**7) + (d + e*x)**4*(5*d**2 - 2*d*e + 3*e**2)*(4*d**4 + 5*d**3*e + 3*d**2*e**2 - d*e**3 + 2*e**4)/(4*e**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_287():
    f = (d + e*x)**2*(5*x**2 + 2*x + 3)*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)
    F = 6*d**2*x + d*x**2*(7*d + 12*e)/2 + 20*e**2*x**9/9 + e*x**8*(5*d - 17*e/8) + x**7*(20*d**2/7 - 34*d*e/7 + 17*e**2/7) - x**6*(17*d**2/6 - 17*d*e/3 + 2*e**2/3) + x**5*(17*d**2/5 - 8*d*e/5 + 21*e**2/5) - x**4*(d**2 - 21*d*e/2 - 7*e**2/4) + x**3*(7*d**2 + 14*d*e/3 + 2*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_288():
    f = (d + e*x)*(5*x**2 + 2*x + 3)*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)
    F = 6*d*x + 5*e*x**8/2 + x**7*(20*d/7 - 17*e/7) - x**6*(17*d/6 - 17*e/6) + x**5*(17*d/5 - 4*e/5) - x**4*(d - 21*e/4) + x**3*(7*d + 7*e/3) + x**2*(7*d/2 + 3*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_289():
    f = (5*x**2 + 2*x + 3)*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)
    F = 20*x**7/7 - 17*x**6/6 + 17*x**5/5 - x**4 + 7*x**3 + 7*x**2/2 + 6*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_290():
    f = (5*x**2 + 2*x + 3)*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)/(d + e*x)
    F = 10*x**6/(3*e) - x**5*(20*d + 17*e)/(5*e**2) + x**4*(20*d**2 + 17*d*e + 17*e**2)/(4*e**3) - x**3*(20*d**3 + 17*d**2*e + 17*d*e**2 + 4*e**3)/(3*e**4) + x**2*(20*d**4 + 17*d**3*e + 17*d**2*e**2 + 4*d*e**3 + 21*e**4)/(2*e**5) - x*(20*d**5 + 17*d**4*e + 17*d**3*e**2 + 4*d**2*e**3 + 21*d*e**4 - 7*e**5)/e**6 + (5*d**2 - 2*d*e + 3*e**2)*(4*d**4 + 5*d**3*e + 3*d**2*e**2 - d*e**3 + 2*e**4)*log(d + e*x)/e**7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_291():
    f = (5*x**2 + 2*x + 3)*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)/(d + e*x)**2
    F = 4*x**5/e**2 - x**4*(40*d + 17*e)/(4*e**3) + x**3*(60*d**2 + 34*d*e + 17*e**2)/(3*e**4) - x**2*(80*d**3 + 51*d**2*e + 34*d*e**2 + 4*e**3)/(2*e**5) + x*(100*d**4 + 68*d**3*e + 51*d**2*e**2 + 8*d*e**3 + 21*e**4)/e**6 - (120*d**5 + 85*d**4*e + 68*d**3*e**2 + 12*d**2*e**3 + 42*d*e**4 - 7*e**5)*log(d + e*x)/e**7 - (5*d**2 - 2*d*e + 3*e**2)*(4*d**4 + 5*d**3*e + 3*d**2*e**2 - d*e**3 + 2*e**4)/(e**7*(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_292():
    f = (5*x**2 + 2*x + 3)*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)/(d + e*x)**3
    F = 5*x**4/e**3 - x**3*(60*d + 17*e)/(3*e**4) + x**2*(120*d**2 + 51*d*e + 17*e**2)/(2*e**5) - x*(200*d**3 + 102*d**2*e + 51*d*e**2 + 4*e**3)/e**6 + (300*d**4 + 170*d**3*e + 102*d**2*e**2 + 12*d*e**3 + 21*e**4)*log(d + e*x)/e**7 + (120*d**5 + 85*d**4*e + 68*d**3*e**2 + 12*d**2*e**3 + 42*d*e**4 - 7*e**5)/(e**7*(d + e*x)) - (5*d**2 - 2*d*e + 3*e**2)*(4*d**4 + 5*d**3*e + 3*d**2*e**2 - d*e**3 + 2*e**4)/(2*e**7*(d + e*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_293():
    f = (d + e*x)**3*(5*x**2 + 2*x + 3)**2*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)
    F = 25*(d + e*x)**12/(3*e**9) - (d + e*x)**11*(800*d + 45*e)/(11*e**9) + (d + e*x)**10*(2800*d**2 + 315*d*e + 111*e**2)/(10*e**9) - (d + e*x)**9*(5600*d**3 + 945*d**2*e + 666*d*e**2 + 37*e**3)/(9*e**9) + (d + e*x)**8*(7000*d**4 + 1575*d**3*e + 1665*d**2*e**2 + 185*d*e**3 + 148*e**4)/(8*e**9) - (d + e*x)**7*(5600*d**5 + 1575*d**4*e + 2220*d**3*e**2 + 370*d**2*e**3 + 592*d*e**4 - 65*e**5)/(7*e**9) + (d + e*x)**6*(2800*d**6 + 945*d**5*e + 1665*d**4*e**2 + 370*d**3*e**3 + 888*d**2*e**4 - 195*d*e**5 + 107*e**6)/(6*e**9) - (d + e*x)**5*(5*d**2 - 2*d*e + 3*e**2)*(160*d**5 + 127*d**4*e + 88*d**3*e**2 - 4*d**2*e**3 + 64*d*e**4 - 11*e**5)/(5*e**9) + (d + e*x)**4*(5*d**2 - 2*d*e + 3*e**2)**2*(4*d**4 + 5*d**3*e + 3*d**2*e**2 - d*e**3 + 2*e**4)/(4*e**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_294():
    f = (d + e*x)**2*(5*x**2 + 2*x + 3)**2*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)
    F = 18*d**2*x + 3*d*x**2*(11*d + 12*e)/2 + 100*e**2*x**11/11 + e*x**10*(20*d - 9*e/2) + x**9*(100*d**2/9 - 10*d*e + 37*e**2/3) - x**8*(45*d**2/8 - 111*d*e/4 + 37*e**2/8) + x**7*(111*d**2/7 - 74*d*e/7 + 148*e**2/7) - x**6*(37*d**2/6 - 148*d*e/3 - 65*e**2/6) + x**5*(148*d**2/5 + 26*d*e + 107*e**2/5) + x**4*(65*d**2/4 + 107*d*e/2 + 33*e**2/4) + x**3*(107*d**2/3 + 22*d*e + 6*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_295():
    f = (d + e*x)*(5*x**2 + 2*x + 3)**2*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)
    F = 18*d*x + 10*e*x**10 + x**9*(100*d/9 - 5*e) - x**8*(45*d/8 - 111*e/8) + x**7*(111*d/7 - 37*e/7) - x**6*(37*d/6 - 74*e/3) + x**5*(148*d/5 + 13*e) + x**4*(65*d/4 + 107*e/4) + x**3*(107*d/3 + 11*e) + x**2*(33*d/2 + 9*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_296():
    f = (5*x**2 + 2*x + 3)**2*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)
    F = 100*x**9/9 - 45*x**8/8 + 111*x**7/7 - 37*x**6/6 + 148*x**5/5 + 65*x**4/4 + 107*x**3/3 + 33*x**2/2 + 18*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_297():
    f = (5*x**2 + 2*x + 3)**2*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)/(d + e*x)
    F = 25*x**8/(2*e) - x**7*(100*d + 45*e)/(7*e**2) + x**6*(100*d**2 + 45*d*e + 111*e**2)/(6*e**3) - x**5*(100*d**3 + 45*d**2*e + 111*d*e**2 + 37*e**3)/(5*e**4) + x**4*(100*d**4 + 45*d**3*e + 111*d**2*e**2 + 37*d*e**3 + 148*e**4)/(4*e**5) - x**3*(100*d**5 + 45*d**4*e + 111*d**3*e**2 + 37*d**2*e**3 + 148*d*e**4 - 65*e**5)/(3*e**6) + x**2*(100*d**6 + 45*d**5*e + 111*d**4*e**2 + 37*d**3*e**3 + 148*d**2*e**4 - 65*d*e**5 + 107*e**6)/(2*e**7) - x*(100*d**7 + 45*d**6*e + 111*d**5*e**2 + 37*d**4*e**3 + 148*d**3*e**4 - 65*d**2*e**5 + 107*d*e**6 - 33*e**7)/e**8 + (5*d**2 - 2*d*e + 3*e**2)**2*(4*d**4 + 5*d**3*e + 3*d**2*e**2 - d*e**3 + 2*e**4)*log(d + e*x)/e**9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_298():
    f = (5*x**2 + 2*x + 3)**2*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)/(d + e*x)**2
    F = 100*x**7/(7*e**2) - x**6*(200*d + 45*e)/(6*e**3) + x**5*(300*d**2 + 90*d*e + 111*e**2)/(5*e**4) - x**4*(400*d**3 + 135*d**2*e + 222*d*e**2 + 37*e**3)/(4*e**5) + x**3*(500*d**4 + 180*d**3*e + 333*d**2*e**2 + 74*d*e**3 + 148*e**4)/(3*e**6) - x**2*(600*d**5 + 225*d**4*e + 444*d**3*e**2 + 111*d**2*e**3 + 296*d*e**4 - 65*e**5)/(2*e**7) + x*(700*d**6 + 270*d**5*e + 555*d**4*e**2 + 148*d**3*e**3 + 444*d**2*e**4 - 130*d*e**5 + 107*e**6)/e**8 - (5*d**2 - 2*d*e + 3*e**2)*(160*d**5 + 127*d**4*e + 88*d**3*e**2 - 4*d**2*e**3 + 64*d*e**4 - 11*e**5)*log(d + e*x)/e**9 - (5*d**2 - 2*d*e + 3*e**2)**2*(4*d**4 + 5*d**3*e + 3*d**2*e**2 - d*e**3 + 2*e**4)/(e**9*(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_299():
    f = (5*x**2 + 2*x + 3)**2*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)/(d + e*x)**3
    F = 50*x**6/(3*e**3) - x**5*(60*d + 9*e)/e**4 + x**4*(600*d**2 + 135*d*e + 111*e**2)/(4*e**5) - x**3*(1000*d**3 + 270*d**2*e + 333*d*e**2 + 37*e**3)/(3*e**6) + x**2*(1500*d**4 + 450*d**3*e + 666*d**2*e**2 + 111*d*e**3 + 148*e**4)/(2*e**7) - x*(2100*d**5 + 675*d**4*e + 1110*d**3*e**2 + 222*d**2*e**3 + 444*d*e**4 - 65*e**5)/e**8 + (2800*d**6 + 945*d**5*e + 1665*d**4*e**2 + 370*d**3*e**3 + 888*d**2*e**4 - 195*d*e**5 + 107*e**6)*log(d + e*x)/e**9 + (5*d**2 - 2*d*e + 3*e**2)*(160*d**5 + 127*d**4*e + 88*d**3*e**2 - 4*d**2*e**3 + 64*d*e**4 - 11*e**5)/(e**9*(d + e*x)) - (5*d**2 - 2*d*e + 3*e**2)**2*(4*d**4 + 5*d**3*e + 3*d**2*e**2 - d*e**3 + 2*e**4)/(2*e**9*(d + e*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_300():
    f = (5*x**2 + 2*x + 3)**2*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)/(d + e*x)**4
    F = 20*x**5/e**4 - x**4*(400*d + 45*e)/(4*e**5) + x**3*(1000*d**2 + 180*d*e + 111*e**2)/(3*e**6) - x**2*(2000*d**3 + 450*d**2*e + 444*d*e**2 + 37*e**3)/(2*e**7) + x*(3500*d**4 + 900*d**3*e + 1110*d**2*e**2 + 148*d*e**3 + 148*e**4)/e**8 - (5600*d**5 + 1575*d**4*e + 2220*d**3*e**2 + 370*d**2*e**3 + 592*d*e**4 - 65*e**5)*log(d + e*x)/e**9 - (2800*d**6 + 945*d**5*e + 1665*d**4*e**2 + 370*d**3*e**3 + 888*d**2*e**4 - 195*d*e**5 + 107*e**6)/(e**9*(d + e*x)) + (5*d**2 - 2*d*e + 3*e**2)*(160*d**5 + 127*d**4*e + 88*d**3*e**2 - 4*d**2*e**3 + 64*d*e**4 - 11*e**5)/(2*e**9*(d + e*x)**2) - (5*d**2 - 2*d*e + 3*e**2)**2*(4*d**4 + 5*d**3*e + 3*d**2*e**2 - d*e**3 + 2*e**4)/(3*e**9*(d + e*x)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_301():
    f = (d + e*x)**3*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)/(5*x**2 + 2*x + 3)
    F = 2*e**3*x**6/15 + e**2*x**5*(12*d/25 - 33*e/125) + 3*e*x**4*(100*d**2 - 165*d*e + 27*e**2)/500 + x**3*(500*d**3 - 2475*d**2*e + 1215*d*e**2 + 458*e**3)/1875 - x**2*(4125*d**3 - 6075*d**2*e - 6870*d*e**2 + 881*e**3)/6250 + x*(10125*d**3 + 34350*d**2*e - 13215*d*e**2 - 5108*e**3)/15625 - sqrt(14)*(52875*d**3 + 449175*d**2*e - 274845*d*e**2 - 53189*e**3)*atan(sqrt(14)*(5*x + 1)/14)/1093750 + (57250*d**3 - 66075*d**2*e - 76620*d*e**2 + 23431*e**3)*log(5*x**2 + 2*x + 3)/156250
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_302():
    f = (d + e*x)**2*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)/(5*x**2 + 2*x + 3)
    F = 4*e**2*x**5/25 + e*x**4*(2*d/5 - 33*e/100) + x**3*(4*d**2/15 - 22*d*e/25 + 27*e**2/125) - x**2*(825*d**2 - 810*d*e - 458*e**2)/1250 + x*(2025*d**2 + 4580*d*e - 881*e**2)/3125 + (5725*d**2 - 4405*d*e - 2554*e**2)*log(5*x**2 + 2*x + 3)/15625 - sqrt(14)*(10575*d**2 + 59890*d*e - 18323*e**2)*atan(sqrt(14)*(5*x + 1)/14)/218750
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_303():
    f = (d + e*x)*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)/(5*x**2 + 2*x + 3)
    F = e*x**4/5 + x**3*(4*d/15 - 11*e/25) - x**2*(33*d/50 - 81*e/250) + x*(81*d/125 + 458*e/625) - sqrt(14)*(2115*d + 5989*e)*atan(sqrt(14)*(5*x + 1)/14)/43750 + (2290*d - 881*e)*log(5*x**2 + 2*x + 3)/6250
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_304():
    f = (4*x**4 - 5*x**3 + 3*x**2 + x + 2)/(5*x**2 + 2*x + 3)
    F = 4*x**3/15 - 33*x**2/50 + 81*x/125 + 229*log(5*x**2 + 2*x + 3)/625 - 423*sqrt(14)*atan(sqrt(14)*(5*x + 1)/14)/8750
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_305():
    f = (4*x**4 - 5*x**3 + 3*x**2 + x + 2)/((d + e*x)*(5*x**2 + 2*x + 3))
    F = -sqrt(14)*(423*d - 1367*e)*atan(sqrt(14)*(5*x + 1)/14)/(1750*(5*d**2 - 2*d*e + 3*e**2)) + (458*d - 7*e)*log(5*x**2 + 2*x + 3)/(1250*d**2 - 500*d*e + 750*e**2) + 2*x**2/(5*e) - x*(20*d + 33*e)/(25*e**2) + (4*d**4 + 5*d**3*e + 3*d**2*e**2 - d*e**3 + 2*e**4)*log(d + e*x)/(e**3*(5*d**2 - 2*d*e + 3*e**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_306():
    f = (4*x**4 - 5*x**3 + 3*x**2 + x + 2)/((d + e*x)**2*(5*x**2 + 2*x + 3))
    F = (229*d**2 - 7*d*e - 136*e**2)*log(5*x**2 + 2*x + 3)/(25*(5*d**2 - 2*d*e + 3*e**2)**2) - sqrt(14)*(423*d**2 - 2734*d*e + 293*e**2)*atan(sqrt(14)*(5*x + 1)/14)/(350*(5*d**2 - 2*d*e + 3*e**2)**2) + 4*x/(5*e**2) - (40*d**5 + d**4*e + 28*d**3*e**2 + 44*d**2*e**3 - 2*d*e**4 + e**5)*log(d + e*x)/(e**3*(5*d**2 - 2*d*e + 3*e**2)**2) - (4*d**4 + 5*d**3*e + 3*d**2*e**2 - d*e**3 + 2*e**4)/(e**3*(d + e*x)*(5*d**2 - 2*d*e + 3*e**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_307():
    f = (4*x**4 - 5*x**3 + 3*x**2 + x + 2)/((d + e*x)**3*(5*x**2 + 2*x + 3))
    F = -sqrt(14)*(423*d**3 - 4101*d**2*e + 879*d*e**2 + 703*e**3)*atan(sqrt(14)*(5*x + 1)/14)/(70*(5*d**2 - 2*d*e + 3*e**2)**3) + (458*d**3 - 21*d**2*e - 816*d*e**2 + 113*e**3)*log(5*x**2 + 2*x + 3)/(10*(5*d**2 - 2*d*e + 3*e**2)**3) + (100*d**6 - 120*d**5*e + 228*d**4*e**2 - 242*d**3*e**3 + 141*d**2*e**4 + 120*d*e**5 - e**6)*log(d + e*x)/(e**3*(5*d**2 - 2*d*e + 3*e**2)**3) + (40*d**5 + d**4*e + 28*d**3*e**2 + 44*d**2*e**3 - 2*d*e**4 + e**5)/(e**3*(d + e*x)*(5*d**2 - 2*d*e + 3*e**2)**2) - (4*d**4 + 5*d**3*e + 3*d**2*e**2 - d*e**3 + 2*e**4)/(2*e**3*(d + e*x)**2*(5*d**2 - 2*d*e + 3*e**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_308():
    f = (d + e*x)**3*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)/(5*x**2 + 2*x + 3)**2
    F = e**3*x**4/25 + e**2*x**3*(4*d/25 - 41*e/375) + e*x**2*(840*d**2 - 1722*d*e + 373*e**2)/3500 + x*(2800*d**3 - 17220*d**2*e + 9921*d*e**2 + 6053*e**3)/17500 - (d + e*x)**3*(423*x + 1367)/(17500*x**2 + 7000*x + 10500) - (1025*d**3 - 1545*d**2*e - 2601*d*e**2 + 832*e**3)*log(5*x**2 + 2*x + 3)/6250 + sqrt(14)*(32825*d**3 + 317565*d**2*e - 221643*d*e**2 - 67499*e**3)*atan(sqrt(14)*(5*x + 1)/14)/1225000
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_309():
    f = (d + e*x)**2*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)/(5*x**2 + 2*x + 3)**2
    F = 4*e**2*x**3/75 + e*x**2*(4*d/25 - 41*e/250) + x*(2800*d**2 - 11480*d*e + 3307*e**2)/17500 - (d + e*x)**2*(423*x + 1367)/(17500*x**2 + 7000*x + 10500) - (1025*d**2 - 1030*d*e - 867*e**2)*log(5*x**2 + 2*x + 3)/6250 + sqrt(14)*(32825*d**2 + 211710*d*e - 73881*e**2)*atan(sqrt(14)*(5*x + 1)/14)/1225000
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_310():
    f = (d + e*x)*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)/(5*x**2 + 2*x + 3)**2
    F = 2*e*x**2/25 + x*(4*d/25 - 41*e/125) - (d + e*x)*(423*x + 1367)/(17500*x**2 + 7000*x + 10500) - (205*d - 103*e)*log(5*x**2 + 2*x + 3)/1250 + sqrt(14)*(6565*d + 21171*e)*atan(sqrt(14)*(5*x + 1)/14)/245000
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_311():
    f = (4*x**4 - 5*x**3 + 3*x**2 + x + 2)/(5*x**2 + 2*x + 3)**2
    F = 4*x/25 - (423*x + 1367)/(17500*x**2 + 7000*x + 10500) - 41*log(5*x**2 + 2*x + 3)/250 + 1313*sqrt(14)*atan(sqrt(14)*(5*x + 1)/14)/49000
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_312():
    f = (4*x**4 - 5*x**3 + 3*x**2 + x + 2)/((d + e*x)*(5*x**2 + 2*x + 3)**2)
    F = -(1367*d - 293*e + x*(423*d - 1367*e))/((3500*d**2 - 1400*d*e + 2100*e**2)*(5*x**2 + 2*x + 3)) - (205*d**3 - 61*d**2*e + 23*d*e**2 + 14*e**3)*log(5*x**2 + 2*x + 3)/(50*(5*d**2 - 2*d*e + 3*e**2)**2) + sqrt(14)*(6565*d**3 - 26423*d**2*e + 11089*d*e**2 - 6623*e**3)*atan(sqrt(14)*(5*x + 1)/14)/(9800*(5*d**2 - 2*d*e + 3*e**2)**2) + (4*d**4 + 5*d**3*e + 3*d**2*e**2 - d*e**3 + 2*e**4)*log(d + e*x)/(e*(5*d**2 - 2*d*e + 3*e**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_313():
    f = (4*x**4 - 5*x**3 + 3*x**2 + x + 2)/((d + e*x)**2*(5*x**2 + 2*x + 3)**2)
    F = -(1367*d**2 - 586*d*e - 703*e**2 + x*(423*d**2 - 2734*d*e + 293*e**2))/(140*(5*d**2 - 2*d*e + 3*e**2)**2*(5*x**2 + 2*x + 3)) + (41*d**4 - 8*d**3*e - 60*d**2*e**2 + 24*d*e**3 - 5*e**4)*log(d + e*x)/(5*d**2 - 2*d*e + 3*e**2)**3 - (41*d**4 - 8*d**3*e - 60*d**2*e**2 + 24*d*e**3 - 5*e**4)*log(5*x**2 + 2*x + 3)/(2*(5*d**2 - 2*d*e + 3*e**2)**3) + sqrt(14)*(1313*d**4 - 10044*d**3*e + 4290*d**2*e**2 + 156*d*e**3 - 271*e**4)*atan(sqrt(14)*(5*x + 1)/14)/(392*(5*d**2 - 2*d*e + 3*e**2)**3) - (4*d**4 + 5*d**3*e + 3*d**2*e**2 - d*e**3 + 2*e**4)/(e*(d + e*x)*(5*d**2 - 2*d*e + 3*e**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_314():
    f = (4*x**4 - 5*x**3 + 3*x**2 + x + 2)/((d + e*x)**3*(5*x**2 + 2*x + 3)**2)
    F = -(1367*d**3 - 879*d**2*e - 2109*d*e**2 + 457*e**3 + x*(423*d**3 - 4101*d**2*e + 879*d*e**2 + 703*e**3))/(28*(5*d**2 - 2*d*e + 3*e**2)**3*(5*x**2 + 2*x + 3)) + (205*d**5 - 19*d**4*e - 846*d**3*e**2 + 396*d**2*e**3 + 57*d*e**4 - 21*e**5)*log(d + e*x)/(5*d**2 - 2*d*e + 3*e**2)**4 - (205*d**5 - 19*d**4*e - 846*d**3*e**2 + 396*d**2*e**3 + 57*d*e**4 - 21*e**5)*log(5*x**2 + 2*x + 3)/(2*(5*d**2 - 2*d*e + 3*e**2)**4) + sqrt(14)*(6565*d**5 - 74017*d**4*e + 35022*d**3*e**2 + 42858*d**2*e**3 - 17247*d*e**4 + 579*e**5)*atan(sqrt(14)*(5*x + 1)/14)/(392*(5*d**2 - 2*d*e + 3*e**2)**4) - (41*d**4 - 8*d**3*e - 60*d**2*e**2 + 24*d*e**3 - 5*e**4)/((d + e*x)*(5*d**2 - 2*d*e + 3*e**2)**3) - (4*d**4 + 5*d**3*e + 3*d**2*e**2 - d*e**3 + 2*e**4)/(2*e*(d + e*x)**2*(5*d**2 - 2*d*e + 3*e**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_315():
    f = (d + e*x)**3*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)/(5*x**2 + 2*x + 3)**3
    F = 2*e**3*x**2/125 + e**2*x*(83065*d - 126009*e)/980000 + 3*e*(100*d**2 - 245*d*e + 47*e**2)*log(5*x**2 + 2*x + 3)/6250 - (d + e*x)**3*(423*x + 1367)/(7000*(5*x**2 + 2*x + 3)**2) + (d + e*x)**2*(34347*d - 6315*e + x*(11015*d + 49177*e))/(980000*x**2 + 392000*x + 588000) + sqrt(14)*(1059375*d**3 - 2565525*d**2*e + 222255*d*e**2 + 1669047*e**3)*atan(sqrt(14)*(5*x + 1)/14)/68600000
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_316():
    f = (d + e*x)**2*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)/(5*x**2 + 2*x + 3)**3
    F = 4*e**2*x/125 + e*(40*d - 49*e)*log(5*x**2 + 2*x + 3)/1250 - (d + e*x)**2*(423*x + 1367)/(7000*(5*x**2 + 2*x + 3)**2) + (d + e*x)*(34347*d - 6413*e + x*(11015*d + 42765*e))/(980000*x**2 + 392000*x + 588000) + sqrt(14)*(211875*d**2 - 342070*d*e + 14817*e**2)*atan(sqrt(14)*(5*x + 1)/14)/13720000
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_317():
    f = (d + e*x)*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)/(5*x**2 + 2*x + 3)**3
    F = 2*e*log(5*x**2 + 2*x + 3)/125 - (d + e*x)*(423*x + 1367)/(7000*(5*x**2 + 2*x + 3)**2) + sqrt(14)*(42375*d - 34207*e)*atan(sqrt(14)*(5*x + 1)/14)/2744000 + (34347*d - 6511*e + x*(11015*d + 36353*e))/(980000*x**2 + 392000*x + 588000)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_318():
    f = (4*x**4 - 5*x**3 + 3*x**2 + x + 2)/(5*x**2 + 2*x + 3)**3
    F = -(423*x + 1367)/(7000*(5*x**2 + 2*x + 3)**2) + (11015*x + 34347)/(980000*x**2 + 392000*x + 588000) + 339*sqrt(14)*atan(sqrt(14)*(5*x + 1)/14)/21952
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_319():
    f = (4*x**4 - 5*x**3 + 3*x**2 + x + 2)/((d + e*x)*(5*x**2 + 2*x + 3)**3)
    F = e*(4*d**4 + 5*d**3*e + 3*d**2*e**2 - d*e**3 + 2*e**4)*log(d + e*x)/(5*d**2 - 2*d*e + 3*e**2)**3 - e*(4*d**4 + 5*d**3*e + 3*d**2*e**2 - d*e**3 + 2*e**4)*log(5*x**2 + 2*x + 3)/(2*(5*d**2 - 2*d*e + 3*e**2)**3) - (1367*d - 293*e + x*(423*d - 1367*e))/((7000*d**2 - 2800*d*e + 4200*e**2)*(5*x**2 + 2*x + 3)**2) + (171735*d**3 - 92989*d**2*e + 36207*d*e**2 + 1831*e**3 + x*(55075*d**3 - 225825*d**2*e + 90875*d*e**2 - 45725*e**3))/(39200*(5*d**2 - 2*d*e + 3*e**2)**2*(5*x**2 + 2*x + 3)) + sqrt(14)*(42375*d**5 - 16643*d**4*e + 58530*d**3*e**2 - 56058*d**2*e**3 + 31811*d*e**4 - 8623*e**5)*atan(sqrt(14)*(5*x + 1)/14)/(21952*(5*d**2 - 2*d*e + 3*e**2)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_320():
    f = (4*x**4 - 5*x**3 + 3*x**2 + x + 2)/((d + e*x)**2*(5*x**2 + 2*x + 3)**3)
    F = e*(40*d**5 + 83*d**4*e + 12*d**3*e**2 - 76*d**2*e**3 + 46*d*e**4 - 9*e**5)*log(d + e*x)/(5*d**2 - 2*d*e + 3*e**2)**4 - e*(40*d**5 + 83*d**4*e + 12*d**3*e**2 - 76*d**2*e**3 + 46*d*e**4 - 9*e**5)*log(5*x**2 + 2*x + 3)/(2*(5*d**2 - 2*d*e + 3*e**2)**4) - e*(4*d**4 + 5*d**3*e + 3*d**2*e**2 - d*e**3 + 2*e**4)/((d + e*x)*(5*d**2 - 2*d*e + 3*e**2)**3) - (1367*d**2 - 586*d*e - 703*e**2 + x*(423*d**2 - 2734*d*e + 293*e**2))/(280*(5*d**2 - 2*d*e + 3*e**2)**2*(5*x**2 + 2*x + 3)**2) + (171735*d**4 - 117284*d**3*e - 200502*d**2*e**2 + 104428*d*e**3 - 23189*e**4 + x*(55075*d**4 - 429620*d**3*e + 173490*d**2*e**2 + 51740*d*e**3 - 17945*e**4))/(7840*(5*d**2 - 2*d*e + 3*e**2)**3*(5*x**2 + 2*x + 3)) + sqrt(14)*(211875*d**6 + 3070*d**5*e + 209039*d**4*e**2 - 921444*d**3*e**3 + 380621*d**2*e**4 - 49586*d*e**5 - 43695*e**6)*atan(sqrt(14)*(5*x + 1)/14)/(21952*(5*d**2 - 2*d*e + 3*e**2)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_321():
    f = (2*x + 5)*sqrt(2*x**2 - x + 3)*(5*x**4 - x**3 + 3*x**2 + x + 2)
    F = 5*(2*x + 5)**4*(2*x**2 - x + 3)**(sympy.S(3)/2)/112 - 823*(2*x + 5)**3*(2*x**2 - x + 3)**(sympy.S(3)/2)/1344 + 11433*(2*x + 5)**2*(2*x**2 - x + 3)**(sympy.S(3)/2)/4480 + (205740*x - 51435)*sqrt(2*x**2 - x + 3)/32768 - (295276*x + 1005757)*(2*x**2 - x + 3)**(sympy.S(3)/2)/71680 - 1183005*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/131072
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_322():
    f = sqrt(2*x**2 - x + 3)*(5*x**4 - x**3 + 3*x**2 + x + 2)
    F = 5*x**3*(2*x**2 - x + 3)**(sympy.S(3)/2)/12 + 7*x**2*(2*x**2 - x + 3)**(sympy.S(3)/2)/80 - 71*x*(2*x**2 - x + 3)**(sympy.S(3)/2)/1280 + (18436*x - 4609)*sqrt(2*x**2 - x + 3)/16384 + 287*(2*x**2 - x + 3)**(sympy.S(3)/2)/5120 - 106007*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/65536
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_323():
    f = sqrt(2*x**2 - x + 3)*(5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x + 5)
    F = (489587 - 80844*x)*sqrt(2*x**2 - x + 3)/4096 + (2*x + 5)**2*(2*x**2 - x + 3)**(sympy.S(3)/2)/16 - (254*x + 635)*(2*x**2 - x + 3)**(sympy.S(3)/2)/128 + 4535*(2*x**2 - x + 3)**(sympy.S(3)/2)/768 + 5627989*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/16384 - 11001*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/32
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_324():
    f = sqrt(2*x**2 - x + 3)*(5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x + 5)**2
    F = -(1996953 - 333380*x)*sqrt(2*x**2 - x + 3)/18432 + (10*x + 25)*(2*x**2 - x + 3)**(sympy.S(3)/2)/64 - 541*(2*x**2 - x + 3)**(sympy.S(3)/2)/384 - 2551847*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/8192 + 239201*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/768 - 3667*(2*x**2 - x + 3)**(sympy.S(3)/2)/(1152*x + 2880)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_325():
    f = sqrt(2*x**2 - x + 3)*(5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x + 5)**3
    F = (3305325 - 550495*x)*sqrt(2*x**2 - x + 3)/82944 + 5*(2*x**2 - x + 3)**(sympy.S(3)/2)/48 + 117315*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/1024 - 12670805*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/110592 + 357391*(2*x**2 - x + 3)**(sympy.S(3)/2)/(165888*x + 414720) - 3667*(2*x**2 - x + 3)**(sympy.S(3)/2)/(1152*(2*x + 5)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_326():
    f = sqrt(2*x**2 - x + 3)*(5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x + 5)**4
    F = -(44378877 - 7400779*x)*sqrt(2*x**2 - x + 3)/5971968 - 10939*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/512 + 170114729*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/7962624 - 6467659*(2*x**2 - x + 3)**(sympy.S(3)/2)/(11943936*x + 29859840) + 158527*(2*x**2 - x + 3)**(sympy.S(3)/2)/(82944*(2*x + 5)**2) - 3667*(2*x**2 - x + 3)**(sympy.S(3)/2)/(1728*(2*x + 5)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_327():
    f = sqrt(2*x**2 - x + 3)*(5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x + 5)**5
    F = (67313372*x + 369856585)*sqrt(2*x**2 - x + 3)/(191102976*x + 477757440) + 259*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/128 - 4640586097*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/2293235712 - 9363383*(2*x**2 - x + 3)**(sympy.S(3)/2)/(23887872*(2*x + 5)**2) + 593771*(2*x**2 - x + 3)**(sympy.S(3)/2)/(497664*(2*x + 5)**3) - 3667*(2*x**2 - x + 3)**(sympy.S(3)/2)/(2304*(2*x + 5)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_328():
    f = sqrt(2*x**2 - x + 3)*(5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x + 5)**6
    F = -5*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/64 + 12895597463*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/165112971264 - (3174439702*x + 4583087983)*sqrt(2*x**2 - x + 3)/(6879707136*(2*x + 5)**2) - 38732321*(2*x**2 - x + 3)**(sympy.S(3)/2)/(179159040*(2*x + 5)**3) + 711961*(2*x**2 - x + 3)**(sympy.S(3)/2)/(829440*(2*x + 5)**4) - 3667*(2*x**2 - x + 3)**(sympy.S(3)/2)/(2880*(2*x + 5)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_329():
    f = sqrt(2*x**2 - x + 3)*(5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x + 5)**7
    F = -26972675*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/7925422620672 + (25799950*x - 19936325)*sqrt(2*x**2 - x + 3)/(330225942528*(2*x + 5)**2) + 87677717*(2*x**2 - x + 3)**(sympy.S(3)/2)/(8599633920*(2*x + 5)**3) - 5703277*(2*x**2 - x + 3)**(sympy.S(3)/2)/(39813120*(2*x + 5)**4) + 92239*(2*x**2 - x + 3)**(sympy.S(3)/2)/(138240*(2*x + 5)**5) - 3667*(2*x**2 - x + 3)**(sympy.S(3)/2)/(3456*(2*x + 5)**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_330():
    f = sqrt(2*x**2 - x + 3)*(5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x + 5)**8
    F = -289071245*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/570630428688384 + (276502930*x - 213661355)*sqrt(2*x**2 - x + 3)/(23776267862016*(2*x + 5)**2) + 246159769*(2*x**2 - x + 3)**(sympy.S(3)/2)/(866843099136*(2*x + 5)**3) + 19414831*(2*x**2 - x + 3)**(sympy.S(3)/2)/(4013162496*(2*x + 5)**4) - 1464037*(2*x**2 - x + 3)**(sympy.S(3)/2)/(13934592*(2*x + 5)**5) + 948341*(2*x**2 - x + 3)**(sympy.S(3)/2)/(1741824*(2*x + 5)**6) - 3667*(2*x**2 - x + 3)**(sympy.S(3)/2)/(4032*(2*x + 5)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_331():
    f = (2*x + 5)*(2*x**2 - x + 3)**(sympy.S(3)/2)*(5*x**4 - x**3 + 3*x**2 + x + 2)
    F = -(92727 - 370908*x)*(2*x**2 - x + 3)**(sympy.S(3)/2)/131072 + 5*(2*x + 5)**4*(2*x**2 - x + 3)**(sympy.S(5)/2)/144 - 1121*(2*x + 5)**3*(2*x**2 - x + 3)**(sympy.S(5)/2)/2304 + 69415*(2*x + 5)**2*(2*x**2 - x + 3)**(sympy.S(5)/2)/32256 - (647700*x + 1984191)*(2*x**2 - x + 3)**(sympy.S(5)/2)/143360 + (25592652*x - 6398163)*sqrt(2*x**2 - x + 3)/2097152 - 147157749*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/8388608
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_332():
    f = (2*x**2 - x + 3)**(sympy.S(3)/2)*(5*x**4 - x**3 + 3*x**2 + x + 2)
    F = 5*x**3*(2*x**2 - x + 3)**(sympy.S(5)/2)/16 + 23*x**2*(2*x**2 - x + 3)**(sympy.S(5)/2)/448 + 125*x*(2*x**2 - x + 3)**(sympy.S(5)/2)/3584 - (8597 - 34388*x)*(2*x**2 - x + 3)**(sympy.S(3)/2)/65536 + (2372772*x - 593193)*sqrt(2*x**2 - x + 3)/1048576 + 1167*(2*x**2 - x + 3)**(sympy.S(5)/2)/14336 - 13643439*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/4194304
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_333():
    f = (2*x**2 - x + 3)**(sympy.S(3)/2)*(5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x + 5)
    F = (500141 - 123060*x)*(2*x**2 - x + 3)**(sympy.S(3)/2)/12288 + (141051019 - 23482924*x)*sqrt(2*x**2 - x + 3)/65536 + 5*(2*x + 5)**2*(2*x**2 - x + 3)**(sympy.S(5)/2)/112 - (622*x + 1555)*(2*x**2 - x + 3)**(sympy.S(5)/2)/448 + 3505*(2*x**2 - x + 3)**(sympy.S(5)/2)/896 + 1622009981*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/262144 - 99009*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_334():
    f = (2*x**2 - x + 3)**(sympy.S(3)/2)*(5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x + 5)**2
    F = -(909513 - 226052*x)*(2*x**2 - x + 3)**(sympy.S(3)/2)/18432 - (85448933 - 14243732*x)*sqrt(2*x**2 - x + 3)/32768 + (10*x + 25)*(2*x**2 - x + 3)**(sympy.S(5)/2)/96 - 839*(2*x**2 - x + 3)**(sympy.S(5)/2)/960 - 982669459*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/131072 + 959625*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/128 - 3667*(2*x**2 - x + 3)**(sympy.S(5)/2)/(1152*x + 2880)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_335():
    f = (2*x**2 - x + 3)**(sympy.S(3)/2)*(5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x + 5)**3
    F = (2154633 - 534617*x)*(2*x**2 - x + 3)**(sympy.S(3)/2)/82944 + (33741483 - 5623292*x)*sqrt(2*x**2 - x + 3)/24576 + (2*x**2 - x + 3)**(sympy.S(5)/2)/16 + 129342063*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/32768 - 8083915*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/2048 + 438065*(2*x**2 - x + 3)**(sympy.S(5)/2)/(165888*x + 414720) - 3667*(2*x**2 - x + 3)**(sympy.S(5)/2)/(1152*(2*x + 5)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_336():
    f = (2*x**2 - x + 3)**(sympy.S(3)/2)*(5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x + 5)**4
    F = -(135068604 - 22512089*x)*sqrt(2*x**2 - x + 3)/331776 - (138006843 - 34265045*x)*(2*x**2 - x + 3)**(sympy.S(3)/2)/17915904 - 19176431*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/16384 + 517762327*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/442368 - 32865365*(2*x**2 - x + 3)**(sympy.S(5)/2)/(35831808*x + 89579520) + 556255*(2*x**2 - x + 3)**(sympy.S(5)/2)/(248832*(2*x + 5)**2) - 3667*(2*x**2 - x + 3)**(sympy.S(5)/2)/(1728*(2*x + 5)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_337():
    f = (2*x**2 - x + 3)**(sympy.S(3)/2)*(5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x + 5)**5
    F = (2339916063 - 389975609*x)*sqrt(2*x**2 - x + 3)/31850496 + (67865260*x + 762984903)*(2*x**2 - x + 3)**(sympy.S(3)/2)/(191102976*x + 477757440) + 432565*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/2048 - 8969688643*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/42467328 - 14477995*(2*x**2 - x + 3)**(sympy.S(5)/2)/(23887872*(2*x + 5)**2) + 224815*(2*x**2 - x + 3)**(sympy.S(5)/2)/(165888*(2*x + 5)**3) - 3667*(2*x**2 - x + 3)**(sympy.S(5)/2)/(2304*(2*x + 5)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_338():
    f = (2*x**2 - x + 3)**(sympy.S(3)/2)*(5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x + 5)**6
    F = -23775*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/1024 + 70991525167*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/3057647616 - (1028823716*x + 5658774871)*sqrt(2*x**2 - x + 3)/(254803968*x + 637009920) + (44773976*x + 246012435)*(2*x**2 - x + 3)**(sympy.S(3)/2)/(95551488*(2*x + 5)**2) - 3730507*(2*x**2 - x + 3)**(sympy.S(5)/2)/(11943936*(2*x + 5)**3) + 158527*(2*x**2 - x + 3)**(sympy.S(5)/2)/(165888*(2*x + 5)**4) - 3667*(2*x**2 - x + 3)**(sympy.S(5)/2)/(2880*(2*x + 5)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_339():
    f = (2*x**2 - x + 3)**(sympy.S(3)/2)*(5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x + 5)**7
    F = (27596573612*x + 151764102421)*sqrt(2*x**2 - x + 3)/(110075314176*x + 275188285440) + 369*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/256 - 1903976002333*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/1320903770112 - (6793718806*x + 9802984711)*(2*x**2 - x + 3)**(sympy.S(3)/2)/(13759414272*(2*x + 5)**3) - 14087245*(2*x**2 - x + 3)**(sympy.S(5)/2)/(71663616*(2*x + 5)**4) + 182165*(2*x**2 - x + 3)**(sympy.S(5)/2)/(248832*(2*x + 5)**5) - 3667*(2*x**2 - x + 3)**(sympy.S(5)/2)/(3456*(2*x + 5)**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_340():
    f = (2*x**2 - x + 3)**(sympy.S(3)/2)*(5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x + 5)**8
    F = -5*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/128 + 412760561351*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/10567230160896 - (101679102454*x + 146583836191)*sqrt(2*x**2 - x + 3)/(440301256704*(2*x + 5)**2) - (411822458*x + 463558457)*(2*x**2 - x + 3)**(sympy.S(3)/2)/(2293235712*(2*x + 5)**4) - 1930441*(2*x**2 - x + 3)**(sympy.S(5)/2)/(13934592*(2*x + 5)**5) + 114335*(2*x**2 - x + 3)**(sympy.S(5)/2)/(193536*(2*x + 5)**6) - 3667*(2*x**2 - x + 3)**(sympy.S(5)/2)/(4032*(2*x + 5)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_341():
    f = (2*x + 5)*(5*x**4 - x**3 + 3*x**2 + x + 2)/sqrt(2*x**2 - x + 3)
    F = (2*x + 5)**4*sqrt(2*x**2 - x + 3)/16 - 105*(2*x + 5)**3*sqrt(2*x**2 - x + 3)/128 + 761*(2*x + 5)**2*sqrt(2*x**2 - x + 3)/256 - (4676*x + 19227)*sqrt(2*x**2 - x + 3)/2048 - 85429*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/8192
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_342():
    f = (5*x**4 - x**3 + 3*x**2 + x + 2)/sqrt(2*x**2 - x + 3)
    F = 5*x**3*sqrt(2*x**2 - x + 3)/8 + 19*x**2*sqrt(2*x**2 - x + 3)/96 - 409*x*sqrt(2*x**2 - x + 3)/768 - 505*sqrt(2*x**2 - x + 3)/1024 - 6863*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/4096
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_343():
    f = (5*x**4 - x**3 + 3*x**2 + x + 2)/((2*x + 5)*sqrt(2*x**2 - x + 3))
    F = 5*(2*x + 5)**2*sqrt(2*x**2 - x + 3)/48 - (674*x + 1685)*sqrt(2*x**2 - x + 3)/192 + 1669*sqrt(2*x**2 - x + 3)/128 + 9657*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/512 - 3667*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/192
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_344():
    f = (5*x**4 - x**3 + 3*x**2 + x + 2)/((2*x + 5)**2*sqrt(2*x**2 - x + 3))
    F = (10*x + 25)*sqrt(2*x**2 - x + 3)/32 - 243*sqrt(2*x**2 - x + 3)/64 - 2943*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/256 + 158527*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/13824 - 3667*sqrt(2*x**2 - x + 3)/(1152*x + 2880)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_345():
    f = (5*x**4 - x**3 + 3*x**2 + x + 2)/((2*x + 5)**3*sqrt(2*x**2 - x + 3))
    F = 5*sqrt(2*x**2 - x + 3)/16 + 149*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/64 - 1546507*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/663552 + 92239*sqrt(2*x**2 - x + 3)/(55296*x + 138240) - 3667*sqrt(2*x**2 - x + 3)/(1152*(2*x + 5)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_346():
    f = (5*x**4 - x**3 + 3*x**2 + x + 2)/((2*x + 5)**4*sqrt(2*x**2 - x + 3))
    F = -5*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/32 + 22389491*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/143327232 - 3163415*sqrt(2*x**2 - x + 3)/(11943936*x + 29859840) + 394907*sqrt(2*x**2 - x + 3)/(248832*(2*x + 5)**2) - 3667*sqrt(2*x**2 - x + 3)/(1728*(2*x + 5)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_347():
    f = (5*x**4 - x**3 + 3*x**2 + x + 2)/((2*x + 5)**5*sqrt(2*x**2 - x + 3))
    F = 2053207*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/41278242816 + 26800085*sqrt(2*x**2 - x + 3)/(3439853568*x + 8599633920) - 16295969*sqrt(2*x**2 - x + 3)/(71663616*(2*x + 5)**2) + 513097*sqrt(2*x**2 - x + 3)/(497664*(2*x + 5)**3) - 3667*sqrt(2*x**2 - x + 3)/(2304*(2*x + 5)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_348():
    f = (2*x + 5)**2*(5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x**2 - x + 3)**(sympy.S(3)/2)
    F = 5*x**3*sqrt(2*x**2 - x + 3)/4 + 153*x**2*sqrt(2*x**2 - x + 3)/16 + 2645*x*sqrt(2*x**2 - x + 3)/128 + (2132*x - 1384)/(23*sqrt(2*x**2 - x + 3)) - 13153*sqrt(2*x**2 - x + 3)/512 + 144217*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/2048
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_349():
    f = (2*x + 5)*(5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x**2 - x + 3)**(sympy.S(3)/2)
    F = 5*x**2*sqrt(2*x**2 - x + 3)/6 + 193*x*sqrt(2*x**2 - x + 3)/48 + (373*x - 53)/(23*sqrt(2*x**2 - x + 3)) + 33*sqrt(2*x**2 - x + 3)/64 + 3111*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/256
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_350():
    f = (5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x**2 - x + 3)**(sympy.S(3)/2)
    F = 5*x*sqrt(2*x**2 - x + 3)/8 + (219*x + 89)/(92*sqrt(2*x**2 - x + 3)) + 27*sqrt(2*x**2 - x + 3)/32 + 213*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/128
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_351():
    f = (5*x**4 - x**3 + 3*x**2 + x + 2)/((2*x + 5)*(2*x**2 - x + 3)**(sympy.S(3)/2))
    F = (917*x + 1191)/(3312*sqrt(2*x**2 - x + 3)) + 5*sqrt(2*x**2 - x + 3)/8 + 39*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/32 - 3667*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/3456
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_352():
    f = (5*x**4 - x**3 + 3*x**2 + x + 2)/((2*x + 5)**2*(2*x**2 - x + 3)**(sympy.S(3)/2))
    F = (2203*x + 9897)/(119232*sqrt(2*x**2 - x + 3)) - 5*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/16 + 25951*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/82944 - 3667*sqrt(2*x**2 - x + 3)/(20736*x + 51840)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_353():
    f = (5*x**4 - x**3 + 3*x**2 + x + 2)/((2*x + 5)**3*(2*x**2 - x + 3)**(sympy.S(3)/2))
    F = (65991 - 8779*x)/(4292352*sqrt(2*x**2 - x + 3)) - 52631*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/11943936 + 115369*sqrt(2*x**2 - x + 3)/(2985984*x + 7464960) - 3667*sqrt(2*x**2 - x + 3)/(20736*(2*x + 5)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_354():
    f = (5*x**4 - x**3 + 3*x**2 + x + 2)/((2*x + 5)**4*(2*x**2 - x + 3)**(sympy.S(3)/2))
    F = (369609 - 175877*x)/(154524672*sqrt(2*x**2 - x + 3)) - 3505819*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/2579890176 + 430799*sqrt(2*x**2 - x + 3)/(214990848*x + 537477120) + 152885*sqrt(2*x**2 - x + 3)/(4478976*(2*x + 5)**2) - 3667*sqrt(2*x**2 - x + 3)/(31104*(2*x + 5)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_355():
    f = (2*x + 5)**2*(5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x**2 - x + 3)**(sympy.S(5)/2)
    F = 5*x*sqrt(2*x**2 - x + 3)/4 + (75928 - 81532*x)/(1587*sqrt(2*x**2 - x + 3)) + (2132*x - 1384)/(69*(2*x**2 - x + 3)**(sympy.S(3)/2)) + 247*sqrt(2*x**2 - x + 3)/16 - 1471*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_356():
    f = (2*x + 5)*(5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x**2 - x + 3)**(sympy.S(5)/2)
    F = (6055 - 28981*x)/(3174*sqrt(2*x**2 - x + 3)) + (373*x - 53)/(69*(2*x**2 - x + 3)**(sympy.S(3)/2)) + 5*sqrt(2*x**2 - x + 3)/4 - 71*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_357():
    f = (5*x**4 - x**3 + 3*x**2 + x + 2)/(2*x**2 - x + 3)**(sympy.S(5)/2)
    F = (219*x + 89)/(276*(2*x**2 - x + 3)**(sympy.S(3)/2)) - (2604*x + 1465)/(2116*sqrt(2*x**2 - x + 3)) - 5*sqrt(2)*asinh(sqrt(23)*(1 - 4*x)/23)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_358():
    f = (5*x**4 - x**3 + 3*x**2 + x + 2)/((2*x + 5)*(2*x**2 - x + 3)**(sympy.S(5)/2))
    F = (917*x + 1191)/(9936*(2*x**2 - x + 3)**(sympy.S(3)/2)) - (146729*x + 335337)/(1371168*sqrt(2*x**2 - x + 3)) - 3667*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/62208
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_359():
    f = (5*x**4 - x**3 + 3*x**2 + x + 2)/((2*x + 5)**2*(2*x**2 - x + 3)**(sympy.S(5)/2))
    F = -(1255878 - 62021*x)/(24681024*sqrt(2*x**2 - x + 3)) + (2203*x + 9897)/(357696*(2*x**2 - x + 3)**(sympy.S(3)/2)) - 2821*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/4478976 - 3667*sqrt(2*x**2 - x + 3)/(373248*x + 933120)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_360():
    f = (5*x**4 - x**3 + 3*x**2 + x + 2)/((2*x + 5)**3*(2*x**2 - x + 3)**(sympy.S(5)/2))
    F = (65991 - 8779*x)/(12877056*(2*x**2 - x + 3)**(sympy.S(3)/2)) - (4679797 - 2148263*x)/(592344576*sqrt(2*x**2 - x + 3)) + 774079*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/644972544 - 45979*sqrt(2*x**2 - x + 3)/(53747712*x + 134369280) - 3667*sqrt(2*x**2 - x + 3)/(373248*(2*x + 5)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_361():
    f = (5*x**4 - x**3 + 3*x**2 + x + 2)/((2*x + 5)**4*(2*x**2 - x + 3)**(sympy.S(5)/2))
    F = (369609 - 175877*x)/(463574016*(2*x**2 - x + 3)**(sympy.S(3)/2)) - (27754539 - 31190998*x)/(31986607104*sqrt(2*x**2 - x + 3)) + 4778789*sqrt(2)*atanh(sqrt(2)*(17 - 22*x)/(24*sqrt(2*x**2 - x + 3)))/15479341056 + 475357*sqrt(2*x**2 - x + 3)/(3869835264*x + 9674588160) - 89137*sqrt(2*x**2 - x + 3)/(80621568*(2*x + 5)**2) - 3667*sqrt(2*x**2 - x + 3)/(559872*(2*x + 5)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_362():
    f = (f + g*x + h*x**2 + i*x**3 + j*x**4)/(a + b*x + c*x**2)**(sympy.S(5)/2)
    F = (-2*a*b**3*j + 2*a*b**2*c*i + 4*a*c**2*(-a*i + c*g) - 2*b*c*(-3*a**2*j + a*c*h + c**2*f) - 2*x*(b**4*j - b**2*c*(4*a*j + b*i) + 2*c**4*f - c**3*(2*a*h + b*g) + c**2*(2*a**2*j + 3*a*b*i + b**2*h)))/(3*c**3*(-4*a*c + b**2)*(a + b*x + c*x**2)**(sympy.S(3)/2)) - (48*a**2*c**3*i - 2*b**5*j + 2*b**4*c*i - 2*b**3*c*(-10*a*j + c*h) + 4*b**2*c**2*(-3*a*i + 2*c*g) - 8*b*c**2*(8*a**2*j + a*c*h + 2*c**2*f) - 2*c*x*(-4*b**4*j + b**2*c*(28*a*j + b*i) + 16*c**4*f - c**3*(-8*a*h + 8*b*g) + 2*c**2*(-16*a**2*j - 6*a*b*i + b**2*h)))/(3*c**3*(-4*a*c + b**2)**2*sqrt(a + b*x + c*x**2)) + j*atanh((b + 2*c*x)/(2*sqrt(c)*sqrt(a + b*x + c*x**2)))/c**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_363():
    f = (f + g*x + h*x**2 + i*x**3 + j*x**4)/(a + b*x - c*x**2)**(sympy.S(5)/2)
    F = (2*a*b**3*j + 2*a*b**2*c*i + 4*a*c**2*(a*i + c*g) - 2*b*c*(-3*a**2*j - a*c*h + c**2*f) + 2*x*(b**4*j + b**2*c*(4*a*j + b*i) + 2*c**4*f + c**3*(2*a*h + b*g) + c**2*(2*a**2*j + 3*a*b*i + b**2*h)))/(3*c**3*(4*a*c + b**2)*(a + b*x - c*x**2)**(sympy.S(3)/2)) - (48*a**2*c**3*i + 2*b**5*j + 2*b**4*c*i + 2*b**3*c*(10*a*j + c*h) + 4*b**2*c**2*(3*a*i + 2*c*g) + 8*b*c**2*(8*a**2*j - a*c*h + 2*c**2*f) - 2*c*x*(-4*b**4*j - b**2*c*(28*a*j + b*i) + 16*c**4*f + 8*c**3*(-a*h + b*g) + 2*c**2*(-16*a**2*j - 6*a*b*i + b**2*h)))/(3*c**3*(4*a*c + b**2)**2*sqrt(a + b*x - c*x**2)) - j*atan((b - 2*c*x)/(2*sqrt(c)*sqrt(a + b*x - c*x**2)))/c**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_364():
    f = (d + e*x)**m*(5*x**2 + 2*x + 3)**3*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)
    F = (d + e*x)**(m + 1)*(5*d**2 - 2*d*e + 3*e**2)**3*(4*d**4 + 5*d**3*e + 3*d**2*e**2 - d*e**3 + 2*e**4)/(e**11*(m + 1)) - (d + e*x)**(m + 2)*(5*d**2 - 2*d*e + 3*e**2)**2*(200*d**5 + 169*d**4*e + 108*d**3*e**2 - 20*d**2*e**3 + 86*d*e**4 - 15*e**5)/(e**11*(m + 2)) + (d + e*x)**(m + 3)*(15*d**2 - 6*d*e + 9*e**2)*(1500*d**6 + 660*d**5*e + 792*d**4*e**2 + 58*d**3*e**3 + 547*d**2*e**4 - 156*d*e**5 + 53*e**6)/(e**11*(m + 3)) - (d + e*x)**(m + 4)*(60000*d**7 + 2100*d**6*e + 42840*d**5*e**2 + 3430*d**4*e**3 + 19980*d**3*e**4 - 5100*d**2*e**5 + 4436*d*e**6 - 574*e**7)/(e**11*(m + 4)) + (d + e*x)**(m + 5)*(105000*d**6 + 3150*d**5*e + 53550*d**4*e**2 + 3430*d**3*e**3 + 14985*d**2*e**4 - 2550*d*e**5 + 1109*e**6)/(e**11*(m + 5)) - (d + e*x)**(m + 6)*(126000*d**5 + 3150*d**4*e + 42840*d**3*e**2 + 2058*d**2*e**3 + 5994*d*e**4 - 510*e**5)/(e**11*(m + 6)) + (d + e*x)**(m + 7)*(105000*d**4 + 2100*d**3*e + 21420*d**2*e**2 + 686*d*e**3 + 999*e**4)/(e**11*(m + 7)) - (d + e*x)**(m + 8)*(60000*d**3 + 900*d**2*e + 6120*d*e**2 + 98*e**3)/(e**11*(m + 8)) + (d + e*x)**(m + 9)*(22500*d**2 + 225*d*e + 765*e**2)/(e**11*(m + 9)) - (d + e*x)**(m + 10)*(5000*d + 25*e)/(e**11*(m + 10)) + 500*(d + e*x)**(m + 11)/(e**11*(m + 11))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_365():
    f = (d + e*x)**m*(5*x**2 + 2*x + 3)**2*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)
    F = (d + e*x)**(m + 1)*(5*d**2 - 2*d*e + 3*e**2)**2*(4*d**4 + 5*d**3*e + 3*d**2*e**2 - d*e**3 + 2*e**4)/(e**9*(m + 1)) - (d + e*x)**(m + 2)*(5*d**2 - 2*d*e + 3*e**2)*(160*d**5 + 127*d**4*e + 88*d**3*e**2 - 4*d**2*e**3 + 64*d*e**4 - 11*e**5)/(e**9*(m + 2)) + (d + e*x)**(m + 3)*(2800*d**6 + 945*d**5*e + 1665*d**4*e**2 + 370*d**3*e**3 + 888*d**2*e**4 - 195*d*e**5 + 107*e**6)/(e**9*(m + 3)) - (d + e*x)**(m + 4)*(5600*d**5 + 1575*d**4*e + 2220*d**3*e**2 + 370*d**2*e**3 + 592*d*e**4 - 65*e**5)/(e**9*(m + 4)) + (d + e*x)**(m + 5)*(7000*d**4 + 1575*d**3*e + 1665*d**2*e**2 + 185*d*e**3 + 148*e**4)/(e**9*(m + 5)) - (d + e*x)**(m + 6)*(5600*d**3 + 945*d**2*e + 666*d*e**2 + 37*e**3)/(e**9*(m + 6)) + (d + e*x)**(m + 7)*(2800*d**2 + 315*d*e + 111*e**2)/(e**9*(m + 7)) - (d + e*x)**(m + 8)*(800*d + 45*e)/(e**9*(m + 8)) + 100*(d + e*x)**(m + 9)/(e**9*(m + 9))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_366():
    f = (d + e*x)**m*(5*x**2 + 2*x + 3)*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)
    F = (d + e*x)**(m + 1)*(5*d**2 - 2*d*e + 3*e**2)*(4*d**4 + 5*d**3*e + 3*d**2*e**2 - d*e**3 + 2*e**4)/(e**7*(m + 1)) - (d + e*x)**(m + 2)*(120*d**5 + 85*d**4*e + 68*d**3*e**2 + 12*d**2*e**3 + 42*d*e**4 - 7*e**5)/(e**7*(m + 2)) + (d + e*x)**(m + 3)*(300*d**4 + 170*d**3*e + 102*d**2*e**2 + 12*d*e**3 + 21*e**4)/(e**7*(m + 3)) - (d + e*x)**(m + 4)*(400*d**3 + 170*d**2*e + 68*d*e**2 + 4*e**3)/(e**7*(m + 4)) + (d + e*x)**(m + 5)*(300*d**2 + 85*d*e + 17*e**2)/(e**7*(m + 5)) - (d + e*x)**(m + 6)*(120*d + 17*e)/(e**7*(m + 6)) + 20*(d + e*x)**(m + 7)/(e**7*(m + 7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_367():
    f = (d + e*x)**m*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)/(5*x**2 + 2*x + 3)
    F = -(-423*sqrt(14) + 6412*I)*(d + e*x)**(m + 1)*hyper((1, m + 1), (m + 2,), (5*d + 5*e*x)/(5*d - e + sqrt(14)*I*e))/((m + 1)*(17500*I*d - 3500*e*(sqrt(14) + I))) - (423*sqrt(14) + 6412*I)*(d + e*x)**(m + 1)*hyper((1, m + 1), (m + 2,), (5*d + 5*e*x)/(5*d - e*(1 + sqrt(14)*I)))/((m + 1)*(17500*I*d - 3500*e*(-sqrt(14) + I))) + (d + e*x)**(m + 1)*(100*d**2 + 165*d*e + 81*e**2)/(125*e**3*(m + 1)) - (d + e*x)**(m + 2)*(40*d + 33*e)/(25*e**3*(m + 2)) + 4*(d + e*x)**(m + 3)/(5*e**3*(m + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_368():
    f = (d + e*x)**m*(4*x**4 - 5*x**3 + 3*x**2 + x + 2)/(5*x**2 + 2*x + 3)**2
    F = -(d + e*x)**(m + 1)*(1367*d - 293*e + x*(423*d - 1367*e))/((3500*d**2 - 1400*d*e + 2100*e**2)*(5*x**2 + 2*x + 3)) + (d + e*x)**(m + 1)*(80360*d**2 - 5922*d*e*m - 32144*d*e + 19138*e**2*m + 48216*e**2 + sqrt(14)*I*(6565*d**2 - 2*d*e*(1313 - 3206*m) + e**2*(3939 - 98*m)))*hyper((1, m + 1), (m + 2,), (5*d + 5*e*x)/(5*d - e + sqrt(14)*I*e))/((98000*d + 19600*I*e*(sqrt(14) + I))*(m + 1)*(5*d**2 - 2*d*e + 3*e**2)) + (d + e*x)**(m + 1)*(80360*d**2 - 5922*d*e*m - 32144*d*e + 19138*e**2*m + 48216*e**2 - sqrt(14)*I*(6565*d**2 - 2*d*e*(1313 - 3206*m) + e**2*(3939 - 98*m)))*hyper((1, m + 1), (m + 2,), (5*d + 5*e*x)/(5*d - e*(1 + sqrt(14)*I)))/((98000*d - 19600*e*(1 + sqrt(14)*I))*(m + 1)*(5*d**2 - 2*d*e + 3*e**2)) + 4*(d + e*x)**(m + 1)/(25*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_369():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4 + i*x**5)/(a + b*x + c*x**2)**3
    F = i*log(a + b*x + c*x**2)/(2*c**3) - (-30*a**2*b*c**2*i + 10*a*b**3*c*i - b**5*i + 12*c**5*d - c**4*(-4*a*f + 6*b*e) + 2*c**3*(6*a**2*h - 3*a*b*g + b**2*f))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**3*(-4*a*c + b**2)**(sympy.S(5)/2)) - (-a*b**4*i + a*b**3*c*h - a*b**2*c*(-4*a*i + c*g) - 2*a*c**2*(a**2*i - a*c*g + c**2*e) + b*c**2*(-3*a**2*h + a*c*f + c**2*d) + x*(-b**5*i + b**3*c*(5*a*i + b*h) - b*c**2*(5*a**2*i + 4*a*b*h + b**2*g) + 2*c**5*d - c**4*(2*a*f + b*e) + c**3*(2*a**2*h + 3*a*b*g + b**2*f)))/(2*c**4*(-4*a*c + b**2)*(a + b*x + c*x**2)**2) + (-16*a**2*c**3*(-2*a*i + c*g) - b**6*i + b**5*c*h - b**4*c*(-11*a*i + c*g) + b**3*c**2*(-8*a*h + c*f) - b**2*c**2*(39*a**2*i - 5*a*c*g + 3*c**2*e) + 2*b*c**3*(11*a**2*h + a*c*f + 3*c**2*d) + 2*c*x*(a*b*c**2*(25*a*i + 8*b*h) + 2*b**5*i - b**3*c*(15*a*i + b*h) + 6*c**5*d - c**4*(-2*a*f + 3*b*e) + c**3*(-10*a**2*h - 3*a*b*g + b**2*f)))/(2*c**4*(-4*a*c + b**2)**2*(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_370():
    f = (d + e*x + f*x**2 + g*x**3 + h*x**4 + j*x**5 + k*x**6 + l*x**7 + m*x**8)/(a + b*x + c*x**2)
    F = m*x**7/(7*c) + x**6*(-b*m + c*l)/(6*c**2) + x**5*(b**2*m + c**2*k - c*(a*m + b*l))/(5*c**3) + x**4*(-b**3*m + b*c*(2*a*m + b*l) + c**3*j - c**2*(a*l + b*k))/(4*c**4) + x**3*(b**4*m - b**2*c*(3*a*m + b*l) + c**4*h - c**3*(a*k + b*j) + c**2*(a**2*m + 2*a*b*l + b**2*k))/(3*c**5) + x**2*(-b**5*m + b**3*c*(4*a*m + b*l) - b*c**2*(3*a**2*m + 3*a*b*l + b**2*k) + c**5*g - c**4*(a*j + b*h) + c**3*(a**2*l + 2*a*b*k + b**2*j))/(2*c**6) + x*(b**6*m - b**4*c*(5*a*m + b*l) + b**2*c**2*(6*a**2*m + 4*a*b*l + b**2*k) + c**6*f - c**5*(a*h + b*g) + c**4*(a**2*k + 2*a*b*j + b**2*h) - c**3*(a**3*m + 3*a**2*b*l + 3*a*b**2*k + b**3*j))/c**7 + (-b**7*m + b**5*c*(6*a*m + b*l) - b**3*c**2*(10*a**2*m + 5*a*b*l + b**2*k) + b*c**3*(4*a**3*m + 6*a**2*b*l + 4*a*b**2*k + b**3*j) + c**7*e - c**6*(a*g + b*f) + c**5*(a**2*j + 2*a*b*h + b**2*g) - c**4*(a**3*l + 3*a**2*b*k + 3*a*b**2*j + b**3*h))*log(a + b*x + c*x**2)/(2*c**8) - (b**8*m - b**6*c*(8*a*m + b*l) + b**4*c**2*(20*a**2*m + 7*a*b*l + b**2*k) - b**2*c**3*(16*a**3*m + 14*a**2*b*l + 6*a*b**2*k + b**3*j) + 2*c**8*d - c**7*(2*a*f + b*e) + c**6*(2*a**2*h + 3*a*b*g + b**2*f) - c**5*(2*a**3*k + 5*a**2*b*j + 4*a*b**2*h + b**3*g) + c**4*(2*a**4*m + 7*a**3*b*l + 9*a**2*b**2*k + 5*a*b**3*j + b**4*h))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**8*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_371():
    f = (-7*x**2 + 4*x + 1)**3*(x**2 + 5*x + 2)*sqrt(5*x**2 + 2*x + 3)
    F = -343*x**7*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/50 - 50519*x**6*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/2250 + 190939*x**5*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/3000 - 888751*x**4*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/105000 - 90960857*x**3*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/1575000 + 98060877*x**2*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/4375000 + 1045360143*x*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/43750000 - (385799915*x + 77159983)*sqrt(5*x**2 + 2*x + 3)/31250000 - 1968340667*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/131250000 - 540119881*sqrt(5)*asinh(sqrt(14)*(5*x + 1)/14)/78125000
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_372():
    f = (-7*x**2 + 4*x + 1)**2*(x**2 + 5*x + 2)*sqrt(5*x**2 + 2*x + 3)
    F = 49*x**5*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/40 + 989*x**4*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/200 - 25277*x**3*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/3000 - 77509*x**2*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/25000 + 1781669*x*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/250000 - (12608615*x + 2521723)*sqrt(5*x**2 + 2*x + 3)/1250000 + 198439*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/750000 - 17652061*sqrt(5)*asinh(sqrt(14)*(5*x + 1)/14)/3125000
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_373():
    f = (-7*x**2 + 4*x + 1)*(x**2 + 5*x + 2)*sqrt(5*x**2 + 2*x + 3)
    F = -7*x**3*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/30 - 289*x**2*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/250 + 2149*x*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/2500 - (23165*x + 4633)*sqrt(5*x**2 + 2*x + 3)/12500 + 7819*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/7500 - 32431*sqrt(5)*asinh(sqrt(14)*(5*x + 1)/14)/31250
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_374():
    f = (x**2 + 5*x + 2)*sqrt(5*x**2 + 2*x + 3)/(-7*x**2 + 4*x + 1)
    F = (-x/14 + sympy.S(-397)/490)*sqrt(5*x**2 + 2*x + 3) - 8233*sqrt(5)*asinh(sqrt(14)*(5*x + 1)/14)/8575 - 3*sqrt(sympy.S(497041)/11 - 146555*sqrt(11)/11)*atanh((x*(17 - 5*sqrt(11)) - sqrt(11) + 23)/(sqrt(250 - 34*sqrt(11))*sqrt(5*x**2 + 2*x + 3)))/343 + 3*sqrt(146555*sqrt(11)/11 + sympy.S(497041)/11)*atanh((x*(5*sqrt(11) + 17) + sqrt(11) + 23)/(sqrt(34*sqrt(11) + 250)*sqrt(5*x**2 + 2*x + 3)))/343
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_375():
    f = (x**2 + 5*x + 2)*sqrt(5*x**2 + 2*x + 3)/(-7*x**2 + 4*x + 1)**2
    F = (183*x + 9)*sqrt(5*x**2 + 2*x + 3)/(-1078*x**2 + 616*x + 154) + sqrt(5)*asinh(sqrt(14)*(5*x + 1)/14)/49 - sqrt(3557521*sqrt(11)/127 + sympy.S(325022311)/1397)*atanh((x*(17 - 5*sqrt(11)) - sqrt(11) + 23)/(sqrt(250 - 34*sqrt(11))*sqrt(5*x**2 + 2*x + 3)))/2156 + sqrt(sympy.S(325022311)/1397 - 3557521*sqrt(11)/127)*atanh((x*(5*sqrt(11) + 17) + sqrt(11) + 23)/(sqrt(34*sqrt(11) + 250)*sqrt(5*x**2 + 2*x + 3)))/2156
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_376():
    f = (x**2 + 5*x + 2)*sqrt(5*x**2 + 2*x + 3)/(-7*x**2 + 4*x + 1)**3
    F = -(272941 - 813113*x)*sqrt(5*x**2 + 2*x + 3)/(-12047728*x**2 + 6884416*x + 1721104) + (183*x + 9)*sqrt(5*x**2 + 2*x + 3)/(308*(-7*x**2 + 4*x + 1)**2) - sqrt(sympy.S(6492253020949)/1397 - 1079924461*sqrt(11)/127)*atanh((x*(17 - 5*sqrt(11)) - sqrt(11) + 23)/(sqrt(250 - 34*sqrt(11))*sqrt(5*x**2 + 2*x + 3)))/491744 + sqrt(1079924461*sqrt(11)/127 + sympy.S(6492253020949)/1397)*atanh((x*(5*sqrt(11) + 17) + sqrt(11) + 23)/(sqrt(34*sqrt(11) + 250)*sqrt(5*x**2 + 2*x + 3)))/491744
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_377():
    f = (-7*x**2 + 4*x + 1)**3*(x**2 + 5*x + 2)*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)
    F = -343*x**7*(5*x**2 + 2*x + 3)**(sympy.S(5)/2)/60 - 61103*x**6*(5*x**2 + 2*x + 3)**(sympy.S(5)/2)/3300 + 1031177*x**5*(5*x**2 + 2*x + 3)**(sympy.S(5)/2)/20625 - 796559*x**4*(5*x**2 + 2*x + 3)**(sympy.S(5)/2)/123750 - 190236913*x**3*(5*x**2 + 2*x + 3)**(sympy.S(5)/2)/4950000 + 2173004363*x**2*(5*x**2 + 2*x + 3)**(sympy.S(5)/2)/173250000 + 837379699*x*(5*x**2 + 2*x + 3)**(sympy.S(5)/2)/72187500 - (114202995*x + 22840599)*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/62500000 - (2398262895*x + 479652579)*sqrt(5*x**2 + 2*x + 3)/312500000 - 6133820867*(5*x**2 + 2*x + 3)**(sympy.S(5)/2)/1203125000 - 3357568053*sqrt(5)*asinh(sqrt(14)*(5*x + 1)/14)/781250000
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_378():
    f = (-7*x**2 + 4*x + 1)**2*(x**2 + 5*x + 2)*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)
    F = 49*x**5*(5*x**2 + 2*x + 3)**(sympy.S(5)/2)/50 + 581*x**4*(5*x**2 + 2*x + 3)**(sympy.S(5)/2)/150 - 18379*x**3*(5*x**2 + 2*x + 3)**(sympy.S(5)/2)/3000 - 219271*x**2*(5*x**2 + 2*x + 3)**(sympy.S(5)/2)/105000 + 86721*x*(5*x**2 + 2*x + 3)**(sympy.S(5)/2)/21875 - (3452805*x + 690561)*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/1250000 - (72508905*x + 14501781)*sqrt(5*x**2 + 2*x + 3)/6250000 + 505667*(5*x**2 + 2*x + 3)**(sympy.S(5)/2)/2187500 - 101512467*sqrt(5)*asinh(sqrt(14)*(5*x + 1)/14)/15625000
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_379():
    f = (-7*x**2 + 4*x + 1)*(x**2 + 5*x + 2)*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)
    F = -7*x**3*(5*x**2 + 2*x + 3)**(sympy.S(5)/2)/40 - 1163*x**2*(5*x**2 + 2*x + 3)**(sympy.S(5)/2)/1400 + 2809*x*(5*x**2 + 2*x + 3)**(sympy.S(5)/2)/5250 - (91985*x + 18397)*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/150000 - (643895*x + 128779)*sqrt(5*x**2 + 2*x + 3)/250000 + 149509*(5*x**2 + 2*x + 3)**(sympy.S(5)/2)/262500 - 901453*sqrt(5)*asinh(sqrt(14)*(5*x + 1)/14)/625000
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_380():
    f = (x**2 + 5*x + 2)*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/(-7*x**2 + 4*x + 1)
    F = -(x/28 + sympy.S(267)/980)*(5*x**2 + 2*x + 3)**(sympy.S(3)/2) - (588315*x + 1714863)*sqrt(5*x**2 + 2*x + 3)/240100 - 34425687*sqrt(5)*asinh(sqrt(14)*(5*x + 1)/14)/4201750 - 6*sqrt(sympy.S(16197805214)/11 - 4868244470*sqrt(11)/11)*atanh((x*(17 - 5*sqrt(11)) - sqrt(11) + 23)/(sqrt(250 - 34*sqrt(11))*sqrt(5*x**2 + 2*x + 3)))/16807 + 6*sqrt(4868244470*sqrt(11)/11 + sympy.S(16197805214)/11)*atanh((x*(5*sqrt(11) + 17) + sqrt(11) + 23)/(sqrt(34*sqrt(11) + 250)*sqrt(5*x**2 + 2*x + 3)))/16807
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_381():
    f = (x**2 + 5*x + 2)*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/(-7*x**2 + 4*x + 1)**2
    F = (183*x + 9)*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/(-1078*x**2 + 616*x + 154) + (3395*x + 5826)*sqrt(5*x**2 + 2*x + 3)/3773 + 16691*sqrt(5)*asinh(sqrt(14)*(5*x + 1)/14)/12005 - sqrt(sympy.S(52175400311)/22 - 1195943321*sqrt(11)/2)*atanh((x*(17 - 5*sqrt(11)) - sqrt(11) + 23)/(sqrt(250 - 34*sqrt(11))*sqrt(5*x**2 + 2*x + 3)))/26411 - sqrt(1195943321*sqrt(11)/2 + sympy.S(52175400311)/22)*atanh((x*(5*sqrt(11) + 17) + sqrt(11) + 23)/(sqrt(34*sqrt(11) + 250)*sqrt(5*x**2 + 2*x + 3)))/26411
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_382():
    f = (x**2 + 5*x + 2)*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/(-7*x**2 + 4*x + 1)**3
    F = -(9495 - 37088*x)*sqrt(5*x**2 + 2*x + 3)/(-166012*x**2 + 94864*x + 23716) + (183*x + 9)*(5*x**2 + 2*x + 3)**(sympy.S(3)/2)/(308*(-7*x**2 + 4*x + 1)**2) - 5*sqrt(5)*asinh(sqrt(14)*(5*x + 1)/14)/343 - sqrt(sympy.S(62294197250171)/2794 - 189585522005*sqrt(11)/254)*atanh((x*(17 - 5*sqrt(11)) - sqrt(11) + 23)/(sqrt(250 - 34*sqrt(11))*sqrt(5*x**2 + 2*x + 3)))/332024 + sqrt(189585522005*sqrt(11)/254 + sympy.S(62294197250171)/2794)*atanh((x*(5*sqrt(11) + 17) + sqrt(11) + 23)/(sqrt(34*sqrt(11) + 250)*sqrt(5*x**2 + 2*x + 3)))/332024
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_383():
    f = (-7*x**2 + 4*x + 1)**3*(x**2 + 5*x + 2)/sqrt(5*x**2 + 2*x + 3)
    F = -343*x**7*sqrt(5*x**2 + 2*x + 3)/40 - 1141*x**6*sqrt(5*x**2 + 2*x + 3)/40 + 26159*x**5*sqrt(5*x**2 + 2*x + 3)/300 - 47807*x**4*sqrt(5*x**2 + 2*x + 3)/3750 - 5160533*x**3*sqrt(5*x**2 + 2*x + 3)/50000 + 40722851*x**2*sqrt(5*x**2 + 2*x + 3)/750000 + 5793077*x*sqrt(5*x**2 + 2*x + 3)/75000 - 16515809*sqrt(5*x**2 + 2*x + 3)/156250 - 77513689*sqrt(5)*asinh(sqrt(14)*(5*x + 1)/14)/3125000
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_384():
    f = (-7*x**2 + 4*x + 1)**2*(x**2 + 5*x + 2)/sqrt(5*x**2 + 2*x + 3)
    F = 49*x**5*sqrt(5*x**2 + 2*x + 3)/30 + 5131*x**4*sqrt(5*x**2 + 2*x + 3)/750 - 33259*x**3*sqrt(5*x**2 + 2*x + 3)/2500 - 207427*x**2*sqrt(5*x**2 + 2*x + 3)/37500 + 36073*x*sqrt(5*x**2 + 2*x + 3)/1875 - 22053*sqrt(5*x**2 + 2*x + 3)/31250 - 1719097*sqrt(5)*asinh(sqrt(14)*(5*x + 1)/14)/156250
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_385():
    f = (-7*x**2 + 4*x + 1)*(x**2 + 5*x + 2)/sqrt(5*x**2 + 2*x + 3)
    F = -7*x**3*sqrt(5*x**2 + 2*x + 3)/20 - 571*x**2*sqrt(5*x**2 + 2*x + 3)/300 + 59*x*sqrt(5*x**2 + 2*x + 3)/30 + 463*sqrt(5*x**2 + 2*x + 3)/125 - 1901*sqrt(5)*asinh(sqrt(14)*(5*x + 1)/14)/1250
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_386():
    f = (x**2 + 5*x + 2)/((-7*x**2 + 4*x + 1)*sqrt(5*x**2 + 2*x + 3))
    F = -sqrt(5)*asinh(sqrt(14)*(5*x + 1)/14)/35 - 3*sqrt(sympy.S(4091)/2794 - 1055*sqrt(11)/2794)*atanh((x*(17 - 5*sqrt(11)) - sqrt(11) + 23)/(sqrt(250 - 34*sqrt(11))*sqrt(5*x**2 + 2*x + 3)))/14 + 3*sqrt(1055*sqrt(11)/2794 + sympy.S(4091)/2794)*atanh((x*(5*sqrt(11) + 17) + sqrt(11) + 23)/(sqrt(34*sqrt(11) + 250)*sqrt(5*x**2 + 2*x + 3)))/14
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_387():
    f = (x**2 + 5*x + 2)/((-7*x**2 + 4*x + 1)**2*sqrt(5*x**2 + 2*x + 3))
    F = -(120 - 1113*x)*sqrt(5*x**2 + 2*x + 3)/(-39116*x**2 + 22352*x + 5588) - sqrt(1275971*sqrt(11)/254 + sympy.S(3027900955)/2794)*atanh((x*(17 - 5*sqrt(11)) - sqrt(11) + 23)/(sqrt(250 - 34*sqrt(11))*sqrt(5*x**2 + 2*x + 3)))/11176 + sqrt(sympy.S(3027900955)/2794 - 1275971*sqrt(11)/254)*atanh((x*(5*sqrt(11) + 17) + sqrt(11) + 23)/(sqrt(34*sqrt(11) + 250)*sqrt(5*x**2 + 2*x + 3)))/11176
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_388():
    f = (x**2 + 5*x + 2)/((-7*x**2 + 4*x + 1)**3*sqrt(5*x**2 + 2*x + 3))
    F = -(120 - 1113*x)*sqrt(5*x**2 + 2*x + 3)/(11176*(-7*x**2 + 4*x + 1)**2) - (2868383 - 8325590*x)*sqrt(5*x**2 + 2*x + 3)/(-437160416*x**2 + 249805952*x + 62451488) - (275591617 - 17771075*sqrt(11))*atanh((x*(17 - 5*sqrt(11)) - sqrt(11) + 23)/(sqrt(250 - 34*sqrt(11))*sqrt(5*x**2 + 2*x + 3)))/(124902976*sqrt(2750 - 374*sqrt(11))) + (17771075*sqrt(11) + 275591617)*atanh((x*(5*sqrt(11) + 17) + sqrt(11) + 23)/(sqrt(34*sqrt(11) + 250)*sqrt(5*x**2 + 2*x + 3)))/(124902976*sqrt(374*sqrt(11) + 2750))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_389():
    f = (-7*x**2 + 4*x + 1)**3*(x**2 + 5*x + 2)/(5*x**2 + 2*x + 3)**(sympy.S(3)/2)
    F = -343*x**5*sqrt(5*x**2 + 2*x + 3)/150 - 25921*x**4*sqrt(5*x**2 + 2*x + 3)/3750 + 393659*x**3*sqrt(5*x**2 + 2*x + 3)/12500 - 2583293*x**2*sqrt(5*x**2 + 2*x + 3)/187500 - 3192602*x*sqrt(5*x**2 + 2*x + 3)/46875 + (97964912 - 85411472*x)/(546875*sqrt(5*x**2 + 2*x + 3)) + 15715799*sqrt(5*x**2 + 2*x + 3)/156250 + 50047657*sqrt(5)*asinh(sqrt(14)*(5*x + 1)/14)/781250
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_390():
    f = (-7*x**2 + 4*x + 1)**2*(x**2 + 5*x + 2)/(5*x**2 + 2*x + 3)**(sympy.S(3)/2)
    F = 49*x**3*sqrt(5*x**2 + 2*x + 3)/100 + 203*x**2*sqrt(5*x**2 + 2*x + 3)/100 - 8749*x*sqrt(5*x**2 + 2*x + 3)/1250 - (1092816*x + 103864)/(21875*sqrt(5*x**2 + 2*x + 3)) - 5086*sqrt(5*x**2 + 2*x + 3)/3125 + 89583*sqrt(5)*asinh(sqrt(14)*(5*x + 1)/14)/6250
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_391():
    f = (-7*x**2 + 4*x + 1)*(x**2 + 5*x + 2)/(5*x**2 + 2*x + 3)**(sympy.S(3)/2)
    F = -7*x*sqrt(5*x**2 + 2*x + 3)/50 - (4898*x + 4642)/(875*sqrt(5*x**2 + 2*x + 3)) - 261*sqrt(5*x**2 + 2*x + 3)/250 + 149*sqrt(5)*asinh(sqrt(14)*(5*x + 1)/14)/125
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_392():
    f = (x**2 + 5*x + 2)/((-7*x**2 + 4*x + 1)*(5*x**2 + 2*x + 3)**(sympy.S(3)/2))
    F = -(131 - 605*x)/(3556*sqrt(5*x**2 + 2*x + 3)) - 3*sqrt(sympy.S(281693)/1397 - 25015*sqrt(11)/1397)*atanh((x*(17 - 5*sqrt(11)) - sqrt(11) + 23)/(sqrt(250 - 34*sqrt(11))*sqrt(5*x**2 + 2*x + 3)))/1016 + 3*sqrt(25015*sqrt(11)/1397 + sympy.S(281693)/1397)*atanh((x*(5*sqrt(11) + 17) + sqrt(11) + 23)/(sqrt(34*sqrt(11) + 250)*sqrt(5*x**2 + 2*x + 3)))/1016
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_393():
    f = (x**2 + 5*x + 2)/((-7*x**2 + 4*x + 1)**2*(5*x**2 + 2*x + 3)**(sympy.S(3)/2))
    F = -(120 - 1113*x)/((-39116*x**2 + 22352*x + 5588)*sqrt(5*x**2 + 2*x + 3)) - (22755*x + 76567)/(19870928*sqrt(5*x**2 + 2*x + 3)) - (3790801 - 36008*sqrt(11))*atanh((x*(17 - 5*sqrt(11)) - sqrt(11) + 23)/(sqrt(250 - 34*sqrt(11))*sqrt(5*x**2 + 2*x + 3)))/(2838704*sqrt(2750 - 374*sqrt(11))) + (36008*sqrt(11) + 3790801)*atanh((x*(5*sqrt(11) + 17) + sqrt(11) + 23)/(sqrt(34*sqrt(11) + 250)*sqrt(5*x**2 + 2*x + 3)))/(2838704*sqrt(374*sqrt(11) + 2750))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_394():
    f = (x**2 + 5*x + 2)/((-7*x**2 + 4*x + 1)**3*(5*x**2 + 2*x + 3)**(sympy.S(3)/2))
    F = -(120 - 1113*x)/(11176*(-7*x**2 + 4*x + 1)**2*sqrt(5*x**2 + 2*x + 3)) - (2701733 - 9148874*x)/((-437160416*x**2 + 249805952*x + 62451488)*sqrt(5*x**2 + 2*x + 3)) - (5593656875*x + 2306853905)/(222077491328*sqrt(5*x**2 + 2*x + 3)) - (19550020168 - 594061265*sqrt(11))*atanh((x*(17 - 5*sqrt(11)) - sqrt(11) + 23)/(sqrt(250 - 34*sqrt(11))*sqrt(5*x**2 + 2*x + 3)))/(31725355904*sqrt(2750 - 374*sqrt(11))) + (594061265*sqrt(11) + 19550020168)*atanh((x*(5*sqrt(11) + 17) + sqrt(11) + 23)/(sqrt(34*sqrt(11) + 250)*sqrt(5*x**2 + 2*x + 3)))/(31725355904*sqrt(374*sqrt(11) + 2750))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_395():
    f = (A + C*x**2)*(a + c*x**2)**p*(d + f*x**2)**q
    F = A*x*(a + c*x**2)**p*(d + f*x**2)**q*appellf1(sympy.S.Half, -p, -q, sympy.S(3)/2, -c*x**2/a, -f*x**2/d)/((1 + c*x**2/a)**p*(1 + f*x**2/d)**q) + C*x**3*(a + c*x**2)**p*(d + f*x**2)**q*appellf1(sympy.S(3)/2, -p, -q, sympy.S(5)/2, -c*x**2/a, -f*x**2/d)/(3*(1 + c*x**2/a)**p*(1 + f*x**2/d)**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_396():
    f = (A + B*x)*(a + c*x**2)**p*(d + f*x**2)**q
    F = A*x*(a + c*x**2)**p*(d + f*x**2)**q*appellf1(sympy.S.Half, -p, -q, sympy.S(3)/2, -c*x**2/a, -f*x**2/d)/((1 + c*x**2/a)**p*(1 + f*x**2/d)**q) + B*(a + c*x**2)**(p + 1)*(d + f*x**2)**q*hyper((-q, p + 1), (p + 2,), -f*(a + c*x**2)/(-a*f + c*d))/(2*c*(c*(d + f*x**2)/(-a*f + c*d))**q*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_1_Quadratic_1_2_1_9_P_x_d_plus_e_x_pow_m_a_plus_b_x_plus_c_x_pow_2_pow_p_397():
    f = (a + c*x**2)**p*(d + f*x**2)**q*(A + B*x + C*x**2)
    F = A*x*(a + c*x**2)**p*(d + f*x**2)**q*appellf1(sympy.S.Half, -p, -q, sympy.S(3)/2, -c*x**2/a, -f*x**2/d)/((1 + c*x**2/a)**p*(1 + f*x**2/d)**q) + B*(a + c*x**2)**(p + 1)*(d + f*x**2)**q*hyper((-q, p + 1), (p + 2,), -f*(a + c*x**2)/(-a*f + c*d))/(2*c*(c*(d + f*x**2)/(-a*f + c*d))**q*(p + 1)) + C*x**3*(a + c*x**2)**p*(d + f*x**2)**q*appellf1(sympy.S(3)/2, -p, -q, sympy.S(5)/2, -c*x**2/a, -f*x**2/d)/(3*(1 + c*x**2/a)**p*(1 + f*x**2/d)**q)
    assert integrate(f, x) == F

