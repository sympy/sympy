"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.2 Trinomial products/1.2.3 General/1.2.3.2 (d x)^m (a+b x^n+c x^(2 n))^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n, p = symbols('a b c d e f m n p')

def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_1():
    f = (a*x**3 + b*x**6)**(sympy.S(5)/3)
    F = -3*a*(a*x**3 + b*x**6)**(sympy.S(8)/3)/(88*b**2*x**8) + (a*x**3 + b*x**6)**(sympy.S(8)/3)/(11*b*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_2():
    f = (a*x**3 + b*x**6)**(sympy.S(2)/3)
    F = (a*x**3 + b*x**6)**(sympy.S(5)/3)/(5*b*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_3():
    f = (a*x**3 + b*x**6)**(sympy.S(-2)/3)
    F = -(a*x**3 + b*x**6)**(sympy.S(1)/3)/(a*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_4():
    f = (a*x**3 + b*x**6)**(sympy.S(-5)/3)
    F = 1/(2*a*x**2*(a*x**3 + b*x**6)**(sympy.S(2)/3)) - 3*(a*x**3 + b*x**6)**(sympy.S(1)/3)/(4*a**2*x**5) + 9*b*(a*x**3 + b*x**6)**(sympy.S(1)/3)/(4*a**3*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_5():
    f = 1/(x**6 - x**3)
    F = log(1 - x)/3 - log(x**2 + x + 1)/6 - sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/3 + 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_6():
    f = x**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)
    F = a*x**6*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(6*a + 6*b*x**3) + b*x**9*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(9*a + 9*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_7():
    f = x**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)
    F = a*x**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(5*a + 5*b*x**3) + b*x**8*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(8*a + 8*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_8():
    f = x**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)
    F = a*x**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*a + 4*b*x**3) + b*x**7*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*a + 7*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_9():
    f = x**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)
    F = (a + b*x**3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(6*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_10():
    f = x*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)
    F = a*x**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*a + 2*b*x**3) + b*x**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(5*a + 5*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_11():
    f = sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)
    F = a*x*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(a + b*x**3) + b*x**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*a + 4*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_12():
    f = sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/x
    F = a*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)*log(x)/(a + b*x**3) + b*x**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(3*a + 3*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_13():
    f = sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/x**2
    F = -a*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x*(a + b*x**3)) + b*x**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*a + 2*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_14():
    f = sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/x**3
    F = -a*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*x**2*(a + b*x**3)) + b*x*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(a + b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_15():
    f = sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/x**4
    F = -a*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(3*x**3*(a + b*x**3)) + b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)*log(x)/(a + b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_16():
    f = sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/x**5
    F = -a*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*x**4*(a + b*x**3)) - b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_17():
    f = sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/x**6
    F = -a*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(5*x**5*(a + b*x**3)) - b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*x**2*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_18():
    f = sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/x**7
    F = -a*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(6*x**6*(a + b*x**3)) - b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(3*x**3*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_19():
    f = sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/x**8
    F = -a*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*x**7*(a + b*x**3)) - b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*x**4*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_20():
    f = sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/x**9
    F = -a*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(8*x**8*(a + b*x**3)) - b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(5*x**5*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_21():
    f = sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/x**10
    F = -a*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(9*x**9*(a + b*x**3)) - b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(6*x**6*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_22():
    f = sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/x**11
    F = -a*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(10*x**10*(a + b*x**3)) - b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*x**7*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_23():
    f = x**9*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)
    F = a**3*x**10*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(10*a + 10*b*x**3) + 3*a**2*b*x**13*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(13*a + 13*b*x**3) + 3*a*b**2*x**16*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(16*a + 16*b*x**3) + b**3*x**19*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(19*a + 19*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_24():
    f = x**8*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)
    F = a**2*(a + b*x**3)**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(12*b**3) - 2*a*(a + b*x**3)**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(15*b**3) + (a + b*x**3)**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(18*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_25():
    f = x**7*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)
    F = a**3*x**8*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(8*a + 8*b*x**3) + 3*a**2*b*x**11*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(11*a + 11*b*x**3) + 3*a*b**2*x**14*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(14*a + 14*b*x**3) + b**3*x**17*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(17*a + 17*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_26():
    f = x**6*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)
    F = a**3*x**7*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*a + 7*b*x**3) + 3*a**2*b*x**10*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(10*a + 10*b*x**3) + 3*a*b**2*x**13*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(13*a + 13*b*x**3) + b**3*x**16*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(16*a + 16*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_27():
    f = x**5*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)
    F = -a*(a + b*x**3)**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(12*b**2) + (a + b*x**3)**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(15*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_28():
    f = x**4*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)
    F = a**3*x**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(5*a + 5*b*x**3) + 3*a**2*b*x**8*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(8*a + 8*b*x**3) + 3*a*b**2*x**11*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(11*a + 11*b*x**3) + b**3*x**14*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(14*a + 14*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_29():
    f = x**3*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)
    F = a**3*x**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*a + 4*b*x**3) + 3*a**2*b*x**7*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*a + 7*b*x**3) + 3*a*b**2*x**10*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(10*a + 10*b*x**3) + b**3*x**13*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(13*a + 13*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_30():
    f = x**2*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)
    F = (a + b*x**3)*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/(12*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_31():
    f = x*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)
    F = a**3*x**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*a + 2*b*x**3) + 3*a**2*b*x**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(5*a + 5*b*x**3) + 3*a*b**2*x**8*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(8*a + 8*b*x**3) + b**3*x**11*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(11*a + 11*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_32():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)
    F = a**3*x*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/(a + b*x**3)**3 + 3*a**2*b*x**4*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/(4*(a + b*x**3)**3) + 3*a*b**2*x**7*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/(7*(a + b*x**3)**3) + b**3*x**10*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/(10*(a + b*x**3)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_33():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/x
    F = a**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)*log(x)/(a + b*x**3) + a**2*b*x**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(a + b*x**3) + a*b**2*x**6*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*a + 2*b*x**3) + b**3*x**9*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(9*a + 9*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_34():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/x**2
    F = -a**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x*(a + b*x**3)) + 3*a**2*b*x**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*a + 2*b*x**3) + 3*a*b**2*x**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(5*a + 5*b*x**3) + b**3*x**8*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(8*a + 8*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_35():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/x**3
    F = -a**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*x**2*(a + b*x**3)) + 3*a**2*b*x*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(a + b*x**3) + 3*a*b**2*x**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*a + 4*b*x**3) + b**3*x**7*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*a + 7*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_36():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/x**4
    F = -a**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(3*x**3*(a + b*x**3)) + 3*a**2*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)*log(x)/(a + b*x**3) + a*b**2*x**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(a + b*x**3) + b**3*x**6*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(6*a + 6*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_37():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/x**5
    F = -a**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*x**4*(a + b*x**3)) - 3*a**2*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x*(a + b*x**3)) + 3*a*b**2*x**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*a + 2*b*x**3) + b**3*x**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(5*a + 5*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_38():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/x**6
    F = -a**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(5*x**5*(a + b*x**3)) - 3*a**2*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*x**2*(a + b*x**3)) + 3*a*b**2*x*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(a + b*x**3) + b**3*x**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*a + 4*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_39():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/x**7
    F = -a**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(6*x**6*(a + b*x**3)) - a**2*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x**3*(a + b*x**3)) + 3*a*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)*log(x)/(a + b*x**3) + b**3*x**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(3*a + 3*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_40():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/x**8
    F = -a**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*x**7*(a + b*x**3)) - 3*a**2*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*x**4*(a + b*x**3)) - 3*a*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x*(a + b*x**3)) + b**3*x**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*a + 2*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_41():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/x**9
    F = -a**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(8*x**8*(a + b*x**3)) - 3*a**2*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(5*x**5*(a + b*x**3)) - 3*a*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*x**2*(a + b*x**3)) + b**3*x*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(a + b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_42():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/x**10
    F = -a**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(9*x**9*(a + b*x**3)) - a**2*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*x**6*(a + b*x**3)) - a*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x**3*(a + b*x**3)) + b**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)*log(x)/(a + b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_43():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/x**11
    F = -a**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(10*x**10*(a + b*x**3)) - 3*a**2*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*x**7*(a + b*x**3)) - 3*a*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*x**4*(a + b*x**3)) - b**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_44():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/x**12
    F = -a**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(11*x**11*(a + b*x**3)) - 3*a**2*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(8*x**8*(a + b*x**3)) - 3*a*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(5*x**5*(a + b*x**3)) - b**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*x**2*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_45():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/x**13
    F = -(a + b*x**3)**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(12*a*x**12)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_46():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/x**14
    F = -a**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(13*x**13*(a + b*x**3)) - 3*a**2*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(10*x**10*(a + b*x**3)) - 3*a*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*x**7*(a + b*x**3)) - b**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*x**4*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_47():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/x**15
    F = -a**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(14*x**14*(a + b*x**3)) - 3*a**2*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(11*x**11*(a + b*x**3)) - 3*a*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(8*x**8*(a + b*x**3)) - b**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(5*x**5*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_48():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/x**16
    F = -(a + b*x**3)**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(15*a*x**15) + b*(a + b*x**3)**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(60*a**2*x**12)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_49():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)/x**17
    F = -a**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(16*x**16*(a + b*x**3)) - 3*a**2*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(13*x**13*(a + b*x**3)) - 3*a*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(10*x**10*(a + b*x**3)) - b**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*x**7*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_50():
    f = x**13*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = a**5*x**14*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(14*a + 14*b*x**3) + 5*a**4*b*x**17*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(17*a + 17*b*x**3) + a**3*b**2*x**20*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*a + 2*b*x**3) + 10*a**2*b**3*x**23*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(23*a + 23*b*x**3) + 5*a*b**4*x**26*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(26*a + 26*b*x**3) + b**5*x**29*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(29*a + 29*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_51():
    f = x**12*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = a**5*x**13*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(13*a + 13*b*x**3) + 5*a**4*b*x**16*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(16*a + 16*b*x**3) + 10*a**3*b**2*x**19*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(19*a + 19*b*x**3) + 5*a**2*b**3*x**22*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(11*a + 11*b*x**3) + a*b**4*x**25*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(5*a + 5*b*x**3) + b**5*x**28*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(28*a + 28*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_52():
    f = x**11*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = -a**3*(a + b*x**3)**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(18*b**4) + a**2*(a + b*x**3)**6*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*b**4) - a*(a + b*x**3)**7*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(8*b**4) + (a + b*x**3)**8*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(27*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_53():
    f = x**10*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = a**5*x**11*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(11*a + 11*b*x**3) + 5*a**4*b*x**14*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(14*a + 14*b*x**3) + 10*a**3*b**2*x**17*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(17*a + 17*b*x**3) + a**2*b**3*x**20*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*a + 2*b*x**3) + 5*a*b**4*x**23*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(23*a + 23*b*x**3) + b**5*x**26*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(26*a + 26*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_54():
    f = x**9*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = a**5*x**10*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(10*a + 10*b*x**3) + 5*a**4*b*x**13*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(13*a + 13*b*x**3) + 5*a**3*b**2*x**16*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(8*a + 8*b*x**3) + 10*a**2*b**3*x**19*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(19*a + 19*b*x**3) + 5*a*b**4*x**22*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(22*a + 22*b*x**3) + b**5*x**25*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(25*a + 25*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_55():
    f = x**8*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = a**2*(a + b*x**3)**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(18*b**3) - 2*a*(a + b*x**3)**6*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(21*b**3) + (a + b*x**3)**7*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(24*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_56():
    f = x**7*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = a**5*x**8*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(8*a + 8*b*x**3) + 5*a**4*b*x**11*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(11*a + 11*b*x**3) + 5*a**3*b**2*x**14*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*a + 7*b*x**3) + 10*a**2*b**3*x**17*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(17*a + 17*b*x**3) + a*b**4*x**20*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*a + 4*b*x**3) + b**5*x**23*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(23*a + 23*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_57():
    f = x**6*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = a**5*x**7*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*a + 7*b*x**3) + a**4*b*x**10*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*a + 2*b*x**3) + 10*a**3*b**2*x**13*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(13*a + 13*b*x**3) + 5*a**2*b**3*x**16*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(8*a + 8*b*x**3) + 5*a*b**4*x**19*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(19*a + 19*b*x**3) + b**5*x**22*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(22*a + 22*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_58():
    f = x**5*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = -a*(a + b*x**3)**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(18*b**2) + (a + b*x**3)**6*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(21*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_59():
    f = x**4*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = a**5*x**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(5*a + 5*b*x**3) + 5*a**4*b*x**8*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(8*a + 8*b*x**3) + 10*a**3*b**2*x**11*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(11*a + 11*b*x**3) + 5*a**2*b**3*x**14*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*a + 7*b*x**3) + 5*a*b**4*x**17*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(17*a + 17*b*x**3) + b**5*x**20*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(20*a + 20*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_60():
    f = x**3*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = a**5*x**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*a + 4*b*x**3) + 5*a**4*b*x**7*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*a + 7*b*x**3) + a**3*b**2*x**10*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(a + b*x**3) + 10*a**2*b**3*x**13*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(13*a + 13*b*x**3) + 5*a*b**4*x**16*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(16*a + 16*b*x**3) + b**5*x**19*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(19*a + 19*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_61():
    f = x**2*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = (a + b*x**3)*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/(18*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_62():
    f = x*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = a**5*x**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*a + 2*b*x**3) + a**4*b*x**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(a + b*x**3) + 5*a**3*b**2*x**8*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*a + 4*b*x**3) + 10*a**2*b**3*x**11*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(11*a + 11*b*x**3) + 5*a*b**4*x**14*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(14*a + 14*b*x**3) + b**5*x**17*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(17*a + 17*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_63():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = a**5*x*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/(a + b*x**3)**5 + 5*a**4*b*x**4*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/(4*(a + b*x**3)**5) + 10*a**3*b**2*x**7*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/(7*(a + b*x**3)**5) + a**2*b**3*x**10*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/(a + b*x**3)**5 + 5*a*b**4*x**13*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/(13*(a + b*x**3)**5) + b**5*x**16*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/(16*(a + b*x**3)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_64():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x
    F = a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)*log(x)/(a + b*x**3) + 5*a**4*b*x**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(3*a + 3*b*x**3) + 5*a**3*b**2*x**6*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(3*a + 3*b*x**3) + 10*a**2*b**3*x**9*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(9*a + 9*b*x**3) + 5*a*b**4*x**12*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(12*a + 12*b*x**3) + b**5*x**15*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(15*a + 15*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_65():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**2
    F = -a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x*(a + b*x**3)) + 5*a**4*b*x**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*a + 2*b*x**3) + 2*a**3*b**2*x**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(a + b*x**3) + 5*a**2*b**3*x**8*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*a + 4*b*x**3) + 5*a*b**4*x**11*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(11*a + 11*b*x**3) + b**5*x**14*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(14*a + 14*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_66():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**3
    F = -a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*x**2*(a + b*x**3)) + 5*a**4*b*x*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(a + b*x**3) + 5*a**3*b**2*x**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*a + 2*b*x**3) + 10*a**2*b**3*x**7*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*a + 7*b*x**3) + a*b**4*x**10*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*a + 2*b*x**3) + b**5*x**13*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(13*a + 13*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_67():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**4
    F = -a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(3*x**3*(a + b*x**3)) + 5*a**4*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)*log(x)/(a + b*x**3) + 10*a**3*b**2*x**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(3*a + 3*b*x**3) + 5*a**2*b**3*x**6*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(3*a + 3*b*x**3) + 5*a*b**4*x**9*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(9*a + 9*b*x**3) + b**5*x**12*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(12*a + 12*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_68():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**5
    F = -a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*x**4*(a + b*x**3)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x*(a + b*x**3)) + 5*a**3*b**2*x**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(a + b*x**3) + 2*a**2*b**3*x**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(a + b*x**3) + 5*a*b**4*x**8*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(8*a + 8*b*x**3) + b**5*x**11*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(11*a + 11*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_69():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**6
    F = -a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(5*x**5*(a + b*x**3)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*x**2*(a + b*x**3)) + 10*a**3*b**2*x*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(a + b*x**3) + 5*a**2*b**3*x**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*a + 2*b*x**3) + 5*a*b**4*x**7*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*a + 7*b*x**3) + b**5*x**10*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(10*a + 10*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_70():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**7
    F = -a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(6*x**6*(a + b*x**3)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(3*x**3*(a + b*x**3)) + 10*a**3*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)*log(x)/(a + b*x**3) + 10*a**2*b**3*x**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(3*a + 3*b*x**3) + 5*a*b**4*x**6*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(6*a + 6*b*x**3) + b**5*x**9*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(9*a + 9*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_71():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**8
    F = -a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*x**7*(a + b*x**3)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*x**4*(a + b*x**3)) - 10*a**3*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x*(a + b*x**3)) + 5*a**2*b**3*x**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(a + b*x**3) + a*b**4*x**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(a + b*x**3) + b**5*x**8*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(8*a + 8*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_72():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**9
    F = -a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(8*x**8*(a + b*x**3)) - a**4*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x**5*(a + b*x**3)) - 5*a**3*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x**2*(a + b*x**3)) + 10*a**2*b**3*x*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(a + b*x**3) + 5*a*b**4*x**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*a + 4*b*x**3) + b**5*x**7*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*a + 7*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_73():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**10
    F = -a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(9*x**9*(a + b*x**3)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(6*x**6*(a + b*x**3)) - 10*a**3*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(3*x**3*(a + b*x**3)) + 10*a**2*b**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)*log(x)/(a + b*x**3) + 5*a*b**4*x**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(3*a + 3*b*x**3) + b**5*x**6*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(6*a + 6*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_74():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**11
    F = -a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(10*x**10*(a + b*x**3)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*x**7*(a + b*x**3)) - 5*a**3*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*x**4*(a + b*x**3)) - 10*a**2*b**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x*(a + b*x**3)) + 5*a*b**4*x**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*a + 2*b*x**3) + b**5*x**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(5*a + 5*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_75():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**12
    F = -a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(11*x**11*(a + b*x**3)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(8*x**8*(a + b*x**3)) - 2*a**3*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x**5*(a + b*x**3)) - 5*a**2*b**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x**2*(a + b*x**3)) + 5*a*b**4*x*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(a + b*x**3) + b**5*x**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*a + 4*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_76():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**13
    F = -a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(12*x**12*(a + b*x**3)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(9*x**9*(a + b*x**3)) - 5*a**3*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(3*x**6*(a + b*x**3)) - 10*a**2*b**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(3*x**3*(a + b*x**3)) + 5*a*b**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)*log(x)/(a + b*x**3) + b**5*x**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(3*a + 3*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_77():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**14
    F = -a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(13*x**13*(a + b*x**3)) - a**4*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*x**10*(a + b*x**3)) - 10*a**3*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*x**7*(a + b*x**3)) - 5*a**2*b**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*x**4*(a + b*x**3)) - 5*a*b**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x*(a + b*x**3)) + b**5*x**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*a + 2*b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_78():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**15
    F = -a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(14*x**14*(a + b*x**3)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(11*x**11*(a + b*x**3)) - 5*a**3*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*x**8*(a + b*x**3)) - 2*a**2*b**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x**5*(a + b*x**3)) - 5*a*b**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*x**2*(a + b*x**3)) + b**5*x*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(a + b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_79():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**16
    F = -a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(15*x**15*(a + b*x**3)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(12*x**12*(a + b*x**3)) - 10*a**3*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(9*x**9*(a + b*x**3)) - 5*a**2*b**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(3*x**6*(a + b*x**3)) - 5*a*b**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(3*x**3*(a + b*x**3)) + b**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)*log(x)/(a + b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_80():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**17
    F = -a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(16*x**16*(a + b*x**3)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(13*x**13*(a + b*x**3)) - a**3*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x**10*(a + b*x**3)) - 10*a**2*b**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*x**7*(a + b*x**3)) - 5*a*b**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*x**4*(a + b*x**3)) - b**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_81():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**18
    F = -a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(17*x**17*(a + b*x**3)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(14*x**14*(a + b*x**3)) - 10*a**3*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(11*x**11*(a + b*x**3)) - 5*a**2*b**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*x**8*(a + b*x**3)) - a*b**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x**5*(a + b*x**3)) - b**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*x**2*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_82():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**19
    F = -(a + b*x**3)**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(18*a*x**18)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_83():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**20
    F = -a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(19*x**19*(a + b*x**3)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(16*x**16*(a + b*x**3)) - 10*a**3*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(13*x**13*(a + b*x**3)) - a**2*b**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(x**10*(a + b*x**3)) - 5*a*b**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*x**7*(a + b*x**3)) - b**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*x**4*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_84():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**21
    F = -a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(20*x**20*(a + b*x**3)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(17*x**17*(a + b*x**3)) - 5*a**3*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*x**14*(a + b*x**3)) - 10*a**2*b**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(11*x**11*(a + b*x**3)) - 5*a*b**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(8*x**8*(a + b*x**3)) - b**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(5*x**5*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_85():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**22
    F = -(a + b*x**3)**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(21*a*x**21) + b*(a + b*x**3)**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(126*a**2*x**18)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_86():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**23
    F = -a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(22*x**22*(a + b*x**3)) - 5*a**4*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(19*x**19*(a + b*x**3)) - 5*a**3*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(8*x**16*(a + b*x**3)) - 10*a**2*b**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(13*x**13*(a + b*x**3)) - a*b**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(2*x**10*(a + b*x**3)) - b**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*x**7*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_87():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**24
    F = -a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(23*x**23*(a + b*x**3)) - a**4*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(4*x**20*(a + b*x**3)) - 10*a**3*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(17*x**17*(a + b*x**3)) - 5*a**2*b**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(7*x**14*(a + b*x**3)) - 5*a*b**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(11*x**11*(a + b*x**3)) - b**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(8*x**8*(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_88():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)/x**25
    F = -(a + b*x**3)**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(24*a*x**24) + b*(a + b*x**3)**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(84*a**2*x**21) - b**2*(a + b*x**3)**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(504*a**3*x**18)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_89():
    f = x**4/sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)
    F = a**(sympy.S(2)/3)*(a + b*x**3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*b**(sympy.S(5)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - a**(sympy.S(2)/3)*(a + b*x**3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*b**(sympy.S(5)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + sqrt(3)*a**(sympy.S(2)/3)*(a + b*x**3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(5)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + x**2*(a + b*x**3)/(2*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_90():
    f = x**3/sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)
    F = -a**(sympy.S(1)/3)*(a + b*x**3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*b**(sympy.S(4)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + a**(sympy.S(1)/3)*(a + b*x**3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*b**(sympy.S(4)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + sqrt(3)*a**(sympy.S(1)/3)*(a + b*x**3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(4)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + x*(a + b*x**3)/(b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_91():
    f = x**2/sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)
    F = (a + b*x**3)*log(a + b*x**3)/(3*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_92():
    f = x/sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)
    F = -(a + b*x**3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + (a + b*x**3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - sqrt(3)*(a + b*x**3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_93():
    f = 1/sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)
    F = (a + b*x**3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - (a + b*x**3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - sqrt(3)*(a + b*x**3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_94():
    f = 1/(x*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    F = (a + b*x**3)*log(x)/(a*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - (a + b*x**3)*log(a + b*x**3)/(3*a*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_95():
    f = 1/(x**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    F = -(a + b*x**3)/(a*x*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + b**(sympy.S(1)/3)*(a + b*x**3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(4)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - b**(sympy.S(1)/3)*(a + b*x**3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(4)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + sqrt(3)*b**(sympy.S(1)/3)*(a + b*x**3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(4)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_96():
    f = 1/(x**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    F = (-a - b*x**3)/(2*a*x**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - b**(sympy.S(2)/3)*(a + b*x**3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(5)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + b**(sympy.S(2)/3)*(a + b*x**3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(5)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + sqrt(3)*b**(sympy.S(2)/3)*(a + b*x**3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(5)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_97():
    f = 1/(x**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    F = (-a - b*x**3)/(3*a*x**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - b*(a + b*x**3)*log(x)/(a**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + b*(a + b*x**3)*log(a + b*x**3)/(3*a**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_98():
    f = x**4/(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)
    F = -x**2/(6*b*(a + b*x**3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + x**2/(9*a*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - (a + b*x**3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(4)/3)*b**(sympy.S(5)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + (a + b*x**3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(4)/3)*b**(sympy.S(5)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - sqrt(3)*(a + b*x**3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(4)/3)*b**(sympy.S(5)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_99():
    f = x**3/(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)
    F = -x/(6*b*(a + b*x**3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + x/(18*a*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + (a + b*x**3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(5)/3)*b**(sympy.S(4)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - (a + b*x**3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(5)/3)*b**(sympy.S(4)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - sqrt(3)*(a + b*x**3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(5)/3)*b**(sympy.S(4)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_100():
    f = x**2/(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)
    F = -1/(6*b*(a + b*x**3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_101():
    f = x/(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)
    F = x**2/(6*a*(a + b*x**3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 2*x**2/(9*a**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + (a + b*x**3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(27*a**(sympy.S(7)/3)*b**(sympy.S(2)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - (2*a + 2*b*x**3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(7)/3)*b**(sympy.S(2)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - sqrt(3)*(2*a + 2*b*x**3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(7)/3)*b**(sympy.S(2)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_102():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(-3)/2)
    F = x*(a + b*x**3)/(6*a*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)) + 5*x*(a + b*x**3)**2/(18*a**2*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)) + 5*(a + b*x**3)**3*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(8)/3)*b**(sympy.S(1)/3)*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)) - 5*(a + b*x**3)**3*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(54*a**(sympy.S(8)/3)*b**(sympy.S(1)/3)*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)) - 5*sqrt(3)*(a + b*x**3)**3*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(8)/3)*b**(sympy.S(1)/3)*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_103():
    f = 1/(x*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2))
    F = 1/(6*a*(a + b*x**3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 1/(3*a**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + (a + b*x**3)*log(x)/(a**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - (a + b*x**3)*log(a + b*x**3)/(3*a**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_104():
    f = 1/(x**2*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2))
    F = 1/(6*a*x*(a + b*x**3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 7/(18*a**2*x*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - (14*a + 14*b*x**3)/(9*a**3*x*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 14*b**(sympy.S(1)/3)*(a + b*x**3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(10)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - 7*b**(sympy.S(1)/3)*(a + b*x**3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(27*a**(sympy.S(10)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 14*sqrt(3)*b**(sympy.S(1)/3)*(a + b*x**3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(10)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_105():
    f = 1/(x**3*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2))
    F = 1/(6*a*x**2*(a + b*x**3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 4/(9*a**2*x**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - (10*a + 10*b*x**3)/(9*a**3*x**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - 20*b**(sympy.S(2)/3)*(a + b*x**3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(27*a**(sympy.S(11)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 10*b**(sympy.S(2)/3)*(a + b*x**3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(27*a**(sympy.S(11)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 20*sqrt(3)*b**(sympy.S(2)/3)*(a + b*x**3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(27*a**(sympy.S(11)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_106():
    f = 1/(x**4*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2))
    F = -b/(6*a**2*(a + b*x**3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - 2*b/(3*a**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - (a + b*x**3)/(3*a**3*x**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - 3*b*(a + b*x**3)*log(x)/(a**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + b*(a + b*x**3)*log(a + b*x**3)/(a**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_107():
    f = x**6/(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = -x**4/(12*b*(a + b*x**3)**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - x/(27*b**2*(a + b*x**3)**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + x/(162*a*b**2*(a + b*x**3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 5*x/(486*a**2*b**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + (5*a + 5*b*x**3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(729*a**(sympy.S(8)/3)*b**(sympy.S(7)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - (5*a + 5*b*x**3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(1458*a**(sympy.S(8)/3)*b**(sympy.S(7)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - sqrt(3)*(5*a + 5*b*x**3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(729*a**(sympy.S(8)/3)*b**(sympy.S(7)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_108():
    f = x**5/(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = a/(12*b**2*(a + b*x**3)**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - 1/(9*b**2*(a + b*x**3)**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_109():
    f = x**4/(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = -x**2/(12*b*(a + b*x**3)**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + x**2/(54*a*b*(a + b*x**3)**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 7*x**2/(324*a**2*b*(a + b*x**3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 7*x**2/(243*a**3*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - (7*a + 7*b*x**3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(729*a**(sympy.S(10)/3)*b**(sympy.S(5)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + (7*a + 7*b*x**3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(1458*a**(sympy.S(10)/3)*b**(sympy.S(5)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - sqrt(3)*(7*a + 7*b*x**3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(729*a**(sympy.S(10)/3)*b**(sympy.S(5)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_110():
    f = x**3/(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = -x/(12*b*(a + b*x**3)**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + x/(108*a*b*(a + b*x**3)**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + x/(81*a**2*b*(a + b*x**3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 5*x/(243*a**3*b*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - (5*a + 5*b*x**3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(729*a**(sympy.S(11)/3)*b**(sympy.S(4)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + (10*a + 10*b*x**3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(729*a**(sympy.S(11)/3)*b**(sympy.S(4)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - sqrt(3)*(10*a + 10*b*x**3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(729*a**(sympy.S(11)/3)*b**(sympy.S(4)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_111():
    f = x**2/(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = -1/(12*b*(a + b*x**3)*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_112():
    f = x/(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = x**2/(12*a*(a + b*x**3)**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 5*x**2/(54*a**2*(a + b*x**3)**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 35*x**2/(324*a**3*(a + b*x**3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 35*x**2/(243*a**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - (35*a + 35*b*x**3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(729*a**(sympy.S(13)/3)*b**(sympy.S(2)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + (35*a + 35*b*x**3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(1458*a**(sympy.S(13)/3)*b**(sympy.S(2)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - sqrt(3)*(35*a + 35*b*x**3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(729*a**(sympy.S(13)/3)*b**(sympy.S(2)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_113():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(-5)/2)
    F = x*(a + b*x**3)/(12*a*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)) + 11*x*(a + b*x**3)**2/(108*a**2*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)) + 11*x*(a + b*x**3)**3/(81*a**3*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)) + 55*x*(a + b*x**3)**4/(243*a**4*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)) + 110*(a + b*x**3)**5*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(729*a**(sympy.S(14)/3)*b**(sympy.S(1)/3)*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)) - 55*(a + b*x**3)**5*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(729*a**(sympy.S(14)/3)*b**(sympy.S(1)/3)*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)) - 110*sqrt(3)*(a + b*x**3)**5*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(729*a**(sympy.S(14)/3)*b**(sympy.S(1)/3)*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_114():
    f = 1/(x*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2))
    F = 1/(12*a*(a + b*x**3)**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 1/(9*a**2*(a + b*x**3)**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 1/(6*a**3*(a + b*x**3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 1/(3*a**4*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + (a + b*x**3)*log(x)/(a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - (a + b*x**3)*log(a + b*x**3)/(3*a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_115():
    f = 1/(x**2*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2))
    F = 1/(12*a*x*(a + b*x**3)**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 13/(108*a**2*x*(a + b*x**3)**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 65/(324*a**3*x*(a + b*x**3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 455/(972*a**4*x*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - (455*a + 455*b*x**3)/(243*a**5*x*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 455*b**(sympy.S(1)/3)*(a + b*x**3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(729*a**(sympy.S(16)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - 455*b**(sympy.S(1)/3)*(a + b*x**3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(1458*a**(sympy.S(16)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 455*sqrt(3)*b**(sympy.S(1)/3)*(a + b*x**3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(729*a**(sympy.S(16)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_116():
    f = 1/(x**3*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2))
    F = 1/(12*a*x**2*(a + b*x**3)**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 7/(54*a**2*x**2*(a + b*x**3)**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 77/(324*a**3*x**2*(a + b*x**3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 154/(243*a**4*x**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - (385*a + 385*b*x**3)/(243*a**5*x**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - 770*b**(sympy.S(2)/3)*(a + b*x**3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(729*a**(sympy.S(17)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 385*b**(sympy.S(2)/3)*(a + b*x**3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(729*a**(sympy.S(17)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 770*sqrt(3)*b**(sympy.S(2)/3)*(a + b*x**3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(729*a**(sympy.S(17)/3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_117():
    f = 1/(x**4*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2))
    F = -b/(12*a**2*(a + b*x**3)**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - 2*b/(9*a**3*(a + b*x**3)**2*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - b/(2*a**4*(a + b*x**3)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - 4*b/(3*a**5*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - (a + b*x**3)/(3*a**5*x**3*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) - 5*b*(a + b*x**3)*log(x)/(a**6*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)) + 5*b*(a + b*x**3)*log(a + b*x**3)/(3*a**6*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_118():
    f = (d*x)**m*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = a**5*(d*x)**(m + 1)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(d*(a + b*x**3)*(m + 1)) + 5*a**4*b*(d*x)**(m + 4)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(d**4*(a + b*x**3)*(m + 4)) + 10*a**3*b**2*(d*x)**(m + 7)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(d**7*(a + b*x**3)*(m + 7)) + 10*a**2*b**3*(d*x)**(m + 10)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(d**10*(a + b*x**3)*(m + 10)) + 5*a*b**4*(d*x)**(m + 13)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(d**13*(a + b*x**3)*(m + 13)) + b**5*(d*x)**(m + 16)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(d**16*(a + b*x**3)*(m + 16))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_119():
    f = (d*x)**m*(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)
    F = a**3*(d*x)**(m + 1)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(d*(a + b*x**3)*(m + 1)) + 3*a**2*b*(d*x)**(m + 4)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(d**4*(a + b*x**3)*(m + 4)) + 3*a*b**2*(d*x)**(m + 7)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(d**7*(a + b*x**3)*(m + 7)) + b**3*(d*x)**(m + 10)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(d**10*(a + b*x**3)*(m + 10))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_120():
    f = (d*x)**m*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)
    F = a*(d*x)**(m + 1)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(d*(a + b*x**3)*(m + 1)) + b*(d*x)**(m + 4)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)/(d**4*(a + b*x**3)*(m + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_121():
    f = (d*x)**m/sqrt(a**2 + 2*a*b*x**3 + b**2*x**6)
    F = (d*x)**(m + 1)*(a + b*x**3)*hyper((1, m/3 + sympy.S(1)/3), (m/3 + sympy.S(4)/3,), -b*x**3/a)/(a*d*(m + 1)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_122():
    f = (d*x)**m/(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(3)/2)
    F = (d*x)**(m + 1)*(a + b*x**3)*hyper((3, m/3 + sympy.S(1)/3), (m/3 + sympy.S(4)/3,), -b*x**3/a)/(a**3*d*(m + 1)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_123():
    f = (d*x)**m/(a**2 + 2*a*b*x**3 + b**2*x**6)**(sympy.S(5)/2)
    F = (d*x)**(m + 1)*(a + b*x**3)*hyper((5, m/3 + sympy.S(1)/3), (m/3 + sympy.S(4)/3,), -b*x**3/a)/(a**5*d*(m + 1)*sqrt(a**2 + 2*a*b*x**3 + b**2*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_124():
    f = (d*x)**m*(a**2 + 2*a*b*x**3 + b**2*x**6)**p
    F = (d*x)**(m + 1)*(a**2 + 2*a*b*x**3 + b**2*x**6)**p*hyper((-2*p, m/3 + sympy.S(1)/3), (m/3 + sympy.S(4)/3,), -b*x**3/a)/(d*(1 + b*x**3/a)**(2*p)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_125():
    f = x**11*(a**2 + 2*a*b*x**3 + b**2*x**6)**p
    F = -a**3*(a + b*x**3)*(a**2 + 2*a*b*x**3 + b**2*x**6)**p/(3*b**4*(2*p + 1)) + a**2*(a + b*x**3)**2*(a**2 + 2*a*b*x**3 + b**2*x**6)**p/(2*b**4*(p + 1)) - a*(a + b*x**3)**3*(a**2 + 2*a*b*x**3 + b**2*x**6)**p/(b**4*(2*p + 3)) + (a + b*x**3)**4*(a**2 + 2*a*b*x**3 + b**2*x**6)**p/(6*b**4*(p + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_126():
    f = x**8*(a**2 + 2*a*b*x**3 + b**2*x**6)**p
    F = a**2*(a + b*x**3)*(a**2 + 2*a*b*x**3 + b**2*x**6)**p/(3*b**3*(2*p + 1)) - a*(a + b*x**3)**2*(a**2 + 2*a*b*x**3 + b**2*x**6)**p/(3*b**3*(p + 1)) + (a + b*x**3)**3*(a**2 + 2*a*b*x**3 + b**2*x**6)**p/(3*b**3*(2*p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_127():
    f = x**5*(a**2 + 2*a*b*x**3 + b**2*x**6)**p
    F = -a*(a + b*x**3)*(a**2 + 2*a*b*x**3 + b**2*x**6)**p/(3*b**2*(2*p + 1)) + (a + b*x**3)**2*(a**2 + 2*a*b*x**3 + b**2*x**6)**p/(6*b**2*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_128():
    f = x**4*(a**2 + 2*a*b*x**3 + b**2*x**6)**p
    F = x**5*(a**2 + 2*a*b*x**3 + b**2*x**6)**p*hyper((sympy.S(5)/3, -2*p), (sympy.S(8)/3,), -b*x**3/a)/(5*(1 + b*x**3/a)**(2*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_129():
    f = x**3*(a**2 + 2*a*b*x**3 + b**2*x**6)**p
    F = x**4*(a**2 + 2*a*b*x**3 + b**2*x**6)**p*hyper((sympy.S(4)/3, -2*p), (sympy.S(7)/3,), -b*x**3/a)/(4*(1 + b*x**3/a)**(2*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_130():
    f = x**2*(a**2 + 2*a*b*x**3 + b**2*x**6)**p
    F = (a + b*x**3)*(a**2 + 2*a*b*x**3 + b**2*x**6)**p/(3*b*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_131():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**p/x
    F = -(a + b*x**3)*(a**2 + 2*a*b*x**3 + b**2*x**6)**p*hyper((1, 2*p + 1), (2*p + 2,), 1 + b*x**3/a)/(3*a*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_132():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**p/x**2
    F = -(a**2 + 2*a*b*x**3 + b**2*x**6)**p*hyper((sympy.S(-1)/3, -2*p), (sympy.S(2)/3,), -b*x**3/a)/(x*(1 + b*x**3/a)**(2*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_133():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**p/x**3
    F = -(a**2 + 2*a*b*x**3 + b**2*x**6)**p*hyper((sympy.S(-2)/3, -2*p), (sympy.S(1)/3,), -b*x**3/a)/(2*x**2*(1 + b*x**3/a)**(2*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_134():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**p/x**4
    F = b*(a + b*x**3)*(a**2 + 2*a*b*x**3 + b**2*x**6)**p*hyper((2, 2*p + 1), (2*p + 2,), 1 + b*x**3/a)/(3*a**2*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_135():
    f = (a**2 + 2*a*b*x**3 + b**2*x**6)**p/x**5
    F = -(a**2 + 2*a*b*x**3 + b**2*x**6)**p*hyper((sympy.S(-4)/3, -2*p), (sympy.S(-1)/3,), -b*x**3/a)/(4*x**4*(1 + b*x**3/a)**(2*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_136():
    f = x**8/(a + b*x**3 + c*x**6)
    F = -b*log(a + b*x**3 + c*x**6)/(6*c**2) + x**3/(3*c) - (-2*a*c + b**2)*atanh((b + 2*c*x**3)/sqrt(-4*a*c + b**2))/(3*c**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_137():
    f = x**5/(a + b*x**3 + c*x**6)
    F = b*atanh((b + 2*c*x**3)/sqrt(-4*a*c + b**2))/(3*c*sqrt(-4*a*c + b**2)) + log(a + b*x**3 + c*x**6)/(6*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_138():
    f = x**2/(a + b*x**3 + c*x**6)
    F = -2*atanh((b + 2*c*x**3)/sqrt(-4*a*c + b**2))/(3*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_139():
    f = 1/(x*(a + b*x**3 + c*x**6))
    F = b*atanh((b + 2*c*x**3)/sqrt(-4*a*c + b**2))/(3*a*sqrt(-4*a*c + b**2)) + log(x)/a - log(a + b*x**3 + c*x**6)/(6*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_140():
    f = 1/(x**4*(a + b*x**3 + c*x**6))
    F = -1/(3*a*x**3) - b*log(x)/a**2 + b*log(a + b*x**3 + c*x**6)/(6*a**2) - (-2*a*c + b**2)*atanh((b + 2*c*x**3)/sqrt(-4*a*c + b**2))/(3*a**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_141():
    f = x**7/(a + b*x**3 + c*x**6)
    F = x**2/(2*c) + 2**(sympy.S(1)/3)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(5)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(5)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*sqrt(3)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(5)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(5)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(5)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*sqrt(3)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(5)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_142():
    f = x**6/(a + b*x**3 + c*x**6)
    F = x/c - 2**(sympy.S(2)/3)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(4)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(4)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*sqrt(3)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(4)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(4)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(4)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*sqrt(3)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(4)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_143():
    f = x**4/(a + b*x**3 + c*x**6)
    F = 2**(sympy.S(1)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(2)/3)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(1)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(2)/3)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(1)/3)*sqrt(3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(2)/3)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(1)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(2)/3)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(1)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(2)/3)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(1)/3)*sqrt(3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(2)/3)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_144():
    f = x**3/(a + b*x**3 + c*x**6)
    F = -2**(sympy.S(2)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(1)/3)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(2)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(1)/3)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(2)/3)*sqrt(3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(1)/3)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(2)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(1)/3)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(2)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(1)/3)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(2)/3)*sqrt(3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(1)/3)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_145():
    f = x/(a + b*x**3 + c*x**6)
    F = 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(3*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(6*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(1)/3)*sqrt(3)*c**(sympy.S(1)/3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(3*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(3*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(6*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(1)/3)*sqrt(3)*c**(sympy.S(1)/3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(3*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_146():
    f = 1/(a + b*x**3 + c*x**6)
    F = -2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(3*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(6*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(2)/3)*sqrt(3)*c**(sympy.S(2)/3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(3*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(3*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(6*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(2)/3)*sqrt(3)*c**(sympy.S(2)/3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(3*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_147():
    f = 1/(x**2*(a + b*x**3 + c*x**6))
    F = 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*(-b/sqrt(-4*a*c + b**2) + 1)*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*a*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*(-b/sqrt(-4*a*c + b**2) + 1)*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*a*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*sqrt(3)*c**(sympy.S(1)/3)*(-b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*a*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*(b/sqrt(-4*a*c + b**2) + 1)*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*a*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*(b/sqrt(-4*a*c + b**2) + 1)*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*a*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*sqrt(3)*c**(sympy.S(1)/3)*(b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*a*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) - 1/(a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_148():
    f = 1/(x**3*(a + b*x**3 + c*x**6))
    F = -2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*(-b/sqrt(-4*a*c + b**2) + 1)*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*a*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*(-b/sqrt(-4*a*c + b**2) + 1)*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*a*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*sqrt(3)*c**(sympy.S(2)/3)*(-b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*a*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*(b/sqrt(-4*a*c + b**2) + 1)*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*a*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*(b/sqrt(-4*a*c + b**2) + 1)*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*a*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*sqrt(3)*c**(sympy.S(2)/3)*(b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*a*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 1/(2*a*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_149():
    f = x**11/(x**6 + 4*x**3 + 3)
    F = x**6/6 - 4*x**3/3 - log(x**3 + 1)/6 + 9*log(x**3 + 3)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_150():
    f = x**8/(x**6 + 4*x**3 + 3)
    F = x**3/3 + log(x**3 + 1)/6 - 3*log(x**3 + 3)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_151():
    f = x**5/(x**6 + 4*x**3 + 3)
    F = -log(x**3 + 1)/6 + log(x**3 + 3)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_152():
    f = 1/(x*(x**6 + 4*x**3 + 3))
    F = log(x)/3 - log(x**3 + 1)/6 + log(x**3 + 3)/18
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_153():
    f = 1/(x**4*(x**6 + 4*x**3 + 3))
    F = -4*log(x)/9 + log(x**3 + 1)/6 - log(x**3 + 3)/54 - 1/(9*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_154():
    f = 1/(x**7*(x**6 + 4*x**3 + 3))
    F = 13*log(x)/27 - log(x**3 + 1)/6 + log(x**3 + 3)/162 + 4/(27*x**3) - 1/(18*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_155():
    f = x**10/(x**6 + 4*x**3 + 3)
    F = x**5/5 - 2*x**2 + log(x + 1)/6 - 3*3**(sympy.S(2)/3)*log(x + 3**(sympy.S(1)/3))/2 - log(x**2 - x + 1)/12 + 3*3**(sympy.S(2)/3)*log(x**2 - 3**(sympy.S(1)/3)*x + 3**(sympy.S(2)/3))/4 - 9*3**(sympy.S(1)/6)*atan(3**(sympy.S(1)/6)*(-2*x + 3**(sympy.S(1)/3))/3)/2 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_156():
    f = x**9/(x**6 + 4*x**3 + 3)
    F = x**4/4 - 4*x - log(x + 1)/6 + 3*3**(sympy.S(1)/3)*log(x + 3**(sympy.S(1)/3))/2 + log(x**2 - x + 1)/12 - 3*3**(sympy.S(1)/3)*log(x**2 - 3**(sympy.S(1)/3)*x + 3**(sympy.S(2)/3))/4 - 3*3**(sympy.S(5)/6)*atan(3**(sympy.S(1)/6)*(-2*x + 3**(sympy.S(1)/3))/3)/2 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_157():
    f = x**7/(x**6 + 4*x**3 + 3)
    F = x**2/2 - log(x + 1)/6 + 3**(sympy.S(2)/3)*log(x + 3**(sympy.S(1)/3))/2 + log(x**2 - x + 1)/12 - 3**(sympy.S(2)/3)*log(x**2 - 3**(sympy.S(1)/3)*x + 3**(sympy.S(2)/3))/4 + 3*3**(sympy.S(1)/6)*atan(3**(sympy.S(1)/6)*(-2*x + 3**(sympy.S(1)/3))/3)/2 - sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_158():
    f = x**6/(x**6 + 4*x**3 + 3)
    F = x + log(x + 1)/6 - 3**(sympy.S(1)/3)*log(x + 3**(sympy.S(1)/3))/2 - log(x**2 - x + 1)/12 + 3**(sympy.S(1)/3)*log(x**2 - 3**(sympy.S(1)/3)*x + 3**(sympy.S(2)/3))/4 + 3**(sympy.S(5)/6)*atan(3**(sympy.S(1)/6)*(-2*x + 3**(sympy.S(1)/3))/3)/2 - sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_159():
    f = x**4/(x**6 + 4*x**3 + 3)
    F = log(x + 1)/6 - 3**(sympy.S(2)/3)*log(x + 3**(sympy.S(1)/3))/6 - log(x**2 - x + 1)/12 + 3**(sympy.S(2)/3)*log(x**2 - 3**(sympy.S(1)/3)*x + 3**(sympy.S(2)/3))/12 - 3**(sympy.S(1)/6)*atan(3**(sympy.S(1)/6)*(-2*x + 3**(sympy.S(1)/3))/3)/2 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_160():
    f = x**3/(x**6 + 4*x**3 + 3)
    F = -log(x + 1)/6 + 3**(sympy.S(1)/3)*log(x + 3**(sympy.S(1)/3))/6 + log(x**2 - x + 1)/12 - 3**(sympy.S(1)/3)*log(x**2 - 3**(sympy.S(1)/3)*x + 3**(sympy.S(2)/3))/12 - 3**(sympy.S(5)/6)*atan(3**(sympy.S(1)/6)*(-2*x + 3**(sympy.S(1)/3))/3)/6 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_161():
    f = x/(x**6 + 4*x**3 + 3)
    F = -log(x + 1)/6 + 3**(sympy.S(2)/3)*log(x + 3**(sympy.S(1)/3))/18 + log(x**2 - x + 1)/12 - 3**(sympy.S(2)/3)*log(x**2 - 3**(sympy.S(1)/3)*x + 3**(sympy.S(2)/3))/36 + 3**(sympy.S(1)/6)*atan(3**(sympy.S(1)/6)*(-2*x + 3**(sympy.S(1)/3))/3)/6 - sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_162():
    f = 1/(x**6 + 4*x**3 + 3)
    F = log(x + 1)/6 - 3**(sympy.S(1)/3)*log(x + 3**(sympy.S(1)/3))/18 - log(x**2 - x + 1)/12 + 3**(sympy.S(1)/3)*log(x**2 - 3**(sympy.S(1)/3)*x + 3**(sympy.S(2)/3))/36 + 3**(sympy.S(5)/6)*atan(3**(sympy.S(1)/6)*(-2*x + 3**(sympy.S(1)/3))/3)/18 - sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_163():
    f = 1/(x**2*(x**6 + 4*x**3 + 3))
    F = log(x + 1)/6 - 3**(sympy.S(2)/3)*log(x + 3**(sympy.S(1)/3))/54 - log(x**2 - x + 1)/12 + 3**(sympy.S(2)/3)*log(x**2 - 3**(sympy.S(1)/3)*x + 3**(sympy.S(2)/3))/108 - 3**(sympy.S(1)/6)*atan(3**(sympy.S(1)/6)*(-2*x + 3**(sympy.S(1)/3))/3)/18 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/6 - 1/(3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_164():
    f = 1/(x**3*(x**6 + 4*x**3 + 3))
    F = -log(x + 1)/6 + 3**(sympy.S(1)/3)*log(x + 3**(sympy.S(1)/3))/54 + log(x**2 - x + 1)/12 - 3**(sympy.S(1)/3)*log(x**2 - 3**(sympy.S(1)/3)*x + 3**(sympy.S(2)/3))/108 - 3**(sympy.S(5)/6)*atan(3**(sympy.S(1)/6)*(-2*x + 3**(sympy.S(1)/3))/3)/54 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/6 - 1/(6*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_165():
    f = 1/(x**5*(x**6 + 4*x**3 + 3))
    F = -log(x + 1)/6 + 3**(sympy.S(2)/3)*log(x + 3**(sympy.S(1)/3))/162 + log(x**2 - x + 1)/12 - 3**(sympy.S(2)/3)*log(x**2 - 3**(sympy.S(1)/3)*x + 3**(sympy.S(2)/3))/324 + 3**(sympy.S(1)/6)*atan(3**(sympy.S(1)/6)*(-2*x + 3**(sympy.S(1)/3))/3)/54 - sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/6 + 4/(9*x) - 1/(12*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_166():
    f = 1/(x**6*(x**6 + 4*x**3 + 3))
    F = log(x + 1)/6 - 3**(sympy.S(1)/3)*log(x + 3**(sympy.S(1)/3))/162 - log(x**2 - x + 1)/12 + 3**(sympy.S(1)/3)*log(x**2 - 3**(sympy.S(1)/3)*x + 3**(sympy.S(2)/3))/324 + 3**(sympy.S(5)/6)*atan(3**(sympy.S(1)/6)*(-2*x + 3**(sympy.S(1)/3))/3)/162 - sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/6 + 2/(9*x**2) - 1/(15*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_167():
    f = x**6/(x**6 - x**3 + 1)
    F = x + 2**(sympy.S(2)/3)*(3 - sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 - sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 - sqrt(3)*I)**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(3 + sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 + sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 + sqrt(3)*I)**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(3 - sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 - 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 - sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 - sqrt(3)*I)**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(3 + sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 + 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 + sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 + sqrt(3)*I)**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(-sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half - sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 - sqrt(3)*I)**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half + sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 + sqrt(3)*I)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_168():
    f = x**5/(x**6 - x**3 + 1)
    F = log(x**6 - x**3 + 1)/6 - sqrt(3)*atan(sqrt(3)*(1 - 2*x**3)/3)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_169():
    f = x**4/(x**6 - x**3 + 1)
    F = 2**(sympy.S(1)/3)*(3 + sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 - sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 - sqrt(3)*I)**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*(3 - sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 + sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 + sqrt(3)*I)**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*(3 + sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 - 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 - sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 - sqrt(3)*I)**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*(3 - sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 + 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 + sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 + sqrt(3)*I)**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*(sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half - sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 - sqrt(3)*I)**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*(-sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half + sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 + sqrt(3)*I)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_170():
    f = x**3/(x**6 - x**3 + 1)
    F = 2**(sympy.S(2)/3)*(3 + sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 - sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 - sqrt(3)*I)**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(3 - sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 + sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 + sqrt(3)*I)**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(3 + sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 - 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 - sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 - sqrt(3)*I)**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(3 - sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 + 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 + sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 + sqrt(3)*I)**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half - sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 - sqrt(3)*I)**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(-sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half + sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 + sqrt(3)*I)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_171():
    f = x**2/(x**6 - x**3 + 1)
    F = -2*sqrt(3)*atan(sqrt(3)*(1 - 2*x**3)/3)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_172():
    f = x/(x**6 - x**3 + 1)
    F = sqrt(3)*I*log(-2**(sympy.S(1)/3)*x + (1 - sqrt(3)*I)**(sympy.S(1)/3))/(9*(sympy.S.Half - sqrt(3)*I/2)**(sympy.S(1)/3)) - sqrt(3)*I*log(-2**(sympy.S(1)/3)*x + (1 + sqrt(3)*I)**(sympy.S(1)/3))/(9*(sympy.S.Half + sqrt(3)*I/2)**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*sqrt(3)*I*log(2**(sympy.S(2)/3)*x**2 + x*(2 - 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 - sqrt(3)*I)**(sympy.S(2)/3))/(18*(1 - sqrt(3)*I)**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*sqrt(3)*I*log(2**(sympy.S(2)/3)*x**2 + x*(2 + 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 + sqrt(3)*I)**(sympy.S(2)/3))/(18*(1 + sqrt(3)*I)**(sympy.S(1)/3)) + I*atan(sqrt(3)*(2*x/(sympy.S.Half - sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(3*(sympy.S.Half - sqrt(3)*I/2)**(sympy.S(1)/3)) - I*atan(sqrt(3)*(2*x/(sympy.S.Half + sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(3*(sympy.S.Half + sqrt(3)*I/2)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_173():
    f = 1/(x*(x**6 - x**3 + 1))
    F = log(x) - log(x**6 - x**3 + 1)/6 - sqrt(3)*atan(sqrt(3)*(1 - 2*x**3)/3)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_174():
    f = 1/(x**2*(x**6 - x**3 + 1))
    F = -2**(sympy.S(1)/3)*(3 - sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 - sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 - sqrt(3)*I)**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*(3 + sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 + sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 + sqrt(3)*I)**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*(3 - sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 - 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 - sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 - sqrt(3)*I)**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*(3 + sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 + 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 + sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 + sqrt(3)*I)**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*(-sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half - sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 - sqrt(3)*I)**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*(sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half + sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 + sqrt(3)*I)**(sympy.S(1)/3)) - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_175():
    f = 1/(x**3*(x**6 - x**3 + 1))
    F = -2**(sympy.S(2)/3)*(3 - sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 - sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 - sqrt(3)*I)**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(3 + sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 + sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 + sqrt(3)*I)**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(3 - sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 - 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 - sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 - sqrt(3)*I)**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(3 + sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 + 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 + sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 + sqrt(3)*I)**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(-sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half - sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 - sqrt(3)*I)**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half + sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 + sqrt(3)*I)**(sympy.S(2)/3)) - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_176():
    f = 1/(x**4*(x**6 - x**3 + 1))
    F = log(x) - log(x**6 - x**3 + 1)/6 + sqrt(3)*atan(sqrt(3)*(1 - 2*x**3)/3)/9 - 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_177():
    f = 1/(x**5*(x**6 - x**3 + 1))
    F = -2**(sympy.S(1)/3)*(3 + sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 - sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 - sqrt(3)*I)**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*(3 - sqrt(3)*I)*log(-2**(sympy.S(1)/3)*x + (1 + sqrt(3)*I)**(sympy.S(1)/3))/(18*(1 + sqrt(3)*I)**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*(3 + sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 - 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 - sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 - sqrt(3)*I)**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*(3 - sqrt(3)*I)*log(2**(sympy.S(2)/3)*x**2 + x*(2 + 2*sqrt(3)*I)**(sympy.S(1)/3) + (1 + sqrt(3)*I)**(sympy.S(2)/3))/(36*(1 + sqrt(3)*I)**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*(sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half - sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 - sqrt(3)*I)**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*(-sqrt(3) + I)*atan(sqrt(3)*(2*x/(sympy.S.Half + sqrt(3)*I/2)**(sympy.S(1)/3) + 1)/3)/(6*(1 + sqrt(3)*I)**(sympy.S(1)/3)) - 1/x - 1/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_178():
    f = 1/(x**6 + x**3 + 2)
    F = -sqrt(7)*I*log(2**(sympy.S(1)/3)*x + (1 - sqrt(7)*I)**(sympy.S(1)/3))/(21*(sympy.S.Half - sqrt(7)*I/2)**(sympy.S(2)/3)) + sqrt(7)*I*log(2**(sympy.S(1)/3)*x + (1 + sqrt(7)*I)**(sympy.S(1)/3))/(21*(sympy.S.Half + sqrt(7)*I/2)**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*sqrt(7)*I*log(2**(sympy.S(2)/3)*x**2 - x*(2 - 2*sqrt(7)*I)**(sympy.S(1)/3) + (1 - sqrt(7)*I)**(sympy.S(2)/3))/(42*(1 - sqrt(7)*I)**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*sqrt(7)*I*log(2**(sympy.S(2)/3)*x**2 - x*(2 + 2*sqrt(7)*I)**(sympy.S(1)/3) + (1 + sqrt(7)*I)**(sympy.S(2)/3))/(42*(1 + sqrt(7)*I)**(sympy.S(2)/3)) + sqrt(21)*I*atan(sqrt(3)*(-2*x/(sympy.S.Half - sqrt(7)*I/2)**(sympy.S(1)/3) + 1)/3)/(21*(sympy.S.Half - sqrt(7)*I/2)**(sympy.S(2)/3)) - sqrt(21)*I*atan(sqrt(3)*(-2*x/(sympy.S.Half + sqrt(7)*I/2)**(sympy.S(1)/3) + 1)/3)/(21*(sympy.S.Half + sqrt(7)*I/2)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_179():
    f = x**2/(x**6 + x**3 + 2)
    F = 2*sqrt(7)*atan(sqrt(7)*(2*x**3 + 1)/7)/21
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_180():
    f = x**3/(x**6 + x**3 + 2)
    F = 2**(sympy.S(2)/3)*(7 + sqrt(7)*I)*log(2**(sympy.S(1)/3)*x + (1 - sqrt(7)*I)**(sympy.S(1)/3))/(42*(1 - sqrt(7)*I)**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(7 - sqrt(7)*I)*log(2**(sympy.S(1)/3)*x + (1 + sqrt(7)*I)**(sympy.S(1)/3))/(42*(1 + sqrt(7)*I)**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(7 + sqrt(7)*I)*log(2**(sympy.S(2)/3)*x**2 - x*(2 - 2*sqrt(7)*I)**(sympy.S(1)/3) + (1 - sqrt(7)*I)**(sympy.S(2)/3))/(84*(1 - sqrt(7)*I)**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(7 - sqrt(7)*I)*log(2**(sympy.S(2)/3)*x**2 - x*(2 + 2*sqrt(7)*I)**(sympy.S(1)/3) + (1 + sqrt(7)*I)**(sympy.S(2)/3))/(84*(1 + sqrt(7)*I)**(sympy.S(2)/3)) - sqrt(21)*I*(sympy.S.Half - sqrt(7)*I/2)**(sympy.S(1)/3)*atan(sqrt(3)*(-2*x/(sympy.S.Half - sqrt(7)*I/2)**(sympy.S(1)/3) + 1)/3)/21 + sqrt(21)*I*(sympy.S.Half + sqrt(7)*I/2)**(sympy.S(1)/3)*atan(sqrt(3)*(-2*x/(sympy.S.Half + sqrt(7)*I/2)**(sympy.S(1)/3) + 1)/3)/21
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_181():
    f = x**14*sqrt(a + b*x**3 + c*x**6)
    F = -b*x**6*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(20*c**2) + x**9*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(18*c) - (7*b*(-28*a*c + 15*b**2) - 6*c*x**3*(-20*a*c + 21*b**2))*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(2880*c**4) + (b + 2*c*x**3)*sqrt(a + b*x**3 + c*x**6)*(16*a**2*c**2 - 56*a*b**2*c + 21*b**4)/(1536*c**5) - (-4*a*c + b**2)*(16*a**2*c**2 - 56*a*b**2*c + 21*b**4)*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/(3072*c**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_182():
    f = x**11*sqrt(a + b*x**3 + c*x**6)
    F = -b*(b + 2*c*x**3)*(-12*a*c + 7*b**2)*sqrt(a + b*x**3 + c*x**6)/(384*c**4) + b*(-12*a*c + 7*b**2)*(-4*a*c + b**2)*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/(768*c**(sympy.S(9)/2)) + x**6*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(15*c) + (a + b*x**3 + c*x**6)**(sympy.S(3)/2)*(-32*a*c + 35*b**2 - 42*b*c*x**3)/(720*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_183():
    f = x**8*sqrt(a + b*x**3 + c*x**6)
    F = -5*b*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(72*c**2) + x**3*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(12*c) + (b + 2*c*x**3)*(-4*a*c + 5*b**2)*sqrt(a + b*x**3 + c*x**6)/(192*c**3) - (-4*a*c + b**2)*(-4*a*c + 5*b**2)*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/(384*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_184():
    f = x**5*sqrt(a + b*x**3 + c*x**6)
    F = -b*(b + 2*c*x**3)*sqrt(a + b*x**3 + c*x**6)/(24*c**2) + b*(-4*a*c + b**2)*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/(48*c**(sympy.S(5)/2)) + (a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(9*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_185():
    f = x**2*sqrt(a + b*x**3 + c*x**6)
    F = (b + 2*c*x**3)*sqrt(a + b*x**3 + c*x**6)/(12*c) - (-4*a*c + b**2)*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/(24*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_186():
    f = sqrt(a + b*x**3 + c*x**6)/x
    F = -sqrt(a)*atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/3 + b*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/(6*sqrt(c)) + sqrt(a + b*x**3 + c*x**6)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_187():
    f = sqrt(a + b*x**3 + c*x**6)/x**4
    F = sqrt(c)*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/3 - sqrt(a + b*x**3 + c*x**6)/(3*x**3) - b*atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/(6*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_188():
    f = sqrt(a + b*x**3 + c*x**6)/x**7
    F = -(2*a + b*x**3)*sqrt(a + b*x**3 + c*x**6)/(12*a*x**6) + (-4*a*c + b**2)*atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/(24*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_189():
    f = sqrt(a + b*x**3 + c*x**6)/x**10
    F = -(a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(9*a*x**9) + b*(2*a + b*x**3)*sqrt(a + b*x**3 + c*x**6)/(24*a**2*x**6) - b*(-4*a*c + b**2)*atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/(48*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_190():
    f = sqrt(a + b*x**3 + c*x**6)/x**13
    F = -(a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(12*a*x**12) + 5*b*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(72*a**2*x**9) - (2*a + b*x**3)*(-4*a*c + 5*b**2)*sqrt(a + b*x**3 + c*x**6)/(192*a**3*x**6) + (-4*a*c + b**2)*(-4*a*c + 5*b**2)*atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/(384*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_191():
    f = sqrt(a + b*x**3 + c*x**6)/x**16
    F = -(a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(15*a*x**15) + 7*b*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(120*a**2*x**12) - (-32*a*c + 35*b**2)*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(720*a**3*x**9) + b*(2*a + b*x**3)*(-12*a*c + 7*b**2)*sqrt(a + b*x**3 + c*x**6)/(384*a**4*x**6) - b*(-12*a*c + 7*b**2)*(-4*a*c + b**2)*atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/(768*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_192():
    f = x**3*sqrt(a + b*x**3 + c*x**6)
    F = x**4*sqrt(a + b*x**3 + c*x**6)*appellf1(sympy.S(4)/3, sympy.S(-1)/2, sympy.S(-1)/2, sympy.S(7)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(4*sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_193():
    f = x*sqrt(a + b*x**3 + c*x**6)
    F = x**2*sqrt(a + b*x**3 + c*x**6)*appellf1(sympy.S(2)/3, sympy.S(-1)/2, sympy.S(-1)/2, sympy.S(5)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(2*sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_194():
    f = sqrt(a + b*x**3 + c*x**6)
    F = x*sqrt(a + b*x**3 + c*x**6)*appellf1(sympy.S(1)/3, sympy.S(-1)/2, sympy.S(-1)/2, sympy.S(4)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_195():
    f = sqrt(a + b*x**3 + c*x**6)/x**2
    F = -sqrt(a + b*x**3 + c*x**6)*appellf1(sympy.S(-1)/3, sympy.S(-1)/2, sympy.S(-1)/2, sympy.S(2)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(x*sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_196():
    f = sqrt(a + b*x**3 + c*x**6)/x**3
    F = -sqrt(a + b*x**3 + c*x**6)*appellf1(sympy.S(-2)/3, sympy.S(-1)/2, sympy.S(-1)/2, sympy.S(1)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(2*x**2*sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_197():
    f = x**14*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)
    F = -11*b*x**6*(a + b*x**3 + c*x**6)**(sympy.S(5)/2)/(336*c**2) + x**9*(a + b*x**3 + c*x**6)**(sympy.S(5)/2)/(24*c) - (3*b*(-124*a*c + 77*b**2) - 10*c*x**3*(-28*a*c + 33*b**2))*(a + b*x**3 + c*x**6)**(sympy.S(5)/2)/(13440*c**4) + (b + 2*c*x**3)*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)*(16*a**2*c**2 - 72*a*b**2*c + 33*b**4)/(6144*c**5) - (b + 2*c*x**3)*(-4*a*c + b**2)*sqrt(a + b*x**3 + c*x**6)*(16*a**2*c**2 - 72*a*b**2*c + 33*b**4)/(16384*c**6) + (-4*a*c + b**2)**2*(16*a**2*c**2 - 72*a*b**2*c + 33*b**4)*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/(32768*c**(sympy.S(13)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_198():
    f = x**11*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)
    F = -b*(b + 2*c*x**3)*(-4*a*c + 3*b**2)*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(384*c**4) + b*(b + 2*c*x**3)*(-4*a*c + b**2)*(-4*a*c + 3*b**2)*sqrt(a + b*x**3 + c*x**6)/(1024*c**5) - b*(-4*a*c + b**2)**2*(-4*a*c + 3*b**2)*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/(2048*c**(sympy.S(11)/2)) + x**6*(a + b*x**3 + c*x**6)**(sympy.S(5)/2)/(21*c) + (a + b*x**3 + c*x**6)**(sympy.S(5)/2)*(-16*a*c + 21*b**2 - 30*b*c*x**3)/(840*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_199():
    f = x**8*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)
    F = -7*b*(a + b*x**3 + c*x**6)**(sympy.S(5)/2)/(180*c**2) + x**3*(a + b*x**3 + c*x**6)**(sympy.S(5)/2)/(18*c) + (b + 2*c*x**3)*(-4*a*c + 7*b**2)*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(576*c**3) - (b + 2*c*x**3)*(-4*a*c + b**2)*(-4*a*c + 7*b**2)*sqrt(a + b*x**3 + c*x**6)/(1536*c**4) + (-4*a*c + b**2)**2*(-4*a*c + 7*b**2)*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/(3072*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_200():
    f = x**5*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)
    F = -b*(b + 2*c*x**3)*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(48*c**2) + b*(b + 2*c*x**3)*(-4*a*c + b**2)*sqrt(a + b*x**3 + c*x**6)/(128*c**3) - b*(-4*a*c + b**2)**2*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/(256*c**(sympy.S(7)/2)) + (a + b*x**3 + c*x**6)**(sympy.S(5)/2)/(15*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_201():
    f = x**2*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)
    F = (b + 2*c*x**3)*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(24*c) - (b + 2*c*x**3)*(-4*a*c + b**2)*sqrt(a + b*x**3 + c*x**6)/(64*c**2) + (-4*a*c + b**2)**2*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/(128*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_202():
    f = (a + b*x**3 + c*x**6)**(sympy.S(3)/2)/x
    F = -a**(sympy.S(3)/2)*atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/3 - b*(-12*a*c + b**2)*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/(48*c**(sympy.S(3)/2)) + (a + b*x**3 + c*x**6)**(sympy.S(3)/2)/9 + sqrt(a + b*x**3 + c*x**6)*(8*a*c + b**2 + 2*b*c*x**3)/(24*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_203():
    f = (a + b*x**3 + c*x**6)**(sympy.S(3)/2)/x**4
    F = -sqrt(a)*b*atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/2 + (3*b/4 + c*x**3/2)*sqrt(a + b*x**3 + c*x**6) - (a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(3*x**3) + (4*a*c + b**2)*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/(8*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_204():
    f = (a + b*x**3 + c*x**6)**(sympy.S(3)/2)/x**7
    F = b*sqrt(c)*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/2 - (b - 2*c*x**3)*sqrt(a + b*x**3 + c*x**6)/(4*x**3) - (a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(6*x**6) - (4*a*c + b**2)*atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/(8*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_205():
    f = (a + b*x**3 + c*x**6)**(sympy.S(3)/2)/x**10
    F = c**(sympy.S(3)/2)*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/3 - (a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(9*x**9) - (2*a*b + x**3*(8*a*c + b**2))*sqrt(a + b*x**3 + c*x**6)/(24*a*x**6) + b*(-12*a*c + b**2)*atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/(48*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_206():
    f = (a + b*x**3 + c*x**6)**(sympy.S(3)/2)/x**13
    F = -(2*a + b*x**3)*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(24*a*x**12) + (2*a + b*x**3)*(-4*a*c + b**2)*sqrt(a + b*x**3 + c*x**6)/(64*a**2*x**6) - (-4*a*c + b**2)**2*atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/(128*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_207():
    f = (a + b*x**3 + c*x**6)**(sympy.S(3)/2)/x**16
    F = -(a + b*x**3 + c*x**6)**(sympy.S(5)/2)/(15*a*x**15) + b*(2*a + b*x**3)*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(48*a**2*x**12) - b*(2*a + b*x**3)*(-4*a*c + b**2)*sqrt(a + b*x**3 + c*x**6)/(128*a**3*x**6) + b*(-4*a*c + b**2)**2*atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/(256*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_208():
    f = (a + b*x**3 + c*x**6)**(sympy.S(3)/2)/x**19
    F = -(a + b*x**3 + c*x**6)**(sympy.S(5)/2)/(18*a*x**18) + 7*b*(a + b*x**3 + c*x**6)**(sympy.S(5)/2)/(180*a**2*x**15) - (2*a + b*x**3)*(-4*a*c + 7*b**2)*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(576*a**3*x**12) + (2*a + b*x**3)*(-4*a*c + b**2)*(-4*a*c + 7*b**2)*sqrt(a + b*x**3 + c*x**6)/(1536*a**4*x**6) - (-4*a*c + b**2)**2*(-4*a*c + 7*b**2)*atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/(3072*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_209():
    f = (a + b*x**3 + c*x**6)**(sympy.S(3)/2)/x**22
    F = -(a + b*x**3 + c*x**6)**(sympy.S(5)/2)/(21*a*x**21) + b*(a + b*x**3 + c*x**6)**(sympy.S(5)/2)/(28*a**2*x**18) - (-16*a*c + 21*b**2)*(a + b*x**3 + c*x**6)**(sympy.S(5)/2)/(840*a**3*x**15) + b*(2*a + b*x**3)*(-4*a*c + 3*b**2)*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)/(384*a**4*x**12) - b*(2*a + b*x**3)*(-4*a*c + b**2)*(-4*a*c + 3*b**2)*sqrt(a + b*x**3 + c*x**6)/(1024*a**5*x**6) + b*(-4*a*c + b**2)**2*(-4*a*c + 3*b**2)*atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/(2048*a**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_210():
    f = x**3*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)
    F = a*x**4*sqrt(a + b*x**3 + c*x**6)*appellf1(sympy.S(4)/3, sympy.S(-3)/2, sympy.S(-3)/2, sympy.S(7)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(4*sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_211():
    f = x*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)
    F = a*x**2*sqrt(a + b*x**3 + c*x**6)*appellf1(sympy.S(2)/3, sympy.S(-3)/2, sympy.S(-3)/2, sympy.S(5)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(2*sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_212():
    f = (a + b*x**3 + c*x**6)**(sympy.S(3)/2)
    F = a*x*sqrt(a + b*x**3 + c*x**6)*appellf1(sympy.S(1)/3, sympy.S(-3)/2, sympy.S(-3)/2, sympy.S(4)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_213():
    f = (a + b*x**3 + c*x**6)**(sympy.S(3)/2)/x**2
    F = -a*sqrt(a + b*x**3 + c*x**6)*appellf1(sympy.S(-1)/3, sympy.S(-3)/2, sympy.S(-3)/2, sympy.S(2)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(x*sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_214():
    f = (a + b*x**3 + c*x**6)**(sympy.S(3)/2)/x**3
    F = -a*sqrt(a + b*x**3 + c*x**6)*appellf1(sympy.S(-2)/3, sympy.S(-3)/2, sympy.S(-3)/2, sympy.S(1)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(2*x**2*sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_215():
    f = x**14/sqrt(a + b*x**3 + c*x**6)
    F = -7*b*x**6*sqrt(a + b*x**3 + c*x**6)/(72*c**2) + x**9*sqrt(a + b*x**3 + c*x**6)/(12*c) - (5*b*(-44*a*c + 21*b**2) - 2*c*x**3*(-36*a*c + 35*b**2))*sqrt(a + b*x**3 + c*x**6)/(576*c**4) + (48*a**2*c**2 - 120*a*b**2*c + 35*b**4)*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/(384*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_216():
    f = x**11/sqrt(a + b*x**3 + c*x**6)
    F = -b*(-12*a*c + 5*b**2)*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/(48*c**(sympy.S(7)/2)) + x**6*sqrt(a + b*x**3 + c*x**6)/(9*c) + sqrt(a + b*x**3 + c*x**6)*(-16*a*c + 15*b**2 - 10*b*c*x**3)/(72*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_217():
    f = x**8/sqrt(a + b*x**3 + c*x**6)
    F = -b*sqrt(a + b*x**3 + c*x**6)/(4*c**2) + x**3*sqrt(a + b*x**3 + c*x**6)/(6*c) + (-4*a*c + 3*b**2)*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/(24*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_218():
    f = x**5/sqrt(a + b*x**3 + c*x**6)
    F = -b*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/(6*c**(sympy.S(3)/2)) + sqrt(a + b*x**3 + c*x**6)/(3*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_219():
    f = x**2/sqrt(a + b*x**3 + c*x**6)
    F = atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/(3*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_220():
    f = 1/(x*sqrt(a + b*x**3 + c*x**6))
    F = -atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/(3*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_221():
    f = 1/(x**4*sqrt(a + b*x**3 + c*x**6))
    F = -sqrt(a + b*x**3 + c*x**6)/(3*a*x**3) + b*atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/(6*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_222():
    f = 1/(x**7*sqrt(a + b*x**3 + c*x**6))
    F = -sqrt(a + b*x**3 + c*x**6)/(6*a*x**6) + b*sqrt(a + b*x**3 + c*x**6)/(4*a**2*x**3) - (-4*a*c + 3*b**2)*atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/(24*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_223():
    f = 1/(x**10*sqrt(a + b*x**3 + c*x**6))
    F = -sqrt(a + b*x**3 + c*x**6)/(9*a*x**9) + 5*b*sqrt(a + b*x**3 + c*x**6)/(36*a**2*x**6) - (-16*a*c + 15*b**2)*sqrt(a + b*x**3 + c*x**6)/(72*a**3*x**3) + b*(-12*a*c + 5*b**2)*atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/(48*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_224():
    f = 1/(x**13*sqrt(a + b*x**3 + c*x**6))
    F = -sqrt(a + b*x**3 + c*x**6)/(12*a*x**12) + 7*b*sqrt(a + b*x**3 + c*x**6)/(72*a**2*x**9) - (-36*a*c + 35*b**2)*sqrt(a + b*x**3 + c*x**6)/(288*a**3*x**6) + 5*b*(-44*a*c + 21*b**2)*sqrt(a + b*x**3 + c*x**6)/(576*a**4*x**3) - (48*a**2*c**2 - 120*a*b**2*c + 35*b**4)*atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/(384*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_225():
    f = x**3/sqrt(a + b*x**3 + c*x**6)
    F = x**4*sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(4)/3, sympy.S.Half, sympy.S.Half, sympy.S(7)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(4*sqrt(a + b*x**3 + c*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_226():
    f = x/sqrt(a + b*x**3 + c*x**6)
    F = x**2*sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(2)/3, sympy.S.Half, sympy.S.Half, sympy.S(5)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(2*sqrt(a + b*x**3 + c*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_227():
    f = 1/sqrt(a + b*x**3 + c*x**6)
    F = x*sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(1)/3, sympy.S.Half, sympy.S.Half, sympy.S(4)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/sqrt(a + b*x**3 + c*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_228():
    f = 1/(x**2*sqrt(a + b*x**3 + c*x**6))
    F = -sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(-1)/3, sympy.S.Half, sympy.S.Half, sympy.S(2)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(x*sqrt(a + b*x**3 + c*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_229():
    f = 1/(x**3*sqrt(a + b*x**3 + c*x**6))
    F = -sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(-2)/3, sympy.S.Half, sympy.S.Half, sympy.S(1)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(2*x**2*sqrt(a + b*x**3 + c*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_230():
    f = x**14/(a + b*x**3 + c*x**6)**(sympy.S(3)/2)
    F = -2*b*x**6*sqrt(a + b*x**3 + c*x**6)/(3*c*(-4*a*c + b**2)) + 2*x**9*(2*a + b*x**3)/((-12*a*c + 3*b**2)*sqrt(a + b*x**3 + c*x**6)) - (b*(-52*a*c + 15*b**2) - 2*c*x**3*(-12*a*c + 5*b**2))*sqrt(a + b*x**3 + c*x**6)/(12*c**3*(-4*a*c + b**2)) + (-4*a*c + 5*b**2)*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/(8*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_231():
    f = x**11/(a + b*x**3 + c*x**6)**(sympy.S(3)/2)
    F = -b*atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/(2*c**(sympy.S(5)/2)) + 2*x**6*(2*a + b*x**3)/((-12*a*c + 3*b**2)*sqrt(a + b*x**3 + c*x**6)) + sqrt(a + b*x**3 + c*x**6)*(-8*a*c + 3*b**2 - 2*b*c*x**3)/(3*c**2*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_232():
    f = x**8/(a + b*x**3 + c*x**6)**(sympy.S(3)/2)
    F = -2*b*sqrt(a + b*x**3 + c*x**6)/(3*c*(-4*a*c + b**2)) + 2*x**3*(2*a + b*x**3)/((-12*a*c + 3*b**2)*sqrt(a + b*x**3 + c*x**6)) + atanh((b + 2*c*x**3)/(2*sqrt(c)*sqrt(a + b*x**3 + c*x**6)))/(3*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_233():
    f = x**5/(a + b*x**3 + c*x**6)**(sympy.S(3)/2)
    F = (4*a + 2*b*x**3)/((-12*a*c + 3*b**2)*sqrt(a + b*x**3 + c*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_234():
    f = x**2/(a + b*x**3 + c*x**6)**(sympy.S(3)/2)
    F = -(2*b + 4*c*x**3)/((-12*a*c + 3*b**2)*sqrt(a + b*x**3 + c*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_235():
    f = 1/(x*(a + b*x**3 + c*x**6)**(sympy.S(3)/2))
    F = (-4*a*c + 2*b**2 + 2*b*c*x**3)/(3*a*(-4*a*c + b**2)*sqrt(a + b*x**3 + c*x**6)) - atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/(3*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_236():
    f = 1/(x**4*(a + b*x**3 + c*x**6)**(sympy.S(3)/2))
    F = (-4*a*c + 2*b**2 + 2*b*c*x**3)/(3*a*x**3*(-4*a*c + b**2)*sqrt(a + b*x**3 + c*x**6)) - (-8*a*c + 3*b**2)*sqrt(a + b*x**3 + c*x**6)/(3*a**2*x**3*(-4*a*c + b**2)) + b*atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/(2*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_237():
    f = 1/(x**7*(a + b*x**3 + c*x**6)**(sympy.S(3)/2))
    F = (-4*a*c + 2*b**2 + 2*b*c*x**3)/(3*a*x**6*(-4*a*c + b**2)*sqrt(a + b*x**3 + c*x**6)) - (-12*a*c + 5*b**2)*sqrt(a + b*x**3 + c*x**6)/(6*a**2*x**6*(-4*a*c + b**2)) + b*(-52*a*c + 15*b**2)*sqrt(a + b*x**3 + c*x**6)/(12*a**3*x**3*(-4*a*c + b**2)) - (-4*a*c + 5*b**2)*atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/(8*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_238():
    f = 1/(x**10*(a + b*x**3 + c*x**6)**(sympy.S(3)/2))
    F = (-4*a*c + 2*b**2 + 2*b*c*x**3)/(3*a*x**9*(-4*a*c + b**2)*sqrt(a + b*x**3 + c*x**6)) - (-16*a*c + 7*b**2)*sqrt(a + b*x**3 + c*x**6)/(9*a**2*x**9*(-4*a*c + b**2)) + b*(-116*a*c + 35*b**2)*sqrt(a + b*x**3 + c*x**6)/(36*a**3*x**6*(-4*a*c + b**2)) - sqrt(a + b*x**3 + c*x**6)*(256*a**2*c**2 - 460*a*b**2*c + 105*b**4)/(72*a**4*x**3*(-4*a*c + b**2)) + 5*b*(-12*a*c + 7*b**2)*atanh((2*a + b*x**3)/(2*sqrt(a)*sqrt(a + b*x**3 + c*x**6)))/(48*a**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_239():
    f = x**3/(a + b*x**3 + c*x**6)**(sympy.S(3)/2)
    F = x**4*sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(4)/3, sympy.S(3)/2, sympy.S(3)/2, sympy.S(7)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(a + b*x**3 + c*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_240():
    f = x/(a + b*x**3 + c*x**6)**(sympy.S(3)/2)
    F = x**2*sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(2)/3, sympy.S(3)/2, sympy.S(3)/2, sympy.S(5)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(2*a*sqrt(a + b*x**3 + c*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_241():
    f = (a + b*x**3 + c*x**6)**(sympy.S(-3)/2)
    F = x*sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(1)/3, sympy.S(3)/2, sympy.S(3)/2, sympy.S(4)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(a*sqrt(a + b*x**3 + c*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_242():
    f = 1/(x**2*(a + b*x**3 + c*x**6)**(sympy.S(3)/2))
    F = -sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(-1)/3, sympy.S(3)/2, sympy.S(3)/2, sympy.S(2)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(a*x*sqrt(a + b*x**3 + c*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_243():
    f = 1/(x**3*(a + b*x**3 + c*x**6)**(sympy.S(3)/2))
    F = -sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(-2)/3, sympy.S(3)/2, sympy.S(3)/2, sympy.S(1)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(2*a*x**2*sqrt(a + b*x**3 + c*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_244():
    f = (d*x)**m*(a + b*x**3 + c*x**6)**2
    F = a**2*(d*x)**(m + 1)/(d*(m + 1)) + 2*a*b*(d*x)**(m + 4)/(d**4*(m + 4)) + 2*b*c*(d*x)**(m + 10)/(d**10*(m + 10)) + c**2*(d*x)**(m + 13)/(d**13*(m + 13)) + (d*x)**(m + 7)*(2*a*c + b**2)/(d**7*(m + 7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_245():
    f = (d*x)**m*(a + b*x**3 + c*x**6)
    F = a*(d*x)**(m + 1)/(d*(m + 1)) + b*(d*x)**(m + 4)/(d**4*(m + 4)) + c*(d*x)**(m + 7)/(d**7*(m + 7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_246():
    f = (d*x)**m/(a + b*x**3 + c*x**6)
    F = -2*c*(d*x)**(m + 1)*hyper((1, m/3 + sympy.S(1)/3), (m/3 + sympy.S(4)/3,), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(d*(b + sqrt(-4*a*c + b**2))*(m + 1)*sqrt(-4*a*c + b**2)) + 2*c*(d*x)**(m + 1)*hyper((1, m/3 + sympy.S(1)/3), (m/3 + sympy.S(4)/3,), -2*c*x**3/(b - sqrt(-4*a*c + b**2)))/(d*(b - sqrt(-4*a*c + b**2))*(m + 1)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_247():
    f = (d*x)**m/(a + b*x**3 + c*x**6)**2
    F = -c*(d*x)**(m + 1)*(-4*a*c*(5 - m) + b**2*(2 - m) - b*(2 - m)*sqrt(-4*a*c + b**2))*hyper((1, m/3 + sympy.S(1)/3), (m/3 + sympy.S(4)/3,), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(3*a*d*(b + sqrt(-4*a*c + b**2))*(m + 1)*(-4*a*c + b**2)**(sympy.S(3)/2)) + c*(d*x)**(m + 1)*(-4*a*c*(5 - m) + b**2*(2 - m) + b*(2 - m)*sqrt(-4*a*c + b**2))*hyper((1, m/3 + sympy.S(1)/3), (m/3 + sympy.S(4)/3,), -2*c*x**3/(b - sqrt(-4*a*c + b**2)))/(3*a*d*(b - sqrt(-4*a*c + b**2))*(m + 1)*(-4*a*c + b**2)**(sympy.S(3)/2)) + (d*x)**(m + 1)*(-2*a*c + b**2 + b*c*x**3)/(3*a*d*(-4*a*c + b**2)*(a + b*x**3 + c*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_248():
    f = (d*x)**m*(a + b*x**3 + c*x**6)**(sympy.S(3)/2)
    F = a*(d*x)**(m + 1)*sqrt(a + b*x**3 + c*x**6)*appellf1(m/3 + sympy.S(1)/3, sympy.S(-3)/2, sympy.S(-3)/2, m/3 + sympy.S(4)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(d*(m + 1)*sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_249():
    f = (d*x)**m*sqrt(a + b*x**3 + c*x**6)
    F = (d*x)**(m + 1)*sqrt(a + b*x**3 + c*x**6)*appellf1(m/3 + sympy.S(1)/3, sympy.S(-1)/2, sympy.S(-1)/2, m/3 + sympy.S(4)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(d*(m + 1)*sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_250():
    f = (d*x)**m/sqrt(a + b*x**3 + c*x**6)
    F = (d*x)**(m + 1)*sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(m/3 + sympy.S(1)/3, sympy.S.Half, sympy.S.Half, m/3 + sympy.S(4)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(d*(m + 1)*sqrt(a + b*x**3 + c*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_251():
    f = (d*x)**m/(a + b*x**3 + c*x**6)**(sympy.S(3)/2)
    F = (d*x)**(m + 1)*sqrt(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(m/3 + sympy.S(1)/3, sympy.S(3)/2, sympy.S(3)/2, m/3 + sympy.S(4)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(a*d*(m + 1)*sqrt(a + b*x**3 + c*x**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_252():
    f = (d*x)**m*(a + b*x**3 + c*x**6)**p
    F = (d*x)**(m + 1)*(a + b*x**3 + c*x**6)**p*appellf1(m/3 + sympy.S(1)/3, -p, -p, m/3 + sympy.S(4)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(d*(m + 1)*(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_253():
    f = x**8*(a + b*x**3 + c*x**6)**p
    F = 2**p*(-(b + 2*c*x**3 - sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))**(-p - 1)*(2*a*c - b**2*(p + 2))*(a + b*x**3 + c*x**6)**(p + 1)*hyper((-p, p + 1), (p + 2,), (b + 2*c*x**3 + sqrt(-4*a*c + b**2))/(2*sqrt(-4*a*c + b**2)))/(3*c**2*(p + 1)*(2*p + 3)*sqrt(-4*a*c + b**2)) - b*(p + 2)*(a + b*x**3 + c*x**6)**(p + 1)/(6*c**2*(p + 1)*(2*p + 3)) + x**3*(a + b*x**3 + c*x**6)**(p + 1)/(3*c*(2*p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_254():
    f = x**5*(a + b*x**3 + c*x**6)**p
    F = 2**p*b*(-(b + 2*c*x**3 - sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))**(-p - 1)*(a + b*x**3 + c*x**6)**(p + 1)*hyper((-p, p + 1), (p + 2,), (b + 2*c*x**3 + sqrt(-4*a*c + b**2))/(2*sqrt(-4*a*c + b**2)))/(3*c*(p + 1)*sqrt(-4*a*c + b**2)) + (a + b*x**3 + c*x**6)**(p + 1)/(6*c*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_255():
    f = x**2*(a + b*x**3 + c*x**6)**p
    F = -2**(p + 1)*(-(b + 2*c*x**3 - sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2))**(-p - 1)*(a + b*x**3 + c*x**6)**(p + 1)*hyper((-p, p + 1), (p + 2,), (b + 2*c*x**3 + sqrt(-4*a*c + b**2))/(2*sqrt(-4*a*c + b**2)))/(3*(p + 1)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_256():
    f = x**4*(a + b*x**3 + c*x**6)**p
    F = x**5*(a + b*x**3 + c*x**6)**p*appellf1(sympy.S(5)/3, -p, -p, sympy.S(8)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(5*(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_257():
    f = x**3*(a + b*x**3 + c*x**6)**p
    F = x**4*(a + b*x**3 + c*x**6)**p*appellf1(sympy.S(4)/3, -p, -p, sympy.S(7)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(4*(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_258():
    f = x*(a + b*x**3 + c*x**6)**p
    F = x**2*(a + b*x**3 + c*x**6)**p*appellf1(sympy.S(2)/3, -p, -p, sympy.S(5)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(2*(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_259():
    f = (a + b*x**3 + c*x**6)**p
    F = x*(a + b*x**3 + c*x**6)**p*appellf1(sympy.S(1)/3, -p, -p, sympy.S(4)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/((2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_260():
    f = (a + b*x**3 + c*x**6)**p/x
    F = 2**(2*p - 1)*(a + b*x**3 + c*x**6)**p*appellf1(-2*p, -p, -p, 1 - 2*p, -(b - sqrt(-4*a*c + b**2))/(2*c*x**3), -(b + sqrt(-4*a*c + b**2))/(2*c*x**3))/(3*p*((b + 2*c*x**3 - sqrt(-4*a*c + b**2))/(c*x**3))**p*((b + 2*c*x**3 + sqrt(-4*a*c + b**2))/(c*x**3))**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_261():
    f = (a + b*x**3 + c*x**6)**p/x**2
    F = -(a + b*x**3 + c*x**6)**p*appellf1(sympy.S(-1)/3, -p, -p, sympy.S(2)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(x*(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_262():
    f = (a + b*x**3 + c*x**6)**p/x**3
    F = -(a + b*x**3 + c*x**6)**p*appellf1(sympy.S(-2)/3, -p, -p, sympy.S(1)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(2*x**2*(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_263():
    f = (a + b*x**3 + c*x**6)**p/x**4
    F = -4**p*(a + b*x**3 + c*x**6)**p*appellf1(1 - 2*p, -p, -p, 2 - 2*p, -(b - sqrt(-4*a*c + b**2))/(2*c*x**3), -(b + sqrt(-4*a*c + b**2))/(2*c*x**3))/(x**3*((b + 2*c*x**3 - sqrt(-4*a*c + b**2))/(c*x**3))**p*((b + 2*c*x**3 + sqrt(-4*a*c + b**2))/(c*x**3))**p*(3 - 6*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_264():
    f = (a + b*x**3 + c*x**6)**p/x**5
    F = -(a + b*x**3 + c*x**6)**p*appellf1(sympy.S(-4)/3, -p, -p, sympy.S(-1)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(4*x**4*(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_265():
    f = (a + b*x**3 + c*x**6)**p/x**6
    F = -(a + b*x**3 + c*x**6)**p*appellf1(sympy.S(-5)/3, -p, -p, sympy.S(-2)/3, -2*c*x**3/(b - sqrt(-4*a*c + b**2)), -2*c*x**3/(b + sqrt(-4*a*c + b**2)))/(5*x**5*(2*c*x**3/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**3/(b + sqrt(-4*a*c + b**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_266():
    f = (a + b*x**3 + c*x**6)**p/x**7
    F = -2**(2*p - 1)*(a + b*x**3 + c*x**6)**p*appellf1(2 - 2*p, -p, -p, 3 - 2*p, -(b - sqrt(-4*a*c + b**2))/(2*c*x**3), -(b + sqrt(-4*a*c + b**2))/(2*c*x**3))/(x**6*((b + 2*c*x**3 - sqrt(-4*a*c + b**2))/(c*x**3))**p*((b + 2*c*x**3 + sqrt(-4*a*c + b**2))/(c*x**3))**p*(3 - 3*p))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_267():
    f = x**m/(x**8 + 2*x**4 + 1)
    F = x**(m + 1)*hyper((2, m/4 + sympy.S(1)/4), (m/4 + sympy.S(5)/4,), -x**4)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_268():
    f = x**9/(x**8 + 2*x**4 + 1)
    F = -x**6/(4*x**4 + 4) + 3*x**2/4 - 3*atan(x**2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_269():
    f = x**7/(x**8 + 2*x**4 + 1)
    F = log(x**4 + 1)/4 + 1/(4*x**4 + 4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_270():
    f = x**5/(x**8 + 2*x**4 + 1)
    F = -x**2/(4*x**4 + 4) + atan(x**2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_271():
    f = x**3/(x**8 + 2*x**4 + 1)
    F = -1/(4*x**4 + 4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_272():
    f = x/(x**8 + 2*x**4 + 1)
    F = x**2/(4*x**4 + 4) + atan(x**2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_273():
    f = 1/(x*(x**8 + 2*x**4 + 1))
    F = log(x) - log(x**4 + 1)/4 + 1/(4*x**4 + 4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_274():
    f = 1/(x**3*(x**8 + 2*x**4 + 1))
    F = -3*atan(x**2)/4 - 3/(4*x**2) + 1/(4*x**2*(x**4 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_275():
    f = 1/(x**5*(x**8 + 2*x**4 + 1))
    F = -2*log(x) + log(x**4 + 1)/2 - 1/(4*x**4 + 4) - 1/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_276():
    f = 1/(x**7*(x**8 + 2*x**4 + 1))
    F = 5*atan(x**2)/4 + 5/(4*x**2) - 5/(12*x**6) + 1/(4*x**6*(x**4 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_277():
    f = x**8/(x**8 + 2*x**4 + 1)
    F = -x**5/(4*x**4 + 4) + 5*x/4 + 5*sqrt(2)*log(x**2 - sqrt(2)*x + 1)/32 - 5*sqrt(2)*log(x**2 + sqrt(2)*x + 1)/32 - 5*sqrt(2)*atan(sqrt(2)*x - 1)/16 - 5*sqrt(2)*atan(sqrt(2)*x + 1)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_278():
    f = x**6/(x**8 + 2*x**4 + 1)
    F = -x**3/(4*x**4 + 4) + 3*sqrt(2)*log(x**2 - sqrt(2)*x + 1)/32 - 3*sqrt(2)*log(x**2 + sqrt(2)*x + 1)/32 + 3*sqrt(2)*atan(sqrt(2)*x - 1)/16 + 3*sqrt(2)*atan(sqrt(2)*x + 1)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_279():
    f = x**4/(x**8 + 2*x**4 + 1)
    F = -x/(4*x**4 + 4) - sqrt(2)*log(x**2 - sqrt(2)*x + 1)/32 + sqrt(2)*log(x**2 + sqrt(2)*x + 1)/32 + sqrt(2)*atan(sqrt(2)*x - 1)/16 + sqrt(2)*atan(sqrt(2)*x + 1)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_280():
    f = x**2/(x**8 + 2*x**4 + 1)
    F = x**3/(4*x**4 + 4) + sqrt(2)*log(x**2 - sqrt(2)*x + 1)/32 - sqrt(2)*log(x**2 + sqrt(2)*x + 1)/32 + sqrt(2)*atan(sqrt(2)*x - 1)/16 + sqrt(2)*atan(sqrt(2)*x + 1)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_281():
    f = 1/(x**8 + 2*x**4 + 1)
    F = x/(4*x**4 + 4) - 3*sqrt(2)*log(x**2 - sqrt(2)*x + 1)/32 + 3*sqrt(2)*log(x**2 + sqrt(2)*x + 1)/32 + 3*sqrt(2)*atan(sqrt(2)*x - 1)/16 + 3*sqrt(2)*atan(sqrt(2)*x + 1)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_282():
    f = 1/(x**2*(x**8 + 2*x**4 + 1))
    F = -5*sqrt(2)*log(x**2 - sqrt(2)*x + 1)/32 + 5*sqrt(2)*log(x**2 + sqrt(2)*x + 1)/32 - 5*sqrt(2)*atan(sqrt(2)*x - 1)/16 - 5*sqrt(2)*atan(sqrt(2)*x + 1)/16 - 5/(4*x) + 1/(4*x*(x**4 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_283():
    f = 1/(x**4*(x**8 + 2*x**4 + 1))
    F = 7*sqrt(2)*log(x**2 - sqrt(2)*x + 1)/32 - 7*sqrt(2)*log(x**2 + sqrt(2)*x + 1)/32 - 7*sqrt(2)*atan(sqrt(2)*x - 1)/16 - 7*sqrt(2)*atan(sqrt(2)*x + 1)/16 - 7/(12*x**3) + 1/(4*x**3*(x**4 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_284():
    f = 1/(x**6*(x**8 + 2*x**4 + 1))
    F = 9*sqrt(2)*log(x**2 - sqrt(2)*x + 1)/32 - 9*sqrt(2)*log(x**2 + sqrt(2)*x + 1)/32 + 9*sqrt(2)*atan(sqrt(2)*x - 1)/16 + 9*sqrt(2)*atan(sqrt(2)*x + 1)/16 + 9/(4*x) - 9/(20*x**5) + 1/(4*x**5*(x**4 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_285():
    f = 1/(x**8*(x**8 + 2*x**4 + 1))
    F = -11*sqrt(2)*log(x**2 - sqrt(2)*x + 1)/32 + 11*sqrt(2)*log(x**2 + sqrt(2)*x + 1)/32 + 11*sqrt(2)*atan(sqrt(2)*x - 1)/16 + 11*sqrt(2)*atan(sqrt(2)*x + 1)/16 + 11/(12*x**3) - 11/(28*x**7) + 1/(4*x**7*(x**4 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_286():
    f = x**m/(x**8 - 2*x**4 + 1)
    F = x**(m + 1)*hyper((2, m/4 + sympy.S(1)/4), (m/4 + sympy.S(5)/4,), x**4)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_287():
    f = x**9/(x**8 - 2*x**4 + 1)
    F = x**6/(4 - 4*x**4) + 3*x**2/4 - 3*atanh(x**2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_288():
    f = x**7/(x**8 - 2*x**4 + 1)
    F = log(1 - x**4)/4 + 1/(4 - 4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_289():
    f = x**5/(x**8 - 2*x**4 + 1)
    F = x**2/(4 - 4*x**4) - atanh(x**2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_290():
    f = x**3/(x**8 - 2*x**4 + 1)
    F = 1/(4 - 4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_291():
    f = x/(x**8 - 2*x**4 + 1)
    F = x**2/(4 - 4*x**4) + atanh(x**2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_292():
    f = 1/(x*(x**8 - 2*x**4 + 1))
    F = log(x) - log(1 - x**4)/4 + 1/(4 - 4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_293():
    f = 1/(x**3*(x**8 - 2*x**4 + 1))
    F = 3*atanh(x**2)/4 - 3/(4*x**2) + 1/(4*x**2*(1 - x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_294():
    f = 1/(x**5*(x**8 - 2*x**4 + 1))
    F = 2*log(x) - log(1 - x**4)/2 + 1/(4 - 4*x**4) - 1/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_295():
    f = 1/(x**7*(x**8 - 2*x**4 + 1))
    F = 5*atanh(x**2)/4 - 5/(4*x**2) - 5/(12*x**6) + 1/(4*x**6*(1 - x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_296():
    f = x**8/(x**8 - 2*x**4 + 1)
    F = x**5/(4 - 4*x**4) + 5*x/4 - 5*atan(x)/8 - 5*atanh(x)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_297():
    f = x**6/(x**8 - 2*x**4 + 1)
    F = x**3/(4 - 4*x**4) + 3*atan(x)/8 - 3*atanh(x)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_298():
    f = x**4/(x**8 - 2*x**4 + 1)
    F = x/(4 - 4*x**4) - atan(x)/8 - atanh(x)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_299():
    f = x**2/(x**8 - 2*x**4 + 1)
    F = x**3/(4 - 4*x**4) - atan(x)/8 + atanh(x)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_300():
    f = 1/(x**8 - 2*x**4 + 1)
    F = x/(4 - 4*x**4) + 3*atan(x)/8 + 3*atanh(x)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_301():
    f = 1/(x**2*(x**8 - 2*x**4 + 1))
    F = -5*atan(x)/8 + 5*atanh(x)/8 - 5/(4*x) + 1/(4*x*(1 - x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_302():
    f = 1/(x**4*(x**8 - 2*x**4 + 1))
    F = 7*atan(x)/8 + 7*atanh(x)/8 - 7/(12*x**3) + 1/(4*x**3*(1 - x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_303():
    f = 1/(x**6*(x**8 - 2*x**4 + 1))
    F = -9*atan(x)/8 + 9*atanh(x)/8 - 9/(4*x) - 9/(20*x**5) + 1/(4*x**5*(1 - x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_304():
    f = 1/(x**8*(x**8 - 2*x**4 + 1))
    F = 11*atan(x)/8 + 11*atanh(x)/8 - 11/(12*x**3) - 11/(28*x**7) + 1/(4*x**7*(1 - x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_305():
    f = x**m/(a + b*x**4 + c*x**8)
    F = -2*c*x**(m + 1)*hyper((1, m/4 + sympy.S(1)/4), (m/4 + sympy.S(5)/4,), -2*c*x**4/(b + sqrt(-4*a*c + b**2)))/((b + sqrt(-4*a*c + b**2))*(m + 1)*sqrt(-4*a*c + b**2)) + 2*c*x**(m + 1)*hyper((1, m/4 + sympy.S(1)/4), (m/4 + sympy.S(5)/4,), -2*c*x**4/(b - sqrt(-4*a*c + b**2)))/((b - sqrt(-4*a*c + b**2))*(m + 1)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_306():
    f = x**11/(a + b*x**4 + c*x**8)
    F = -b*log(a + b*x**4 + c*x**8)/(8*c**2) + x**4/(4*c) - (-2*a*c + b**2)*atanh((b + 2*c*x**4)/sqrt(-4*a*c + b**2))/(4*c**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_307():
    f = x**9/(a + b*x**4 + c*x**8)
    F = x**2/(2*c) - sqrt(2)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x**2/sqrt(b - sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2))) - sqrt(2)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x**2/sqrt(b + sqrt(-4*a*c + b**2)))/(4*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_308():
    f = x**7/(a + b*x**4 + c*x**8)
    F = b*atanh((b + 2*c*x**4)/sqrt(-4*a*c + b**2))/(4*c*sqrt(-4*a*c + b**2)) + log(a + b*x**4 + c*x**8)/(8*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_309():
    f = x**5/(a + b*x**4 + c*x**8)
    F = -sqrt(2)*sqrt(b - sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x**2/sqrt(b - sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(-4*a*c + b**2)) + sqrt(2)*sqrt(b + sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x**2/sqrt(b + sqrt(-4*a*c + b**2)))/(4*sqrt(c)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_310():
    f = x**3/(a + b*x**4 + c*x**8)
    F = -atanh((b + 2*c*x**4)/sqrt(-4*a*c + b**2))/(2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_311():
    f = x/(a + b*x**4 + c*x**8)
    F = -sqrt(2)*sqrt(c)*atan(sqrt(2)*sqrt(c)*x**2/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)) + sqrt(2)*sqrt(c)*atan(sqrt(2)*sqrt(c)*x**2/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(b - sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_312():
    f = 1/(x*(a + b*x**4 + c*x**8))
    F = b*atanh((b + 2*c*x**4)/sqrt(-4*a*c + b**2))/(4*a*sqrt(-4*a*c + b**2)) + log(x)/a - log(a + b*x**4 + c*x**8)/(8*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_313():
    f = 1/(x**3*(a + b*x**4 + c*x**8))
    F = -sqrt(2)*sqrt(c)*(-b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(2)*sqrt(c)*x**2/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(b + sqrt(-4*a*c + b**2))) - sqrt(2)*sqrt(c)*(b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(2)*sqrt(c)*x**2/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*sqrt(b - sqrt(-4*a*c + b**2))) - 1/(2*a*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_314():
    f = 1/(x**5*(a + b*x**4 + c*x**8))
    F = -1/(4*a*x**4) - b*log(x)/a**2 + b*log(a + b*x**4 + c*x**8)/(8*a**2) - (-2*a*c + b**2)*atanh((b + 2*c*x**4)/sqrt(-4*a*c + b**2))/(4*a**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_315():
    f = x**10/(a + b*x**4 + c*x**8)
    F = x**3/(3*c) - 2**(sympy.S(1)/4)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(7)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) + 2**(sympy.S(1)/4)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(7)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) - 2**(sympy.S(1)/4)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(7)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) + 2**(sympy.S(1)/4)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(7)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_316():
    f = x**8/(a + b*x**4 + c*x**8)
    F = x/c + 2**(sympy.S(3)/4)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(5)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) + 2**(sympy.S(3)/4)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(5)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) + 2**(sympy.S(3)/4)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(5)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) + 2**(sympy.S(3)/4)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(5)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_317():
    f = x**6/(a + b*x**4 + c*x**8)
    F = -2**(sympy.S(1)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(1)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(1)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(1)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(3)/4)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_318():
    f = x**4/(a + b*x**4 + c*x**8)
    F = 2**(sympy.S(3)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(1)/4)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(3)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(1)/4)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(3)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(1)/4)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(3)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(1)/4)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_319():
    f = x**2/(a + b*x**4 + c*x**8)
    F = 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_320():
    f = 1/(a + b*x**4 + c*x**8)
    F = -2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(2*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_321():
    f = 1/(x**2*(a + b*x**4 + c*x**8))
    F = -2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(b/sqrt(-4*a*c + b**2) + 1)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*a*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) + 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(b/sqrt(-4*a*c + b**2) + 1)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*a*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) - 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(-b/sqrt(-4*a*c + b**2) + 1)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*a*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) + 2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*(-b/sqrt(-4*a*c + b**2) + 1)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*a*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)) - 1/(a*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_322():
    f = 1/(x**4*(a + b*x**4 + c*x**8))
    F = 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(b/sqrt(-4*a*c + b**2) + 1)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*a*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(b/sqrt(-4*a*c + b**2) + 1)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*a*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(-b/sqrt(-4*a*c + b**2) + 1)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*a*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) + 2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*(-b/sqrt(-4*a*c + b**2) + 1)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*a*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) - 1/(3*a*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_323():
    f = x**m/(x**8 + x**4 + 1)
    F = 2*sqrt(3)*x**(m + 1)*hyper((1, m/4 + sympy.S(1)/4), (m/4 + sympy.S(5)/4,), -2*x**4/(1 - sqrt(3)*I))/(3*(sqrt(3) + I)*(m + 1)) - 2*sqrt(3)*x**(m + 1)*hyper((1, m/4 + sympy.S(1)/4), (m/4 + sympy.S(5)/4,), -2*x**4/(1 + sqrt(3)*I))/(3*(-sqrt(3) + I)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_324():
    f = x**11/(x**8 + x**4 + 1)
    F = x**4/4 - log(x**8 + x**4 + 1)/8 - sqrt(3)*atan(sqrt(3)*(2*x**4 + 1)/3)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_325():
    f = x**9/(x**8 + x**4 + 1)
    F = x**2/2 + sqrt(3)*atan(sqrt(3)*(1 - 2*x**2)/3)/6 - sqrt(3)*atan(sqrt(3)*(2*x**2 + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_326():
    f = x**7/(x**8 + x**4 + 1)
    F = log(x**8 + x**4 + 1)/8 - sqrt(3)*atan(sqrt(3)*(2*x**4 + 1)/3)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_327():
    f = x**5/(x**8 + x**4 + 1)
    F = log(x**4 - x**2 + 1)/8 - log(x**4 + x**2 + 1)/8 - sqrt(3)*atan(sqrt(3)*(1 - 2*x**2)/3)/12 + sqrt(3)*atan(sqrt(3)*(2*x**2 + 1)/3)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_328():
    f = x**3/(x**8 + x**4 + 1)
    F = sqrt(3)*atan(sqrt(3)*(2*x**4 + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_329():
    f = x/(x**8 + x**4 + 1)
    F = -log(x**4 - x**2 + 1)/8 + log(x**4 + x**2 + 1)/8 - sqrt(3)*atan(sqrt(3)*(1 - 2*x**2)/3)/12 + sqrt(3)*atan(sqrt(3)*(2*x**2 + 1)/3)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_330():
    f = 1/(x*(x**8 + x**4 + 1))
    F = log(x) - log(x**8 + x**4 + 1)/8 - sqrt(3)*atan(sqrt(3)*(2*x**4 + 1)/3)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_331():
    f = 1/(x**3*(x**8 + x**4 + 1))
    F = sqrt(3)*atan(sqrt(3)*(1 - 2*x**2)/3)/6 - sqrt(3)*atan(sqrt(3)*(2*x**2 + 1)/3)/6 - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_332():
    f = 1/(x**5*(x**8 + x**4 + 1))
    F = -log(x) + log(x**8 + x**4 + 1)/8 - sqrt(3)*atan(sqrt(3)*(2*x**4 + 1)/3)/12 - 1/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_333():
    f = 1/(x**7*(x**8 + x**4 + 1))
    F = log(x**4 - x**2 + 1)/8 - log(x**4 + x**2 + 1)/8 - sqrt(3)*atan(sqrt(3)*(1 - 2*x**2)/3)/12 + sqrt(3)*atan(sqrt(3)*(2*x**2 + 1)/3)/12 + 1/(2*x**2) - 1/(6*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_334():
    f = x**8/(x**8 + x**4 + 1)
    F = x + log(x**2 - x + 1)/8 - log(x**2 + x + 1)/8 + sqrt(3)*log(x**2 - sqrt(3)*x + 1)/24 - sqrt(3)*log(x**2 + sqrt(3)*x + 1)/24 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/12 - sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/12 - atan(2*x - sqrt(3))/4 - atan(2*x + sqrt(3))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_335():
    f = x**6/(x**8 + x**4 + 1)
    F = sqrt(3)*log(x**2 - sqrt(3)*x + 1)/12 - sqrt(3)*log(x**2 + sqrt(3)*x + 1)/12 - sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/6 + sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_336():
    f = x**4/(x**8 + x**4 + 1)
    F = -log(x**2 - x + 1)/8 + log(x**2 + x + 1)/8 + sqrt(3)*log(x**2 - sqrt(3)*x + 1)/24 - sqrt(3)*log(x**2 + sqrt(3)*x + 1)/24 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/12 - sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/12 + atan(2*x - sqrt(3))/4 + atan(2*x + sqrt(3))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_337():
    f = x**2/(x**8 + x**4 + 1)
    F = log(x**2 - x + 1)/8 - log(x**2 + x + 1)/8 - sqrt(3)*log(x**2 - sqrt(3)*x + 1)/24 + sqrt(3)*log(x**2 + sqrt(3)*x + 1)/24 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/12 - sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/12 + atan(2*x - sqrt(3))/4 + atan(2*x + sqrt(3))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_338():
    f = 1/(x**8 + x**4 + 1)
    F = -sqrt(3)*log(x**2 - sqrt(3)*x + 1)/12 + sqrt(3)*log(x**2 + sqrt(3)*x + 1)/12 - sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/6 + sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_339():
    f = 1/(x**2*(x**8 + x**4 + 1))
    F = -log(x**2 - x + 1)/8 + log(x**2 + x + 1)/8 - sqrt(3)*log(x**2 - sqrt(3)*x + 1)/24 + sqrt(3)*log(x**2 + sqrt(3)*x + 1)/24 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/12 - sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/12 - atan(2*x - sqrt(3))/4 - atan(2*x + sqrt(3))/4 - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_340():
    f = 1/(x**4*(x**8 + x**4 + 1))
    F = log(x**2 - x + 1)/8 - log(x**2 + x + 1)/8 + sqrt(3)*log(x**2 - sqrt(3)*x + 1)/24 - sqrt(3)*log(x**2 + sqrt(3)*x + 1)/24 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/12 - sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/12 - atan(2*x - sqrt(3))/4 - atan(2*x + sqrt(3))/4 - 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_341():
    f = 1/(x**6*(x**8 + x**4 + 1))
    F = sqrt(3)*log(x**2 - sqrt(3)*x + 1)/12 - sqrt(3)*log(x**2 + sqrt(3)*x + 1)/12 - sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/6 + sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/6 + 1/x - 1/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_342():
    f = 1/(x**8*(x**8 + x**4 + 1))
    F = -log(x**2 - x + 1)/8 + log(x**2 + x + 1)/8 + sqrt(3)*log(x**2 - sqrt(3)*x + 1)/24 - sqrt(3)*log(x**2 + sqrt(3)*x + 1)/24 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/12 - sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/12 + atan(2*x - sqrt(3))/4 + atan(2*x + sqrt(3))/4 + 1/(3*x**3) - 1/(7*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_343():
    f = x**m/(x**8 - x**4 + 1)
    F = 2*sqrt(3)*x**(m + 1)*hyper((1, m/4 + sympy.S(1)/4), (m/4 + sympy.S(5)/4,), 2*x**4/(1 - sqrt(3)*I))/(3*(sqrt(3) + I)*(m + 1)) - 2*sqrt(3)*x**(m + 1)*hyper((1, m/4 + sympy.S(1)/4), (m/4 + sympy.S(5)/4,), 2*x**4/(1 + sqrt(3)*I))/(3*(-sqrt(3) + I)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_344():
    f = x**11/(x**8 - x**4 + 1)
    F = x**4/4 + log(x**8 - x**4 + 1)/8 + sqrt(3)*atan(sqrt(3)*(1 - 2*x**4)/3)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_345():
    f = x**9/(x**8 - x**4 + 1)
    F = x**2/2 + sqrt(3)*log(x**4 - sqrt(3)*x**2 + 1)/12 - sqrt(3)*log(x**4 + sqrt(3)*x**2 + 1)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_346():
    f = x**7/(x**8 - x**4 + 1)
    F = log(x**8 - x**4 + 1)/8 - sqrt(3)*atan(sqrt(3)*(1 - 2*x**4)/3)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_347():
    f = x**5/(x**8 - x**4 + 1)
    F = sqrt(3)*log(x**4 - sqrt(3)*x**2 + 1)/24 - sqrt(3)*log(x**4 + sqrt(3)*x**2 + 1)/24 + atan(2*x**2 - sqrt(3))/4 + atan(2*x**2 + sqrt(3))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_348():
    f = x**3/(x**8 - x**4 + 1)
    F = -sqrt(3)*atan(sqrt(3)*(1 - 2*x**4)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_349():
    f = x/(x**8 - x**4 + 1)
    F = -sqrt(3)*log(x**4 - sqrt(3)*x**2 + 1)/24 + sqrt(3)*log(x**4 + sqrt(3)*x**2 + 1)/24 + atan(2*x**2 - sqrt(3))/4 + atan(2*x**2 + sqrt(3))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_350():
    f = 1/(x*(x**8 - x**4 + 1))
    F = log(x) - log(x**8 - x**4 + 1)/8 - sqrt(3)*atan(sqrt(3)*(1 - 2*x**4)/3)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_351():
    f = 1/(x**3*(x**8 - x**4 + 1))
    F = -sqrt(3)*log(x**4 - sqrt(3)*x**2 + 1)/12 + sqrt(3)*log(x**4 + sqrt(3)*x**2 + 1)/12 - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_352():
    f = 1/(x**5*(x**8 - x**4 + 1))
    F = log(x) - log(x**8 - x**4 + 1)/8 + sqrt(3)*atan(sqrt(3)*(1 - 2*x**4)/3)/12 - 1/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_353():
    f = 1/(x**7*(x**8 - x**4 + 1))
    F = -sqrt(3)*log(x**4 - sqrt(3)*x**2 + 1)/24 + sqrt(3)*log(x**4 + sqrt(3)*x**2 + 1)/24 - atan(2*x**2 - sqrt(3))/4 - atan(2*x**2 + sqrt(3))/4 - 1/(2*x**2) - 1/(6*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_354():
    f = x**8/(x**8 - x**4 + 1)
    F = x - sqrt(sympy.S(2)/3 - sqrt(3)/3)*log(x**2 - x*sqrt(2 - sqrt(3)) + 1)/8 + sqrt(sympy.S(2)/3 - sqrt(3)/3)*log(x**2 + x*sqrt(2 - sqrt(3)) + 1)/8 + sqrt(sqrt(3)/3 + sympy.S(2)/3)*log(x**2 - x*sqrt(sqrt(3) + 2) + 1)/8 - sqrt(sqrt(3)/3 + sympy.S(2)/3)*log(x**2 + x*sqrt(sqrt(3) + 2) + 1)/8 - atan((-2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/(4*sqrt(3*sqrt(3) + 6)) + atan((2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/(4*sqrt(3*sqrt(3) + 6)) + atan((-2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/(4*sqrt(6 - 3*sqrt(3))) - atan((2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/(4*sqrt(6 - 3*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_355():
    f = x**6/(x**8 - x**4 + 1)
    F = sqrt(6)*log(x**2 - x*sqrt(2 - sqrt(3)) + 1)/24 - sqrt(6)*log(x**2 + x*sqrt(2 - sqrt(3)) + 1)/24 + sqrt(6)*log(x**2 - x*sqrt(sqrt(3) + 2) + 1)/24 - sqrt(6)*log(x**2 + x*sqrt(sqrt(3) + 2) + 1)/24 - sqrt(6)*atan((-2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/12 + sqrt(6)*atan((2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/12 - sqrt(6)*atan((-2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/12 + sqrt(6)*atan((2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_356():
    f = x**4/(x**8 - x**4 + 1)
    F = -log(x**2 - x*sqrt(2 - sqrt(3)) + 1)/(8*sqrt(6 - 3*sqrt(3))) + log(x**2 + x*sqrt(2 - sqrt(3)) + 1)/(8*sqrt(6 - 3*sqrt(3))) + log(x**2 - x*sqrt(sqrt(3) + 2) + 1)/(8*sqrt(3*sqrt(3) + 6)) - log(x**2 + x*sqrt(sqrt(3) + 2) + 1)/(8*sqrt(3*sqrt(3) + 6)) - atan((-2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/(4*sqrt(6 - 3*sqrt(3))) + atan((2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/(4*sqrt(6 - 3*sqrt(3))) + atan((-2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/(4*sqrt(3*sqrt(3) + 6)) - atan((2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/(4*sqrt(3*sqrt(3) + 6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_357():
    f = x**2/(x**8 - x**4 + 1)
    F = log(x**2 - x*sqrt(2 - sqrt(3)) + 1)/(8*sqrt(6 - 3*sqrt(3))) - log(x**2 + x*sqrt(2 - sqrt(3)) + 1)/(8*sqrt(6 - 3*sqrt(3))) - log(x**2 - x*sqrt(sqrt(3) + 2) + 1)/(8*sqrt(3*sqrt(3) + 6)) + log(x**2 + x*sqrt(sqrt(3) + 2) + 1)/(8*sqrt(3*sqrt(3) + 6)) - sqrt(sqrt(3)/3 + sympy.S(2)/3)*atan((-2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/4 + sqrt(sqrt(3)/3 + sympy.S(2)/3)*atan((2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/4 + sqrt(sympy.S(2)/3 - sqrt(3)/3)*atan((-2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/4 - sqrt(sympy.S(2)/3 - sqrt(3)/3)*atan((2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_358():
    f = 1/(x**8 - x**4 + 1)
    F = -sqrt(6)*log(x**2 - x*sqrt(2 - sqrt(3)) + 1)/24 + sqrt(6)*log(x**2 + x*sqrt(2 - sqrt(3)) + 1)/24 - sqrt(6)*log(x**2 - x*sqrt(sqrt(3) + 2) + 1)/24 + sqrt(6)*log(x**2 + x*sqrt(sqrt(3) + 2) + 1)/24 - sqrt(6)*atan((-2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/12 + sqrt(6)*atan((2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/12 - sqrt(6)*atan((-2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/12 + sqrt(6)*atan((2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_359():
    f = 1/(x**2*(x**8 - x**4 + 1))
    F = sqrt(sympy.S(2)/3 - sqrt(3)/3)*log(x**2 - x*sqrt(2 - sqrt(3)) + 1)/8 - sqrt(sympy.S(2)/3 - sqrt(3)/3)*log(x**2 + x*sqrt(2 - sqrt(3)) + 1)/8 - sqrt(sqrt(3)/3 + sympy.S(2)/3)*log(x**2 - x*sqrt(sqrt(3) + 2) + 1)/8 + sqrt(sqrt(3)/3 + sympy.S(2)/3)*log(x**2 + x*sqrt(sqrt(3) + 2) + 1)/8 - atan((-2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/(4*sqrt(3*sqrt(3) + 6)) + atan((2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/(4*sqrt(3*sqrt(3) + 6)) + atan((-2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/(4*sqrt(6 - 3*sqrt(3))) - atan((2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/(4*sqrt(6 - 3*sqrt(3))) - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_360():
    f = 1/(x**4*(x**8 - x**4 + 1))
    F = sqrt(sympy.S(2)/3 - sqrt(3)/3)*log(x**2 - x*sqrt(2 - sqrt(3)) + 1)/8 - sqrt(sympy.S(2)/3 - sqrt(3)/3)*log(x**2 + x*sqrt(2 - sqrt(3)) + 1)/8 - sqrt(sqrt(3)/3 + sympy.S(2)/3)*log(x**2 - x*sqrt(sqrt(3) + 2) + 1)/8 + sqrt(sqrt(3)/3 + sympy.S(2)/3)*log(x**2 + x*sqrt(sqrt(3) + 2) + 1)/8 + sqrt(sympy.S(2)/3 - sqrt(3)/3)*atan((-2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/4 - sqrt(sympy.S(2)/3 - sqrt(3)/3)*atan((2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/4 - sqrt(sqrt(3)/3 + sympy.S(2)/3)*atan((-2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/4 + sqrt(sqrt(3)/3 + sympy.S(2)/3)*atan((2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/4 - 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_361():
    f = 1/(x**6*(x**8 - x**4 + 1))
    F = -sqrt(6)*log(x**2 - x*sqrt(2 - sqrt(3)) + 1)/24 + sqrt(6)*log(x**2 + x*sqrt(2 - sqrt(3)) + 1)/24 - sqrt(6)*log(x**2 - x*sqrt(sqrt(3) + 2) + 1)/24 + sqrt(6)*log(x**2 + x*sqrt(sqrt(3) + 2) + 1)/24 + sqrt(6)*atan((-2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/12 - sqrt(6)*atan((2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/12 + sqrt(6)*atan((-2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/12 - sqrt(6)*atan((2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/12 - 1/x - 1/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_362():
    f = 1/(x**8*(x**8 - x**4 + 1))
    F = sqrt(sqrt(3)/3 + sympy.S(2)/3)*log(x**2 - x*sqrt(2 - sqrt(3)) + 1)/8 - sqrt(sqrt(3)/3 + sympy.S(2)/3)*log(x**2 + x*sqrt(2 - sqrt(3)) + 1)/8 - sqrt(sympy.S(2)/3 - sqrt(3)/3)*log(x**2 - x*sqrt(sqrt(3) + 2) + 1)/8 + sqrt(sympy.S(2)/3 - sqrt(3)/3)*log(x**2 + x*sqrt(sqrt(3) + 2) + 1)/8 + sqrt(sqrt(3)/3 + sympy.S(2)/3)*atan((-2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/4 - sqrt(sqrt(3)/3 + sympy.S(2)/3)*atan((2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/4 - sqrt(sympy.S(2)/3 - sqrt(3)/3)*atan((-2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/4 + sqrt(sympy.S(2)/3 - sqrt(3)/3)*atan((2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/4 - 1/(3*x**3) - 1/(7*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_363():
    f = x**m/(x**8 + 3*x**4 + 1)
    F = 2*sqrt(5)*x**(m + 1)*hyper((1, m/4 + sympy.S(1)/4), (m/4 + sympy.S(5)/4,), -2*x**4/(3 - sqrt(5)))/(5*(3 - sqrt(5))*(m + 1)) - 2*sqrt(5)*x**(m + 1)*hyper((1, m/4 + sympy.S(1)/4), (m/4 + sympy.S(5)/4,), -2*x**4/(sqrt(5) + 3))/(5*(sqrt(5) + 3)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_364():
    f = x**11/(x**8 + 3*x**4 + 1)
    F = x**4/4 - (sympy.S(3)/8 - 7*sqrt(5)/40)*log(2*x**4 - sqrt(5) + 3) - (sympy.S(3)/8 + 7*sqrt(5)/40)*log(2*x**4 + sqrt(5) + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_365():
    f = x**9/(x**8 + 3*x**4 + 1)
    F = x**2/2 + sqrt(sympy.S(9)/5 - 4*sqrt(5)/5)*atan(x**2*sqrt(sqrt(5)/2 + sympy.S(3)/2))/2 - sqrt(4*sqrt(5)/5 + sympy.S(9)/5)*atan(sqrt(2)*x**2/sqrt(sqrt(5) + 3))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_366():
    f = x**7/(x**8 + 3*x**4 + 1)
    F = (sympy.S(1)/8 - 3*sqrt(5)/40)*log(2*x**4 - sqrt(5) + 3) + (sympy.S(1)/8 + 3*sqrt(5)/40)*log(2*x**4 + sqrt(5) + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_367():
    f = x**5/(x**8 + 3*x**4 + 1)
    F = -sqrt(sympy.S(3)/10 - sqrt(5)/10)*atan(x**2*sqrt(sqrt(5)/2 + sympy.S(3)/2))/2 + sqrt(sqrt(5)/10 + sympy.S(3)/10)*atan(sqrt(2)*x**2/sqrt(sqrt(5) + 3))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_368():
    f = x**3/(x**8 + 3*x**4 + 1)
    F = -sqrt(5)*atanh(sqrt(5)*(2*x**4 + 3)/5)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_369():
    f = x/(x**8 + 3*x**4 + 1)
    F = sqrt(sqrt(5)/10 + sympy.S(3)/10)*atan(x**2*sqrt(sqrt(5)/2 + sympy.S(3)/2))/2 - atan(sqrt(2)*x**2/sqrt(sqrt(5) + 3))/sqrt(10*sqrt(5) + 30)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_370():
    f = 1/(x*(x**8 + 3*x**4 + 1))
    F = log(x) - (sympy.S(1)/8 + 3*sqrt(5)/40)*log(2*x**4 - sqrt(5) + 3) - (sympy.S(1)/8 - 3*sqrt(5)/40)*log(2*x**4 + sqrt(5) + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_371():
    f = 1/(x**3*(x**8 + 3*x**4 + 1))
    F = -sqrt(10)*(sqrt(5) + 3)**(sympy.S(3)/2)*atan(x**2*sqrt(sqrt(5)/2 + sympy.S(3)/2))/40 + sqrt(sympy.S(9)/5 - 4*sqrt(5)/5)*atan(sqrt(2)*x**2/sqrt(sqrt(5) + 3))/2 - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_372():
    f = 1/(x**5*(x**8 + 3*x**4 + 1))
    F = -3*log(x) + (sympy.S(3)/8 + 7*sqrt(5)/40)*log(2*x**4 - sqrt(5) + 3) + (sympy.S(3)/8 - 7*sqrt(5)/40)*log(2*x**4 + sqrt(5) + 3) - 1/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_373():
    f = 1/(x**7*(x**8 + 3*x**4 + 1))
    F = sqrt(11*sqrt(5)/2 + sympy.S(123)/10)*atan(x**2*sqrt(sqrt(5)/2 + sympy.S(3)/2))/2 - sqrt(sympy.S(123)/10 - 11*sqrt(5)/2)*atan(sqrt(2)*x**2/sqrt(sqrt(5) + 3))/2 + 3/(2*x**2) - 1/(6*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_374():
    f = x**8/(x**8 + 3*x**4 + 1)
    F = x - 2**(sympy.S(1)/4)*sqrt(5)*(123 - 55*sqrt(5))**(sympy.S(1)/4)*log(2*x**2 - 2*x*(6 - 2*sqrt(5))**(sympy.S(1)/4) + sqrt(6 - 2*sqrt(5)))/40 + 2**(sympy.S(1)/4)*sqrt(5)*(123 - 55*sqrt(5))**(sympy.S(1)/4)*log(2*x**2 + 2*x*(6 - 2*sqrt(5))**(sympy.S(1)/4) + sqrt(6 - 2*sqrt(5)))/40 + 2**(sympy.S(1)/4)*sqrt(5)*(55*sqrt(5) + 123)**(sympy.S(1)/4)*log(2*x**2 - 2*x*(2*sqrt(5) + 6)**(sympy.S(1)/4) + sqrt(2*sqrt(5) + 6))/40 - 2**(sympy.S(1)/4)*sqrt(5)*(55*sqrt(5) + 123)**(sympy.S(1)/4)*log(2*x**2 + 2*x*(2*sqrt(5) + 6)**(sympy.S(1)/4) + sqrt(2*sqrt(5) + 6))/40 + 2**(sympy.S(1)/4)*sqrt(5)*(123 - 55*sqrt(5))**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(3 - sqrt(5))**(sympy.S(1)/4) - 1)/20 + 2**(sympy.S(1)/4)*sqrt(5)*(123 - 55*sqrt(5))**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(3 - sqrt(5))**(sympy.S(1)/4) + 1)/20 - 2**(sympy.S(1)/4)*sqrt(5)*(55*sqrt(5) + 123)**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4) - 1)/20 - 2**(sympy.S(1)/4)*sqrt(5)*(55*sqrt(5) + 123)**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4) + 1)/20
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_375():
    f = x**6/(x**8 + 3*x**4 + 1)
    F = -2**(sympy.S(3)/4)*sqrt(5)*(3 - sqrt(5))**(sympy.S(3)/4)*log(2*x**2 - 2*x*(6 - 2*sqrt(5))**(sympy.S(1)/4) + sqrt(6 - 2*sqrt(5)))/80 + 2**(sympy.S(3)/4)*sqrt(5)*(3 - sqrt(5))**(sympy.S(3)/4)*log(2*x**2 + 2*x*(6 - 2*sqrt(5))**(sympy.S(1)/4) + sqrt(6 - 2*sqrt(5)))/80 + 2**(sympy.S(3)/4)*sqrt(5)*(sqrt(5) + 3)**(sympy.S(3)/4)*log(2*x**2 - 2*x*(2*sqrt(5) + 6)**(sympy.S(1)/4) + sqrt(2*sqrt(5) + 6))/80 - 2**(sympy.S(3)/4)*sqrt(5)*(sqrt(5) + 3)**(sympy.S(3)/4)*log(2*x**2 + 2*x*(2*sqrt(5) + 6)**(sympy.S(1)/4) + sqrt(2*sqrt(5) + 6))/80 - 2**(sympy.S(3)/4)*sqrt(5)*(3 - sqrt(5))**(sympy.S(3)/4)*atan(2**(sympy.S(3)/4)*x/(3 - sqrt(5))**(sympy.S(1)/4) - 1)/40 - 2**(sympy.S(3)/4)*sqrt(5)*(3 - sqrt(5))**(sympy.S(3)/4)*atan(2**(sympy.S(3)/4)*x/(3 - sqrt(5))**(sympy.S(1)/4) + 1)/40 + 2**(sympy.S(3)/4)*sqrt(5)*(sqrt(5) + 3)**(sympy.S(3)/4)*atan(2**(sympy.S(3)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4) - 1)/40 + 2**(sympy.S(3)/4)*sqrt(5)*(sqrt(5) + 3)**(sympy.S(3)/4)*atan(2**(sympy.S(3)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4) + 1)/40
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_376():
    f = x**4/(x**8 + 3*x**4 + 1)
    F = 2**(sympy.S(1)/4)*sqrt(5)*(3 - sqrt(5))**(sympy.S(1)/4)*log(2*x**2 - 2*x*(6 - 2*sqrt(5))**(sympy.S(1)/4) + sqrt(6 - 2*sqrt(5)))/40 - 2**(sympy.S(1)/4)*sqrt(5)*(3 - sqrt(5))**(sympy.S(1)/4)*log(2*x**2 + 2*x*(6 - 2*sqrt(5))**(sympy.S(1)/4) + sqrt(6 - 2*sqrt(5)))/40 - 2**(sympy.S(1)/4)*sqrt(5)*(sqrt(5) + 3)**(sympy.S(1)/4)*log(2*x**2 - 2*x*(2*sqrt(5) + 6)**(sympy.S(1)/4) + sqrt(2*sqrt(5) + 6))/40 + 2**(sympy.S(1)/4)*sqrt(5)*(sqrt(5) + 3)**(sympy.S(1)/4)*log(2*x**2 + 2*x*(2*sqrt(5) + 6)**(sympy.S(1)/4) + sqrt(2*sqrt(5) + 6))/40 - 2**(sympy.S(1)/4)*sqrt(5)*(3 - sqrt(5))**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(3 - sqrt(5))**(sympy.S(1)/4) - 1)/20 - 2**(sympy.S(1)/4)*sqrt(5)*(3 - sqrt(5))**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(3 - sqrt(5))**(sympy.S(1)/4) + 1)/20 + 2**(sympy.S(1)/4)*sqrt(5)*(sqrt(5) + 3)**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4) - 1)/20 + 2**(sympy.S(1)/4)*sqrt(5)*(sqrt(5) + 3)**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4) + 1)/20
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_377():
    f = x**2/(x**8 + 3*x**4 + 1)
    F = sqrt(5)*log(2*x**2 - 2*x*(6 - 2*sqrt(5))**(sympy.S(1)/4) + sqrt(6 - 2*sqrt(5)))/(20*(6 - 2*sqrt(5))**(sympy.S(1)/4)) - sqrt(5)*log(2*x**2 + 2*x*(6 - 2*sqrt(5))**(sympy.S(1)/4) + sqrt(6 - 2*sqrt(5)))/(20*(6 - 2*sqrt(5))**(sympy.S(1)/4)) - sqrt(5)*log(2*x**2 - 2*x*(2*sqrt(5) + 6)**(sympy.S(1)/4) + sqrt(2*sqrt(5) + 6))/(20*(2*sqrt(5) + 6)**(sympy.S(1)/4)) + sqrt(5)*log(2*x**2 + 2*x*(2*sqrt(5) + 6)**(sympy.S(1)/4) + sqrt(2*sqrt(5) + 6))/(20*(2*sqrt(5) + 6)**(sympy.S(1)/4)) + sqrt(5)*atan(2**(sympy.S(3)/4)*x/(3 - sqrt(5))**(sympy.S(1)/4) - 1)/(10*(6 - 2*sqrt(5))**(sympy.S(1)/4)) + sqrt(5)*atan(2**(sympy.S(3)/4)*x/(3 - sqrt(5))**(sympy.S(1)/4) + 1)/(10*(6 - 2*sqrt(5))**(sympy.S(1)/4)) - sqrt(5)*atan(2**(sympy.S(3)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4) - 1)/(10*(2*sqrt(5) + 6)**(sympy.S(1)/4)) - sqrt(5)*atan(2**(sympy.S(3)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4) + 1)/(10*(2*sqrt(5) + 6)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_378():
    f = 1/(x**2*(x**8 + 3*x**4 + 1))
    F = -sqrt(5)*(110*sqrt(5) + 246)**(sympy.S(1)/4)*log(2*x**2 - 2*x*(6 - 2*sqrt(5))**(sympy.S(1)/4) + sqrt(6 - 2*sqrt(5)))/40 + sqrt(5)*(110*sqrt(5) + 246)**(sympy.S(1)/4)*log(2*x**2 + 2*x*(6 - 2*sqrt(5))**(sympy.S(1)/4) + sqrt(6 - 2*sqrt(5)))/40 + (6150 - 2750*sqrt(5))**(sympy.S(1)/4)*log(2*x**2 - 2*x*(2*sqrt(5) + 6)**(sympy.S(1)/4) + sqrt(2*sqrt(5) + 6))/40 - (6150 - 2750*sqrt(5))**(sympy.S(1)/4)*log(2*x**2 + 2*x*(2*sqrt(5) + 6)**(sympy.S(1)/4) + sqrt(2*sqrt(5) + 6))/40 - sqrt(5)*(110*sqrt(5) + 246)**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(3 - sqrt(5))**(sympy.S(1)/4) - 1)/20 - sqrt(5)*(110*sqrt(5) + 246)**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(3 - sqrt(5))**(sympy.S(1)/4) + 1)/20 + (6150 - 2750*sqrt(5))**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4) - 1)/20 + (6150 - 2750*sqrt(5))**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4) + 1)/20 - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_379():
    f = 1/(x**4*(x**8 + 3*x**4 + 1))
    F = 2**(sympy.S(1)/4)*sqrt(5)*(377*sqrt(5) + 843)**(sympy.S(1)/4)*log(2*x**2 - 2*x*(6 - 2*sqrt(5))**(sympy.S(1)/4) + sqrt(6 - 2*sqrt(5)))/40 - 2**(sympy.S(1)/4)*sqrt(5)*(377*sqrt(5) + 843)**(sympy.S(1)/4)*log(2*x**2 + 2*x*(6 - 2*sqrt(5))**(sympy.S(1)/4) + sqrt(6 - 2*sqrt(5)))/40 - 2**(sympy.S(1)/4)*sqrt(5)*(843 - 377*sqrt(5))**(sympy.S(1)/4)*log(2*x**2 - 2*x*(2*sqrt(5) + 6)**(sympy.S(1)/4) + sqrt(2*sqrt(5) + 6))/40 + 2**(sympy.S(1)/4)*sqrt(5)*(843 - 377*sqrt(5))**(sympy.S(1)/4)*log(2*x**2 + 2*x*(2*sqrt(5) + 6)**(sympy.S(1)/4) + sqrt(2*sqrt(5) + 6))/40 - 2**(sympy.S(1)/4)*sqrt(5)*(377*sqrt(5) + 843)**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(3 - sqrt(5))**(sympy.S(1)/4) - 1)/20 - 2**(sympy.S(1)/4)*sqrt(5)*(377*sqrt(5) + 843)**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(3 - sqrt(5))**(sympy.S(1)/4) + 1)/20 + 2**(sympy.S(1)/4)*sqrt(5)*(843 - 377*sqrt(5))**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4) - 1)/20 + 2**(sympy.S(1)/4)*sqrt(5)*(843 - 377*sqrt(5))**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4) + 1)/20 - 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_380():
    f = x**m/(x**8 - 3*x**4 + 1)
    F = 2*sqrt(5)*x**(m + 1)*hyper((1, m/4 + sympy.S(1)/4), (m/4 + sympy.S(5)/4,), 2*x**4/(3 - sqrt(5)))/(5*(3 - sqrt(5))*(m + 1)) - 2*sqrt(5)*x**(m + 1)*hyper((1, m/4 + sympy.S(1)/4), (m/4 + sympy.S(5)/4,), 2*x**4/(sqrt(5) + 3))/(5*(sqrt(5) + 3)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_381():
    f = x**11/(x**8 - 3*x**4 + 1)
    F = x**4/4 + (sympy.S(3)/8 - 7*sqrt(5)/40)*log(-2*x**4 - sqrt(5) + 3) + (sympy.S(3)/8 + 7*sqrt(5)/40)*log(-2*x**4 + sqrt(5) + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_382():
    f = x**9/(x**8 - 3*x**4 + 1)
    F = x**2/2 + sqrt(sympy.S(9)/5 - 4*sqrt(5)/5)*atanh(x**2*sqrt(sqrt(5)/2 + sympy.S(3)/2))/2 - sqrt(4*sqrt(5)/5 + sympy.S(9)/5)*atanh(sqrt(2)*x**2/sqrt(sqrt(5) + 3))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_383():
    f = x**7/(x**8 - 3*x**4 + 1)
    F = (sympy.S(1)/8 - 3*sqrt(5)/40)*log(-2*x**4 - sqrt(5) + 3) + (sympy.S(1)/8 + 3*sqrt(5)/40)*log(-2*x**4 + sqrt(5) + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_384():
    f = x**5/(x**8 - 3*x**4 + 1)
    F = sqrt(sympy.S(3)/10 - sqrt(5)/10)*atanh(x**2*sqrt(sqrt(5)/2 + sympy.S(3)/2))/2 - sqrt(sqrt(5)/10 + sympy.S(3)/10)*atanh(sqrt(2)*x**2/sqrt(sqrt(5) + 3))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_385():
    f = x**3/(x**8 - 3*x**4 + 1)
    F = sqrt(5)*atanh(sqrt(5)*(3 - 2*x**4)/5)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_386():
    f = x/(x**8 - 3*x**4 + 1)
    F = sqrt(sqrt(5)/10 + sympy.S(3)/10)*atanh(x**2*sqrt(sqrt(5)/2 + sympy.S(3)/2))/2 - atanh(sqrt(2)*x**2/sqrt(sqrt(5) + 3))/sqrt(10*sqrt(5) + 30)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_387():
    f = 1/(x*(x**8 - 3*x**4 + 1))
    F = log(x) - (sympy.S(1)/8 + 3*sqrt(5)/40)*log(-2*x**4 - sqrt(5) + 3) - (sympy.S(1)/8 - 3*sqrt(5)/40)*log(-2*x**4 + sqrt(5) + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_388():
    f = 1/(x**3*(x**8 - 3*x**4 + 1))
    F = sqrt(10)*(sqrt(5) + 3)**(sympy.S(3)/2)*atanh(x**2*sqrt(sqrt(5)/2 + sympy.S(3)/2))/40 - sqrt(sympy.S(9)/5 - 4*sqrt(5)/5)*atanh(sqrt(2)*x**2/sqrt(sqrt(5) + 3))/2 - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_389():
    f = 1/(x**5*(x**8 - 3*x**4 + 1))
    F = 3*log(x) - (sympy.S(3)/8 + 7*sqrt(5)/40)*log(-2*x**4 - sqrt(5) + 3) - (sympy.S(3)/8 - 7*sqrt(5)/40)*log(-2*x**4 + sqrt(5) + 3) - 1/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_390():
    f = 1/(x**7*(x**8 - 3*x**4 + 1))
    F = sqrt(11*sqrt(5)/2 + sympy.S(123)/10)*atanh(x**2*sqrt(sqrt(5)/2 + sympy.S(3)/2))/2 - sqrt(sympy.S(123)/10 - 11*sqrt(5)/2)*atanh(sqrt(2)*x**2/sqrt(sqrt(5) + 3))/2 - 3/(2*x**2) - 1/(6*x**6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_391():
    f = x**8/(x**8 - 3*x**4 + 1)
    F = x + sqrt(5)*(984 - 440*sqrt(5))**(sympy.S(1)/4)*atan(x*(sqrt(5)/2 + sympy.S(3)/2)**(sympy.S(1)/4))/20 - sqrt(5)*(55*sqrt(5)/2 + sympy.S(123)/2)**(sympy.S(1)/4)*atan(2**(sympy.S(1)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4))/10 + sqrt(5)*(984 - 440*sqrt(5))**(sympy.S(1)/4)*atanh(x*(sqrt(5)/2 + sympy.S(3)/2)**(sympy.S(1)/4))/20 - sqrt(5)*(55*sqrt(5)/2 + sympy.S(123)/2)**(sympy.S(1)/4)*atanh(2**(sympy.S(1)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4))/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_392():
    f = x**6/(x**8 - 3*x**4 + 1)
    F = -sqrt(5)*(144 - 64*sqrt(5))**(sympy.S(1)/4)*atan(x*(sqrt(5)/2 + sympy.S(3)/2)**(sympy.S(1)/4))/20 + 2**(sympy.S(1)/4)*sqrt(5)*(sqrt(5) + 3)**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4))/20 + sqrt(5)*(144 - 64*sqrt(5))**(sympy.S(1)/4)*atanh(x*(sqrt(5)/2 + sympy.S(3)/2)**(sympy.S(1)/4))/20 - 2**(sympy.S(1)/4)*sqrt(5)*(sqrt(5) + 3)**(sympy.S(3)/4)*atanh(2**(sympy.S(1)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4))/20
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_393():
    f = x**4/(x**8 - 3*x**4 + 1)
    F = sqrt(5)*(sympy.S(3)/2 - sqrt(5)/2)**(sympy.S(1)/4)*atan(x*(sqrt(5)/2 + sympy.S(3)/2)**(sympy.S(1)/4))/10 - sqrt(5)*(sqrt(5)/2 + sympy.S(3)/2)**(sympy.S(1)/4)*atan(2**(sympy.S(1)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4))/10 + sqrt(5)*(sympy.S(3)/2 - sqrt(5)/2)**(sympy.S(1)/4)*atanh(x*(sqrt(5)/2 + sympy.S(3)/2)**(sympy.S(1)/4))/10 - sqrt(5)*(sqrt(5)/2 + sympy.S(3)/2)**(sympy.S(1)/4)*atanh(2**(sympy.S(1)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4))/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_394():
    f = 1/(x**8 - 3*x**4 + 1)
    F = sqrt(5)*(4*sqrt(5) + 9)**(sympy.S(1)/4)*atan(x*(sqrt(5)/2 + sympy.S(3)/2)**(sympy.S(1)/4))/10 - 2**(sympy.S(3)/4)*sqrt(5)*atan(2**(sympy.S(1)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4))/(10*(sqrt(5) + 3)**(sympy.S(3)/4)) + sqrt(5)*(4*sqrt(5) + 9)**(sympy.S(1)/4)*atanh(x*(sqrt(5)/2 + sympy.S(3)/2)**(sympy.S(1)/4))/10 - 2**(sympy.S(3)/4)*sqrt(5)*atanh(2**(sympy.S(1)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4))/(10*(sqrt(5) + 3)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_395():
    f = 1/(x**2*(x**8 - 3*x**4 + 1))
    F = -2**(sympy.S(3)/4)*sqrt(5)*(sqrt(5) + 3)**(sympy.S(5)/4)*atan(x*(sqrt(5)/2 + sympy.S(3)/2)**(sympy.S(1)/4))/40 + sqrt(5)*(984 - 440*sqrt(5))**(sympy.S(1)/4)*atan(2**(sympy.S(1)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4))/20 + 2**(sympy.S(3)/4)*sqrt(5)*(sqrt(5) + 3)**(sympy.S(5)/4)*atanh(x*(sqrt(5)/2 + sympy.S(3)/2)**(sympy.S(1)/4))/40 - sqrt(5)*(984 - 440*sqrt(5))**(sympy.S(1)/4)*atanh(2**(sympy.S(1)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4))/20 - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_396():
    f = 1/(x**4*(x**8 - 3*x**4 + 1))
    F = sqrt(5)*(377*sqrt(5)/2 + sympy.S(843)/2)**(sympy.S(1)/4)*atan(x*(sqrt(5)/2 + sympy.S(3)/2)**(sympy.S(1)/4))/10 - sqrt(5)*(sympy.S(843)/2 - 377*sqrt(5)/2)**(sympy.S(1)/4)*atan(2**(sympy.S(1)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4))/10 + sqrt(5)*(377*sqrt(5)/2 + sympy.S(843)/2)**(sympy.S(1)/4)*atanh(x*(sqrt(5)/2 + sympy.S(3)/2)**(sympy.S(1)/4))/10 - sqrt(5)*(sympy.S(843)/2 - 377*sqrt(5)/2)**(sympy.S(1)/4)*atanh(2**(sympy.S(1)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4))/10 - 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_397():
    f = 1/(x**6*(x**8 - 3*x**4 + 1))
    F = -sqrt(5)*(1292*sqrt(5) + 2889)**(sympy.S(1)/4)*atan(x*(sqrt(5)/2 + sympy.S(3)/2)**(sympy.S(1)/4))/10 + sqrt(5)*(2889 - 1292*sqrt(5))**(sympy.S(1)/4)*atan(2**(sympy.S(1)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4))/10 + sqrt(5)*(1292*sqrt(5) + 2889)**(sympy.S(1)/4)*atanh(x*(sqrt(5)/2 + sympy.S(3)/2)**(sympy.S(1)/4))/10 - sqrt(5)*(2889 - 1292*sqrt(5))**(sympy.S(1)/4)*atanh(2**(sympy.S(1)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4))/10 - 3/x - 1/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_398():
    f = 1/(x**8*(x**8 - 3*x**4 + 1))
    F = sqrt(5)*(17711*sqrt(5)/2 + sympy.S(39603)/2)**(sympy.S(1)/4)*atan(x*(sqrt(5)/2 + sympy.S(3)/2)**(sympy.S(1)/4))/10 - sqrt(5)*(sympy.S(39603)/2 - 17711*sqrt(5)/2)**(sympy.S(1)/4)*atan(2**(sympy.S(1)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4))/10 + sqrt(5)*(17711*sqrt(5)/2 + sympy.S(39603)/2)**(sympy.S(1)/4)*atanh(x*(sqrt(5)/2 + sympy.S(3)/2)**(sympy.S(1)/4))/10 - sqrt(5)*(sympy.S(39603)/2 - 17711*sqrt(5)/2)**(sympy.S(1)/4)*atanh(2**(sympy.S(1)/4)*x/(sqrt(5) + 3)**(sympy.S(1)/4))/10 - 1/x**3 - 1/(7*x**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_399():
    f = x**3/(x**8 + 3*x**4 + 2)
    F = log(x**4 + 1)/4 - log(x**4 + 2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_400():
    f = x**11/(x**8 + 3*x**4 + 2)
    F = x**4/4 + log(x**4 + 1)/4 - log(x**4 + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_401():
    f = x**9/(x**10 + x**5 + 2)
    F = log(x**10 + x**5 + 2)/10 - sqrt(7)*atan(sqrt(7)*(2*x**5 + 1)/7)/35
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_402():
    f = x**4/(x**10 + x**5 + 2)
    F = 2*sqrt(7)*atan(sqrt(7)*(2*x**5 + 1)/7)/35
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_403():
    f = 1/(x*(x**10 + x**5 + 1))
    F = log(x) - log(x**10 + x**5 + 1)/10 - sqrt(3)*atan(sqrt(3)*(2*x**5 + 1)/3)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_404():
    f = 1/(x**6*(x**10 + x**5 + 1))
    F = -log(x) + log(x**10 + x**5 + 1)/10 - sqrt(3)*atan(sqrt(3)*(2*x**5 + 1)/3)/15 - 1/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_405():
    f = 1/(x**11 + x**6 + x)
    F = log(x) - log(x**10 + x**5 + 1)/10 - sqrt(3)*atan(sqrt(3)*(2*x**5 + 1)/3)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_406():
    f = x**3/(a/x**2 + b/x + c)
    F = -b*x**3/(3*c**2) - b*x*(-2*a*c + b**2)/c**4 + b*(5*a**2*c**2 - 5*a*b**2*c + b**4)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**5*sqrt(-4*a*c + b**2)) + x**4/(4*c) + x**2*(-a*c + b**2)/(2*c**3) + (a**2*c**2 - 3*a*b**2*c + b**4)*log(a + b*x + c*x**2)/(2*c**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_407():
    f = x**2/(a/x**2 + b/x + c)
    F = -b*x**2/(2*c**2) - b*(-2*a*c + b**2)*log(a + b*x + c*x**2)/(2*c**4) + x**3/(3*c) + x*(-a*c + b**2)/c**3 - (2*a**2*c**2 - 4*a*b**2*c + b**4)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**4*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_408():
    f = x/(a/x**2 + b/x + c)
    F = -b*x/c**2 + b*(-3*a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**3*sqrt(-4*a*c + b**2)) + x**2/(2*c) + (-a*c + b**2)*log(a + b*x + c*x**2)/(2*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_409():
    f = 1/(a/x**2 + b/x + c)
    F = -b*log(a + b*x + c*x**2)/(2*c**2) + x/c - (-2*a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_410():
    f = 1/(x*(a/x**2 + b/x + c))
    F = b*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c*sqrt(-4*a*c + b**2)) + log(a + b*x + c*x**2)/(2*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_411():
    f = 1/(x**2*(a/x**2 + b/x + c))
    F = 2*atanh((2*a/x + b)/sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_412():
    f = 1/(x**3*(a/x**2 + b/x + c))
    F = b*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(a*sqrt(-4*a*c + b**2)) + log(x)/a - log(a + b*x + c*x**2)/(2*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_413():
    f = 1/(x**4*(a/x**2 + b/x + c))
    F = -1/(a*x) - b*log(x)/a**2 + b*log(a + b*x + c*x**2)/(2*a**2) - (-2*a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(a**2*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_414():
    f = 1/(x**5*(a/x**2 + b/x + c))
    F = -1/(2*a*x**2) + b/(a**2*x) + b*(-3*a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(a**3*sqrt(-4*a*c + b**2)) + (-a*c + b**2)*log(x)/a**3 - (-a*c + b**2)*log(a + b*x + c*x**2)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_415():
    f = 1/(x**6*(a/x**2 + b/x + c))
    F = -1/(3*a*x**3) + b/(2*a**2*x**2) - (-a*c + b**2)/(a**3*x) - b*(-2*a*c + b**2)*log(x)/a**4 + b*(-2*a*c + b**2)*log(a + b*x + c*x**2)/(2*a**4) - (2*a**2*c**2 - 4*a*b**2*c + b**4)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(a**4*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_416():
    f = x/(a/x**2 + b/x + c)**2
    F = -b*x**3/(c*(-4*a*c + b**2)) - b*x*(-11*a*c + 3*b**2)/(c**3*(-4*a*c + b**2)) + b*(30*a**2*c**2 - 20*a*b**2*c + 3*b**4)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**4*(-4*a*c + b**2)**(sympy.S(3)/2)) + x**4*(2*a + b*x)/((-4*a*c + b**2)*(a + b*x + c*x**2)) + x**2*(-8*a*c + 3*b**2)/(2*c**2*(-4*a*c + b**2)) + (-2*a*c + 3*b**2)*log(a + b*x + c*x**2)/(2*c**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_417():
    f = (a/x**2 + b/x + c)**(-2)
    F = -b*x**2/(c*(-4*a*c + b**2)) - b*log(a + b*x + c*x**2)/c**3 + x**3*(2*a + b*x)/((-4*a*c + b**2)*(a + b*x + c*x**2)) + x*(-6*a*c + 2*b**2)/(c**2*(-4*a*c + b**2)) - (12*a**2*c**2 - 12*a*b**2*c + 2*b**4)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**3*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_418():
    f = 1/(x*(a/x**2 + b/x + c)**2)
    F = -b*x/(c*(-4*a*c + b**2)) + b*(-6*a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**2*(-4*a*c + b**2)**(sympy.S(3)/2)) + x**2*(2*a + b*x)/((-4*a*c + b**2)*(a + b*x + c*x**2)) + log(a + b*x + c*x**2)/(2*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_419():
    f = 1/(x**2*(a/x**2 + b/x + c)**2)
    F = -4*a*atanh((2*a/x + b)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) + (2*a/x + b)/((-4*a*c + b**2)*(a/x**2 + b/x + c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_420():
    f = 1/(x**3*(a/x**2 + b/x + c)**2)
    F = -2*b*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) + (2*a + b*x)/((-4*a*c + b**2)*(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_421():
    f = 1/(x**4*(a/x**2 + b/x + c)**2)
    F = 4*c*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(3)/2) - (b + 2*c*x)/((-4*a*c + b**2)*(a + b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_422():
    f = 1/(x**5*(a/x**2 + b/x + c)**2)
    F = (-2*a*c + b**2 + b*c*x)/(a*(-4*a*c + b**2)*(a + b*x + c*x**2)) + b*(-6*a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(a**2*(-4*a*c + b**2)**(sympy.S(3)/2)) + log(x)/a**2 - log(a + b*x + c*x**2)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_423():
    f = 1/(x**6*(a/x**2 + b/x + c)**2)
    F = (-2*a*c + b**2 + b*c*x)/(a*x*(-4*a*c + b**2)*(a + b*x + c*x**2)) + (6*a*c - 2*b**2)/(a**2*x*(-4*a*c + b**2)) - 2*b*log(x)/a**3 + b*log(a + b*x + c*x**2)/a**3 - (12*a**2*c**2 - 12*a*b**2*c + 2*b**4)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(a**3*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_424():
    f = 1/(x**7*(a/x**2 + b/x + c)**2)
    F = (-2*a*c + b**2 + b*c*x)/(a*x**2*(-4*a*c + b**2)*(a + b*x + c*x**2)) - (-8*a*c + 3*b**2)/(2*a**2*x**2*(-4*a*c + b**2)) + b*(-11*a*c + 3*b**2)/(a**3*x*(-4*a*c + b**2)) + b*(30*a**2*c**2 - 20*a*b**2*c + 3*b**4)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(a**4*(-4*a*c + b**2)**(sympy.S(3)/2)) + (-2*a*c + 3*b**2)*log(x)/a**4 - (-2*a*c + 3*b**2)*log(a + b*x + c*x**2)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_425():
    f = (a/x**2 + b/x + c)**(-3)
    F = -3*b*x**2*(-6*a*c + b**2)/(2*c**2*(-4*a*c + b**2)**2) - 3*b*log(a + b*x + c*x**2)/(2*c**4) + x**5*(2*a + b*x)/((-8*a*c + 2*b**2)*(a + b*x + c*x**2)**2) + x**3*(a*(-10*a*c + b**2) + b*x*(-7*a*c + b**2))/(c*(-4*a*c + b**2)**2*(a + b*x + c*x**2)) + x*(30*a**2*c**2 - 21*a*b**2*c + 3*b**4)/(c**3*(-4*a*c + b**2)**2) - (-60*a**3*c**3 + 90*a**2*b**2*c**2 - 30*a*b**4*c + 3*b**6)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**4*(-4*a*c + b**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_426():
    f = 1/(x*(a/x**2 + b/x + c)**3)
    F = -b*x*(-7*a*c + b**2)/(c**2*(-4*a*c + b**2)**2) + b*(30*a**2*c**2 - 10*a*b**2*c + b**4)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(c**3*(-4*a*c + b**2)**(sympy.S(5)/2)) + x**4*(2*a + b*x)/((-8*a*c + 2*b**2)*(a + b*x + c*x**2)**2) + x**2*(a*(-16*a*c + b**2) + b*x*(-10*a*c + b**2))/(2*c*(-4*a*c + b**2)**2*(a + b*x + c*x**2)) + log(a + b*x + c*x**2)/(2*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_427():
    f = 1/(x**2*(a/x**2 + b/x + c)**3)
    F = 12*a**2*atanh((2*a/x + b)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2) - 3*a*(2*a/x + b)/((-4*a*c + b**2)**2*(a/x**2 + b/x + c)) + (2*a/x + b)/((-8*a*c + 2*b**2)*(a/x**2 + b/x + c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_428():
    f = 1/(x**3*(a/x**2 + b/x + c)**3)
    F = 6*a*b*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2) + 3*b*x*(2*a + b*x)/(2*(-4*a*c + b**2)**2*(a + b*x + c*x**2)) - x**3*(b + 2*c*x)/((-8*a*c + 2*b**2)*(a + b*x + c*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_429():
    f = 1/(x**4*(a/x**2 + b/x + c)**3)
    F = x*(2*a + b*x)/((-8*a*c + 2*b**2)*(a + b*x + c*x**2)**2) + (3*a*b + x*(2*a*c + b**2))/((-4*a*c + b**2)**2*(a + b*x + c*x**2)) - (4*a*c + 2*b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_430():
    f = 1/(x**5*(a/x**2 + b/x + c)**3)
    F = 6*b*c*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2) - 3*b*(b + 2*c*x)/(2*(-4*a*c + b**2)**2*(a + b*x + c*x**2)) + (2*a + b*x)/((-8*a*c + 2*b**2)*(a + b*x + c*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_431():
    f = 1/(x**6*(a/x**2 + b/x + c)**3)
    F = -12*c**2*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(-4*a*c + b**2)**(sympy.S(5)/2) + 3*c*(b + 2*c*x)/((-4*a*c + b**2)**2*(a + b*x + c*x**2)) + (-b - 2*c*x)/((-8*a*c + 2*b**2)*(a + b*x + c*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_432():
    f = 1/(x**7*(a/x**2 + b/x + c)**3)
    F = (-2*a*c + b**2 + b*c*x)/(2*a*(-4*a*c + b**2)*(a + b*x + c*x**2)**2) + (16*a**2*c**2 - 15*a*b**2*c + 2*b**4 + 2*b*c*x*(-7*a*c + b**2))/(2*a**2*(-4*a*c + b**2)**2*(a + b*x + c*x**2)) + b*(30*a**2*c**2 - 10*a*b**2*c + b**4)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(a**3*(-4*a*c + b**2)**(sympy.S(5)/2)) + log(x)/a**3 - log(a + b*x + c*x**2)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_433():
    f = 1/(x**8*(a/x**2 + b/x + c)**3)
    F = (-2*a*c + b**2 + b*c*x)/(2*a*x*(-4*a*c + b**2)*(a + b*x + c*x**2)**2) + (20*a**2*c**2 - 20*a*b**2*c + 3*b**4 + 3*b*c*x*(-6*a*c + b**2))/(2*a**2*x*(-4*a*c + b**2)**2*(a + b*x + c*x**2)) - (-15*a*c + 3*b**2)*(-2*a*c + b**2)/(a**3*x*(-4*a*c + b**2)**2) - 3*b*log(x)/a**4 + 3*b*log(a + b*x + c*x**2)/(2*a**4) - (-60*a**3*c**3 + 90*a**2*b**2*c**2 - 30*a*b**4*c + 3*b**6)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(a**4*(-4*a*c + b**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_434():
    f = x**2/(15 + 13/x + 2/x**2)
    F = x**3/45 - 13*x**2/450 + 139*x/3375 - 16*log(3*x + 2)/567 + log(5*x + 1)/4375
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_435():
    f = x/(15 + 13/x + 2/x**2)
    F = x**2/30 - 13*x/225 + 8*log(3*x + 2)/189 - log(5*x + 1)/875
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_436():
    f = 1/(15 + 13/x + 2/x**2)
    F = x/15 - 4*log(3*x + 2)/63 + log(5*x + 1)/175
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_437():
    f = 1/(x*(15 + 13/x + 2/x**2))
    F = 2*log(3*x + 2)/21 - log(5*x + 1)/35
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_438():
    f = 1/(x**2*(15 + 13/x + 2/x**2))
    F = -log(3 + 2/x)/7 + log(5 + 1/x)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_439():
    f = 1/(x**3*(15 + 13/x + 2/x**2))
    F = log(x)/2 + 3*log(3*x + 2)/14 - 5*log(5*x + 1)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_440():
    f = 1/(x**4*(15 + 13/x + 2/x**2))
    F = -13*log(x)/4 - 9*log(3*x + 2)/28 + 25*log(5*x + 1)/7 - 1/(2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_441():
    f = 1/(x**5*(15 + 13/x + 2/x**2))
    F = 139*log(x)/8 + 27*log(3*x + 2)/56 - 125*log(5*x + 1)/7 + 13/(4*x) - 1/(4*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_442():
    f = 1/(x**6*(15 + 13/x + 2/x**2))
    F = -1417*log(x)/16 - 81*log(3*x + 2)/112 + 625*log(5*x + 1)/7 - 139/(8*x) + 13/(8*x**2) - 1/(6*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_443():
    f = (a + b/x + c/x**2)**(sympy.S(5)/2)
    F = 5*a**(sympy.S(3)/2)*b*atanh((2*a + b/x)/(2*sqrt(a)*sqrt(a + b/x + c/x**2)))/2 + x*(a + b/x + c/x**2)**(sympy.S(5)/2) - 5*(7*b + 6*c/x)*(a + b/x + c/x**2)**(sympy.S(3)/2)/24 - 5*(b*(44*a*c + b**2) + 2*c*(12*a*c + b**2)/x)*sqrt(a + b/x + c/x**2)/(64*c) + (-240*a**2*c**2 - 120*a*b**2*c + 5*b**4)*atanh((b + 2*c/x)/(2*sqrt(c)*sqrt(a + b/x + c/x**2)))/(128*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_444():
    f = (a + b/x + c/x**2)**(sympy.S(3)/2)
    F = 3*sqrt(a)*b*atanh((2*a + b/x)/(2*sqrt(a)*sqrt(a + b/x + c/x**2)))/2 + x*(a + b/x + c/x**2)**(sympy.S(3)/2) - 3*(3*b + 2*c/x)*sqrt(a + b/x + c/x**2)/4 - (12*a*c + 3*b**2)*atanh((b + 2*c/x)/(2*sqrt(c)*sqrt(a + b/x + c/x**2)))/(8*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_445():
    f = sqrt(a + b/x + c/x**2)
    F = -sqrt(c)*atanh((b + 2*c/x)/(2*sqrt(c)*sqrt(a + b/x + c/x**2))) + x*sqrt(a + b/x + c/x**2) + b*atanh((2*a + b/x)/(2*sqrt(a)*sqrt(a + b/x + c/x**2)))/(2*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_446():
    f = 1/sqrt(a + b/x + c/x**2)
    F = x*sqrt(a + b/x + c/x**2)/a - b*atanh((2*a + b/x)/(2*sqrt(a)*sqrt(a + b/x + c/x**2)))/(2*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_447():
    f = (a + b/x + c/x**2)**(sympy.S(-3)/2)
    F = -x*(-4*a*c + 2*b**2 + 2*b*c/x)/(a*(-4*a*c + b**2)*sqrt(a + b/x + c/x**2)) + x*(-8*a*c + 3*b**2)*sqrt(a + b/x + c/x**2)/(a**2*(-4*a*c + b**2)) - 3*b*atanh((2*a + b/x)/(2*sqrt(a)*sqrt(a + b/x + c/x**2)))/(2*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_448():
    f = (a + b/x + c/x**2)**(sympy.S(-5)/2)
    F = -x*(-4*a*c + 2*b**2 + 2*b*c/x)/(3*a*(-4*a*c + b**2)*(a + b/x + c/x**2)**(sympy.S(3)/2)) - x*(64*a**2*c**2 - 64*a*b**2*c + 10*b**4 + 2*b*c*(-28*a*c + 5*b**2)/x)/(3*a**2*(-4*a*c + b**2)**2*sqrt(a + b/x + c/x**2)) + x*sqrt(a + b/x + c/x**2)*(128*a**2*c**2 - 100*a*b**2*c + 15*b**4)/(3*a**3*(-4*a*c + b**2)**2) - 5*b*atanh((2*a + b/x)/(2*sqrt(a)*sqrt(a + b/x + c/x**2)))/(2*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_449():
    f = sqrt(a**2 + 2*a*b/x + b**2/x**2)
    F = a*x*sqrt(a**2 + 2*a*b/x + b**2/x**2)/(a + b/x) - b*sqrt(a**2 + 2*a*b/x + b**2/x**2)*log(1/x)/(a + b/x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_450():
    f = 1/(a/x**4 + b/x**2 + c)
    F = x/c - sqrt(2)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2))) - sqrt(2)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_451():
    f = 1/(a/x**6 + b/x**3 + c)
    F = x/c - 2**(sympy.S(2)/3)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(4)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(4)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*sqrt(3)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(4)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x + (b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(6*c**(sympy.S(4)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**2 - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(12*c**(sympy.S(4)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*sqrt(3)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x/(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(6*c**(sympy.S(4)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_452():
    f = 1/(a/x**8 + b/x**4 + c)
    F = x/c + 2**(sympy.S(3)/4)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(5)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) + 2**(sympy.S(3)/4)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(5)/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) + 2**(sympy.S(3)/4)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(5)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) + 2**(sympy.S(3)/4)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(4*c**(sympy.S(5)/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_453():
    f = sqrt(a + b*sqrt(x) + c*x)/x
    F = -2*sqrt(a)*atanh((2*a + b*sqrt(x))/(2*sqrt(a)*sqrt(a + b*sqrt(x) + c*x))) + b*atanh((b + 2*c*sqrt(x))/(2*sqrt(c)*sqrt(a + b*sqrt(x) + c*x)))/sqrt(c) + 2*sqrt(a + b*sqrt(x) + c*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_454():
    f = (b**2/(4*c) + b*sqrt(x) + c*x)**2
    F = -b*(b + 2*c*sqrt(x))**5/(160*c**4) + (b + 2*c*sqrt(x))**6/(192*c**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_455():
    f = 1/sqrt(a**2 + 2*a*b*sqrt(x) + b**2*x)
    F = -2*a*(a + b*sqrt(x))*log(a + b*sqrt(x))/(b**2*sqrt(a**2 + 2*a*b*sqrt(x) + b**2*x)) + 2*sqrt(a**2 + 2*a*b*sqrt(x) + b**2*x)/b**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_456():
    f = (a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**(sympy.S(7)/2)
    F = 3*a**2*(a + b*x**(sympy.S(1)/3))**7*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))/(8*b**3) - 2*a*(a + b*x**(sympy.S(1)/3))**8*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))/(3*b**3) + 3*(a + b*x**(sympy.S(1)/3))**9*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))/(10*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_457():
    f = (a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**(sympy.S(5)/2)
    F = a**2*(a + b*x**(sympy.S(1)/3))**5*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))/(2*b**3) - 6*a*(a + b*x**(sympy.S(1)/3))**6*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))/(7*b**3) + 3*(a + b*x**(sympy.S(1)/3))**7*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))/(8*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_458():
    f = (a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**(sympy.S(3)/2)
    F = 3*a**2*(a + b*x**(sympy.S(1)/3))**3*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))/(4*b**3) - 6*a*(a + b*x**(sympy.S(1)/3))**4*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))/(5*b**3) + (a + b*x**(sympy.S(1)/3))**5*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_459():
    f = sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))
    F = a*x*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))/(a + b*x**(sympy.S(1)/3)) + 3*b*x**(sympy.S(4)/3)*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))/(4*a + 4*b*x**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_460():
    f = 1/sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))
    F = 3*a**2*(a + b*x**(sympy.S(1)/3))*log(a + b*x**(sympy.S(1)/3))/(b**3*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))) - 3*a*x**(sympy.S(1)/3)*(a + b*x**(sympy.S(1)/3))/(b**2*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))) + x**(sympy.S(2)/3)*(3*a + 3*b*x**(sympy.S(1)/3))/(2*b*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_461():
    f = (a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**(sympy.S(-3)/2)
    F = -3*a**2/(2*b**3*(a + b*x**(sympy.S(1)/3))*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))) + 6*a/(b**3*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))) + (3*a + 3*b*x**(sympy.S(1)/3))*log(a + b*x**(sympy.S(1)/3))/(b**3*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_462():
    f = (a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**(sympy.S(-5)/2)
    F = -3*a**2/(4*b**3*(a + b*x**(sympy.S(1)/3))**3*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))) + 2*a/(b**3*(a + b*x**(sympy.S(1)/3))**2*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))) - 3/(2*b**3*(a + b*x**(sympy.S(1)/3))*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_463():
    f = (a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**(sympy.S(-7)/2)
    F = -a**2/(2*b**3*(a + b*x**(sympy.S(1)/3))**5*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))) + 6*a/(5*b**3*(a + b*x**(sympy.S(1)/3))**4*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))) - 3/(4*b**3*(a + b*x**(sympy.S(1)/3))**3*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_464():
    f = (a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**(sympy.S(-9)/2)
    F = -3*a**2/(8*b**3*(a + b*x**(sympy.S(1)/3))**7*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))) + 6*a/(7*b**3*(a + b*x**(sympy.S(1)/3))**6*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))) - 1/(2*b**3*(a + b*x**(sympy.S(1)/3))**5*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_465():
    f = (a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**(sympy.S(-11)/2)
    F = -3*a**2/(10*b**3*(a + b*x**(sympy.S(1)/3))**9*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))) + 2*a/(3*b**3*(a + b*x**(sympy.S(1)/3))**8*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))) - 3/(8*b**3*(a + b*x**(sympy.S(1)/3))**7*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_466():
    f = (d*x)**m*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p
    F = x*(d*x)**m*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p*hyper((-2*p, 3*m + 3), (3*m + 4,), -b*x**(sympy.S(1)/3)/a)/((1 + b*x**(sympy.S(1)/3)/a)**(2*p)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_467():
    f = x**2*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p
    F = 3*a**9*(1 + b*x**(sympy.S(1)/3)/a)**9*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(b**9*(2*p + 9)) - 12*a**9*(1 + b*x**(sympy.S(1)/3)/a)**8*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(b**9*(p + 4)) + 84*a**9*(1 + b*x**(sympy.S(1)/3)/a)**7*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(b**9*(2*p + 7)) - 84*a**9*(1 + b*x**(sympy.S(1)/3)/a)**6*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(b**9*(p + 3)) + 210*a**9*(1 + b*x**(sympy.S(1)/3)/a)**5*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(b**9*(2*p + 5)) - 84*a**9*(1 + b*x**(sympy.S(1)/3)/a)**4*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(b**9*(p + 2)) + 84*a**9*(1 + b*x**(sympy.S(1)/3)/a)**3*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(b**9*(2*p + 3)) - 12*a**9*(1 + b*x**(sympy.S(1)/3)/a)**2*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(b**9*(p + 1)) + 3*a**9*(1 + b*x**(sympy.S(1)/3)/a)*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(b**9*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_468():
    f = x*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p
    F = 3*a**6*(1 + b*x**(sympy.S(1)/3)/a)**6*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(2*b**6*(p + 3)) - 15*a**6*(1 + b*x**(sympy.S(1)/3)/a)**5*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(b**6*(2*p + 5)) + 15*a**6*(1 + b*x**(sympy.S(1)/3)/a)**4*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(b**6*(p + 2)) - 30*a**6*(1 + b*x**(sympy.S(1)/3)/a)**3*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(b**6*(2*p + 3)) + 15*a**6*(1 + b*x**(sympy.S(1)/3)/a)**2*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(2*b**6*(p + 1)) - 3*a**6*(1 + b*x**(sympy.S(1)/3)/a)*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(b**6*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_469():
    f = (a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p
    F = 3*a**2*(a + b*x**(sympy.S(1)/3))*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(b**3*(2*p + 1)) - 3*a*(a + b*x**(sympy.S(1)/3))**2*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(b**3*(p + 1)) + 3*(a + b*x**(sympy.S(1)/3))**3*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(b**3*(2*p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_470():
    f = (a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/x
    F = -(3 + 3*b*x**(sympy.S(1)/3)/a)*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p*hyper((1, 2*p + 1), (2*p + 2,), 1 + b*x**(sympy.S(1)/3)/a)/(2*p + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_471():
    f = (a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/x**2
    F = 3*b**3*(1 + b*x**(sympy.S(1)/3)/a)*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p*hyper((4, 2*p + 1), (2*p + 2,), 1 + b*x**(sympy.S(1)/3)/a)/(a**3*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_472():
    f = (a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/x**2 - 2*b**3*p*(1 - 2*p)*(1 - p)*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(3*a**3*x)
    F = -(a + b*x**(sympy.S(1)/3))*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(a*x) + b*(1 - p)*(a + b*x**(sympy.S(1)/3))*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(a**2*x**(sympy.S(2)/3)) - b**2*(1 - 2*p)*(1 - p)*(a + b*x**(sympy.S(1)/3))*(a**2 + 2*a*b*x**(sympy.S(1)/3) + b**2*x**(sympy.S(2)/3))**p/(a**3*x**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_473():
    f = (a**2 + 2*a*b*x**(sympy.S(1)/4) + b**2*sqrt(x))**(sympy.S(-3)/2)
    F = 2*a**3/(b**4*(a + b*x**(sympy.S(1)/4))*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/4) + b**2*sqrt(x))) - 12*a**2/(b**4*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/4) + b**2*sqrt(x))) - 12*a*(a + b*x**(sympy.S(1)/4))*log(a + b*x**(sympy.S(1)/4))/(b**4*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/4) + b**2*sqrt(x))) + x**(sympy.S(1)/4)*(4*a + 4*b*x**(sympy.S(1)/4))/(b**3*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/4) + b**2*sqrt(x)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_474():
    f = (a**2 + 2*a*b*x**(sympy.S(1)/6) + b**2*x**(sympy.S(1)/3))**(sympy.S(-5)/2)
    F = 3*a**5/(2*b**6*(a + b*x**(sympy.S(1)/6))**3*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/6) + b**2*x**(sympy.S(1)/3))) - 10*a**4/(b**6*(a + b*x**(sympy.S(1)/6))**2*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/6) + b**2*x**(sympy.S(1)/3))) + 30*a**3/(b**6*(a + b*x**(sympy.S(1)/6))*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/6) + b**2*x**(sympy.S(1)/3))) - 60*a**2/(b**6*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/6) + b**2*x**(sympy.S(1)/3))) - 30*a*(a + b*x**(sympy.S(1)/6))*log(a + b*x**(sympy.S(1)/6))/(b**6*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/6) + b**2*x**(sympy.S(1)/3))) + x**(sympy.S(1)/6)*(6*a + 6*b*x**(sympy.S(1)/6))/(b**5*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/6) + b**2*x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_475():
    f = (a**2 + 2*a*b/sqrt(x) + b**2/x)**(sympy.S(3)/2)
    F = a**3*x*sqrt(a**2 + 2*a*b/sqrt(x) + b**2/x)/(a + b/sqrt(x)) + 6*a**2*b*sqrt(x)*sqrt(a**2 + 2*a*b/sqrt(x) + b**2/x)/(a + b/sqrt(x)) + 6*a*b**2*sqrt(a**2 + 2*a*b/sqrt(x) + b**2/x)*log(sqrt(x))/(a + b/sqrt(x)) - 2*b**3*sqrt(a**2 + 2*a*b/sqrt(x) + b**2/x)/(sqrt(x)*(a + b/sqrt(x)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_476():
    f = (a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))**(sympy.S(7)/2)
    F = a**7*x*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))/(a + b/x**(sympy.S(1)/3)) + 21*a**6*b*x**(sympy.S(2)/3)*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))/(2*a + 2*b/x**(sympy.S(1)/3)) + 63*a**5*b**2*x**(sympy.S(1)/3)*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))/(a + b/x**(sympy.S(1)/3)) + 105*a**4*b**3*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))*log(x**(sympy.S(1)/3))/(a + b/x**(sympy.S(1)/3)) - 105*a**3*b**4*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))/(x**(sympy.S(1)/3)*(a + b/x**(sympy.S(1)/3))) - 63*a**2*b**5*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))/(x**(sympy.S(2)/3)*(2*a + 2*b/x**(sympy.S(1)/3))) - 7*a*b**6*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))/(x*(a + b/x**(sympy.S(1)/3))) - 3*b**7*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))/(x**(sympy.S(4)/3)*(4*a + 4*b/x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_477():
    f = (a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))**(sympy.S(5)/2)
    F = a**5*x*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))/(a + b/x**(sympy.S(1)/3)) + 15*a**4*b*x**(sympy.S(2)/3)*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))/(2*a + 2*b/x**(sympy.S(1)/3)) + 30*a**3*b**2*x**(sympy.S(1)/3)*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))/(a + b/x**(sympy.S(1)/3)) + 30*a**2*b**3*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))*log(x**(sympy.S(1)/3))/(a + b/x**(sympy.S(1)/3)) - 15*a*b**4*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))/(x**(sympy.S(1)/3)*(a + b/x**(sympy.S(1)/3))) - 3*b**5*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))/(x**(sympy.S(2)/3)*(2*a + 2*b/x**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_478():
    f = (a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))**(sympy.S(3)/2)
    F = a**3*x*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))/(a + b/x**(sympy.S(1)/3)) + 9*a**2*b*x**(sympy.S(2)/3)*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))/(2*a + 2*b/x**(sympy.S(1)/3)) + 9*a*b**2*x**(sympy.S(1)/3)*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))/(a + b/x**(sympy.S(1)/3)) + 3*b**3*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))*log(x**(sympy.S(1)/3))/(a + b/x**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_479():
    f = sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))
    F = a*x*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))/(a + b/x**(sympy.S(1)/3)) + 3*b*x**(sympy.S(2)/3)*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))/(2*a + 2*b/x**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_480():
    f = 1/sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))
    F = x*(a + b/x**(sympy.S(1)/3))/(a*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))) - 3*b*x**(sympy.S(2)/3)*(a + b/x**(sympy.S(1)/3))/(2*a**2*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))) + 3*b**2*x**(sympy.S(1)/3)*(a + b/x**(sympy.S(1)/3))/(a**3*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))) - 3*b**3*(a + b/x**(sympy.S(1)/3))*log(a*x**(sympy.S(1)/3) + b)/(a**4*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_481():
    f = (a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))**(sympy.S(-3)/2)
    F = x*(a + b/x**(sympy.S(1)/3))/(a**3*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))) - 9*b*x**(sympy.S(2)/3)*(a + b/x**(sympy.S(1)/3))/(2*a**4*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))) + 18*b**2*x**(sympy.S(1)/3)*(a + b/x**(sympy.S(1)/3))/(a**5*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))) + 3*b**5*(a + b/x**(sympy.S(1)/3))/(2*a**6*(a*x**(sympy.S(1)/3) + b)**2*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))) - 15*b**4*(a + b/x**(sympy.S(1)/3))/(a**6*(a*x**(sympy.S(1)/3) + b)*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))) - 30*b**3*(a + b/x**(sympy.S(1)/3))*log(a*x**(sympy.S(1)/3) + b)/(a**6*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_482():
    f = (a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))**(sympy.S(-5)/2)
    F = x*(a + b/x**(sympy.S(1)/3))/(a**5*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))) - 15*b*x**(sympy.S(2)/3)*(a + b/x**(sympy.S(1)/3))/(2*a**6*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))) + 45*b**2*x**(sympy.S(1)/3)*(a + b/x**(sympy.S(1)/3))/(a**7*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))) + 3*b**7*(a + b/x**(sympy.S(1)/3))/(4*a**8*(a*x**(sympy.S(1)/3) + b)**4*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))) - 7*b**6*(a + b/x**(sympy.S(1)/3))/(a**8*(a*x**(sympy.S(1)/3) + b)**3*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))) + 63*b**5*(a + b/x**(sympy.S(1)/3))/(2*a**8*(a*x**(sympy.S(1)/3) + b)**2*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))) - 105*b**4*(a + b/x**(sympy.S(1)/3))/(a**8*(a*x**(sympy.S(1)/3) + b)*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3))) - 105*b**3*(a + b/x**(sympy.S(1)/3))*log(a*x**(sympy.S(1)/3) + b)/(a**8*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/3) + b**2/x**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_483():
    f = (a**2 + 2*a*b/x**(sympy.S(1)/4) + b**2/sqrt(x))**(sympy.S(5)/2)
    F = a**5*x*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/4) + b**2/sqrt(x))/(a + b/x**(sympy.S(1)/4)) + 20*a**4*b*x**(sympy.S(3)/4)*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/4) + b**2/sqrt(x))/(3*a + 3*b/x**(sympy.S(1)/4)) + 20*a**3*b**2*sqrt(x)*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/4) + b**2/sqrt(x))/(a + b/x**(sympy.S(1)/4)) + 40*a**2*b**3*x**(sympy.S(1)/4)*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/4) + b**2/sqrt(x))/(a + b/x**(sympy.S(1)/4)) + 20*a*b**4*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/4) + b**2/sqrt(x))*log(x**(sympy.S(1)/4))/(a + b/x**(sympy.S(1)/4)) - 4*b**5*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/4) + b**2/sqrt(x))/(x**(sympy.S(1)/4)*(a + b/x**(sympy.S(1)/4)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_484():
    f = (a**2 + 2*a*b/x**(sympy.S(1)/5) + b**2/x**(sympy.S(2)/5))**(sympy.S(5)/2)
    F = a**5*x*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/5) + b**2/x**(sympy.S(2)/5))/(a + b/x**(sympy.S(1)/5)) + 25*a**4*b*x**(sympy.S(4)/5)*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/5) + b**2/x**(sympy.S(2)/5))/(4*a + 4*b/x**(sympy.S(1)/5)) + 50*a**3*b**2*x**(sympy.S(3)/5)*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/5) + b**2/x**(sympy.S(2)/5))/(3*a + 3*b/x**(sympy.S(1)/5)) + 25*a**2*b**3*x**(sympy.S(2)/5)*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/5) + b**2/x**(sympy.S(2)/5))/(a + b/x**(sympy.S(1)/5)) + 25*a*b**4*x**(sympy.S(1)/5)*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/5) + b**2/x**(sympy.S(2)/5))/(a + b/x**(sympy.S(1)/5)) + 5*b**5*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/5) + b**2/x**(sympy.S(2)/5))*log(x**(sympy.S(1)/5))/(a + b/x**(sympy.S(1)/5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_485():
    f = (a**2 + 2*a*b*x**(sympy.S(1)/5) + b**2*x**(sympy.S(2)/5))**(sympy.S(-5)/2)
    F = -5*a**4/(4*b**5*(a + b*x**(sympy.S(1)/5))**3*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/5) + b**2*x**(sympy.S(2)/5))) + 20*a**3/(3*b**5*(a + b*x**(sympy.S(1)/5))**2*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/5) + b**2*x**(sympy.S(2)/5))) - 15*a**2/(b**5*(a + b*x**(sympy.S(1)/5))*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/5) + b**2*x**(sympy.S(2)/5))) + 20*a/(b**5*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/5) + b**2*x**(sympy.S(2)/5))) + (5*a + 5*b*x**(sympy.S(1)/5))*log(a + b*x**(sympy.S(1)/5))/(b**5*sqrt(a**2 + 2*a*b*x**(sympy.S(1)/5) + b**2*x**(sympy.S(2)/5)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_486():
    f = (a**2 + 2*a*b/x**(sympy.S(1)/6) + b**2/x**(sympy.S(1)/3))**(sympy.S(7)/2)
    F = a**7*x*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/6) + b**2/x**(sympy.S(1)/3))/(a + b/x**(sympy.S(1)/6)) + 42*a**6*b*x**(sympy.S(5)/6)*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/6) + b**2/x**(sympy.S(1)/3))/(5*a + 5*b/x**(sympy.S(1)/6)) + 63*a**5*b**2*x**(sympy.S(2)/3)*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/6) + b**2/x**(sympy.S(1)/3))/(2*a + 2*b/x**(sympy.S(1)/6)) + 70*a**4*b**3*sqrt(x)*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/6) + b**2/x**(sympy.S(1)/3))/(a + b/x**(sympy.S(1)/6)) + 105*a**3*b**4*x**(sympy.S(1)/3)*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/6) + b**2/x**(sympy.S(1)/3))/(a + b/x**(sympy.S(1)/6)) + 126*a**2*b**5*x**(sympy.S(1)/6)*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/6) + b**2/x**(sympy.S(1)/3))/(a + b/x**(sympy.S(1)/6)) + 42*a*b**6*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/6) + b**2/x**(sympy.S(1)/3))*log(x**(sympy.S(1)/6))/(a + b/x**(sympy.S(1)/6)) - 6*b**7*sqrt(a**2 + 2*a*b/x**(sympy.S(1)/6) + b**2/x**(sympy.S(1)/3))/(x**(sympy.S(1)/6)*(a + b/x**(sympy.S(1)/6)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_487():
    f = x**(4*n - 1)/(b*x**n + c*x**(2*n))
    F = b**2*log(b + c*x**n)/(c**3*n) - b*x**n/(c**2*n) + x**(2*n)/(2*c*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_488():
    f = x**(3*n - 1)/(b*x**n + c*x**(2*n))
    F = -b*log(b + c*x**n)/(c**2*n) + x**n/(c*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_489():
    f = x**(2*n - 1)/(b*x**n + c*x**(2*n))
    F = log(b + c*x**n)/(c*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_490():
    f = x**(n - 1)/(b*x**n + c*x**(2*n))
    F = log(x)/b - log(b + c*x**n)/(b*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_491():
    f = x**(-n - 1)/(b*x**n + c*x**(2*n))
    F = -1/(2*b*n*x**(2*n)) + c/(b**2*n*x**n) + c**2*log(x)/b**3 - c**2*log(b + c*x**n)/(b**3*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_492():
    f = x**(-2*n - 1)/(b*x**n + c*x**(2*n))
    F = -1/(3*b*n*x**(3*n)) + c/(2*b**2*n*x**(2*n)) - c**2/(b**3*n*x**n) - c**3*log(x)/b**4 + c**3*log(b + c*x**n)/(b**4*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_493():
    f = x**(-3*n - 1)/(b*x**n + c*x**(2*n))
    F = -1/(4*b*n*x**(4*n)) + c/(3*b**2*n*x**(3*n)) - c**2/(2*b**3*n*x**(2*n)) + c**3/(b**4*n*x**n) + c**4*log(x)/b**5 - c**4*log(b + c*x**n)/(b**5*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_494():
    f = x**(n/4 - 1)/(b*x**n + c*x**(2*n))
    F = -4/(3*b*n*x**(3*n/4)) + sqrt(2)*c**(sympy.S(3)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x**(n/4) + sqrt(b) + sqrt(c)*x**(n/2))/(2*b**(sympy.S(7)/4)*n) - sqrt(2)*c**(sympy.S(3)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x**(n/4) + sqrt(b) + sqrt(c)*x**(n/2))/(2*b**(sympy.S(7)/4)*n) + sqrt(2)*c**(sympy.S(3)/4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x**(n/4)/b**(sympy.S(1)/4))/(b**(sympy.S(7)/4)*n) - sqrt(2)*c**(sympy.S(3)/4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x**(n/4)/b**(sympy.S(1)/4))/(b**(sympy.S(7)/4)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_495():
    f = x**(n/3 - 1)/(b*x**n + c*x**(2*n))
    F = -3/(2*b*n*x**(2*n/3)) - c**(sympy.S(2)/3)*log(b**(sympy.S(1)/3) + c**(sympy.S(1)/3)*x**(n/3))/(b**(sympy.S(5)/3)*n) + c**(sympy.S(2)/3)*log(b**(sympy.S(2)/3) - b**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x**(n/3) + c**(sympy.S(2)/3)*x**(2*n/3))/(2*b**(sympy.S(5)/3)*n) + sqrt(3)*c**(sympy.S(2)/3)*atan(sqrt(3)*(b**(sympy.S(1)/3) - 2*c**(sympy.S(1)/3)*x**(n/3))/(3*b**(sympy.S(1)/3)))/(b**(sympy.S(5)/3)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_496():
    f = x**(n/2 - 1)/(b*x**n + c*x**(2*n))
    F = -2/(b*n*x**(n/2)) + 2*sqrt(c)*atan(sqrt(b)/(sqrt(c)*x**(n/2)))/(b**(sympy.S(3)/2)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_497():
    f = x**(-n/2 - 1)/(b*x**n + c*x**(2*n))
    F = -2/(3*b*n*x**(3*n/2)) + 2*c/(b**2*n*x**(n/2)) - 2*c**(sympy.S(3)/2)*atan(sqrt(b)/(sqrt(c)*x**(n/2)))/(b**(sympy.S(5)/2)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_498():
    f = x**(-n/3 - 1)/(b*x**n + c*x**(2*n))
    F = -3/(4*b*n*x**(4*n/3)) + 3*c/(b**2*n*x**(n/3)) - c**(sympy.S(4)/3)*log(b**(sympy.S(1)/3)/x**(n/3) + c**(sympy.S(1)/3))/(b**(sympy.S(7)/3)*n) + c**(sympy.S(4)/3)*log(b**(sympy.S(2)/3)/x**(2*n/3) - b**(sympy.S(1)/3)*c**(sympy.S(1)/3)/x**(n/3) + c**(sympy.S(2)/3))/(2*b**(sympy.S(7)/3)*n) + sqrt(3)*c**(sympy.S(4)/3)*atan(sqrt(3)*(-2*b**(sympy.S(1)/3)/x**(n/3) + c**(sympy.S(1)/3))/(3*c**(sympy.S(1)/3)))/(b**(sympy.S(7)/3)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_499():
    f = x**(-n/4 - 1)/(b*x**n + c*x**(2*n))
    F = -4/(5*b*n*x**(5*n/4)) + 4*c/(b**2*n*x**(n/4)) + sqrt(2)*c**(sympy.S(5)/4)*log(-sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)/x**(n/4) + sqrt(b)/x**(n/2) + sqrt(c))/(2*b**(sympy.S(9)/4)*n) - sqrt(2)*c**(sympy.S(5)/4)*log(sqrt(2)*b**(sympy.S(1)/4)*c**(sympy.S(1)/4)/x**(n/4) + sqrt(b)/x**(n/2) + sqrt(c))/(2*b**(sympy.S(9)/4)*n) - sqrt(2)*c**(sympy.S(5)/4)*atan(sqrt(2)*b**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*x**(n/4)) - 1)/(b**(sympy.S(9)/4)*n) - sqrt(2)*c**(sympy.S(5)/4)*atan(sqrt(2)*b**(sympy.S(1)/4)/(c**(sympy.S(1)/4)*x**(n/4)) + 1)/(b**(sympy.S(9)/4)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_500():
    f = x**(-n*(p - 1) - 1)*(b*x**n + c*x**(2*n))**p
    F = (b*x**n + c*x**(2*n))**(p + 1)/(c*n*x**(n*(p + 1))*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_501():
    f = x**(-n*(2*p + 1) - 1)*(b*x**n + c*x**(2*n))**p
    F = -(b*x**n + c*x**(2*n))**(p + 1)/(b*n*x**(2*n*(p + 1))*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_502():
    f = x**(2*n - 1)*(a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(5)/2)
    F = -a*(a + b*x**n)**6*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/(6*n*(a*b**2 + b**3*x**n)) + (a + b*x**n)**7*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/(7*n*(a*b**2 + b**3*x**n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_503():
    f = x**(2*n - 1)*(a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(3)/2)
    F = -a*(a + b*x**n)**4*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/(4*n*(a*b**2 + b**3*x**n)) + (a + b*x**n)**5*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/(5*n*(a*b**2 + b**3*x**n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_504():
    f = x**(2*n - 1)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))
    F = a*x**(2*n)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/(2*n*(a + b*x**n)) + b**2*x**(3*n)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/(3*n*(a*b + b**2*x**n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_505():
    f = x**(2*n - 1)/sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))
    F = -a*(a + b*x**n)*log(a + b*x**n)/(b**2*n*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))) + x**n*(a + b*x**n)/(b*n*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_506():
    f = x**(2*n - 1)/(a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(3)/2)
    F = x**(2*n)/(2*a*n*(a + b*x**n)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_507():
    f = x**(2*n - 1)/(a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(5)/2)
    F = a/(4*b**2*n*(a + b*x**n)**3*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))) - 1/(3*b**2*n*(a + b*x**n)**2*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_508():
    f = x**(2*n - 1)/(a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(7)/2)
    F = a/(6*b**2*n*(a + b*x**n)**5*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))) - 1/(5*b**2*n*(a + b*x**n)**4*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_509():
    f = (d*x)**m*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))
    F = a*(d*x)**(m + 1)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/(d*(a + b*x**n)*(m + 1)) + b**2*x**(n + 1)*(d*x)**m*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/((a*b + b**2*x**n)*(m + n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_510():
    f = x**2*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))
    F = a*x**3*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/(3*a + 3*b*x**n) + b**2*x**(n + 3)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/((n + 3)*(a*b + b**2*x**n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_511():
    f = x*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))
    F = a*x**2*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/(2*a + 2*b*x**n) + b**2*x**(n + 2)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/((n + 2)*(a*b + b**2*x**n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_512():
    f = sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))
    F = a*x*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/(a + b*x**n) + b**2*x**(n + 1)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/((n + 1)*(a*b + b**2*x**n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_513():
    f = sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/x
    F = a*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))*log(x)/(a + b*x**n) + b**2*x**n*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/(n*(a*b + b**2*x**n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_514():
    f = sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/x**2
    F = -a*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/(x*(a + b*x**n)) - b**2*x**(n - 1)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/((1 - n)*(a*b + b**2*x**n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_515():
    f = sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/x**3
    F = -a*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/(2*x**2*(a + b*x**n)) - b**2*x**(n - 2)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/((2 - n)*(a*b + b**2*x**n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_516():
    f = (d*x)**m*(a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(3)/2)
    F = a**3*(d*x)**(m + 1)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/(d*(a + b*x**n)*(m + 1)) + 3*a**2*b**2*x**(n + 1)*(d*x)**m*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/((a*b + b**2*x**n)*(m + n + 1)) + 3*a*b**3*x**(2*n + 1)*(d*x)**m*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/((a*b + b**2*x**n)*(m + 2*n + 1)) + b**4*x**(3*n + 1)*(d*x)**m*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/((a*b + b**2*x**n)*(m + 3*n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_517():
    f = x**2*(a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(3)/2)
    F = a**3*x**3*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/(3*a + 3*b*x**n) + 3*a**2*b**2*x**(n + 3)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/((n + 3)*(a*b + b**2*x**n)) + 3*a*b**3*x**(2*n + 3)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/((2*n + 3)*(a*b + b**2*x**n)) + b**4*x**(3*n + 3)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/((3*n + 3)*(a*b + b**2*x**n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_518():
    f = x*(a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(3)/2)
    F = a**3*x**2*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/(2*a + 2*b*x**n) + 3*a**2*b**2*x**(n + 2)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/((n + 2)*(a*b + b**2*x**n)) + 3*a*b**3*x**(2*n + 2)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/((2*n + 2)*(a*b + b**2*x**n)) + b**4*x**(3*n + 2)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/((3*n + 2)*(a*b + b**2*x**n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_519():
    f = (a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(3)/2)
    F = a**3*x*(a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(3)/2)/(a + b*x**n)**3 + 3*a**2*b**4*x**(n + 1)*(a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(3)/2)/((n + 1)*(a*b + b**2*x**n)**3) + 3*a*b**5*x**(2*n + 1)*(a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(3)/2)/((2*n + 1)*(a*b + b**2*x**n)**3) + b**6*x**(3*n + 1)*(a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(3)/2)/((3*n + 1)*(a*b + b**2*x**n)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_520():
    f = (a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(3)/2)/x
    F = a**3*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))*log(x)/(a + b*x**n) + 3*a**2*b**2*x**n*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/(n*(a*b + b**2*x**n)) + 3*a*b**3*x**(2*n)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/(2*n*(a*b + b**2*x**n)) + b**4*x**(3*n)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/(3*n*(a*b + b**2*x**n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_521():
    f = (a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(3)/2)/x**2
    F = -a**3*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/(x*(a + b*x**n)) - 3*a**2*b**2*x**(n - 1)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/((1 - n)*(a*b + b**2*x**n)) - 3*a*b**3*x**(2*n - 1)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/((1 - 2*n)*(a*b + b**2*x**n)) - b**4*x**(3*n - 1)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/((1 - 3*n)*(a*b + b**2*x**n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_522():
    f = (a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(3)/2)/x**3
    F = -a**3*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/(2*x**2*(a + b*x**n)) - 3*a**2*b**2*x**(n - 2)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/((2 - n)*(a*b + b**2*x**n)) - 3*a*b**3*x**(2*n - 2)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/((2 - 2*n)*(a*b + b**2*x**n)) - b**4*x**(3*n - 2)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))/((2 - 3*n)*(a*b + b**2*x**n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_523():
    f = (d*x)**m/sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))
    F = (d*x)**(m + 1)*(a + b*x**n)*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(a*d*(m + 1)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_524():
    f = x**2/sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))
    F = x**3*(a + b*x**n)*hyper((1, 3/n), ((n + 3)/n,), -b*x**n/a)/(3*a*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_525():
    f = x/sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))
    F = x**2*(a + b*x**n)*hyper((1, 2/n), ((n + 2)/n,), -b*x**n/a)/(2*a*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_526():
    f = 1/sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))
    F = x*(a + b*x**n)*hyper((1, 1/n), (1 + 1/n,), -b*x**n/a)/(a*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_527():
    f = 1/(x*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n)))
    F = (a + b*x**n)*log(x)/(a*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))) - (a + b*x**n)*log(a + b*x**n)/(a*n*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_528():
    f = 1/(x**2*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n)))
    F = -(a + b*x**n)*hyper((1, -1/n), (-(1 - n)/n,), -b*x**n/a)/(a*x*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_529():
    f = 1/(x**3*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n)))
    F = -(a + b*x**n)*hyper((1, -2/n), (-(2 - n)/n,), -b*x**n/a)/(2*a*x**2*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_530():
    f = (d*x)**m/(a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(3)/2)
    F = (d*x)**(m + 1)*(a + b*x**n)*hyper((3, (m + 1)/n), ((m + n + 1)/n,), -b*x**n/a)/(a**3*d*(m + 1)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_531():
    f = x**2/(a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(3)/2)
    F = x**3*(a + b*x**n)*hyper((3, 3/n), ((n + 3)/n,), -b*x**n/a)/(3*a**3*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_532():
    f = x/(a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(3)/2)
    F = x**2*(a + b*x**n)*hyper((3, 2/n), ((n + 2)/n,), -b*x**n/a)/(2*a**3*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_533():
    f = (a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(-3)/2)
    F = x*(a + b*x**n)**3*hyper((3, 1/n), (1 + 1/n,), -b*x**n/a)/(a**3*(a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_534():
    f = 1/(x*(a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(3)/2))
    F = 1/(2*a*n*(a + b*x**n)*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))) + 1/(a**2*n*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))) + (a + b*x**n)*log(x)/(a**3*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n))) - (a + b*x**n)*log(a + b*x**n)/(a**3*n*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_535():
    f = 1/(x**2*(a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(3)/2))
    F = -(a + b*x**n)*hyper((3, -1/n), (-(1 - n)/n,), -b*x**n/a)/(a**3*x*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_536():
    f = 1/(x**3*(a**2 + 2*a*b*x**n + b**2*x**(2*n))**(sympy.S(3)/2))
    F = -(a + b*x**n)*hyper((3, -2/n), (-(2 - n)/n,), -b*x**n/a)/(2*a**3*x**2*sqrt(a**2 + 2*a*b*x**n + b**2*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_537():
    f = (a**2 + 2*a*b/x**(1/(2*p + 1)) + b**2/x**(2/(2*p + 1)))**p
    F = x*(a + b*x**(1/(-2*p - 1)))*(a**2 + 2*a*b*x**(1/(-2*p - 1)) + b**2/x**(2/(2*p + 1)))**p/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_538():
    f = (a**2 + 2*a*b*x**n + b**2*x**(2*n))**((-n - 1)/(2*n))
    F = x*(a + b*x**n)/(a*(a**2 + 2*a*b*x**n + b**2*x**(2*n))**((n + 1)/(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_539():
    f = (a**2 + 2*a*b/x**(1/(2*p + 2)) + b**2/x**(1/(p + 1)))**p
    F = x*(a + b/x**(1/(2*p + 2)))*(2*p + 2)*(a**2 + 2*a*b/x**(1/(2*p + 2)) + b**2/x**(1/(p + 1)))**p/(a*(2*p + 1)) - x*(a + b/x**(1/(2*p + 2)))**2*(a**2 + 2*a*b/x**(1/(2*p + 2)) + b**2/x**(1/(p + 1)))**p/(a**2*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_540():
    f = (a**2 + 2*a*b*x**n + b**2*x**(2*n))**((-2*n - 1)/(2*n))
    F = x*(a + b*x**n)*(a**2 + 2*a*b*x**n + b**2*x**(2*n))**(-1 - 1/(2*n))/(a*(n + 1)) + n*x*(a + b*x**n)**2*(a**2 + 2*a*b*x**n + b**2*x**(2*n))**(-1 - 1/(2*n))/(a**2*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_541():
    f = x**(2*n - 1)*(a**2 + 2*a*b*x**n + b**2*x**(2*n))**p
    F = a**2*(1 + b*x**n/a)**2*(a**2 + 2*a*b*x**n + b**2*x**(2*n))**p/(2*b**2*n*(p + 1)) - a**2*(1 + b*x**n/a)*(a**2 + 2*a*b*x**n + b**2*x**(2*n))**p/(b**2*n*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_542():
    f = x**(4*n - 1)/(a + b*x**n + c*x**(2*n))
    F = -b*x**n/(c**2*n) + b*(-3*a*c + b**2)*atanh((b + 2*c*x**n)/sqrt(-4*a*c + b**2))/(c**3*n*sqrt(-4*a*c + b**2)) + x**(2*n)/(2*c*n) + (-a*c + b**2)*log(a + b*x**n + c*x**(2*n))/(2*c**3*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_543():
    f = x**(3*n - 1)/(a + b*x**n + c*x**(2*n))
    F = -b*log(a + b*x**n + c*x**(2*n))/(2*c**2*n) + x**n/(c*n) - (-2*a*c + b**2)*atanh((b + 2*c*x**n)/sqrt(-4*a*c + b**2))/(c**2*n*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_544():
    f = x**(2*n - 1)/(a + b*x**n + c*x**(2*n))
    F = b*atanh((b + 2*c*x**n)/sqrt(-4*a*c + b**2))/(c*n*sqrt(-4*a*c + b**2)) + log(a + b*x**n + c*x**(2*n))/(2*c*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_545():
    f = x**(n - 1)/(a + b*x**n + c*x**(2*n))
    F = -2*atanh((b + 2*c*x**n)/sqrt(-4*a*c + b**2))/(n*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_546():
    f = x**(-n - 1)/(a + b*x**n + c*x**(2*n))
    F = -1/(a*n*x**n) - b*log(x)/a**2 + b*log(a + b*x**n + c*x**(2*n))/(2*a**2*n) - (-2*a*c + b**2)*atanh((b + 2*c*x**n)/sqrt(-4*a*c + b**2))/(a**2*n*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_547():
    f = x**(-2*n - 1)/(a + b*x**n + c*x**(2*n))
    F = -1/(2*a*n*x**(2*n)) + b/(a**2*n*x**n) + b*(-3*a*c + b**2)*atanh((b + 2*c*x**n)/sqrt(-4*a*c + b**2))/(a**3*n*sqrt(-4*a*c + b**2)) + (-a*c + b**2)*log(x)/a**3 - (-a*c + b**2)*log(a + b*x**n + c*x**(2*n))/(2*a**3*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_548():
    f = x**(-3*n - 1)/(a + b*x**n + c*x**(2*n))
    F = -1/(3*a*n*x**(3*n)) + b/(2*a**2*n*x**(2*n)) - (-a*c + b**2)/(a**3*n*x**n) - b*(-2*a*c + b**2)*log(x)/a**4 + b*(-2*a*c + b**2)*log(a + b*x**n + c*x**(2*n))/(2*a**4*n) - (2*a**2*c**2 - 4*a*b**2*c + b**4)*atanh((b + 2*c*x**n)/sqrt(-4*a*c + b**2))/(a**4*n*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_549():
    f = x**(n/4 - 1)/(a + b*x**n + c*x**(2*n))
    F = -2*2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x**(n/4)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(n*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*sqrt(-4*a*c + b**2)) - 2*2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x**(n/4)/(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(n*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*sqrt(-4*a*c + b**2)) + 2*2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*atan(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x**(n/4)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(n*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*sqrt(-4*a*c + b**2)) + 2*2**(sympy.S(3)/4)*c**(sympy.S(3)/4)*atanh(2**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x**(n/4)/(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4))/(n*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_550():
    f = x**(n/3 - 1)/(a + b*x**n + c*x**(2*n))
    F = -2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x**(n/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(n*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**(2*n/3) - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x**(n/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(2*n*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(2)/3)*sqrt(3)*c**(sympy.S(2)/3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x**(n/3)/(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(n*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)*sqrt(-4*a*c + b**2)) + 2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x**(n/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(n*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*log(2**(sympy.S(2)/3)*c**(sympy.S(2)/3)*x**(2*n/3) - 2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x**(n/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(2*n*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)*sqrt(-4*a*c + b**2)) - 2**(sympy.S(2)/3)*sqrt(3)*c**(sympy.S(2)/3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*c**(sympy.S(1)/3)*x**(n/3)/(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3) + 1)/3)/(n*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_551():
    f = x**(n/2 - 1)/(a + b*x**n + c*x**(2*n))
    F = -2*sqrt(2)*sqrt(c)*atan(sqrt(2)*sqrt(c)*x**(n/2)/sqrt(b + sqrt(-4*a*c + b**2)))/(n*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)) + 2*sqrt(2)*sqrt(c)*atan(sqrt(2)*sqrt(c)*x**(n/2)/sqrt(b - sqrt(-4*a*c + b**2)))/(n*sqrt(b - sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_552():
    f = x**(-n/2 - 1)/(a + b*x**n + c*x**(2*n))
    F = -2/(a*n*x**(n/2)) + sqrt(2)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(a)/(x**(n/2)*sqrt(b - sqrt(-4*a*c + b**2))))/(a**(sympy.S(3)/2)*n*sqrt(b - sqrt(-4*a*c + b**2))) + sqrt(2)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(a)/(x**(n/2)*sqrt(b + sqrt(-4*a*c + b**2))))/(a**(sympy.S(3)/2)*n*sqrt(b + sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_553():
    f = x**(-n/3 - 1)/(a + b*x**n + c*x**(2*n))
    F = -3/(a*n*x**(n/3)) + 2**(sympy.S(2)/3)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3)/x**(n/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(2*a**(sympy.S(4)/3)*n*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*a**(sympy.S(2)/3)/x**(2*n/3) - 2**(sympy.S(1)/3)*a**(sympy.S(1)/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)/x**(n/3) + (b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(4*a**(sympy.S(4)/3)*n*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*sqrt(3)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*a**(sympy.S(1)/3)/(x**(n/3)*(b - sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) + 1)/3)/(2*a**(sympy.S(4)/3)*n*(b - sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) + 2**(sympy.S(2)/3)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3)/x**(n/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3))/(2*a**(sympy.S(4)/3)*n*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*log(2**(sympy.S(2)/3)*a**(sympy.S(2)/3)/x**(2*n/3) - 2**(sympy.S(1)/3)*a**(sympy.S(1)/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)/x**(n/3) + (b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))/(4*a**(sympy.S(4)/3)*n*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3)) - 2**(sympy.S(2)/3)*sqrt(3)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*a**(sympy.S(1)/3)/(x**(n/3)*(b + sqrt(-4*a*c + b**2))**(sympy.S(1)/3)) + 1)/3)/(2*a**(sympy.S(4)/3)*n*(b + sqrt(-4*a*c + b**2))**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_554():
    f = x**(-n/4 - 1)/(a + b*x**n + c*x**(2*n))
    F = -4/(a*n*x**(n/4)) - 2**(sympy.S(3)/4)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*a**(sympy.S(1)/4)/(x**(n/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)))/(a**(sympy.S(5)/4)*n*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) - 2**(sympy.S(3)/4)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*a**(sympy.S(1)/4)/(x**(n/4)*(-b + sqrt(-4*a*c + b**2))**(sympy.S(1)/4)))/(a**(sympy.S(5)/4)*n*(-b + sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) - 2**(sympy.S(3)/4)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(2**(sympy.S(1)/4)*a**(sympy.S(1)/4)/(x**(n/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)))/(a**(sympy.S(5)/4)*n*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4)) - 2**(sympy.S(3)/4)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(2**(sympy.S(1)/4)*a**(sympy.S(1)/4)/(x**(n/4)*(-b - sqrt(-4*a*c + b**2))**(sympy.S(1)/4)))/(a**(sympy.S(5)/4)*n*(-b - sqrt(-4*a*c + b**2))**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_555():
    f = x**2/(a + b*x**n + c*x**(2*n))
    F = -2*c*x**3*hyper((1, 3/n), ((n + 3)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(-12*a*c + 3*b**2 + 3*b*sqrt(-4*a*c + b**2)) - 2*c*x**3*hyper((1, 3/n), ((n + 3)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(-12*a*c + 3*b**2 - 3*b*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_556():
    f = x/(a + b*x**n + c*x**(2*n))
    F = -c*x**2*hyper((1, 2/n), ((n + 2)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2)) - c*x**2*hyper((1, 2/n), ((n + 2)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_557():
    f = 1/(a + b*x**n + c*x**(2*n))
    F = -2*c*x*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2)) - 2*c*x*hyper((1, 1/n), (1 + 1/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_558():
    f = 1/(x*(a + b*x**n + c*x**(2*n)))
    F = b*atanh((b + 2*c*x**n)/sqrt(-4*a*c + b**2))/(a*n*sqrt(-4*a*c + b**2)) + log(x)/a - log(a + b*x**n + c*x**(2*n))/(2*a*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_559():
    f = 1/(x**2*(a + b*x**n + c*x**(2*n)))
    F = 2*c*hyper((1, -1/n), (-(1 - n)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(x*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))) + 2*c*hyper((1, -1/n), (-(1 - n)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(x*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_560():
    f = 1/(x**3*(a + b*x**n + c*x**(2*n)))
    F = c*hyper((1, -2/n), (-(2 - n)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(x**2*(-4*a*c + b**2 + b*sqrt(-4*a*c + b**2))) + c*hyper((1, -2/n), (-(2 - n)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(x**2*(-4*a*c + b**2 - b*sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_561():
    f = x**3*sqrt(a + b*x**n + c*x**(2*n))
    F = x**4*sqrt(a + b*x**n + c*x**(2*n))*appellf1(4/n, sympy.S(-1)/2, sympy.S(-1)/2, (n + 4)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(4*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_562():
    f = x**2*sqrt(a + b*x**n + c*x**(2*n))
    F = x**3*sqrt(a + b*x**n + c*x**(2*n))*appellf1(3/n, sympy.S(-1)/2, sympy.S(-1)/2, (n + 3)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(3*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_563():
    f = x*sqrt(a + b*x**n + c*x**(2*n))
    F = x**2*sqrt(a + b*x**n + c*x**(2*n))*appellf1(2/n, sympy.S(-1)/2, sympy.S(-1)/2, (n + 2)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(2*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_564():
    f = sqrt(a + b*x**n + c*x**(2*n))
    F = x*sqrt(a + b*x**n + c*x**(2*n))*appellf1(1/n, sympy.S(-1)/2, sympy.S(-1)/2, 1 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_565():
    f = sqrt(a + b*x**n + c*x**(2*n))/x
    F = -sqrt(a)*atanh((2*a + b*x**n)/(2*sqrt(a)*sqrt(a + b*x**n + c*x**(2*n))))/n + b*atanh((b + 2*c*x**n)/(2*sqrt(c)*sqrt(a + b*x**n + c*x**(2*n))))/(2*sqrt(c)*n) + sqrt(a + b*x**n + c*x**(2*n))/n
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_566():
    f = sqrt(a + b*x**n + c*x**(2*n))/x**2
    F = -sqrt(a + b*x**n + c*x**(2*n))*appellf1(-1/n, sympy.S(-1)/2, sympy.S(-1)/2, -(1 - n)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(x*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_567():
    f = sqrt(a + b*x**n + c*x**(2*n))/x**3
    F = -sqrt(a + b*x**n + c*x**(2*n))*appellf1(-2/n, sympy.S(-1)/2, sympy.S(-1)/2, -(2 - n)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(2*x**2*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_568():
    f = x**3*(a + b*x**n + c*x**(2*n))**(sympy.S(3)/2)
    F = a*x**4*sqrt(a + b*x**n + c*x**(2*n))*appellf1(4/n, sympy.S(-3)/2, sympy.S(-3)/2, (n + 4)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(4*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_569():
    f = x**2*(a + b*x**n + c*x**(2*n))**(sympy.S(3)/2)
    F = a*x**3*sqrt(a + b*x**n + c*x**(2*n))*appellf1(3/n, sympy.S(-3)/2, sympy.S(-3)/2, (n + 3)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(3*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_570():
    f = x*(a + b*x**n + c*x**(2*n))**(sympy.S(3)/2)
    F = a*x**2*sqrt(a + b*x**n + c*x**(2*n))*appellf1(2/n, sympy.S(-3)/2, sympy.S(-3)/2, (n + 2)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(2*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_571():
    f = (a + b*x**n + c*x**(2*n))**(sympy.S(3)/2)
    F = a*x*sqrt(a + b*x**n + c*x**(2*n))*appellf1(1/n, sympy.S(-3)/2, sympy.S(-3)/2, 1 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_572():
    f = (a + b*x**n + c*x**(2*n))**(sympy.S(3)/2)/x
    F = -a**(sympy.S(3)/2)*atanh((2*a + b*x**n)/(2*sqrt(a)*sqrt(a + b*x**n + c*x**(2*n))))/n - b*(-12*a*c + b**2)*atanh((b + 2*c*x**n)/(2*sqrt(c)*sqrt(a + b*x**n + c*x**(2*n))))/(16*c**(sympy.S(3)/2)*n) + (a + b*x**n + c*x**(2*n))**(sympy.S(3)/2)/(3*n) + sqrt(a + b*x**n + c*x**(2*n))*(8*a*c + b**2 + 2*b*c*x**n)/(8*c*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_573():
    f = (a + b*x**n + c*x**(2*n))**(sympy.S(3)/2)/x**2
    F = -a*sqrt(a + b*x**n + c*x**(2*n))*appellf1(-1/n, sympy.S(-3)/2, sympy.S(-3)/2, -(1 - n)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(x*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_574():
    f = (a + b*x**n + c*x**(2*n))**(sympy.S(3)/2)/x**3
    F = -a*sqrt(a + b*x**n + c*x**(2*n))*appellf1(-2/n, sympy.S(-3)/2, sympy.S(-3)/2, -(2 - n)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(2*x**2*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_575():
    f = x**3/sqrt(a + b*x**n + c*x**(2*n))
    F = x**4*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(4/n, sympy.S.Half, sympy.S.Half, (n + 4)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(4*sqrt(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_576():
    f = x**2/sqrt(a + b*x**n + c*x**(2*n))
    F = x**3*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(3/n, sympy.S.Half, sympy.S.Half, (n + 3)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(3*sqrt(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_577():
    f = x/sqrt(a + b*x**n + c*x**(2*n))
    F = x**2*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(2/n, sympy.S.Half, sympy.S.Half, (n + 2)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(2*sqrt(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_578():
    f = 1/sqrt(a + b*x**n + c*x**(2*n))
    F = x*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(1/n, sympy.S.Half, sympy.S.Half, 1 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/sqrt(a + b*x**n + c*x**(2*n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_579():
    f = 1/(x*sqrt(a + b*x**n + c*x**(2*n)))
    F = -atanh((2*a + b*x**n)/(2*sqrt(a)*sqrt(a + b*x**n + c*x**(2*n))))/(sqrt(a)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_580():
    f = 1/(x**2*sqrt(a + b*x**n + c*x**(2*n)))
    F = -sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(-1/n, sympy.S.Half, sympy.S.Half, -(1 - n)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(x*sqrt(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_581():
    f = 1/(x**3*sqrt(a + b*x**n + c*x**(2*n)))
    F = -sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(-2/n, sympy.S.Half, sympy.S.Half, -(2 - n)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(2*x**2*sqrt(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_582():
    f = x**3/(a + b*x**n + c*x**(2*n))**(sympy.S(3)/2)
    F = x**4*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(4/n, sympy.S(3)/2, sympy.S(3)/2, (n + 4)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_583():
    f = x**2/(a + b*x**n + c*x**(2*n))**(sympy.S(3)/2)
    F = x**3*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(3/n, sympy.S(3)/2, sympy.S(3)/2, (n + 3)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(3*a*sqrt(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_584():
    f = x/(a + b*x**n + c*x**(2*n))**(sympy.S(3)/2)
    F = x**2*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(2/n, sympy.S(3)/2, sympy.S(3)/2, (n + 2)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(2*a*sqrt(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_585():
    f = (a + b*x**n + c*x**(2*n))**(sympy.S(-3)/2)
    F = x*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(1/n, sympy.S(3)/2, sympy.S(3)/2, 1 + 1/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*sqrt(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_586():
    f = 1/(x*(a + b*x**n + c*x**(2*n))**(sympy.S(3)/2))
    F = (-4*a*c + 2*b**2 + 2*b*c*x**n)/(a*n*(-4*a*c + b**2)*sqrt(a + b*x**n + c*x**(2*n))) - atanh((2*a + b*x**n)/(2*sqrt(a)*sqrt(a + b*x**n + c*x**(2*n))))/(a**(sympy.S(3)/2)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_587():
    f = 1/(x**2*(a + b*x**n + c*x**(2*n))**(sympy.S(3)/2))
    F = -sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(-1/n, sympy.S(3)/2, sympy.S(3)/2, -(1 - n)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*x*sqrt(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_588():
    f = 1/(x**3*(a + b*x**n + c*x**(2*n))**(sympy.S(3)/2))
    F = -sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(-2/n, sympy.S(3)/2, sympy.S(3)/2, -(2 - n)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(2*a*x**2*sqrt(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_589():
    f = (d*x)**m*(a + b*x**n + c*x**(2*n))**3
    F = a**3*(d*x)**(m + 1)/(d*(m + 1)) + 3*a**2*b*x**(n + 1)*(d*x)**m/(m + n + 1) + 3*a*x**(2*n + 1)*(d*x)**m*(a*c + b**2)/(m + 2*n + 1) + 3*b*c**2*x**(5*n + 1)*(d*x)**m/(m + 5*n + 1) + b*x**(3*n + 1)*(d*x)**m*(6*a*c + b**2)/(m + 3*n + 1) + c**3*x**(6*n + 1)*(d*x)**m/(m + 6*n + 1) + 3*c*x**(4*n + 1)*(d*x)**m*(a*c + b**2)/(m + 4*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_590():
    f = (d*x)**m*(a + b*x**n + c*x**(2*n))**2
    F = a**2*(d*x)**(m + 1)/(d*(m + 1)) + 2*a*b*x**(n + 1)*(d*x)**m/(m + n + 1) + 2*b*c*x**(3*n + 1)*(d*x)**m/(m + 3*n + 1) + c**2*x**(4*n + 1)*(d*x)**m/(m + 4*n + 1) + x**(2*n + 1)*(d*x)**m*(2*a*c + b**2)/(m + 2*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_591():
    f = (d*x)**m*(a + b*x**n + c*x**(2*n))
    F = a*(d*x)**(m + 1)/(d*(m + 1)) + b*x**(n + 1)*(d*x)**m/(m + n + 1) + c*x**(2*n + 1)*(d*x)**m/(m + 2*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_592():
    f = (d*x)**m/(a + b*x**n + c*x**(2*n))
    F = -2*c*(d*x)**(m + 1)*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(d*(b + sqrt(-4*a*c + b**2))*(m + 1)*sqrt(-4*a*c + b**2)) + 2*c*(d*x)**(m + 1)*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(d*(b - sqrt(-4*a*c + b**2))*(m + 1)*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_593():
    f = (d*x)**m/(a + b*x**n + c*x**(2*n))**2
    F = -c*(d*x)**(m + 1)*(4*a*c*(m - 2*n + 1) - b**2*(m - n + 1) + b*sqrt(-4*a*c + b**2)*(m - n + 1))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*d*n*(b + sqrt(-4*a*c + b**2))*(m + 1)*(-4*a*c + b**2)**(sympy.S(3)/2)) + c*(d*x)**(m + 1)*(-b*(m - n + 1) + (4*a*c*(m - 2*n + 1) - b**2*(m - n + 1))/sqrt(-4*a*c + b**2))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(a*d*n*(b - sqrt(-4*a*c + b**2))*(m + 1)*(-4*a*c + b**2)) + (d*x)**(m + 1)*(-2*a*c + b**2 + b*c*x**n)/(a*d*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_594():
    f = (d*x)**m/(a + b*x**n + c*x**(2*n))**3
    F = (d*x)**(m + 1)*(-2*a*c + b**2 + b*c*x**n)/(2*a*d*n*(-4*a*c + b**2)*(a + b*x**n + c*x**(2*n))**2) - c*(d*x)**(m + 1)*(8*a**2*c**2*(m**2 + m*(2 - 6*n) + 8*n**2 - 6*n + 1) - 6*a*b**2*c*(m**2 + m*(2 - 4*n) + 3*n**2 - 4*n + 1) + b**4*(m**2 + m*(2 - 3*n) + 2*n**2 - 3*n + 1) + b*sqrt(-4*a*c + b**2)*(2*a*c*(2*m - 7*n + 2) - b**2*(m - 2*n + 1))*(m - n + 1))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(2*a**2*d*n**2*(b + sqrt(-4*a*c + b**2))*(m + 1)*(-4*a*c + b**2)**(sympy.S(5)/2)) - c*(d*x)**(m + 1)*(-8*a**2*c**2*(m**2 + m*(2 - 6*n) + 8*n**2 - 6*n + 1) + 6*a*b**2*c*(m**2 + m*(2 - 4*n) + 3*n**2 - 4*n + 1) - b**4*(m**2 + m*(2 - 3*n) + 2*n**2 - 3*n + 1) + b*sqrt(-4*a*c + b**2)*(2*a*c*(2*m - 7*n + 2) - b**2*(m - 2*n + 1))*(m - n + 1))*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -2*c*x**n/(b - sqrt(-4*a*c + b**2)))/(2*a**2*d*n**2*(b - sqrt(-4*a*c + b**2))*(m + 1)*(-4*a*c + b**2)**(sympy.S(5)/2)) - (d*x)**(m + 1)*(4*a**2*c**2*(m - 4*n + 1) - 5*a*b**2*c*(m - 3*n + 1) + b**4*(m - 2*n + 1) - b*c*x**n*(2*a*c*(2*m - 7*n + 2) - b**2*(m - 2*n + 1)))/(2*a**2*d*n**2*(-4*a*c + b**2)**2*(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_595():
    f = (d*x)**m*(a + b*x**n + c*x**(2*n))**(sympy.S(3)/2)
    F = a*(d*x)**(m + 1)*sqrt(a + b*x**n + c*x**(2*n))*appellf1((m + 1)/n, sympy.S(-3)/2, sympy.S(-3)/2, (m + n + 1)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(d*(m + 1)*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_596():
    f = (d*x)**m*sqrt(a + b*x**n + c*x**(2*n))
    F = (d*x)**(m + 1)*sqrt(a + b*x**n + c*x**(2*n))*appellf1((m + 1)/n, sympy.S(-1)/2, sympy.S(-1)/2, (m + n + 1)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(d*(m + 1)*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_597():
    f = (d*x)**m/sqrt(a + b*x**n + c*x**(2*n))
    F = (d*x)**(m + 1)*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1((m + 1)/n, sympy.S.Half, sympy.S.Half, (m + n + 1)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(d*(m + 1)*sqrt(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_598():
    f = (d*x)**m/(a + b*x**n + c*x**(2*n))**(sympy.S(3)/2)
    F = (d*x)**(m + 1)*sqrt(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1((m + 1)/n, sympy.S(3)/2, sympy.S(3)/2, (m + n + 1)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(a*d*(m + 1)*sqrt(a + b*x**n + c*x**(2*n)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_599():
    f = (d*x)**m*(a + b*x**n + c*x**(2*n))**p
    F = (d*x)**(m + 1)*(a + b*x**n + c*x**(2*n))**p*appellf1((m + 1)/n, -p, -p, (m + n + 1)/n, -2*c*x**n/(b - sqrt(-4*a*c + b**2)), -2*c*x**n/(b + sqrt(-4*a*c + b**2)))/(d*(m + 1)*(2*c*x**n/(b - sqrt(-4*a*c + b**2)) + 1)**p*(2*c*x**n/(b + sqrt(-4*a*c + b**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_600():
    f = (d + e*x)**3*(a + b*(d + e*x)**2 + c*(d + e*x)**4)
    F = a*(d + e*x)**4/(4*e) + b*(d + e*x)**6/(6*e) + c*(d + e*x)**8/(8*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_601():
    f = (d + e*x)**3*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2
    F = a**2*(d + e*x)**4/(4*e) + a*b*(d + e*x)**6/(3*e) + b*c*(d + e*x)**10/(5*e) + c**2*(d + e*x)**12/(12*e) + (d + e*x)**8*(2*a*c + b**2)/(8*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_602():
    f = (d + e*x)**3*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**3
    F = a**3*(d + e*x)**4/(4*e) + a**2*b*(d + e*x)**6/(2*e) + 3*a*(d + e*x)**8*(a*c + b**2)/(8*e) + 3*b*c**2*(d + e*x)**14/(14*e) + b*(d + e*x)**10*(6*a*c + b**2)/(10*e) + c**3*(d + e*x)**16/(16*e) + c*(d + e*x)**12*(a*c + b**2)/(4*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_603():
    f = (d*f + e*f*x)**3*(a + b*(d + e*x)**2 + c*(d + e*x)**4)
    F = a*f**3*(d + e*x)**4/(4*e) + b*f**3*(d + e*x)**6/(6*e) + c*f**3*(d + e*x)**8/(8*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_604():
    f = (d*f + e*f*x)**3*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2
    F = a**2*f**3*(d + e*x)**4/(4*e) + a*b*f**3*(d + e*x)**6/(3*e) + b*c*f**3*(d + e*x)**10/(5*e) + c**2*f**3*(d + e*x)**12/(12*e) + f**3*(d + e*x)**8*(2*a*c + b**2)/(8*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_605():
    f = (d*f + e*f*x)**3*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**3
    F = a**3*f**3*(d + e*x)**4/(4*e) + a**2*b*f**3*(d + e*x)**6/(2*e) + 3*a*f**3*(d + e*x)**8*(a*c + b**2)/(8*e) + 3*b*c**2*f**3*(d + e*x)**14/(14*e) + b*f**3*(d + e*x)**10*(6*a*c + b**2)/(10*e) + c**3*f**3*(d + e*x)**16/(16*e) + c*f**3*(d + e*x)**12*(a*c + b**2)/(4*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_606():
    f = (d + e*x)**4/(a + b*(d + e*x)**2 + c*(d + e*x)**4)
    F = x/c - sqrt(2)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*e*sqrt(b - sqrt(-4*a*c + b**2))) - sqrt(2)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*e*sqrt(b + sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_607():
    f = (d + e*x)**3/(a + b*(d + e*x)**2 + c*(d + e*x)**4)
    F = b*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(2*c*e*sqrt(-4*a*c + b**2)) + log(a + b*(d + e*x)**2 + c*(d + e*x)**4)/(4*c*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_608():
    f = (d + e*x)**2/(a + b*(d + e*x)**2 + c*(d + e*x)**4)
    F = -sqrt(2)*sqrt(b - sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(c)*e*sqrt(-4*a*c + b**2)) + sqrt(2)*sqrt(b + sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(c)*e*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_609():
    f = (d + e*x)/(a + b*(d + e*x)**2 + c*(d + e*x)**4)
    F = -atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(e*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_610():
    f = 1/((d + e*x)*(a + b*(d + e*x)**2 + c*(d + e*x)**4))
    F = b*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(2*a*e*sqrt(-4*a*c + b**2)) + log(d + e*x)/(a*e) - log(a + b*(d + e*x)**2 + c*(d + e*x)**4)/(4*a*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_611():
    f = 1/((d + e*x)**2*(a + b*(d + e*x)**2 + c*(d + e*x)**4))
    F = -sqrt(2)*sqrt(c)*(-b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(2*a*e*sqrt(b + sqrt(-4*a*c + b**2))) - sqrt(2)*sqrt(c)*(b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(2*a*e*sqrt(b - sqrt(-4*a*c + b**2))) - 1/(a*e*(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_612():
    f = 1/((d + e*x)**3*(a + b*(d + e*x)**2 + c*(d + e*x)**4))
    F = -1/(2*a*e*(d + e*x)**2) - b*log(d + e*x)/(a**2*e) + b*log(a + b*(d + e*x)**2 + c*(d + e*x)**4)/(4*a**2*e) - (-2*a*c + b**2)*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(2*a**2*e*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_613():
    f = 1/((d + e*x)**4*(a + b*(d + e*x)**2 + c*(d + e*x)**4))
    F = -1/(3*a*e*(d + e*x)**3) + b/(a**2*e*(d + e*x)) + sqrt(2)*sqrt(c)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(2*a**2*e*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*sqrt(c)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(2*a**2*e*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_614():
    f = (d + e*x)**4/(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2
    F = (2*a + b*(d + e*x)**2)*(d + e*x)/(e*(-8*a*c + 2*b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)) + sqrt(2)*(b - (4*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(4*sqrt(c)*e*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*(4*a*c + b**2 + b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(4*sqrt(c)*e*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_615():
    f = (d + e*x)**3/(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2
    F = -b*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(e*(-4*a*c + b**2)**(sympy.S(3)/2)) + (2*a + b*(d + e*x)**2)/(e*(-8*a*c + 2*b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_616():
    f = (d + e*x)**2/(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2
    F = -sqrt(2)*sqrt(c)*(2*b + sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(2*e*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + sqrt(2)*sqrt(c)*(2*b - sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(2*e*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) - (b + 2*c*(d + e*x)**2)*(d + e*x)/(e*(-8*a*c + 2*b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_617():
    f = (d + e*x)/(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2
    F = 2*c*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(e*(-4*a*c + b**2)**(sympy.S(3)/2)) + (-b - 2*c*(d + e*x)**2)/(e*(-8*a*c + 2*b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_618():
    f = (a + b*(d + e*x)**2 + c*(d + e*x)**4)**(-2)
    F = -sqrt(2)*sqrt(c)*(-12*a*c + b**2 - b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*e*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + sqrt(2)*sqrt(c)*(-12*a*c + b**2 + b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*e*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + (d/e + x)*(-2*a*c + b**2 + b*c*e**2*(d/e + x)**2)/(2*a*(-4*a*c + b**2)*(a + b*e**2*(d/e + x)**2 + c*e**4*(d/e + x)**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_619():
    f = 1/((d + e*x)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2)
    F = (-2*a*c + b**2 + b*c*(d + e*x)**2)/(2*a*e*(-4*a*c + b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)) + b*(-6*a*c + b**2)*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(2*a**2*e*(-4*a*c + b**2)**(sympy.S(3)/2)) + log(d + e*x)/(a**2*e) - log(a + b*(d + e*x)**2 + c*(d + e*x)**4)/(4*a**2*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_620():
    f = 1/((d + e*x)**2*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2)
    F = (-2*a*c + b**2 + b*c*(d + e*x)**2)/(2*a*e*(d + e*x)*(-4*a*c + b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)) + sqrt(2)*sqrt(c)*(-16*a*b*c + 3*b**3 - (-10*a*c + 3*b**2)*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a**2*e*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) - sqrt(2)*sqrt(c)*(-16*a*b*c + 3*b**3 + (-10*a*c + 3*b**2)*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a**2*e*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) - (-10*a*c + 3*b**2)/(2*a**2*e*(d + e*x)*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_621():
    f = 1/((d + e*x)**3*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2)
    F = (-2*a*c + b**2 + b*c*(d + e*x)**2)/(2*a*e*(d + e*x)**2*(-4*a*c + b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)) - (-3*a*c + b**2)/(a**2*e*(d + e*x)**2*(-4*a*c + b**2)) - 2*b*log(d + e*x)/(a**3*e) + b*log(a + b*(d + e*x)**2 + c*(d + e*x)**4)/(2*a**3*e) - (6*a**2*c**2 - 6*a*b**2*c + b**4)*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(a**3*e*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_622():
    f = 1/((d + e*x)**4*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2)
    F = (-2*a*c + b**2 + b*c*(d + e*x)**2)/(2*a*e*(d + e*x)**3*(-4*a*c + b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)) - (-14*a*c + 5*b**2)/(6*a**2*e*(d + e*x)**3*(-4*a*c + b**2)) + b*(-19*a*c + 5*b**2)/(2*a**3*e*(d + e*x)*(-4*a*c + b**2)) - sqrt(2)*sqrt(c)*(28*a**2*c**2 - 29*a*b**2*c + 5*b**4 - b*(-19*a*c + 5*b**2)*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a**3*e*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + sqrt(2)*sqrt(c)*(28*a**2*c**2 - 29*a*b**2*c + 5*b**4 + b*(-19*a*c + 5*b**2)*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a**3*e*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_623():
    f = (d + e*x)**4/(a + b*(d + e*x)**2 + c*(d + e*x)**4)**3
    F = -3*sqrt(2)*sqrt(c)*(4*a*c + 3*b**2 + 2*b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(8*e*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(5)/2)) + 3*sqrt(2)*sqrt(c)*(4*a*c + 3*b**2 - 2*b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(8*e*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(5)/2)) + (2*a + b*(d + e*x)**2)*(d + e*x)/(e*(-16*a*c + 4*b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2) - (d + e*x)*(-4*a*c + 7*b**2 + 12*b*c*(d + e*x)**2)/(8*e*(-4*a*c + b**2)**2*(a + b*(d + e*x)**2 + c*(d + e*x)**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_624():
    f = (d + e*x)**3/(a + b*(d + e*x)**2 + c*(d + e*x)**4)**3
    F = 3*b*c*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(e*(-4*a*c + b**2)**(sympy.S(5)/2)) - 3*b*(b + 2*c*(d + e*x)**2)/(4*e*(-4*a*c + b**2)**2*(a + b*(d + e*x)**2 + c*(d + e*x)**4)) + (2*a + b*(d + e*x)**2)/(e*(-16*a*c + 4*b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_625():
    f = (d + e*x)**2/(a + b*(d + e*x)**2 + c*(d + e*x)**4)**3
    F = -(b + 2*c*(d + e*x)**2)*(d + e*x)/(e*(-16*a*c + 4*b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2) + sqrt(2)*sqrt(c)*(20*a*c + b**2 - b*(-52*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(16*a*e*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + sqrt(2)*sqrt(c)*(20*a*c + b**2 + b*(-52*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(16*a*e*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + (d + e*x)*(b*(8*a*c + b**2) + c*(d + e*x)**2*(20*a*c + b**2))/(8*a*e*(-4*a*c + b**2)**2*(a + b*(d + e*x)**2 + c*(d + e*x)**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_626():
    f = (d + e*x)/(a + b*(d + e*x)**2 + c*(d + e*x)**4)**3
    F = -6*c**2*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(e*(-4*a*c + b**2)**(sympy.S(5)/2)) + 3*c*(b + 2*c*(d + e*x)**2)/(2*e*(-4*a*c + b**2)**2*(a + b*(d + e*x)**2 + c*(d + e*x)**4)) + (-b - 2*c*(d + e*x)**2)/(e*(-16*a*c + 4*b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_627():
    f = (a + b*(d + e*x)**2 + c*(d + e*x)**4)**(-3)
    F = (d/e + x)*(-2*a*c + b**2 + b*c*e**2*(d/e + x)**2)/(4*a*(-4*a*c + b**2)*(a + b*e**2*(d/e + x)**2 + c*e**4*(d/e + x)**4)**2) - 3*sqrt(2)*sqrt(c)*(56*a**2*c**2 - 10*a*b**2*c + b**4 - b*(-8*a*c + b**2)*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(16*a**2*e*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(5)/2)) + 3*sqrt(2)*sqrt(c)*(56*a**2*c**2 - 10*a*b**2*c + b**4 + b*(-8*a*c + b**2)*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(16*a**2*e*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(5)/2)) + (d/e + x)*(3*b*c*e**2*(-8*a*c + b**2)*(d/e + x)**2 + (-7*a*c + b**2)*(-4*a*c + 3*b**2))/(8*a**2*(-4*a*c + b**2)**2*(a + b*e**2*(d/e + x)**2 + c*e**4*(d/e + x)**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_628():
    f = 1/((d + e*x)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**3)
    F = (-2*a*c + b**2 + b*c*(d + e*x)**2)/(4*a*e*(-4*a*c + b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2) + (16*a**2*c**2 - 15*a*b**2*c + 2*b**4 + 2*b*c*(d + e*x)**2*(-7*a*c + b**2))/(4*a**2*e*(-4*a*c + b**2)**2*(a + b*(d + e*x)**2 + c*(d + e*x)**4)) + b*(30*a**2*c**2 - 10*a*b**2*c + b**4)*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(2*a**3*e*(-4*a*c + b**2)**(sympy.S(5)/2)) + log(d + e*x)/(a**3*e) - log(a + b*(d + e*x)**2 + c*(d + e*x)**4)/(4*a**3*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_629():
    f = 1/((d + e*x)**2*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**3)
    F = (-2*a*c + b**2 + b*c*(d + e*x)**2)/(4*a*e*(d + e*x)*(-4*a*c + b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2) + (36*a**2*c**2 - 35*a*b**2*c + 5*b**4 + b*c*(d + e*x)**2*(-32*a*c + 5*b**2))/(8*a**2*e*(d + e*x)*(-4*a*c + b**2)**2*(a + b*(d + e*x)**2 + c*(d + e*x)**4)) - 3*sqrt(2)*sqrt(c)*((-12*a*c + 5*b**2)*(-5*a*c + b**2) - (124*a**2*b*c**2 - 47*a*b**3*c + 5*b**5)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(16*a**3*e*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) - 3*sqrt(2)*sqrt(c)*(b*(124*a**2*c**2 - 47*a*b**2*c + 5*b**4)/sqrt(-4*a*c + b**2) + (-12*a*c + 5*b**2)*(-5*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(16*a**3*e*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) - (-36*a*c + 15*b**2)*(-5*a*c + b**2)/(8*a**3*e*(d + e*x)*(-4*a*c + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_630():
    f = 1/((d + e*x)**3*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**3)
    F = (-2*a*c + b**2 + b*c*(d + e*x)**2)/(4*a*e*(d + e*x)**2*(-4*a*c + b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2) + (20*a**2*c**2 - 20*a*b**2*c + 3*b**4 + 3*b*c*(d + e*x)**2*(-6*a*c + b**2))/(4*a**2*e*(d + e*x)**2*(-4*a*c + b**2)**2*(a + b*(d + e*x)**2 + c*(d + e*x)**4)) - (-15*a*c + 3*b**2)*(-2*a*c + b**2)/(2*a**3*e*(d + e*x)**2*(-4*a*c + b**2)**2) - 3*b*log(d + e*x)/(a**4*e) + 3*b*log(a + b*(d + e*x)**2 + c*(d + e*x)**4)/(4*a**4*e) - (-60*a**3*c**3 + 90*a**2*b**2*c**2 - 30*a*b**4*c + 3*b**6)*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(2*a**4*e*(-4*a*c + b**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_631():
    f = (d*f + e*f*x)**4/(a + b*(d + e*x)**2 + c*(d + e*x)**4)
    F = f**4*x/c - sqrt(2)*f**4*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*e*sqrt(b - sqrt(-4*a*c + b**2))) - sqrt(2)*f**4*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*e*sqrt(b + sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_632():
    f = (d*f + e*f*x)**3/(a + b*(d + e*x)**2 + c*(d + e*x)**4)
    F = b*f**3*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(2*c*e*sqrt(-4*a*c + b**2)) + f**3*log(a + b*(d + e*x)**2 + c*(d + e*x)**4)/(4*c*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_633():
    f = (d*f + e*f*x)**2/(a + b*(d + e*x)**2 + c*(d + e*x)**4)
    F = -sqrt(2)*f**2*sqrt(b - sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(c)*e*sqrt(-4*a*c + b**2)) + sqrt(2)*f**2*sqrt(b + sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(c)*e*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_634():
    f = (d*f + e*f*x)/(a + b*(d + e*x)**2 + c*(d + e*x)**4)
    F = -f*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(e*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_635():
    f = 1/((d*f + e*f*x)*(a + b*(d + e*x)**2 + c*(d + e*x)**4))
    F = b*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(2*a*e*f*sqrt(-4*a*c + b**2)) + log(d + e*x)/(a*e*f) - log(a + b*(d + e*x)**2 + c*(d + e*x)**4)/(4*a*e*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_636():
    f = 1/((d*f + e*f*x)**2*(a + b*(d + e*x)**2 + c*(d + e*x)**4))
    F = -sqrt(2)*sqrt(c)*(-b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(2*a*e*f**2*sqrt(b + sqrt(-4*a*c + b**2))) - sqrt(2)*sqrt(c)*(b/sqrt(-4*a*c + b**2) + 1)*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(2*a*e*f**2*sqrt(b - sqrt(-4*a*c + b**2))) - 1/(a*e*f**2*(d + e*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_637():
    f = 1/((d*f + e*f*x)**3*(a + b*(d + e*x)**2 + c*(d + e*x)**4))
    F = -1/(2*a*e*f**3*(d + e*x)**2) - b*log(d + e*x)/(a**2*e*f**3) + b*log(a + b*(d + e*x)**2 + c*(d + e*x)**4)/(4*a**2*e*f**3) - (-2*a*c + b**2)*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(2*a**2*e*f**3*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_638():
    f = 1/((d*f + e*f*x)**4*(a + b*(d + e*x)**2 + c*(d + e*x)**4))
    F = -1/(3*a*e*f**4*(d + e*x)**3) + b/(a**2*e*f**4*(d + e*x)) + sqrt(2)*sqrt(c)*(b - (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(2*a**2*e*f**4*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*sqrt(c)*(b + (-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(2*a**2*e*f**4*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_639():
    f = (d*f + e*f*x)**4/(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2
    F = f**4*(2*a + b*(d + e*x)**2)*(d + e*x)/(e*(-8*a*c + 2*b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)) + sqrt(2)*f**4*(b - (4*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(4*sqrt(c)*e*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*f**4*(4*a*c + b**2 + b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(4*sqrt(c)*e*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_640():
    f = (d*f + e*f*x)**3/(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2
    F = -b*f**3*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(e*(-4*a*c + b**2)**(sympy.S(3)/2)) + f**3*(2*a + b*(d + e*x)**2)/(e*(-8*a*c + 2*b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_641():
    f = (d*f + e*f*x)**2/(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2
    F = -sqrt(2)*sqrt(c)*f**2*(2*b + sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(2*e*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + sqrt(2)*sqrt(c)*f**2*(2*b - sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(2*e*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) - f**2*(b + 2*c*(d + e*x)**2)*(d + e*x)/(e*(-8*a*c + 2*b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_642():
    f = (d*f + e*f*x)/(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2
    F = 2*c*f*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(e*(-4*a*c + b**2)**(sympy.S(3)/2)) - f*(b + 2*c*(d + e*x)**2)/(e*(-8*a*c + 2*b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_643():
    f = 1/((d*f + e*f*x)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2)
    F = (-2*a*c + b**2 + b*c*(d + e*x)**2)/(2*a*e*f*(-4*a*c + b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)) + b*(-6*a*c + b**2)*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(2*a**2*e*f*(-4*a*c + b**2)**(sympy.S(3)/2)) + log(d + e*x)/(a**2*e*f) - log(a + b*(d + e*x)**2 + c*(d + e*x)**4)/(4*a**2*e*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_644():
    f = 1/((d*f + e*f*x)**2*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2)
    F = (-2*a*c + b**2 + b*c*(d + e*x)**2)/(2*a*e*f**2*(d + e*x)*(-4*a*c + b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)) + sqrt(2)*sqrt(c)*(-16*a*b*c + 3*b**3 - (-10*a*c + 3*b**2)*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a**2*e*f**2*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) - sqrt(2)*sqrt(c)*(-16*a*b*c + 3*b**3 + (-10*a*c + 3*b**2)*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a**2*e*f**2*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) - (-10*a*c + 3*b**2)/(2*a**2*e*f**2*(d + e*x)*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_645():
    f = 1/((d*f + e*f*x)**3*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2)
    F = (-2*a*c + b**2 + b*c*(d + e*x)**2)/(2*a*e*f**3*(d + e*x)**2*(-4*a*c + b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)) - (-3*a*c + b**2)/(a**2*e*f**3*(d + e*x)**2*(-4*a*c + b**2)) - 2*b*log(d + e*x)/(a**3*e*f**3) + b*log(a + b*(d + e*x)**2 + c*(d + e*x)**4)/(2*a**3*e*f**3) - (6*a**2*c**2 - 6*a*b**2*c + b**4)*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(a**3*e*f**3*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_646():
    f = 1/((d*f + e*f*x)**4*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2)
    F = (-2*a*c + b**2 + b*c*(d + e*x)**2)/(2*a*e*f**4*(d + e*x)**3*(-4*a*c + b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)) - (-14*a*c + 5*b**2)/(6*a**2*e*f**4*(d + e*x)**3*(-4*a*c + b**2)) + b*(-19*a*c + 5*b**2)/(2*a**3*e*f**4*(d + e*x)*(-4*a*c + b**2)) - sqrt(2)*sqrt(c)*(28*a**2*c**2 - 29*a*b**2*c + 5*b**4 - b*(-19*a*c + 5*b**2)*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a**3*e*f**4*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + sqrt(2)*sqrt(c)*(28*a**2*c**2 - 29*a*b**2*c + 5*b**4 + b*(-19*a*c + 5*b**2)*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a**3*e*f**4*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_647():
    f = (d*f + e*f*x)**4/(a + b*(d + e*x)**2 + c*(d + e*x)**4)**3
    F = -3*sqrt(2)*sqrt(c)*f**4*(4*a*c + 3*b**2 + 2*b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(8*e*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(5)/2)) + 3*sqrt(2)*sqrt(c)*f**4*(4*a*c + 3*b**2 - 2*b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(8*e*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(5)/2)) + f**4*(2*a + b*(d + e*x)**2)*(d + e*x)/(e*(-16*a*c + 4*b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2) - f**4*(d + e*x)*(-4*a*c + 7*b**2 + 12*b*c*(d + e*x)**2)/(8*e*(-4*a*c + b**2)**2*(a + b*(d + e*x)**2 + c*(d + e*x)**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_648():
    f = (d*f + e*f*x)**3/(a + b*(d + e*x)**2 + c*(d + e*x)**4)**3
    F = 3*b*c*f**3*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(e*(-4*a*c + b**2)**(sympy.S(5)/2)) - 3*b*f**3*(b + 2*c*(d + e*x)**2)/(4*e*(-4*a*c + b**2)**2*(a + b*(d + e*x)**2 + c*(d + e*x)**4)) + f**3*(2*a + b*(d + e*x)**2)/(e*(-16*a*c + 4*b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_649():
    f = (d*f + e*f*x)**2/(a + b*(d + e*x)**2 + c*(d + e*x)**4)**3
    F = -f**2*(b + 2*c*(d + e*x)**2)*(d + e*x)/(e*(-16*a*c + 4*b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2) + sqrt(2)*sqrt(c)*f**2*(20*a*c + b**2 - b*(-52*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(16*a*e*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + sqrt(2)*sqrt(c)*f**2*(20*a*c + b**2 + b*(-52*a*c + b**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(16*a*e*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) + f**2*(d + e*x)*(b*(8*a*c + b**2) + c*(d + e*x)**2*(20*a*c + b**2))/(8*a*e*(-4*a*c + b**2)**2*(a + b*(d + e*x)**2 + c*(d + e*x)**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_650():
    f = (d*f + e*f*x)/(a + b*(d + e*x)**2 + c*(d + e*x)**4)**3
    F = -6*c**2*f*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(e*(-4*a*c + b**2)**(sympy.S(5)/2)) + 3*c*f*(b + 2*c*(d + e*x)**2)/(2*e*(-4*a*c + b**2)**2*(a + b*(d + e*x)**2 + c*(d + e*x)**4)) - f*(b + 2*c*(d + e*x)**2)/(e*(-16*a*c + 4*b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_651():
    f = 1/((d*f + e*f*x)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**3)
    F = (-2*a*c + b**2 + b*c*(d + e*x)**2)/(4*a*e*f*(-4*a*c + b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2) + (16*a**2*c**2 - 15*a*b**2*c + 2*b**4 + 2*b*c*(d + e*x)**2*(-7*a*c + b**2))/(4*a**2*e*f*(-4*a*c + b**2)**2*(a + b*(d + e*x)**2 + c*(d + e*x)**4)) + b*(30*a**2*c**2 - 10*a*b**2*c + b**4)*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(2*a**3*e*f*(-4*a*c + b**2)**(sympy.S(5)/2)) + log(d + e*x)/(a**3*e*f) - log(a + b*(d + e*x)**2 + c*(d + e*x)**4)/(4*a**3*e*f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_652():
    f = 1/((d*f + e*f*x)**2*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**3)
    F = (-2*a*c + b**2 + b*c*(d + e*x)**2)/(4*a*e*f**2*(d + e*x)*(-4*a*c + b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2) + (36*a**2*c**2 - 35*a*b**2*c + 5*b**4 + b*c*(d + e*x)**2*(-32*a*c + 5*b**2))/(8*a**2*e*f**2*(d + e*x)*(-4*a*c + b**2)**2*(a + b*(d + e*x)**2 + c*(d + e*x)**4)) - 3*sqrt(2)*sqrt(c)*((-12*a*c + 5*b**2)*(-5*a*c + b**2) - (124*a**2*b*c**2 - 47*a*b**3*c + 5*b**5)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b + sqrt(-4*a*c + b**2)))/(16*a**3*e*f**2*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) - 3*sqrt(2)*sqrt(c)*(b*(124*a**2*c**2 - 47*a*b**2*c + 5*b**4)/sqrt(-4*a*c + b**2) + (-12*a*c + 5*b**2)*(-5*a*c + b**2))*atan(sqrt(2)*sqrt(c)*(d + e*x)/sqrt(b - sqrt(-4*a*c + b**2)))/(16*a**3*e*f**2*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**2) - (-36*a*c + 15*b**2)*(-5*a*c + b**2)/(8*a**3*e*f**2*(d + e*x)*(-4*a*c + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_653():
    f = 1/((d*f + e*f*x)**3*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**3)
    F = (-2*a*c + b**2 + b*c*(d + e*x)**2)/(4*a*e*f**3*(d + e*x)**2*(-4*a*c + b**2)*(a + b*(d + e*x)**2 + c*(d + e*x)**4)**2) + (20*a**2*c**2 - 20*a*b**2*c + 3*b**4 + 3*b*c*(d + e*x)**2*(-6*a*c + b**2))/(4*a**2*e*f**3*(d + e*x)**2*(-4*a*c + b**2)**2*(a + b*(d + e*x)**2 + c*(d + e*x)**4)) - (-15*a*c + 3*b**2)*(-2*a*c + b**2)/(2*a**3*e*f**3*(d + e*x)**2*(-4*a*c + b**2)**2) - 3*b*log(d + e*x)/(a**4*e*f**3) + 3*b*log(a + b*(d + e*x)**2 + c*(d + e*x)**4)/(4*a**4*e*f**3) - (-60*a**3*c**3 + 90*a**2*b**2*c**2 - 30*a*b**4*c + 3*b**6)*atanh((b + 2*c*(d + e*x)**2)/sqrt(-4*a*c + b**2))/(2*a**4*e*f**3*(-4*a*c + b**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_654():
    f = x/sqrt(a + b*(d + e*x)**3 + c*(d + e*x)**6)
    F = -d*(d + e*x)*sqrt(2*c*(d + e*x)**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*(d + e*x)**3/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(1)/3, sympy.S.Half, sympy.S.Half, sympy.S(4)/3, -2*c*(d + e*x)**3/(b - sqrt(-4*a*c + b**2)), -2*c*(d + e*x)**3/(b + sqrt(-4*a*c + b**2)))/(e**2*sqrt(a + b*(d + e*x)**3 + c*(d + e*x)**6)) + (d + e*x)**2*sqrt(2*c*(d + e*x)**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*(d + e*x)**3/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(2)/3, sympy.S.Half, sympy.S.Half, sympy.S(5)/3, -2*c*(d + e*x)**3/(b - sqrt(-4*a*c + b**2)), -2*c*(d + e*x)**3/(b + sqrt(-4*a*c + b**2)))/(2*e**2*sqrt(a + b*(d + e*x)**3 + c*(d + e*x)**6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_655():
    f = x**2/sqrt(a + b*(d + e*x)**3 + c*(d + e*x)**6)
    F = d**2*(d + e*x)*sqrt(2*c*(d + e*x)**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*(d + e*x)**3/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(1)/3, sympy.S.Half, sympy.S.Half, sympy.S(4)/3, -2*c*(d + e*x)**3/(b - sqrt(-4*a*c + b**2)), -2*c*(d + e*x)**3/(b + sqrt(-4*a*c + b**2)))/(e**3*sqrt(a + b*(d + e*x)**3 + c*(d + e*x)**6)) - d*(d + e*x)**2*sqrt(2*c*(d + e*x)**3/(b - sqrt(-4*a*c + b**2)) + 1)*sqrt(2*c*(d + e*x)**3/(b + sqrt(-4*a*c + b**2)) + 1)*appellf1(sympy.S(2)/3, sympy.S.Half, sympy.S.Half, sympy.S(5)/3, -2*c*(d + e*x)**3/(b - sqrt(-4*a*c + b**2)), -2*c*(d + e*x)**3/(b + sqrt(-4*a*c + b**2)))/(e**3*sqrt(a + b*(d + e*x)**3 + c*(d + e*x)**6)) + atanh((b + 2*c*(d + e*x)**3)/(2*sqrt(c)*sqrt(a + b*(d + e*x)**3 + c*(d + e*x)**6)))/(3*sqrt(c)*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_656():
    f = (3*x + 2)**6*((3*x + 2)**14 + (3*x + 2)**7 + 1)
    F = (3*x + 2)**21/63 + (3*x + 2)**14/42 + (3*x + 2)**7/21
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_3_General_1_2_3_2_d_x_pow_m_a_plus_b_x_pow_n_plus_c_x_pow_2_n_pow_p_657():
    f = (3*x + 2)**6*((3*x + 2)**14 + (3*x + 2)**7 + 1)**2
    F = (3*x + 2)**35/105 + (3*x + 2)**28/42 + (3*x + 2)**21/21 + (3*x + 2)**14/21 + (3*x + 2)**7/21
    assert integrate(f, x) == F

