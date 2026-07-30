"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.1 Binomial products/1.1.2 Quadratic/1.1.2.3 (a+b x^2)^p (c+d x^2)^q.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, m, p, q = symbols('a b c d m p q')

def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_1():
    f = (a + b*x**2)*(c + d*x**2)**4
    F = a*c**4*x + b*d**4*x**11/11 + c**3*x**3*(4*a*d + b*c)/3 + 2*c**2*d*x**5*(3*a*d + 2*b*c)/5 + 2*c*d**2*x**7*(2*a*d + 3*b*c)/7 + d**3*x**9*(a*d + 4*b*c)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_2():
    f = (a + b*x**2)*(c + d*x**2)**3
    F = a*c**3*x + b*d**3*x**9/9 + c**2*x**3*(3*a*d + b*c)/3 + 3*c*d*x**5*(a*d + b*c)/5 + d**2*x**7*(a*d + 3*b*c)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_3():
    f = (a + b*x**2)*(c + d*x**2)**2
    F = a*c**2*x + b*d**2*x**7/7 + c*x**3*(2*a*d + b*c)/3 + d*x**5*(a*d + 2*b*c)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_4():
    f = (a + b*x**2)*(c + d*x**2)
    F = a*c*x + b*d*x**5/5 + x**3*(a*d/3 + b*c/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_5():
    f = (a + b*x**2)/(c + d*x**2)
    F = b*x/d - (-a*d + b*c)*atan(sqrt(d)*x/sqrt(c))/(sqrt(c)*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_6():
    f = (a + b*x**2)/(c + d*x**2)**2
    F = -x*(-a*d + b*c)/(2*c*d*(c + d*x**2)) + (a*d + b*c)*atan(sqrt(d)*x/sqrt(c))/(2*c**(sympy.S(3)/2)*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_7():
    f = (a + b*x**2)/(c + d*x**2)**3
    F = -x*(-a*d + b*c)/(4*c*d*(c + d*x**2)**2) + x*(3*a*d + b*c)/(8*c**2*d*(c + d*x**2)) + (3*a*d + b*c)*atan(sqrt(d)*x/sqrt(c))/(8*c**(sympy.S(5)/2)*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_8():
    f = (a + b*x**2)**2*(c + d*x**2)**3
    F = a**2*c**3*x + a*c**2*x**3*(3*a*d + 2*b*c)/3 + b**2*d**3*x**11/11 + b*d**2*x**9*(2*a*d + 3*b*c)/9 + c*x**5*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2)/5 + d*x**7*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_9():
    f = (a + b*x**2)**2*(c + d*x**2)**2
    F = a**2*c**2*x + 2*a*c*x**3*(a*d + b*c)/3 + b**2*d**2*x**9/9 + 2*b*d*x**7*(a*d + b*c)/7 + x**5*(a**2*d**2/5 + 4*a*b*c*d/5 + b**2*c**2/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_10():
    f = (a + b*x**2)**2*(c + d*x**2)
    F = a**2*c*x + a*x**3*(a*d + 2*b*c)/3 + b**2*d*x**7/7 + b*x**5*(2*a*d + b*c)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_11():
    f = (a + b*x**2)**2/(c + d*x**2)
    F = b**2*x**3/(3*d) - b*x*(-2*a*d + b*c)/d**2 + (-a*d + b*c)**2*atan(sqrt(d)*x/sqrt(c))/(sqrt(c)*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_12():
    f = (a + b*x**2)**2/(c + d*x**2)**2
    F = b**2*x/d**2 + x*(-a*d + b*c)**2/(2*c*d**2*(c + d*x**2)) - (-a*d + b*c)*(a*d + 3*b*c)*atan(sqrt(d)*x/sqrt(c))/(2*c**(sympy.S(3)/2)*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_13():
    f = (a + b*x**2)**2/(c + d*x**2)**3
    F = x*(3*a**2/c**2 - 3*b**2/d**2)/(8*c + 8*d*x**2) - x*(a + b*x**2)*(-a*d + b*c)/(4*c*d*(c + d*x**2)**2) + (3*a**2*d**2 + 2*a*b*c*d + 3*b**2*c**2)*atan(sqrt(d)*x/sqrt(c))/(8*c**(sympy.S(5)/2)*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_14():
    f = (a + b*x**2)**3*(c + d*x**2)**3
    F = a**3*c**3*x + a**2*c**2*x**3*(a*d + b*c) + 3*a*c*x**5*(a**2*d**2 + 3*a*b*c*d + b**2*c**2)/5 + b**3*d**3*x**13/13 + 3*b**2*d**2*x**11*(a*d + b*c)/11 + b*d*x**9*(a**2*d**2 + 3*a*b*c*d + b**2*c**2)/3 + x**7*(a*d/7 + b*c/7)*(a**2*d**2 + 8*a*b*c*d + b**2*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_15():
    f = (a + b*x**2)**3*(c + d*x**2)**2
    F = a**3*c**2*x + a**2*c*x**3*(2*a*d + 3*b*c)/3 + a*x**5*(a**2*d**2 + 6*a*b*c*d + 3*b**2*c**2)/5 + b**3*d**2*x**11/11 + b**2*d*x**9*(3*a*d + 2*b*c)/9 + b*x**7*(3*a**2*d**2 + 6*a*b*c*d + b**2*c**2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_16():
    f = (a + b*x**2)**3*(c + d*x**2)
    F = a**3*c*x + a**2*x**3*(a*d + 3*b*c)/3 + 3*a*b*x**5*(a*d + b*c)/5 + b**3*d*x**9/9 + b**2*x**7*(3*a*d + b*c)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_17():
    f = (a + b*x**2)**3/(c + d*x**2)
    F = b**3*x**5/(5*d) - b**2*x**3*(-3*a*d + b*c)/(3*d**2) + b*x*(3*a**2*d**2 - 3*a*b*c*d + b**2*c**2)/d**3 - (-a*d + b*c)**3*atan(sqrt(d)*x/sqrt(c))/(sqrt(c)*d**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_18():
    f = (a + b*x**2)**3/(c + d*x**2)**2
    F = b**3*x**3/(3*d**2) - b**2*x*(-3*a*d + 2*b*c)/d**3 - x*(-a*d + b*c)**3/(2*c*d**3*(c + d*x**2)) + (-a*d + b*c)**2*(a*d + 5*b*c)*atan(sqrt(d)*x/sqrt(c))/(2*c**(sympy.S(3)/2)*d**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_19():
    f = (a + b*x**2)**3/(c + d*x**2)**3
    F = b**3*x/d**3 - x*(-a*d + b*c)**3/(4*c*d**3*(c + d*x**2)**2) + 3*x*(-a*d + b*c)**2*(a*d + 3*b*c)/(8*c**2*d**3*(c + d*x**2)) - (-3*a*d + 3*b*c)*(4*b**2*c**2 + (a*d + b*c)**2)*atan(sqrt(d)*x/sqrt(c))/(8*c**(sympy.S(5)/2)*d**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_20():
    f = (c + d*x**2)**4/(a + b*x**2)
    F = d**4*x**7/(7*b) + d**3*x**5*(-a*d + 4*b*c)/(5*b**2) + d**2*x**3*(a**2*d**2 - 4*a*b*c*d + 6*b**2*c**2)/(3*b**3) + d*x*(-a*d + 2*b*c)*(a**2*d**2 - 2*a*b*c*d + 2*b**2*c**2)/b**4 + (-a*d + b*c)**4*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*b**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_21():
    f = (c + d*x**2)**3/(a + b*x**2)
    F = d**3*x**5/(5*b) + d**2*x**3*(-a*d + 3*b*c)/(3*b**2) + d*x*(a**2*d**2 - 3*a*b*c*d + 3*b**2*c**2)/b**3 + (-a*d + b*c)**3*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_22():
    f = (c + d*x**2)**2/(a + b*x**2)
    F = d**2*x**3/(3*b) + d*x*(-a*d + 2*b*c)/b**2 + (-a*d + b*c)**2*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_23():
    f = (c + d*x**2)/(a + b*x**2)
    F = d*x/b + (-a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_24():
    f = 1/((a + b*x**2)*(c + d*x**2))
    F = -sqrt(d)*atan(sqrt(d)*x/sqrt(c))/(sqrt(c)*(-a*d + b*c)) + sqrt(b)*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_25():
    f = 1/((a + b*x**2)*(c + d*x**2)**2)
    F = -d*x/(2*c*(c + d*x**2)*(-a*d + b*c)) - sqrt(d)*(-a*d + 3*b*c)*atan(sqrt(d)*x/sqrt(c))/(2*c**(sympy.S(3)/2)*(-a*d + b*c)**2) + b**(sympy.S(3)/2)*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_26():
    f = 1/((a + b*x**2)*(c + d*x**2)**3)
    F = -d*x/(4*c*(c + d*x**2)**2*(-a*d + b*c)) - d*x*(-3*a*d + 7*b*c)/(8*c**2*(c + d*x**2)*(-a*d + b*c)**2) - sqrt(d)*(3*a**2*d**2 - 10*a*b*c*d + 15*b**2*c**2)*atan(sqrt(d)*x/sqrt(c))/(8*c**(sympy.S(5)/2)*(-a*d + b*c)**3) + b**(sympy.S(5)/2)*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_27():
    f = (c + d*x**2)**5/(a + b*x**2)**2
    F = d**5*x**7/(7*b**2) + d**4*x**5*(-2*a*d + 5*b*c)/(5*b**3) + d**3*x**3*(3*a**2*d**2 - 10*a*b*c*d + 10*b**2*c**2)/(3*b**4) + d**2*x*(-4*a**3*d**3 + 15*a**2*b*c*d**2 - 20*a*b**2*c**2*d + 10*b**3*c**3)/b**5 + x*(-a*d + b*c)**5/(2*a*b**5*(a + b*x**2)) + (-a*d + b*c)**4*(9*a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*b**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_28():
    f = (c + d*x**2)**4/(a + b*x**2)**2
    F = d**4*x**5/(5*b**2) + 2*d**3*x**3*(-a*d + 2*b*c)/(3*b**3) + d**2*x*(3*a**2*d**2 - 8*a*b*c*d + 6*b**2*c**2)/b**4 + x*(-a*d + b*c)**4/(2*a*b**4*(a + b*x**2)) + (-a*d + b*c)**3*(7*a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*b**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_29():
    f = (c + d*x**2)**3/(a + b*x**2)**2
    F = d**3*x**3/(3*b**2) + d**2*x*(-2*a*d + 3*b*c)/b**3 + x*(-a*d + b*c)**3/(2*a*b**3*(a + b*x**2)) + (-a*d + b*c)**2*(5*a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_30():
    f = (c + d*x**2)**2/(a + b*x**2)**2
    F = d**2*x/b**2 + x*(-a*d + b*c)**2/(2*a*b**2*(a + b*x**2)) + (-a*d + b*c)*(3*a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_31():
    f = (c + d*x**2)/(a + b*x**2)**2
    F = x*(-a*d + b*c)/(2*a*b*(a + b*x**2)) + (a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_32():
    f = 1/((a + b*x**2)**2*(c + d*x**2))
    F = d**(sympy.S(3)/2)*atan(sqrt(d)*x/sqrt(c))/(sqrt(c)*(-a*d + b*c)**2) + b*x/(2*a*(a + b*x**2)*(-a*d + b*c)) + sqrt(b)*(-3*a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_33():
    f = 1/((a + b*x**2)**2*(c + d*x**2)**2)
    F = d**(sympy.S(3)/2)*(-a*d + 5*b*c)*atan(sqrt(d)*x/sqrt(c))/(2*c**(sympy.S(3)/2)*(-a*d + b*c)**3) + b*x/(2*a*(a + b*x**2)*(c + d*x**2)*(-a*d + b*c)) + d*x*(a*d + b*c)/(2*a*c*(c + d*x**2)*(-a*d + b*c)**2) + b**(sympy.S(3)/2)*(-5*a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_34():
    f = 1/((a + b*x**2)**2*(c + d*x**2)**3)
    F = d**(sympy.S(3)/2)*(3*a**2*d**2 - 14*a*b*c*d + 35*b**2*c**2)*atan(sqrt(d)*x/sqrt(c))/(8*c**(sympy.S(5)/2)*(-a*d + b*c)**4) + b*x/(2*a*(a + b*x**2)*(c + d*x**2)**2*(-a*d + b*c)) + d*x*(a*d + 2*b*c)/(4*a*c*(c + d*x**2)**2*(-a*d + b*c)**2) + d*x*(-a*d + 4*b*c)*(3*a*d + b*c)/(8*a*c**2*(c + d*x**2)*(-a*d + b*c)**3) + b**(sympy.S(5)/2)*(-7*a*d + b*c)*atan(sqrt(b)*x/sqrt(a))/(2*a**(sympy.S(3)/2)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_35():
    f = (c + d*x**2)**5/(a + b*x**2)**3
    F = d**5*x**5/(5*b**3) + d**4*x**3*(-3*a*d + 5*b*c)/(3*b**4) + d**3*x*(6*a**2*d**2 - 15*a*b*c*d + 10*b**2*c**2)/b**5 + x*(-a*d + b*c)**5/(4*a*b**5*(a + b*x**2)**2) + x*(-a*d + b*c)**4*(17*a*d + 3*b*c)/(8*a**2*b**5*(a + b*x**2)) + (-a*d + b*c)**3*(63*a**2*d**2 + 14*a*b*c*d + 3*b**2*c**2)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*b**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_36():
    f = (c + d*x**2)**4/(a + b*x**2)**3
    F = d**4*x**3/(3*b**3) + d**3*x*(-3*a*d + 4*b*c)/b**4 + x*(-a*d + b*c)**4/(4*a*b**4*(a + b*x**2)**2) + x*(-a*d + b*c)**3*(13*a*d + 3*b*c)/(8*a**2*b**4*(a + b*x**2)) + (-a*d + b*c)**2*(35*a**2*d**2 + 10*a*b*c*d + 3*b**2*c**2)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*b**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_37():
    f = (c + d*x**2)**3/(a + b*x**2)**3
    F = d**3*x/b**3 + x*(-a*d + b*c)**3/(4*a*b**3*(a + b*x**2)**2) + 3*x*(-a*d + b*c)**2*(3*a*d + b*c)/(8*a**2*b**3*(a + b*x**2)) + (-3*a*d + 3*b*c)*(4*a**2*d**2 + (a*d + b*c)**2)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_38():
    f = (c + d*x**2)**2/(a + b*x**2)**3
    F = x*(-3*d**2/b**2 + 3*c**2/a**2)/(8*a + 8*b*x**2) + x*(c + d*x**2)*(-a*d + b*c)/(4*a*b*(a + b*x**2)**2) + (3*a**2*d**2 + 2*a*b*c*d + 3*b**2*c**2)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_39():
    f = (c + d*x**2)/(a + b*x**2)**3
    F = x*(-a*d + b*c)/(4*a*b*(a + b*x**2)**2) + x*(a*d + 3*b*c)/(8*a**2*b*(a + b*x**2)) + (a*d + 3*b*c)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_40():
    f = 1/((a + b*x**2)**3*(c + d*x**2))
    F = -d**(sympy.S(5)/2)*atan(sqrt(d)*x/sqrt(c))/(sqrt(c)*(-a*d + b*c)**3) + b*x/(4*a*(a + b*x**2)**2*(-a*d + b*c)) + b*x*(-7*a*d + 3*b*c)/(8*a**2*(a + b*x**2)*(-a*d + b*c)**2) + sqrt(b)*(15*a**2*d**2 - 10*a*b*c*d + 3*b**2*c**2)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_41():
    f = 1/((a + b*x**2)**3*(c + d*x**2)**2)
    F = -d**(sympy.S(5)/2)*(-a*d + 7*b*c)*atan(sqrt(d)*x/sqrt(c))/(2*c**(sympy.S(3)/2)*(-a*d + b*c)**4) + b*x/(4*a*(a + b*x**2)**2*(c + d*x**2)*(-a*d + b*c)) + 3*b*x*(-3*a*d + b*c)/(8*a**2*(a + b*x**2)*(c + d*x**2)*(-a*d + b*c)**2) + d*x*(-4*a*d + b*c)*(a*d + 3*b*c)/(8*a**2*c*(c + d*x**2)*(-a*d + b*c)**3) + b**(sympy.S(3)/2)*(35*a**2*d**2 - 14*a*b*c*d + 3*b**2*c**2)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_42():
    f = 1/((a + b*x**2)**3*(c + d*x**2)**3)
    F = -3*d**(sympy.S(5)/2)*(a**2*d**2 - 6*a*b*c*d + 21*b**2*c**2)*atan(sqrt(d)*x/sqrt(c))/(8*c**(sympy.S(5)/2)*(-a*d + b*c)**5) + b*x/(4*a*(a + b*x**2)**2*(c + d*x**2)**2*(-a*d + b*c)) + b*x*(-11*a*d + 3*b*c)/(8*a**2*(a + b*x**2)*(c + d*x**2)**2*(-a*d + b*c)**2) + d*x*(-2*a**2*d**2 - 13*a*b*c*d + 3*b**2*c**2)/(8*a**2*c*(c + d*x**2)**2*(-a*d + b*c)**3) + 3*d*x*(a*d + b*c)*(a**2*d**2 - 6*a*b*c*d + b**2*c**2)/(8*a**2*c**2*(c + d*x**2)*(-a*d + b*c)**4) + 3*b**(sympy.S(5)/2)*(21*a**2*d**2 - 6*a*b*c*d + b**2*c**2)*atan(sqrt(b)*x/sqrt(a))/(8*a**(sympy.S(5)/2)*(-a*d + b*c)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_43():
    f = (x**2 - 1)**3/(x**2 + 1)**4
    F = -x*(1 - x**2)**2/(3*(x**2 + 1)**3) - 2*x/(3*x**2 + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_44():
    f = (x**2 - 1)**4/(x**2 + 1)**5
    F = x*(1 - x**2)**3/(4*(x**2 + 1)**4) + 3*x*(1 - x**2)/(8*(x**2 + 1)**2) + 3*atan(x)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_45():
    f = sqrt(a + b*x**2)*(c + d*x**2)**3
    F = a*(-5*a**3*d**3 + 24*a**2*b*c*d**2 - 48*a*b**2*c**2*d + 64*b**3*c**3)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(128*b**(sympy.S(7)/2)) + d*x*(a + b*x**2)**(sympy.S(3)/2)*(c + d*x**2)**2/(8*b) + d*x*(a + b*x**2)**(sympy.S(3)/2)*(c + d*x**2)*(-5*a*d + 12*b*c)/(48*b**2) + d*x*(a + b*x**2)**(sympy.S(3)/2)*(15*a**2*d**2 - 52*a*b*c*d + 72*b**2*c**2)/(192*b**3) + x*sqrt(a + b*x**2)*(-5*a**3*d**3 + 24*a**2*b*c*d**2 - 48*a*b**2*c**2*d + 64*b**3*c**3)/(128*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_46():
    f = sqrt(a + b*x**2)*(c + d*x**2)**2
    F = a*(a**2*d**2 - 4*a*b*c*d + 8*b**2*c**2)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(16*b**(sympy.S(5)/2)) + d*x*(a + b*x**2)**(sympy.S(3)/2)*(c + d*x**2)/(6*b) + d*x*(a + b*x**2)**(sympy.S(3)/2)*(-3*a*d + 8*b*c)/(24*b**2) + x*sqrt(a + b*x**2)*(a**2*d**2 - 4*a*b*c*d + 8*b**2*c**2)/(16*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_47():
    f = sqrt(a + b*x**2)*(c + d*x**2)
    F = a*(-a*d + 4*b*c)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*b**(sympy.S(3)/2)) + d*x*(a + b*x**2)**(sympy.S(3)/2)/(4*b) + x*sqrt(a + b*x**2)*(-a*d + 4*b*c)/(8*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_48():
    f = sqrt(a + b*x**2)
    F = a*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*sqrt(b)) + x*sqrt(a + b*x**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_49():
    f = sqrt(a + b*x**2)/(c + d*x**2)
    F = sqrt(b)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/d - sqrt(-a*d + b*c)*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(sqrt(c)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_50():
    f = sqrt(a + b*x**2)/(c + d*x**2)**2
    F = a*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(2*c**(sympy.S(3)/2)*sqrt(-a*d + b*c)) + x*sqrt(a + b*x**2)/(2*c*(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_51():
    f = sqrt(a + b*x**2)/(c + d*x**2)**3
    F = a*(-3*a*d + 4*b*c)*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(8*c**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(3)/2)) - d*x*(a + b*x**2)**(sympy.S(3)/2)/(4*c*(c + d*x**2)**2*(-a*d + b*c)) + x*sqrt(a + b*x**2)*(-3*a*d + 4*b*c)/(8*c**2*(c + d*x**2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_52():
    f = sqrt(a + b*x**2)/(c + d*x**2)**4
    F = a*(5*a**2*d**2 - 12*a*b*c*d + 8*b**2*c**2)*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(16*c**(sympy.S(7)/2)*(-a*d + b*c)**(sympy.S(5)/2)) + x*sqrt(a + b*x**2)/(6*c*(c + d*x**2)**3) + x*sqrt(a + b*x**2)*(-5*a*d + 4*b*c)/(24*c**2*(c + d*x**2)**2*(-a*d + b*c)) + x*sqrt(a + b*x**2)*(-5*a*d + 2*b*c)*(-3*a*d + 4*b*c)/(48*c**3*(c + d*x**2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_53():
    f = (a + b*x**2)**(sympy.S(3)/2)*(c + d*x**2)**3
    F = 3*a**2*(-a*d + 4*b*c)*(a**2*d**2 - 2*a*b*c*d + 8*b**2*c**2)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(256*b**(sympy.S(7)/2)) + 3*a*x*sqrt(a + b*x**2)*(-a*d + 4*b*c)*(a**2*d**2 - 2*a*b*c*d + 8*b**2*c**2)/(256*b**3) + d*x*(a + b*x**2)**(sympy.S(5)/2)*(c + d*x**2)**2/(10*b) + d*x*(a + b*x**2)**(sympy.S(5)/2)*(c + d*x**2)*(-5*a*d + 14*b*c)/(80*b**2) + d*x*(a + b*x**2)**(sympy.S(5)/2)*(5*a**2*d**2 - 20*a*b*c*d + 36*b**2*c**2)/(160*b**3) + x*(a + b*x**2)**(sympy.S(3)/2)*(-a*d + 4*b*c)*(a**2*d**2 - 2*a*b*c*d + 8*b**2*c**2)/(128*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_54():
    f = (a + b*x**2)**(sympy.S(3)/2)*(c + d*x**2)**2
    F = a**2*(3*a**2*d**2 - 16*a*b*c*d + 48*b**2*c**2)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(128*b**(sympy.S(5)/2)) + a*x*sqrt(a + b*x**2)*(3*a**2*d**2 - 16*a*b*c*d + 48*b**2*c**2)/(128*b**2) + d*x*(a + b*x**2)**(sympy.S(5)/2)*(c + d*x**2)/(8*b) + d*x*(a + b*x**2)**(sympy.S(5)/2)*(-3*a*d + 10*b*c)/(48*b**2) + x*(a + b*x**2)**(sympy.S(3)/2)*(3*a**2*d**2 - 16*a*b*c*d + 48*b**2*c**2)/(192*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_55():
    f = (a + b*x**2)**(sympy.S(3)/2)*(c + d*x**2)
    F = a**2*(-a*d + 6*b*c)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(16*b**(sympy.S(3)/2)) + a*x*sqrt(a + b*x**2)*(-a*d + 6*b*c)/(16*b) + d*x*(a + b*x**2)**(sympy.S(5)/2)/(6*b) + x*(a + b*x**2)**(sympy.S(3)/2)*(-a*d + 6*b*c)/(24*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_56():
    f = (a + b*x**2)**(sympy.S(3)/2)
    F = 3*a**2*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*sqrt(b)) + 3*a*x*sqrt(a + b*x**2)/8 + x*(a + b*x**2)**(sympy.S(3)/2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_57():
    f = (a + b*x**2)**(sympy.S(3)/2)/(c + d*x**2)
    F = -sqrt(b)*(-3*a*d + 2*b*c)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*d**2) + b*x*sqrt(a + b*x**2)/(2*d) + (-a*d + b*c)**(sympy.S(3)/2)*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(sqrt(c)*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_58():
    f = (a + b*x**2)**(sympy.S(3)/2)/(c + d*x**2)**2
    F = b**(sympy.S(3)/2)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/d**2 - x*sqrt(a + b*x**2)*(-a*d + b*c)/(2*c*d*(c + d*x**2)) - sqrt(-a*d + b*c)*(a*d + 2*b*c)*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(2*c**(sympy.S(3)/2)*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_59():
    f = (a + b*x**2)**(sympy.S(3)/2)/(c + d*x**2)**3
    F = 3*a**2*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(8*c**(sympy.S(5)/2)*sqrt(-a*d + b*c)) + 3*a*x*sqrt(a + b*x**2)/(8*c**2*(c + d*x**2)) + x*(a + b*x**2)**(sympy.S(3)/2)/(4*c*(c + d*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_60():
    f = (a + b*x**2)**(sympy.S(3)/2)/(c + d*x**2)**4
    F = a**2*(-5*a*d + 6*b*c)*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(16*c**(sympy.S(7)/2)*(-a*d + b*c)**(sympy.S(3)/2)) + a*x*sqrt(a + b*x**2)*(-5*a*d + 6*b*c)/(16*c**3*(c + d*x**2)*(-a*d + b*c)) - d*x*(a + b*x**2)**(sympy.S(5)/2)/(6*c*(c + d*x**2)**3*(-a*d + b*c)) + x*(a + b*x**2)**(sympy.S(3)/2)*(-5*a*d + 6*b*c)/(24*c**2*(c + d*x**2)**2*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_61():
    f = (a + b*x**2)**(sympy.S(3)/2)/(c + d*x**2)**5
    F = a**2*(35*a**2*d**2 - 80*a*b*c*d + 48*b**2*c**2)*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(128*c**(sympy.S(9)/2)*(-a*d + b*c)**(sympy.S(5)/2)) - x*sqrt(a + b*x**2)*(-a*d + b*c)/(8*c*d*(c + d*x**2)**4) + x*sqrt(a + b*x**2)*(7*a*d + 2*b*c)/(48*c**2*d*(c + d*x**2)**3) + x*sqrt(a + b*x**2)*(-35*a**2*d**2 + 24*a*b*c*d + 8*b**2*c**2)/(192*c**3*d*(c + d*x**2)**2*(-a*d + b*c)) + x*sqrt(a + b*x**2)*(105*a**3*d**3 - 170*a**2*b*c*d**2 + 40*a*b**2*c**2*d + 16*b**3*c**3)/(384*c**4*d*(c + d*x**2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_62():
    f = (a + b*x**2)**(sympy.S(5)/2)*(c + d*x**2)**3
    F = a**3*(-5*a**3*d**3 + 36*a**2*b*c*d**2 - 120*a*b**2*c**2*d + 320*b**3*c**3)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(1024*b**(sympy.S(7)/2)) + a**2*x*sqrt(a + b*x**2)*(-5*a**3*d**3 + 36*a**2*b*c*d**2 - 120*a*b**2*c**2*d + 320*b**3*c**3)/(1024*b**3) + a*x*(a + b*x**2)**(sympy.S(3)/2)*(-5*a**3*d**3 + 36*a**2*b*c*d**2 - 120*a*b**2*c**2*d + 320*b**3*c**3)/(1536*b**3) + d*x*(a + b*x**2)**(sympy.S(7)/2)*(c + d*x**2)**2/(12*b) + d*x*(a + b*x**2)**(sympy.S(7)/2)*(c + d*x**2)*(-5*a*d + 16*b*c)/(120*b**2) + d*x*(a + b*x**2)**(sympy.S(7)/2)*(15*a**2*d**2 - 68*a*b*c*d + 152*b**2*c**2)/(960*b**3) + x*(a + b*x**2)**(sympy.S(5)/2)*(-5*a**3*d**3 + 36*a**2*b*c*d**2 - 120*a*b**2*c**2*d + 320*b**3*c**3)/(1920*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_63():
    f = (a + b*x**2)**(sympy.S(5)/2)*(c + d*x**2)**2
    F = a**3*(3*a**2*d**2 - 20*a*b*c*d + 80*b**2*c**2)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(256*b**(sympy.S(5)/2)) + a**2*x*sqrt(a + b*x**2)*(3*a**2*d**2 - 20*a*b*c*d + 80*b**2*c**2)/(256*b**2) + a*x*(a + b*x**2)**(sympy.S(3)/2)*(3*a**2*d**2 - 20*a*b*c*d + 80*b**2*c**2)/(384*b**2) + d*x*(a + b*x**2)**(sympy.S(7)/2)*(c + d*x**2)/(10*b) + 3*d*x*(a + b*x**2)**(sympy.S(7)/2)*(-a*d + 4*b*c)/(80*b**2) + x*(a + b*x**2)**(sympy.S(5)/2)*(3*a**2*d**2 - 20*a*b*c*d + 80*b**2*c**2)/(480*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_64():
    f = (a + b*x**2)**(sympy.S(5)/2)*(c + d*x**2)
    F = 5*a**3*(-a*d + 8*b*c)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(128*b**(sympy.S(3)/2)) + 5*a**2*x*sqrt(a + b*x**2)*(-a*d + 8*b*c)/(128*b) + 5*a*x*(a + b*x**2)**(sympy.S(3)/2)*(-a*d + 8*b*c)/(192*b) + d*x*(a + b*x**2)**(sympy.S(7)/2)/(8*b) + x*(a + b*x**2)**(sympy.S(5)/2)*(-a*d + 8*b*c)/(48*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_65():
    f = (a + b*x**2)**(sympy.S(5)/2)
    F = 5*a**3*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(16*sqrt(b)) + 5*a**2*x*sqrt(a + b*x**2)/16 + 5*a*x*(a + b*x**2)**(sympy.S(3)/2)/24 + x*(a + b*x**2)**(sympy.S(5)/2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_66():
    f = (a + b*x**2)**(sympy.S(5)/2)/(c + d*x**2)
    F = sqrt(b)*(15*a**2*d**2 - 20*a*b*c*d + 8*b**2*c**2)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*d**3) + b*x*(a + b*x**2)**(sympy.S(3)/2)/(4*d) - b*x*sqrt(a + b*x**2)*(-7*a*d + 4*b*c)/(8*d**2) - (-a*d + b*c)**(sympy.S(5)/2)*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(sqrt(c)*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_67():
    f = (a + b*x**2)**(sympy.S(5)/2)/(c + d*x**2)**2
    F = -b**(sympy.S(3)/2)*(-5*a*d + 4*b*c)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*d**3) + b*x*sqrt(a + b*x**2)*(-a*d + 2*b*c)/(2*c*d**2) - x*(a + b*x**2)**(sympy.S(3)/2)*(-a*d + b*c)/(2*c*d*(c + d*x**2)) + (-a*d + b*c)**(sympy.S(3)/2)*(a*d + 4*b*c)*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(2*c**(sympy.S(3)/2)*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_68():
    f = (a + b*x**2)**(sympy.S(5)/2)/(c + d*x**2)**3
    F = b**(sympy.S(5)/2)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/d**3 - x*(a + b*x**2)**(sympy.S(3)/2)*(-a*d + b*c)/(4*c*d*(c + d*x**2)**2) - x*sqrt(a + b*x**2)*(-a*d + b*c)*(3*a*d + 4*b*c)/(8*c**2*d**2*(c + d*x**2)) - sqrt(-a*d + b*c)*(3*a**2*d**2 + 4*a*b*c*d + 8*b**2*c**2)*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(8*c**(sympy.S(5)/2)*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_69():
    f = (a + b*x**2)**(sympy.S(5)/2)/(c + d*x**2)**4
    F = 5*a**3*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(16*c**(sympy.S(7)/2)*sqrt(-a*d + b*c)) + 5*a**2*x*sqrt(a + b*x**2)/(16*c**3*(c + d*x**2)) + 5*a*x*(a + b*x**2)**(sympy.S(3)/2)/(24*c**2*(c + d*x**2)**2) + x*(a + b*x**2)**(sympy.S(5)/2)/(6*c*(c + d*x**2)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_70():
    f = (a + b*x**2)**(sympy.S(5)/2)/(c + d*x**2)**5
    F = 5*a**3*(-7*a*d + 8*b*c)*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(128*c**(sympy.S(9)/2)*(-a*d + b*c)**(sympy.S(3)/2)) + 5*a**2*x*sqrt(a + b*x**2)*(-7*a*d + 8*b*c)/(128*c**4*(c + d*x**2)*(-a*d + b*c)) + 5*a*x*(a + b*x**2)**(sympy.S(3)/2)*(-7*a*d + 8*b*c)/(192*c**3*(c + d*x**2)**2*(-a*d + b*c)) - d*x*(a + b*x**2)**(sympy.S(7)/2)/(8*c*(c + d*x**2)**4*(-a*d + b*c)) + x*(a + b*x**2)**(sympy.S(5)/2)*(-7*a*d + 8*b*c)/(48*c**2*(c + d*x**2)**3*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_71():
    f = sqrt(1 - x**2)/(x**2 + 1)
    F = -asin(x) + sqrt(2)*atan(sqrt(2)*x/sqrt(1 - x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_72():
    f = sqrt(x**2 + 1)/(x**2 - 1)
    F = asinh(x) - sqrt(2)*atanh(sqrt(2)*x/sqrt(x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_73():
    f = sqrt(1 - x**2)/(2*x**2 - 1)
    F = -asin(x)/2 - atanh(x/sqrt(1 - x**2))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_74():
    f = (c + d*x**2)**3/sqrt(a + b*x**2)
    F = d*x*sqrt(a + b*x**2)*(c + d*x**2)**2/(6*b) + 5*d*x*sqrt(a + b*x**2)*(c + d*x**2)*(-a*d + 2*b*c)/(24*b**2) + d*x*sqrt(a + b*x**2)*(15*a**2*d**2 - 44*a*b*c*d + 44*b**2*c**2)/(48*b**3) + (-a*d + 2*b*c)*(5*a**2*d**2 - 8*a*b*c*d + 8*b**2*c**2)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(16*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_75():
    f = (c + d*x**2)**2/sqrt(a + b*x**2)
    F = d*x*sqrt(a + b*x**2)*(c + d*x**2)/(4*b) + 3*d*x*sqrt(a + b*x**2)*(-a*d + 2*b*c)/(8*b**2) + (3*a**2*d**2 - 8*a*b*c*d + 8*b**2*c**2)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_76():
    f = (c + d*x**2)/sqrt(a + b*x**2)
    F = d*x*sqrt(a + b*x**2)/(2*b) + (-a*d + 2*b*c)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_77():
    f = 1/sqrt(a + b*x**2)
    F = atanh(sqrt(b)*x/sqrt(a + b*x**2))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_78():
    f = 1/(sqrt(a + b*x**2)*(c + d*x**2))
    F = atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(sqrt(c)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_79():
    f = 1/(sqrt(a + b*x**2)*(c + d*x**2)**2)
    F = -d*x*sqrt(a + b*x**2)/(2*c*(c + d*x**2)*(-a*d + b*c)) + (-a*d + 2*b*c)*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(2*c**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_80():
    f = 1/(sqrt(a + b*x**2)*(c + d*x**2)**3)
    F = -d*x*sqrt(a + b*x**2)/(4*c*(c + d*x**2)**2*(-a*d + b*c)) - 3*d*x*sqrt(a + b*x**2)*(-a*d + 2*b*c)/(8*c**2*(c + d*x**2)*(-a*d + b*c)**2) + (3*a**2*d**2 - 8*a*b*c*d + 8*b**2*c**2)*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(8*c**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_81():
    f = (c + d*x**2)**4/(a + b*x**2)**(sympy.S(3)/2)
    F = d*(-35*a**3*d**3 + 120*a**2*b*c*d**2 - 144*a*b**2*c**2*d + 64*b**3*c**3)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(16*b**(sympy.S(9)/2)) + x*(c + d*x**2)**3*(-a*d + b*c)/(a*b*sqrt(a + b*x**2)) - d*x*sqrt(a + b*x**2)*(c + d*x**2)**2*(-7*a*d + 6*b*c)/(6*a*b**2) - d*x*sqrt(a + b*x**2)*(c + d*x**2)*(35*a**2*d**2 - 64*a*b*c*d + 24*b**2*c**2)/(24*a*b**3) - d*x*sqrt(a + b*x**2)*(-105*a**3*d**3 + 290*a**2*b*c*d**2 - 248*a*b**2*c**2*d + 48*b**3*c**3)/(48*a*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_82():
    f = (c + d*x**2)**3/(a + b*x**2)**(sympy.S(3)/2)
    F = 3*d*(5*a**2*d**2 - 12*a*b*c*d + 8*b**2*c**2)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*b**(sympy.S(7)/2)) + x*(c + d*x**2)**2*(-a*d + b*c)/(a*b*sqrt(a + b*x**2)) - d*x*sqrt(a + b*x**2)*(c + d*x**2)*(-5*a*d + 4*b*c)/(4*a*b**2) - d*x*sqrt(a + b*x**2)*(-5*a*d + 2*b*c)*(-3*a*d + 4*b*c)/(8*a*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_83():
    f = (c + d*x**2)/(a + b*x**2)**(sympy.S(3)/2)
    F = d*atanh(sqrt(b)*x/sqrt(a + b*x**2))/b**(sympy.S(3)/2) + x*(-a*d + b*c)/(a*b*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_84():
    f = (a + b*x**2)**(sympy.S(-3)/2)
    F = x/(a*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_85():
    f = 1/((a + b*x**2)**(sympy.S(3)/2)*(c + d*x**2))
    F = -d*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(sqrt(c)*(-a*d + b*c)**(sympy.S(3)/2)) + b*x/(a*sqrt(a + b*x**2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_86():
    f = 1/((a + b*x**2)**(sympy.S(3)/2)*(c + d*x**2)**2)
    F = -d*x/(2*c*sqrt(a + b*x**2)*(c + d*x**2)*(-a*d + b*c)) - d*(-a*d + 4*b*c)*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(2*c**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(5)/2)) + b*x*(a*d + 2*b*c)/(2*a*c*sqrt(a + b*x**2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_87():
    f = 1/((a + b*x**2)**(sympy.S(3)/2)*(c + d*x**2)**3)
    F = -d*x/(4*c*sqrt(a + b*x**2)*(c + d*x**2)**2*(-a*d + b*c)) - 3*d*(a**2*d**2 - 4*a*b*c*d + 8*b**2*c**2)*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(8*c**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(7)/2)) + b*x*(a*d + 4*b*c)/(4*a*c*sqrt(a + b*x**2)*(c + d*x**2)*(-a*d + b*c)**2) + d*x*sqrt(a + b*x**2)*(-a*d + 4*b*c)*(3*a*d + 2*b*c)/(8*a*c**2*(c + d*x**2)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_88():
    f = (c + d*x**2)**4/(a + b*x**2)**(sympy.S(5)/2)
    F = d**2*(35*a**2*d**2 - 80*a*b*c*d + 48*b**2*c**2)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*b**(sympy.S(9)/2)) + x*(c + d*x**2)**3*(-a*d + b*c)/(3*a*b*(a + b*x**2)**(sympy.S(3)/2)) + x*(c + d*x**2)**2*(-a*d + b*c)*(7*a*d + 2*b*c)/(3*a**2*b**2*sqrt(a + b*x**2)) - d*x*sqrt(a + b*x**2)*(c + d*x**2)*(-35*a**2*d**2 + 24*a*b*c*d + 8*b**2*c**2)/(12*a**2*b**3) - d*x*sqrt(a + b*x**2)*(105*a**3*d**3 - 170*a**2*b*c*d**2 + 40*a*b**2*c**2*d + 16*b**3*c**3)/(24*a**2*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_89():
    f = (c + d*x**2)**3/(a + b*x**2)**(sympy.S(5)/2)
    F = d**2*(-5*a*d + 6*b*c)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*b**(sympy.S(7)/2)) + x*(c + d*x**2)**2*(-a*d + b*c)/(3*a*b*(a + b*x**2)**(sympy.S(3)/2)) + x*(c + d*x**2)*(-a*d + b*c)*(5*a*d + 2*b*c)/(3*a**2*b**2*sqrt(a + b*x**2)) - d*x*sqrt(a + b*x**2)*(-15*a**2*d**2 + 8*a*b*c*d + 4*b**2*c**2)/(6*a**2*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_90():
    f = (c + d*x**2)**2/(a + b*x**2)**(sympy.S(5)/2)
    F = d**2*atanh(sqrt(b)*x/sqrt(a + b*x**2))/b**(sympy.S(5)/2) + x*(c + d*x**2)*(-a*d + b*c)/(3*a*b*(a + b*x**2)**(sympy.S(3)/2)) + x*(-a*d + b*c)*(3*a*d + 2*b*c)/(3*a**2*b**2*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_91():
    f = (c + d*x**2)/(a + b*x**2)**(sympy.S(5)/2)
    F = x*(c + d*x**2)/(3*a*(a + b*x**2)**(sympy.S(3)/2)) + 2*c*x/(3*a**2*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_92():
    f = (a + b*x**2)**(sympy.S(-5)/2)
    F = x/(3*a*(a + b*x**2)**(sympy.S(3)/2)) + 2*x/(3*a**2*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_93():
    f = 1/((a + b*x**2)**(sympy.S(5)/2)*(c + d*x**2))
    F = d**2*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(sqrt(c)*(-a*d + b*c)**(sympy.S(5)/2)) + b*x/(3*a*(a + b*x**2)**(sympy.S(3)/2)*(-a*d + b*c)) + b*x*(-5*a*d + 2*b*c)/(3*a**2*sqrt(a + b*x**2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_94():
    f = 1/((a + b*x**2)**(sympy.S(5)/2)*(c + d*x**2)**2)
    F = -d*x/(2*c*(a + b*x**2)**(sympy.S(3)/2)*(c + d*x**2)*(-a*d + b*c)) + d**2*(-a*d + 6*b*c)*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(2*c**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(7)/2)) + b*x*(3*a*d + 2*b*c)/(6*a*c*(a + b*x**2)**(sympy.S(3)/2)*(-a*d + b*c)**2) + b*x*(-3*a**2*d**2 - 16*a*b*c*d + 4*b**2*c**2)/(6*a**2*c*sqrt(a + b*x**2)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_95():
    f = 1/((a + b*x**2)**(sympy.S(5)/2)*(c + d*x**2)**3)
    F = -d*x/(4*c*(a + b*x**2)**(sympy.S(3)/2)*(c + d*x**2)**2*(-a*d + b*c)) + d**2*(3*a**2*d**2 - 16*a*b*c*d + 48*b**2*c**2)*atanh(x*sqrt(-a*d + b*c)/(sqrt(c)*sqrt(a + b*x**2)))/(8*c**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(9)/2)) + b*x*(3*a*d + 4*b*c)/(12*a*c*(a + b*x**2)**(sympy.S(3)/2)*(c + d*x**2)*(-a*d + b*c)**2) + b*x*(-3*a**2*d**2 - 40*a*b*c*d + 8*b**2*c**2)/(12*a**2*c*sqrt(a + b*x**2)*(c + d*x**2)*(-a*d + b*c)**3) + d*x*sqrt(a + b*x**2)*(9*a**3*d**3 - 42*a**2*b*c*d**2 - 88*a*b**2*c**2*d + 16*b**3*c**3)/(24*a**2*c**2*(c + d*x**2)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_96():
    f = (a + b*x**2)**3/(c + d*x**2)**(sympy.S(11)/2)
    F = 16*a**3*x*(-8*a*d + 9*b*c)/(315*c**5*sqrt(c + d*x**2)*(-a*d + b*c)) + 8*a**2*x*(a + b*x**2)*(-8*a*d + 9*b*c)/(315*c**4*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)) + 2*a*x*(a + b*x**2)**2*(-8*a*d + 9*b*c)/(105*c**3*(c + d*x**2)**(sympy.S(5)/2)*(-a*d + b*c)) - d*x*(a + b*x**2)**4/(9*c*(c + d*x**2)**(sympy.S(9)/2)*(-a*d + b*c)) + x*(a + b*x**2)**3*(-8*a*d + 9*b*c)/(63*c**2*(c + d*x**2)**(sympy.S(7)/2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_97():
    f = (a + b*x**2)**2/(c + d*x**2)**(sympy.S(9)/2)
    F = 8*a**2*x*(-6*a*d + 7*b*c)/(105*c**4*sqrt(c + d*x**2)*(-a*d + b*c)) + 4*a*x*(a + b*x**2)*(-6*a*d + 7*b*c)/(105*c**3*(c + d*x**2)**(sympy.S(3)/2)*(-a*d + b*c)) - d*x*(a + b*x**2)**3/(7*c*(c + d*x**2)**(sympy.S(7)/2)*(-a*d + b*c)) + x*(a + b*x**2)**2*(-6*a*d + 7*b*c)/(35*c**2*(c + d*x**2)**(sympy.S(5)/2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_98():
    f = (a + b*x**2)/(c + d*x**2)**(sympy.S(7)/2)
    F = -x*(-a*d + b*c)/(5*c*d*(c + d*x**2)**(sympy.S(5)/2)) + x*(4*a*d + b*c)/(15*c**2*d*(c + d*x**2)**(sympy.S(3)/2)) + x*(8*a*d + 2*b*c)/(15*c**3*d*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_99():
    f = (c + d*x**2)**(sympy.S(-5)/2)
    F = x/(3*c*(c + d*x**2)**(sympy.S(3)/2)) + 2*x/(3*c**2*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_100():
    f = 1/((a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2))
    F = -d*x/(c*sqrt(c + d*x**2)*(-a*d + b*c)) + b*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(sqrt(a)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_101():
    f = 1/((a + b*x**2)**2*sqrt(c + d*x**2))
    F = b*x*sqrt(c + d*x**2)/(2*a*(a + b*x**2)*(-a*d + b*c)) + (-2*a*d + b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(2*a**(sympy.S(3)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_102():
    f = sqrt(c + d*x**2)/(a + b*x**2)**3
    F = b*x*(c + d*x**2)**(sympy.S(3)/2)/(4*a*(a + b*x**2)**2*(-a*d + b*c)) + x*sqrt(c + d*x**2)*(-4*a*d + 3*b*c)/(8*a**2*(a + b*x**2)*(-a*d + b*c)) + c*(-4*a*d + 3*b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(8*a**(sympy.S(5)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_103():
    f = (c + d*x**2)**(sympy.S(3)/2)/(a + b*x**2)**4
    F = b*x*(c + d*x**2)**(sympy.S(5)/2)/(6*a*(a + b*x**2)**3*(-a*d + b*c)) + x*(c + d*x**2)**(sympy.S(3)/2)*(-6*a*d + 5*b*c)/(24*a**2*(a + b*x**2)**2*(-a*d + b*c)) + c*x*sqrt(c + d*x**2)*(-6*a*d + 5*b*c)/(16*a**3*(a + b*x**2)*(-a*d + b*c)) + c**2*(-6*a*d + 5*b*c)*atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(16*a**(sympy.S(7)/2)*(-a*d + b*c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_104():
    f = 1/(sqrt(c + d*x**2)*(b*c/d + b*x**2))
    F = d*x/(b*c*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_105():
    f = 1/(sqrt(1 - x**2)*(x**2 + 1))
    F = sqrt(2)*atan(sqrt(2)*x/sqrt(1 - x**2))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_106():
    f = 1/((a + b*x**2)*sqrt(c + d*x**2))
    F = atan(x*sqrt(-a*d + b*c)/(sqrt(a)*sqrt(c + d*x**2)))/(sqrt(a)*sqrt(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_107():
    f = (x**2 - 1)/(x**2 + 1)**(sympy.S(3)/2)
    F = -2*x/sqrt(x**2 + 1) + asinh(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_108():
    f = (a - b*x**2)**(sympy.S(2)/3)*(3*a + b*x**2)**3
    F = -36288*3**(sympy.S(1)/4)*a**(sympy.S(13)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(1235*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 24192*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(13)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(1235*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - 72576*a**4*x/(1235*a**(sympy.S(1)/3)*(1 - sqrt(3)) - 1235*(a - b*x**2)**(sympy.S(1)/3)) + 18144*a**3*x*(a - b*x**2)**(sympy.S(2)/3)/1235 - 23544*a**2*x*(a - b*x**2)**(sympy.S(5)/3)/6175 - 378*a*x*(a - b*x**2)**(sympy.S(5)/3)*(3*a + b*x**2)/475 - 3*x*(a - b*x**2)**(sympy.S(5)/3)*(3*a + b*x**2)**2/25
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_109():
    f = (a - b*x**2)**(sympy.S(2)/3)*(3*a + b*x**2)**2
    F = -15552*3**(sympy.S(1)/4)*a**(sympy.S(10)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(1729*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 10368*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(10)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(1729*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - 31104*a**3*x/(1729*a**(sympy.S(1)/3)*(1 - sqrt(3)) - 1729*(a - b*x**2)**(sympy.S(1)/3)) + 7776*a**2*x*(a - b*x**2)**(sympy.S(2)/3)/1729 - 252*a*x*(a - b*x**2)**(sympy.S(5)/3)/247 - 3*x*(a - b*x**2)**(sympy.S(5)/3)*(3*a + b*x**2)/19
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_110():
    f = (a - b*x**2)**(sympy.S(2)/3)*(3*a + b*x**2)
    F = -36*3**(sympy.S(1)/4)*a**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(13*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 24*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(13*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - 72*a**2*x/(13*a**(sympy.S(1)/3)*(1 - sqrt(3)) - 13*(a - b*x**2)**(sympy.S(1)/3)) + 18*a*x*(a - b*x**2)**(sympy.S(2)/3)/13 - 3*x*(a - b*x**2)**(sympy.S(5)/3)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_111():
    f = (a - b*x**2)**(sympy.S(2)/3)/(3*a + b*x**2)
    F = 2**(sympy.S(1)/3)*sqrt(3)*a**(sympy.S(1)/6)*atan(sqrt(3)*sqrt(a)/(sqrt(b)*x))/(3*sqrt(b)) + 2**(sympy.S(1)/3)*sqrt(3)*a**(sympy.S(1)/6)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))/(sqrt(b)*x))/(3*sqrt(b)) - 2**(sympy.S(1)/3)*a**(sympy.S(1)/6)*atanh(sqrt(b)*x/sqrt(a))/(3*sqrt(b)) + 2**(sympy.S(1)/3)*a**(sympy.S(1)/6)*atanh(sqrt(b)*x/(a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))))/sqrt(b) + 3*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(2*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 3*x/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_112():
    f = (a - b*x**2)**(sympy.S(2)/3)/(3*a + b*x**2)**2
    F = x*(a - b*x**2)**(sympy.S(2)/3)/(6*a*(3*a + b*x**2)) - x/(6*a*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))) - 3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(12*a**(sympy.S(2)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + sqrt(2)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(18*a**(sympy.S(2)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_113():
    f = (a - b*x**2)**(sympy.S(2)/3)/(3*a + b*x**2)**3
    F = x*(a - b*x**2)**(sympy.S(2)/3)/(12*a*(3*a + b*x**2)**2) + x*(a - b*x**2)**(sympy.S(2)/3)/(36*a**2*(3*a + b*x**2)) - x/(36*a**2*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))) - 3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(72*a**(sympy.S(5)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + sqrt(2)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(108*a**(sympy.S(5)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*sqrt(a)/(sqrt(b)*x))/(432*a**(sympy.S(11)/6)*sqrt(b)) + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))/(sqrt(b)*x))/(432*a**(sympy.S(11)/6)*sqrt(b)) - 2**(sympy.S(1)/3)*atanh(sqrt(b)*x/sqrt(a))/(432*a**(sympy.S(11)/6)*sqrt(b)) + 2**(sympy.S(1)/3)*atanh(sqrt(b)*x/(a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))))/(144*a**(sympy.S(11)/6)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_114():
    f = (a - b*x**2)**(sympy.S(2)/3)/(3*a + b*x**2)**4
    F = x*(a - b*x**2)**(sympy.S(2)/3)/(18*a*(3*a + b*x**2)**3) + x*(a - b*x**2)**(sympy.S(2)/3)/(54*a**2*(3*a + b*x**2)**2) + x*(a - b*x**2)**(sympy.S(2)/3)/(144*a**3*(3*a + b*x**2)) - x/(144*a**3*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))) - 3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(288*a**(sympy.S(8)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + sqrt(2)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(432*a**(sympy.S(8)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 7*2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*sqrt(a)/(sqrt(b)*x))/(7776*a**(sympy.S(17)/6)*sqrt(b)) + 7*2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))/(sqrt(b)*x))/(7776*a**(sympy.S(17)/6)*sqrt(b)) - 7*2**(sympy.S(1)/3)*atanh(sqrt(b)*x/sqrt(a))/(7776*a**(sympy.S(17)/6)*sqrt(b)) + 7*2**(sympy.S(1)/3)*atanh(sqrt(b)*x/(a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))))/(2592*a**(sympy.S(17)/6)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_115():
    f = (a - b*x**2)**(sympy.S(5)/3)*(3*a + b*x**2)**3
    F = -5619456*3**(sympy.S(1)/4)*a**(sympy.S(16)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(267995*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 3746304*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(16)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(267995*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - 11238912*a**5*x/(267995*a**(sympy.S(1)/3)*(1 - sqrt(3)) - 267995*(a - b*x**2)**(sympy.S(1)/3)) + 2809728*a**4*x*(a - b*x**2)**(sympy.S(2)/3)/267995 + 1404864*a**3*x*(a - b*x**2)**(sympy.S(5)/3)/191425 - 33264*a**2*x*(a - b*x**2)**(sympy.S(8)/3)/14725 - 432*a*x*(a - b*x**2)**(sympy.S(8)/3)*(3*a + b*x**2)/775 - 3*x*(a - b*x**2)**(sympy.S(8)/3)*(3*a + b*x**2)**2/31
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_116():
    f = (a - b*x**2)**(sympy.S(5)/3)*(3*a + b*x**2)**2
    F = -57024*3**(sympy.S(1)/4)*a**(sympy.S(13)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(8645*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 38016*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(13)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(8645*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - 114048*a**4*x/(8645*a**(sympy.S(1)/3)*(1 - sqrt(3)) - 8645*(a - b*x**2)**(sympy.S(1)/3)) + 28512*a**3*x*(a - b*x**2)**(sympy.S(2)/3)/8645 + 14256*a**2*x*(a - b*x**2)**(sympy.S(5)/3)/6175 - 306*a*x*(a - b*x**2)**(sympy.S(8)/3)/475 - 3*x*(a - b*x**2)**(sympy.S(8)/3)*(3*a + b*x**2)/25
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_117():
    f = (a - b*x**2)**(sympy.S(5)/3)*(3*a + b*x**2)
    F = -3600*3**(sympy.S(1)/4)*a**(sympy.S(10)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(1729*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 2400*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(10)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(1729*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - 7200*a**3*x/(1729*a**(sympy.S(1)/3)*(1 - sqrt(3)) - 1729*(a - b*x**2)**(sympy.S(1)/3)) + 1800*a**2*x*(a - b*x**2)**(sympy.S(2)/3)/1729 + 180*a*x*(a - b*x**2)**(sympy.S(5)/3)/247 - 3*x*(a - b*x**2)**(sympy.S(8)/3)/19
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_118():
    f = (a - b*x**2)**(sympy.S(5)/3)/(3*a + b*x**2)
    F = 4*2**(sympy.S(1)/3)*sqrt(3)*a**(sympy.S(7)/6)*atan(sqrt(3)*sqrt(a)/(sqrt(b)*x))/(3*sqrt(b)) + 4*2**(sympy.S(1)/3)*sqrt(3)*a**(sympy.S(7)/6)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))/(sqrt(b)*x))/(3*sqrt(b)) - 4*2**(sympy.S(1)/3)*a**(sympy.S(7)/6)*atanh(sqrt(b)*x/sqrt(a))/(3*sqrt(b)) + 4*2**(sympy.S(1)/3)*a**(sympy.S(7)/6)*atanh(sqrt(b)*x/(a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))))/sqrt(b) + 48*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(7*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - 32*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(7*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 96*a*x/(7*a**(sympy.S(1)/3)*(1 - sqrt(3)) - 7*(a - b*x**2)**(sympy.S(1)/3)) - 3*x*(a - b*x**2)**(sympy.S(2)/3)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_119():
    f = (a - b*x**2)**(sympy.S(5)/3)/(3*a + b*x**2)**2
    F = -2**(sympy.S(1)/3)*sqrt(3)*a**(sympy.S(1)/6)*atan(sqrt(3)*sqrt(a)/(sqrt(b)*x))/(3*sqrt(b)) - 2**(sympy.S(1)/3)*sqrt(3)*a**(sympy.S(1)/6)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))/(sqrt(b)*x))/(3*sqrt(b)) + 2**(sympy.S(1)/3)*a**(sympy.S(1)/6)*atanh(sqrt(b)*x/sqrt(a))/(3*sqrt(b)) - 2**(sympy.S(1)/3)*a**(sympy.S(1)/6)*atanh(sqrt(b)*x/(a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))))/sqrt(b) - 11*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(6*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 11*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(9*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 2*x*(a - b*x**2)**(sympy.S(2)/3)/(9*a + 3*b*x**2) - 11*x/(3*a**(sympy.S(1)/3)*(1 - sqrt(3)) - 3*(a - b*x**2)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_120():
    f = (a - b*x**2)**(sympy.S(5)/3)/(3*a + b*x**2)**3
    F = x*(a - b*x**2)**(sympy.S(2)/3)/(3*(3*a + b*x**2)**2) - x*(a - b*x**2)**(sympy.S(2)/3)/(18*a*(3*a + b*x**2)) + x/(18*a*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))) + 3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(36*a**(sympy.S(2)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - sqrt(2)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(54*a**(sympy.S(2)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*sqrt(a)/(sqrt(b)*x))/(108*a**(sympy.S(5)/6)*sqrt(b)) + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))/(sqrt(b)*x))/(108*a**(sympy.S(5)/6)*sqrt(b)) - 2**(sympy.S(1)/3)*atanh(sqrt(b)*x/sqrt(a))/(108*a**(sympy.S(5)/6)*sqrt(b)) + 2**(sympy.S(1)/3)*atanh(sqrt(b)*x/(a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))))/(36*a**(sympy.S(5)/6)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_121():
    f = (3*a + b*x**2)**4/(a - b*x**2)**(sympy.S(1)/3)
    F = -1897344*3**(sympy.S(1)/4)*a**(sympy.S(13)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(8645*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 1264896*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(13)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(8645*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - 3794688*a**4*x/(8645*a**(sympy.S(1)/3)*(1 - sqrt(3)) - 8645*(a - b*x**2)**(sympy.S(1)/3)) - 1552608*a**3*x*(a - b*x**2)**(sympy.S(2)/3)/43225 - 36288*a**2*x*(a - b*x**2)**(sympy.S(2)/3)*(3*a + b*x**2)/6175 - 18*a*x*(a - b*x**2)**(sympy.S(2)/3)*(3*a + b*x**2)**2/19 - 3*x*(a - b*x**2)**(sympy.S(2)/3)*(3*a + b*x**2)**3/25
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_122():
    f = (3*a + b*x**2)**3/(a - b*x**2)**(sympy.S(1)/3)
    F = -107568*3**(sympy.S(1)/4)*a**(sympy.S(10)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(1729*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 71712*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(10)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(1729*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - 215136*a**3*x/(1729*a**(sympy.S(1)/3)*(1 - sqrt(3)) - 1729*(a - b*x**2)**(sympy.S(1)/3)) - 15768*a**2*x*(a - b*x**2)**(sympy.S(2)/3)/1729 - 324*a*x*(a - b*x**2)**(sympy.S(2)/3)*(3*a + b*x**2)/247 - 3*x*(a - b*x**2)**(sympy.S(2)/3)*(3*a + b*x**2)**2/19
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_123():
    f = (3*a + b*x**2)**2/(a - b*x**2)**(sympy.S(1)/3)
    F = -1620*3**(sympy.S(1)/4)*a**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(91*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 1080*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(91*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - 3240*a**2*x/(91*a**(sympy.S(1)/3)*(1 - sqrt(3)) - 91*(a - b*x**2)**(sympy.S(1)/3)) - 198*a*x*(a - b*x**2)**(sympy.S(2)/3)/91 - 3*x*(a - b*x**2)**(sympy.S(2)/3)*(3*a + b*x**2)/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_124():
    f = (3*a + b*x**2)/(a - b*x**2)**(sympy.S(1)/3)
    F = -36*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(7*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 24*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(7*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - 72*a*x/(7*a**(sympy.S(1)/3)*(1 - sqrt(3)) - 7*(a - b*x**2)**(sympy.S(1)/3)) - 3*x*(a - b*x**2)**(sympy.S(2)/3)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_125():
    f = 1/((a - b*x**2)**(sympy.S(1)/3)*(3*a + b*x**2))
    F = 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*sqrt(a)/(sqrt(b)*x))/(12*a**(sympy.S(5)/6)*sqrt(b)) + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))/(sqrt(b)*x))/(12*a**(sympy.S(5)/6)*sqrt(b)) - 2**(sympy.S(1)/3)*atanh(sqrt(b)*x/sqrt(a))/(12*a**(sympy.S(5)/6)*sqrt(b)) + 2**(sympy.S(1)/3)*atanh(sqrt(b)*x/(a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))))/(4*a**(sympy.S(5)/6)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_126():
    f = 1/((a - b*x**2)**(sympy.S(1)/3)*(3*a + b*x**2)**2)
    F = x*(a - b*x**2)**(sympy.S(2)/3)/(24*a**2*(3*a + b*x**2)) - x/(24*a**2*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))) - 3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(48*a**(sympy.S(5)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + sqrt(2)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(72*a**(sympy.S(5)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*sqrt(a)/(sqrt(b)*x))/(48*a**(sympy.S(11)/6)*sqrt(b)) + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))/(sqrt(b)*x))/(48*a**(sympy.S(11)/6)*sqrt(b)) - 2**(sympy.S(1)/3)*atanh(sqrt(b)*x/sqrt(a))/(48*a**(sympy.S(11)/6)*sqrt(b)) + 2**(sympy.S(1)/3)*atanh(sqrt(b)*x/(a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))))/(16*a**(sympy.S(11)/6)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_127():
    f = 1/((a - b*x**2)**(sympy.S(1)/3)*(3*a + b*x**2)**3)
    F = x*(a - b*x**2)**(sympy.S(2)/3)/(48*a**2*(3*a + b*x**2)**2) + 5*x*(a - b*x**2)**(sympy.S(2)/3)/(288*a**3*(3*a + b*x**2)) - 5*x/(288*a**3*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))) - 5*3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(576*a**(sympy.S(8)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + sqrt(2)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(5*a**(sympy.S(1)/3) - 5*(a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(864*a**(sympy.S(8)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 5*2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*sqrt(a)/(sqrt(b)*x))/(864*a**(sympy.S(17)/6)*sqrt(b)) + 5*2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))/(sqrt(b)*x))/(864*a**(sympy.S(17)/6)*sqrt(b)) - 5*2**(sympy.S(1)/3)*atanh(sqrt(b)*x/sqrt(a))/(864*a**(sympy.S(17)/6)*sqrt(b)) + 5*2**(sympy.S(1)/3)*atanh(sqrt(b)*x/(a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))))/(288*a**(sympy.S(17)/6)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_128():
    f = (3*a + b*x**2)**3/(a - b*x**2)**(sympy.S(4)/3)
    F = 10044*3**(sympy.S(1)/4)*a**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(91*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - 6696*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(91*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 20088*a**2*x/(91*a**(sympy.S(1)/3)*(1 - sqrt(3)) - 91*(a - b*x**2)**(sympy.S(1)/3)) + 2538*a*x*(a - b*x**2)**(sympy.S(2)/3)/91 + 81*x*(a - b*x**2)**(sympy.S(2)/3)*(3*a + b*x**2)/13 + 6*x*(3*a + b*x**2)**2/(a - b*x**2)**(sympy.S(1)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_129():
    f = (3*a + b*x**2)**2/(a - b*x**2)**(sympy.S(4)/3)
    F = 162*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(7*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - 108*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(7*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 324*a*x/(7*a**(sympy.S(1)/3)*(1 - sqrt(3)) - 7*(a - b*x**2)**(sympy.S(1)/3)) + 45*x*(a - b*x**2)**(sympy.S(2)/3)/7 + 6*x*(3*a + b*x**2)/(a - b*x**2)**(sympy.S(1)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_130():
    f = (3*a + b*x**2)/(a - b*x**2)**(sympy.S(4)/3)
    F = 9*3**(sympy.S(1)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(2*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - 3*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(1)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 9*x/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3)) + 6*x/(a - b*x**2)**(sympy.S(1)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_131():
    f = 1/((a - b*x**2)**(sympy.S(4)/3)*(3*a + b*x**2))
    F = 3*x/(8*a**2*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))) + 3*x/(8*a**2*(a - b*x**2)**(sympy.S(1)/3)) + 3*3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(16*a**(sympy.S(5)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - sqrt(2)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(8*a**(sympy.S(5)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*sqrt(a)/(sqrt(b)*x))/(48*a**(sympy.S(11)/6)*sqrt(b)) + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))/(sqrt(b)*x))/(48*a**(sympy.S(11)/6)*sqrt(b)) - 2**(sympy.S(1)/3)*atanh(sqrt(b)*x/sqrt(a))/(48*a**(sympy.S(11)/6)*sqrt(b)) + 2**(sympy.S(1)/3)*atanh(sqrt(b)*x/(a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))))/(16*a**(sympy.S(11)/6)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_132():
    f = 1/((a - b*x**2)**(sympy.S(4)/3)*(3*a + b*x**2)**2)
    F = x/(24*a**2*(a - b*x**2)**(sympy.S(1)/3)*(3*a + b*x**2)) + x/(12*a**3*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))) + x/(12*a**3*(a - b*x**2)**(sympy.S(1)/3)) + 3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(24*a**(sympy.S(8)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - sqrt(2)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(36*a**(sympy.S(8)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*sqrt(a)/(sqrt(b)*x))/(96*a**(sympy.S(17)/6)*sqrt(b)) + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))/(sqrt(b)*x))/(96*a**(sympy.S(17)/6)*sqrt(b)) - 2**(sympy.S(1)/3)*atanh(sqrt(b)*x/sqrt(a))/(96*a**(sympy.S(17)/6)*sqrt(b)) + 2**(sympy.S(1)/3)*atanh(sqrt(b)*x/(a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))))/(32*a**(sympy.S(17)/6)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_133():
    f = 1/((a - b*x**2)**(sympy.S(4)/3)*(3*a + b*x**2)**3)
    F = x/(48*a**2*(a - b*x**2)**(sympy.S(1)/3)*(3*a + b*x**2)**2) + 17*x/(192*a**3*(a - b*x**2)**(sympy.S(1)/3)*(3*a + b*x**2)) - 19*x*(a - b*x**2)**(sympy.S(2)/3)/(1152*a**4*(3*a + b*x**2)) + 19*x/(1152*a**4*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))) + 19*3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(2304*a**(sympy.S(11)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - sqrt(2)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(19*a**(sympy.S(1)/3) - 19*(a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(3456*a**(sympy.S(11)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 7*2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*sqrt(a)/(sqrt(b)*x))/(1728*a**(sympy.S(23)/6)*sqrt(b)) + 7*2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))/(sqrt(b)*x))/(1728*a**(sympy.S(23)/6)*sqrt(b)) - 7*2**(sympy.S(1)/3)*atanh(sqrt(b)*x/sqrt(a))/(1728*a**(sympy.S(23)/6)*sqrt(b)) + 7*2**(sympy.S(1)/3)*atanh(sqrt(b)*x/(a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))))/(576*a**(sympy.S(23)/6)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_134():
    f = (3*a + b*x**2)**4/(a - b*x**2)**(sympy.S(7)/3)
    F = -18468*3**(sympy.S(1)/4)*a**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(91*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 12312*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(7)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(91*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - 36936*a**2*x/(91*a**(sympy.S(1)/3)*(1 - sqrt(3)) - 91*(a - b*x**2)**(sympy.S(1)/3)) - 3240*a*x*(a - b*x**2)**(sympy.S(2)/3)/91 - 81*x*(a - b*x**2)**(sympy.S(2)/3)*(3*a + b*x**2)/13 - 9*x*(3*a + b*x**2)**2/(2*(a - b*x**2)**(sympy.S(1)/3)) + 3*x*(3*a + b*x**2)**3/(2*(a - b*x**2)**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_135():
    f = (3*a + b*x**2)**3/(a - b*x**2)**(sympy.S(7)/3)
    F = -162*3**(sympy.S(1)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(7*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 108*sqrt(2)*3**(sympy.S(3)/4)*a**(sympy.S(4)/3)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(7*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - 324*a*x/(7*a**(sympy.S(1)/3)*(1 - sqrt(3)) - 7*(a - b*x**2)**(sympy.S(1)/3)) - 27*x*(a - b*x**2)**(sympy.S(2)/3)/14 + 3*x*(3*a + b*x**2)**2/(2*(a - b*x**2)**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_136():
    f = (3*a + b*x**2)**2/(a - b*x**2)**(sympy.S(7)/3)
    F = 9*x/(2*(a - b*x**2)**(sympy.S(1)/3)) + 3*x*(3*a + b*x**2)/(2*(a - b*x**2)**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_137():
    f = (3*a + b*x**2)/(a - b*x**2)**(sympy.S(7)/3)
    F = 3*x/(2*(a - b*x**2)**(sympy.S(4)/3)) + 9*x/(4*a*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))) + 9*x/(4*a*(a - b*x**2)**(sympy.S(1)/3)) + 9*3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(8*a**(sympy.S(2)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - 3*sqrt(2)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(4*a**(sympy.S(2)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_138():
    f = 1/((a - b*x**2)**(sympy.S(7)/3)*(3*a + b*x**2))
    F = 3*x/(32*a**2*(a - b*x**2)**(sympy.S(4)/3)) + 21*x/(64*a**3*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))) + 21*x/(64*a**3*(a - b*x**2)**(sympy.S(1)/3)) + 21*3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(128*a**(sympy.S(8)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - 7*sqrt(2)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(64*a**(sympy.S(8)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*sqrt(a)/(sqrt(b)*x))/(192*a**(sympy.S(17)/6)*sqrt(b)) + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))/(sqrt(b)*x))/(192*a**(sympy.S(17)/6)*sqrt(b)) - 2**(sympy.S(1)/3)*atanh(sqrt(b)*x/sqrt(a))/(192*a**(sympy.S(17)/6)*sqrt(b)) + 2**(sympy.S(1)/3)*atanh(sqrt(b)*x/(a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))))/(64*a**(sympy.S(17)/6)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_139():
    f = 1/((a - b*x**2)**(sympy.S(7)/3)*(3*a + b*x**2)**2)
    F = x/(24*a**2*(a - b*x**2)**(sympy.S(4)/3)*(3*a + b*x**2)) + 5*x/(384*a**3*(a - b*x**2)**(sympy.S(4)/3)) + 79*x/(768*a**4*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))) + 79*x/(768*a**4*(a - b*x**2)**(sympy.S(1)/3)) + 79*3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))*elliptic_e(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(1536*a**(sympy.S(11)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) - sqrt(2)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3) + (a - b*x**2)**(sympy.S(2)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)*(79*a**(sympy.S(1)/3) - 79*(a - b*x**2)**(sympy.S(1)/3))*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))), -7 + 4*sqrt(3))/(2304*a**(sympy.S(11)/3)*b*x*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - (a - b*x**2)**(sympy.S(1)/3))**2)) + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*sqrt(a)/(sqrt(b)*x))/(256*a**(sympy.S(23)/6)*sqrt(b)) + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))/(sqrt(b)*x))/(256*a**(sympy.S(23)/6)*sqrt(b)) - 2**(sympy.S(1)/3)*atanh(sqrt(b)*x/sqrt(a))/(256*a**(sympy.S(23)/6)*sqrt(b)) + 3*2**(sympy.S(1)/3)*atanh(sqrt(b)*x/(a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))))/(256*a**(sympy.S(23)/6)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_140():
    f = 1/((-3*a - b*x**2)*(-a + b*x**2)**(sympy.S(1)/3))
    F = -2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*sqrt(a)/(sqrt(b)*x))/(12*sqrt(a)*sqrt(b)*(-a)**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*sqrt(a)*((-a)**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*(-a + b*x**2)**(sympy.S(1)/3))/(sqrt(b)*x*(-a)**(sympy.S(1)/3)))/(12*sqrt(a)*sqrt(b)*(-a)**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*atanh(sqrt(b)*x/sqrt(a))/(12*sqrt(a)*sqrt(b)*(-a)**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*atanh(sqrt(b)*x*(-a)**(sympy.S(1)/3)/(sqrt(a)*((-a)**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*(-a + b*x**2)**(sympy.S(1)/3))))/(4*sqrt(a)*sqrt(b)*(-a)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_141():
    f = 1/((a + b*x**2)**(sympy.S(1)/3)*(3*a - b*x**2))
    F = -2**(sympy.S(1)/3)*atan(sqrt(b)*x/sqrt(a))/(12*a**(sympy.S(5)/6)*sqrt(b)) + 2**(sympy.S(1)/3)*atan(sqrt(b)*x/(a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3))))/(4*a**(sympy.S(5)/6)*sqrt(b)) - 2**(sympy.S(1)/3)*sqrt(3)*atanh(sqrt(3)*sqrt(a)/(sqrt(b)*x))/(12*a**(sympy.S(5)/6)*sqrt(b)) - 2**(sympy.S(1)/3)*sqrt(3)*atanh(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*(a + b*x**2)**(sympy.S(1)/3))/(sqrt(b)*x))/(12*a**(sympy.S(5)/6)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_142():
    f = 1/((c - d*x**2)*(c + 3*d*x**2)**(sympy.S(1)/3))
    F = -2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*sqrt(d)*x/sqrt(c))/(12*c**(sympy.S(5)/6)*sqrt(d)) + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*sqrt(d)*x/(c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*(c + 3*d*x**2)**(sympy.S(1)/3))))/(4*c**(sympy.S(5)/6)*sqrt(d)) - 2**(sympy.S(1)/3)*atanh(sqrt(c)/(sqrt(d)*x))/(4*c**(sympy.S(5)/6)*sqrt(d)) - 2**(sympy.S(1)/3)*atanh(c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*(c + 3*d*x**2)**(sympy.S(1)/3))/(sqrt(d)*x))/(4*c**(sympy.S(5)/6)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_143():
    f = 1/((a - b*x**2)**(sympy.S(1)/3)*(3*a + b*x**2))
    F = 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*sqrt(a)/(sqrt(b)*x))/(12*a**(sympy.S(5)/6)*sqrt(b)) + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))/(sqrt(b)*x))/(12*a**(sympy.S(5)/6)*sqrt(b)) - 2**(sympy.S(1)/3)*atanh(sqrt(b)*x/sqrt(a))/(12*a**(sympy.S(5)/6)*sqrt(b)) + 2**(sympy.S(1)/3)*atanh(sqrt(b)*x/(a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*(a - b*x**2)**(sympy.S(1)/3))))/(4*a**(sympy.S(5)/6)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_144():
    f = 1/((c - 3*d*x**2)**(sympy.S(1)/3)*(c + d*x**2))
    F = 2**(sympy.S(1)/3)*atan(sqrt(c)/(sqrt(d)*x))/(4*c**(sympy.S(5)/6)*sqrt(d)) + 2**(sympy.S(1)/3)*atan(c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*(c - 3*d*x**2)**(sympy.S(1)/3))/(sqrt(d)*x))/(4*c**(sympy.S(5)/6)*sqrt(d)) - 2**(sympy.S(1)/3)*sqrt(3)*atanh(sqrt(3)*sqrt(d)*x/sqrt(c))/(12*c**(sympy.S(5)/6)*sqrt(d)) + 2**(sympy.S(1)/3)*sqrt(3)*atanh(sqrt(3)*sqrt(d)*x/(c**(sympy.S(1)/6)*(c**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*(c - 3*d*x**2)**(sympy.S(1)/3))))/(4*c**(sympy.S(5)/6)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_145():
    f = 1/((1 - x**2)**(sympy.S(1)/3)*(x**2 + 3))
    F = 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)/x)/12 + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1)/x)/12 - 2**(sympy.S(1)/3)*atanh(x)/12 + 2**(sympy.S(1)/3)*atanh(x/(2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_146():
    f = 1/((3 - x**2)*(x**2 + 1)**(sympy.S(1)/3))
    F = -2**(sympy.S(1)/3)*atan(x)/12 + 2**(sympy.S(1)/3)*atan(x/(2**(sympy.S(1)/3)*(x**2 + 1)**(sympy.S(1)/3) + 1))/4 - 2**(sympy.S(1)/3)*sqrt(3)*atanh(sqrt(3)/x)/12 - 2**(sympy.S(1)/3)*sqrt(3)*atanh(sqrt(3)*(-2**(sympy.S(1)/3)*(x**2 + 1)**(sympy.S(1)/3) + 1)/x)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_147():
    f = (3 - x)/((1 - x**2)**(sympy.S(1)/3)*(x**2 + 3))
    F = -2**(sympy.S(1)/3)*log(x**2 + 3)/4 + 3*2**(sympy.S(1)/3)*log(2**(sympy.S(1)/3)*(1 - x)**(sympy.S(1)/3) + (x + 1)**(sympy.S(2)/3))/4 - 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)/3 - 2**(sympy.S(2)/3)*sqrt(3)*(x + 1)**(sympy.S(2)/3)/(3*(1 - x)**(sympy.S(1)/3)))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_148():
    f = (x + 3)/((1 - x**2)**(sympy.S(1)/3)*(x**2 + 3))
    F = 2**(sympy.S(1)/3)*log(x**2 + 3)/4 - 3*2**(sympy.S(1)/3)*log((1 - x)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(x + 1)**(sympy.S(1)/3))/4 - 2**(sympy.S(1)/3)*sqrt(3)*atan(2**(sympy.S(2)/3)*sqrt(3)*(1 - x)**(sympy.S(2)/3)/(3*(x + 1)**(sympy.S(1)/3)) - sqrt(3)/3)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_149():
    f = 1/((a + b*x**2)**(sympy.S(1)/3)*(9*a*d/b + d*x**2))
    F = sqrt(b)*atan(sqrt(b)*x/(3*sqrt(a)))/(12*a**(sympy.S(5)/6)*d) + sqrt(b)*atan((a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))**2/(3*a**(sympy.S(1)/6)*sqrt(b)*x))/(12*a**(sympy.S(5)/6)*d) - sqrt(3)*sqrt(b)*atanh(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - (a + b*x**2)**(sympy.S(1)/3))/(sqrt(b)*x))/(12*a**(sympy.S(5)/6)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_150():
    f = 1/((a - b*x**2)**(sympy.S(1)/3)*(-9*a*d/b + d*x**2))
    F = -sqrt(3)*sqrt(b)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))/(sqrt(b)*x))/(12*a**(sympy.S(5)/6)*d) - sqrt(b)*atanh(sqrt(b)*x/(3*sqrt(a)))/(12*a**(sympy.S(5)/6)*d) + sqrt(b)*atanh((a**(sympy.S(1)/3) - (a - b*x**2)**(sympy.S(1)/3))**2/(3*a**(sympy.S(1)/6)*sqrt(b)*x))/(12*a**(sympy.S(5)/6)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_151():
    f = 1/((-a + b*x**2)**(sympy.S(1)/3)*(-9*a*d/b + d*x**2))
    F = sqrt(3)*sqrt(b)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + (-a + b*x**2)**(sympy.S(1)/3))/(sqrt(b)*x))/(12*a**(sympy.S(5)/6)*d) + sqrt(b)*atanh(sqrt(b)*x/(3*sqrt(a)))/(12*a**(sympy.S(5)/6)*d) - sqrt(b)*atanh((a**(sympy.S(1)/3) + (-a + b*x**2)**(sympy.S(1)/3))**2/(3*a**(sympy.S(1)/6)*sqrt(b)*x))/(12*a**(sympy.S(5)/6)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_152():
    f = 1/((-a - b*x**2)**(sympy.S(1)/3)*(9*a*d/b + d*x**2))
    F = -sqrt(b)*atan(sqrt(b)*x/(3*sqrt(a)))/(12*a**(sympy.S(5)/6)*d) - sqrt(b)*atan((a**(sympy.S(1)/3) + (-a - b*x**2)**(sympy.S(1)/3))**2/(3*a**(sympy.S(1)/6)*sqrt(b)*x))/(12*a**(sympy.S(5)/6)*d) + sqrt(3)*sqrt(b)*atanh(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + (-a - b*x**2)**(sympy.S(1)/3))/(sqrt(b)*x))/(12*a**(sympy.S(5)/6)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_153():
    f = 1/((b*x**2 + 2)**(sympy.S(1)/3)*(d*x**2 + 18*d/b))
    F = 2**(sympy.S(1)/6)*sqrt(b)*atan(sqrt(2)*sqrt(b)*x/6)/(24*d) + 2**(sympy.S(1)/6)*sqrt(b)*atan(2**(sympy.S(5)/6)*(-(b*x**2 + 2)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2/(6*sqrt(b)*x))/(24*d) - 2**(sympy.S(1)/6)*sqrt(3)*sqrt(b)*atanh(2**(sympy.S(1)/6)*sqrt(3)*(-(b*x**2 + 2)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(sqrt(b)*x))/(24*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_154():
    f = 1/((b*x**2 - 2)**(sympy.S(1)/3)*(d*x**2 - 18*d/b))
    F = 2**(sympy.S(1)/6)*sqrt(3)*sqrt(b)*atan(2**(sympy.S(1)/6)*sqrt(3)*((b*x**2 - 2)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(sqrt(b)*x))/(24*d) + 2**(sympy.S(1)/6)*sqrt(b)*atanh(sqrt(2)*sqrt(b)*x/6)/(24*d) - 2**(sympy.S(1)/6)*sqrt(b)*atanh(2**(sympy.S(5)/6)*((b*x**2 - 2)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2/(6*sqrt(b)*x))/(24*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_155():
    f = 1/((3*x**2 + 2)**(sympy.S(1)/3)*(d*x**2 + 6*d))
    F = 2**(sympy.S(1)/6)*sqrt(3)*atan(sqrt(6)*x/6)/(24*d) + 2**(sympy.S(1)/6)*sqrt(3)*atan(2**(sympy.S(5)/6)*sqrt(3)*(-(3*x**2 + 2)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2/(18*x))/(24*d) - 2**(sympy.S(1)/6)*atanh(2**(sympy.S(1)/6)*(-(3*x**2 + 2)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_156():
    f = 1/((2 - 3*x**2)**(sympy.S(1)/3)*(d*x**2 - 6*d))
    F = -2**(sympy.S(1)/6)*atan(2**(sympy.S(1)/6)*(-(2 - 3*x**2)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/x)/(8*d) - 2**(sympy.S(1)/6)*sqrt(3)*atanh(sqrt(6)*x/6)/(24*d) + 2**(sympy.S(1)/6)*sqrt(3)*atanh(2**(sympy.S(5)/6)*sqrt(3)*(-(2 - 3*x**2)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2/(18*x))/(24*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_157():
    f = 1/((3*x**2 - 2)**(sympy.S(1)/3)*(d*x**2 - 6*d))
    F = 2**(sympy.S(1)/6)*atan(2**(sympy.S(1)/6)*((3*x**2 - 2)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/x)/(8*d) + 2**(sympy.S(1)/6)*sqrt(3)*atanh(sqrt(6)*x/6)/(24*d) - 2**(sympy.S(1)/6)*sqrt(3)*atanh(2**(sympy.S(5)/6)*sqrt(3)*((3*x**2 - 2)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2/(18*x))/(24*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_158():
    f = 1/((-3*x**2 - 2)**(sympy.S(1)/3)*(d*x**2 + 6*d))
    F = -2**(sympy.S(1)/6)*sqrt(3)*atan(sqrt(6)*x/6)/(24*d) - 2**(sympy.S(1)/6)*sqrt(3)*atan(2**(sympy.S(5)/6)*sqrt(3)*((-3*x**2 - 2)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2/(18*x))/(24*d) + 2**(sympy.S(1)/6)*atanh(2**(sympy.S(1)/6)*((-3*x**2 - 2)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_159():
    f = 1/((x**2 + 1)**(sympy.S(1)/3)*(x**2 + 9))
    F = atan(x/3)/12 + atan((1 - (x**2 + 1)**(sympy.S(1)/3))**2/(3*x))/12 - sqrt(3)*atanh(sqrt(3)*(1 - (x**2 + 1)**(sympy.S(1)/3))/x)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_160():
    f = 1/((b*x**2 + 1)**(sympy.S(1)/3)*(b*x**2 + 9))
    F = atan(sqrt(b)*x/3)/(12*sqrt(b)) + atan((1 - (b*x**2 + 1)**(sympy.S(1)/3))**2/(3*sqrt(b)*x))/(12*sqrt(b)) - sqrt(3)*atanh(sqrt(3)*(1 - (b*x**2 + 1)**(sympy.S(1)/3))/(sqrt(b)*x))/(12*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_161():
    f = 1/((1 - x**2)**(sympy.S(1)/3)*(9 - x**2))
    F = sqrt(3)*atan(sqrt(3)*(1 - (1 - x**2)**(sympy.S(1)/3))/x)/12 + atanh(x/3)/12 - atanh((1 - (1 - x**2)**(sympy.S(1)/3))**2/(3*x))/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_162():
    f = (a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)
    F = b*x*sqrt(a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2)/(5*d) - c**(sympy.S(3)/2)*sqrt(a + b*x**2)*(-9*a*d + b*c)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(15*d**(sympy.S(3)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + x*sqrt(a + b*x**2)*(3*a**2*d/b + 7*a*c - 2*b*c**2/d)/(15*sqrt(c + d*x**2)) - x*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-6*a*d + 2*b*c)/(15*d) + sqrt(c)*sqrt(a + b*x**2)*(-3*a**2*d**2 - 7*a*b*c*d + 2*b**2*c**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(15*b*d**(sympy.S(3)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_163():
    f = sqrt(a + b*x**2)*sqrt(c + d*x**2)
    F = 2*c**(sympy.S(3)/2)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + x*sqrt(a + b*x**2)*sqrt(c + d*x**2)/3 - sqrt(c)*sqrt(a + b*x**2)*(a*d + b*c)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*b*sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + x*sqrt(a + b*x**2)*(a*d + b*c)/(3*b*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_164():
    f = sqrt(c + d*x**2)/sqrt(a + b*x**2)
    F = -sqrt(c)*sqrt(d)*sqrt(a + b*x**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(b*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + d*x*sqrt(a + b*x**2)/(b*sqrt(c + d*x**2)) + c**(sympy.S(3)/2)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(a*sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_165():
    f = sqrt(c + d*x**2)/(a + b*x**2)**(sympy.S(3)/2)
    F = sqrt(c + d*x**2)*elliptic_e(atan(sqrt(b)*x/sqrt(a)), -a*d/(b*c) + 1)/(sqrt(a)*sqrt(b)*sqrt(a*(c + d*x**2)/(c*(a + b*x**2)))*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_166():
    f = sqrt(c + d*x**2)/(a + b*x**2)**(sympy.S(5)/2)
    F = x*sqrt(c + d*x**2)/(3*a*(a + b*x**2)**(sympy.S(3)/2)) - c**(sympy.S(3)/2)*sqrt(d)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*a**2*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)*(-a*d + b*c)) + sqrt(c + d*x**2)*(-a*d + 2*b*c)*elliptic_e(atan(sqrt(b)*x/sqrt(a)), -a*d/(b*c) + 1)/(3*a**(sympy.S(3)/2)*sqrt(b)*sqrt(a*(c + d*x**2)/(c*(a + b*x**2)))*sqrt(a + b*x**2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_167():
    f = sqrt(c + d*x**2)/(a + b*x**2)**(sympy.S(7)/2)
    F = x*sqrt(c + d*x**2)/(5*a*(a + b*x**2)**(sympy.S(5)/2)) + x*sqrt(c + d*x**2)*(-3*a*d + 4*b*c)/(15*a**2*(a + b*x**2)**(sympy.S(3)/2)*(-a*d + b*c)) - 2*c**(sympy.S(3)/2)*sqrt(d)*sqrt(a + b*x**2)*(-3*a*d + 2*b*c)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(15*a**3*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)*(-a*d + b*c)**2) + sqrt(c + d*x**2)*(3*a**2*d**2 - 13*a*b*c*d + 8*b**2*c**2)*elliptic_e(atan(sqrt(b)*x/sqrt(a)), -a*d/(b*c) + 1)/(15*a**(sympy.S(5)/2)*sqrt(b)*sqrt(a*(c + d*x**2)/(c*(a + b*x**2)))*sqrt(a + b*x**2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_168():
    f = (a + b*x**2)**(sympy.S(3)/2)*(c + d*x**2)**(sympy.S(3)/2)
    F = x*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-2*a**2*d/(35*b) + 9*a*c/35 + b*c**2/(35*d)) - c**(sympy.S(3)/2)*sqrt(a + b*x**2)*(a**2*d**2 - 18*a*b*c*d + b**2*c**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(35*b*d**(sympy.S(3)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + d*x*(a + b*x**2)**(sympy.S(5)/2)*sqrt(c + d*x**2)/(7*b) + x*(a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)*(-2*a*d + 8*b*c)/(35*b) + 2*sqrt(c)*sqrt(a + b*x**2)*(a*d + b*c)*(a**2*d**2 - 6*a*b*c*d + b**2*c**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(35*b**2*d**(sympy.S(3)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) - x*sqrt(a + b*x**2)*(2*a*d + 2*b*c)*(a**2*d**2 - 6*a*b*c*d + b**2*c**2)/(35*b**2*d*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_169():
    f = sqrt(a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2)
    F = c**(sympy.S(3)/2)*sqrt(a + b*x**2)*(-a*d + 9*b*c)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(15*b*sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + d*x*(a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)/(5*b) + x*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-2*a*d + 6*b*c)/(15*b) - sqrt(c)*sqrt(a + b*x**2)*(-2*a**2*d**2 + 7*a*b*c*d + 3*b**2*c**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(15*b**2*sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + x*sqrt(a + b*x**2)*(-2*a**2*d**2 + 7*a*b*c*d + 3*b**2*c**2)/(15*b**2*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_170():
    f = (c + d*x**2)**(sympy.S(3)/2)/sqrt(a + b*x**2)
    F = d*x*sqrt(a + b*x**2)*sqrt(c + d*x**2)/(3*b) - 2*sqrt(c)*sqrt(d)*sqrt(a + b*x**2)*(-a*d + 2*b*c)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*b**2*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + 2*d*x*sqrt(a + b*x**2)*(-a*d + 2*b*c)/(3*b**2*sqrt(c + d*x**2)) + c**(sympy.S(3)/2)*sqrt(a + b*x**2)*(-a*d + 3*b*c)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*a*b*sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_171():
    f = (c + d*x**2)**(sympy.S(3)/2)/(a + b*x**2)**(sympy.S(3)/2)
    F = c**(sympy.S(3)/2)*sqrt(d)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(a*b*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + x*sqrt(c + d*x**2)*(-a*d + b*c)/(a*b*sqrt(a + b*x**2)) + sqrt(c)*sqrt(d)*sqrt(a + b*x**2)*(-2*a*d + b*c)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(a*b**2*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) - d*x*sqrt(a + b*x**2)*(-2*a*d + b*c)/(a*b**2*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_172():
    f = (c + d*x**2)**(sympy.S(3)/2)/(a + b*x**2)**(sympy.S(5)/2)
    F = x*sqrt(c + d*x**2)*(-a*d + b*c)/(3*a*b*(a + b*x**2)**(sympy.S(3)/2)) - c**(sympy.S(3)/2)*sqrt(d)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*a**2*b*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + sqrt(c + d*x**2)*(2*a*d + 2*b*c)*elliptic_e(atan(sqrt(b)*x/sqrt(a)), -a*d/(b*c) + 1)/(3*a**(sympy.S(3)/2)*b**(sympy.S(3)/2)*sqrt(a*(c + d*x**2)/(c*(a + b*x**2)))*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_173():
    f = (c + d*x**2)**(sympy.S(3)/2)/(a + b*x**2)**(sympy.S(7)/2)
    F = x*sqrt(c + d*x**2)*(-a*d + b*c)/(5*a*b*(a + b*x**2)**(sympy.S(5)/2)) + x*sqrt(c + d*x**2)*(2*a*d + 4*b*c)/(15*a**2*b*(a + b*x**2)**(sympy.S(3)/2)) - c**(sympy.S(3)/2)*sqrt(d)*sqrt(a + b*x**2)*(-a*d + 4*b*c)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(15*a**3*b*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)*(-a*d + b*c)) + sqrt(c + d*x**2)*(-2*a**2*d**2 - 3*a*b*c*d + 8*b**2*c**2)*elliptic_e(atan(sqrt(b)*x/sqrt(a)), -a*d/(b*c) + 1)/(15*a**(sympy.S(5)/2)*b**(sympy.S(3)/2)*sqrt(a*(c + d*x**2)/(c*(a + b*x**2)))*sqrt(a + b*x**2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_174():
    f = sqrt(b*x**2 + 2)*sqrt(d*x**2 + 3)
    F = x*sqrt(b*x**2 + 2)*sqrt(d*x**2 + 3)/3 + 2*sqrt(2)*sqrt(b*x**2 + 2)*elliptic_f(atan(sqrt(3)*sqrt(d)*x/3), -3*b/(2*d) + 1)/(sqrt(d)*sqrt((b*x**2 + 2)/(d*x**2 + 3))*sqrt(d*x**2 + 3)) + x*(3*b + 2*d)*sqrt(b*x**2 + 2)/(3*b*sqrt(d*x**2 + 3)) - sqrt(2)*(3*b + 2*d)*sqrt(b*x**2 + 2)*elliptic_e(atan(sqrt(3)*sqrt(d)*x/3), -3*b/(2*d) + 1)/(3*b*sqrt(d)*sqrt((b*x**2 + 2)/(d*x**2 + 3))*sqrt(d*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_175():
    f = sqrt(3 - 6*x**2)*sqrt(4*x**2 + 2)
    F = sqrt(6)*x*sqrt(1 - 4*x**4)/3 + 2*sqrt(3)*elliptic_f(asin(sqrt(2)*x), -1)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_176():
    f = sqrt(4*x**2 + 2)*sqrt(6*x**2 + 3)
    F = 2*sqrt(6)*x**3/3 + sqrt(6)*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_177():
    f = sqrt(b*x**2 + 2)/sqrt(d*x**2 + 3)
    F = x*sqrt(b*x**2 + 2)/sqrt(d*x**2 + 3) - sqrt(2)*sqrt(b*x**2 + 2)*elliptic_e(atan(sqrt(3)*sqrt(d)*x/3), -3*b/(2*d) + 1)/(sqrt(d)*sqrt((b*x**2 + 2)/(d*x**2 + 3))*sqrt(d*x**2 + 3)) + sqrt(2)*sqrt(b*x**2 + 2)*elliptic_f(atan(sqrt(3)*sqrt(d)*x/3), -3*b/(2*d) + 1)/(sqrt(d)*sqrt((b*x**2 + 2)/(d*x**2 + 3))*sqrt(d*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_178():
    f = sqrt(4 - x**2)/sqrt(c + d*x**2)
    F = sqrt(1 + d*x**2/c)*(c + 4*d)*elliptic_f(asin(x/2), -4*d/c)/(d*sqrt(c + d*x**2)) - sqrt(c + d*x**2)*elliptic_e(asin(x/2), -4*d/c)/(d*sqrt(1 + d*x**2/c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_179():
    f = sqrt(x**2 + 4)/sqrt(c + d*x**2)
    F = x*sqrt(c + d*x**2)/(d*sqrt(x**2 + 4)) - sqrt(c + d*x**2)*elliptic_e(atan(x/2), 1 - 4*d/c)/(d*sqrt((c + d*x**2)/(c*(x**2 + 4)))*sqrt(x**2 + 4)) + 4*sqrt(c + d*x**2)*elliptic_f(atan(x/2), 1 - 4*d/c)/(c*sqrt((c + d*x**2)/(c*(x**2 + 4)))*sqrt(x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_180():
    f = sqrt(1 - x**2)/sqrt(2 - 3*x**2)
    F = sqrt(3)*elliptic_e(asin(sqrt(6)*x/2), sympy.S(2)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_181():
    f = sqrt(4 - x**2)/sqrt(2 - 3*x**2)
    F = 2*sqrt(3)*elliptic_e(asin(sqrt(6)*x/2), sympy.S(1)/6)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_182():
    f = sqrt(1 - 4*x**2)/sqrt(2 - 3*x**2)
    F = sqrt(3)*elliptic_e(asin(sqrt(6)*x/2), sympy.S(8)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_183():
    f = sqrt(x**2 + 1)/sqrt(1 - x**2)
    F = elliptic_e(asin(x), -1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_184():
    f = sqrt(x**2 + 1)/sqrt(2 - 3*x**2)
    F = sqrt(3)*elliptic_e(asin(sqrt(6)*x/2), sympy.S(-2)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_185():
    f = sqrt(x**2 + 4)/sqrt(2 - 3*x**2)
    F = 2*sqrt(3)*elliptic_e(asin(sqrt(6)*x/2), sympy.S(-1)/6)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_186():
    f = sqrt(4*x**2 + 1)/sqrt(2 - 3*x**2)
    F = sqrt(3)*elliptic_e(asin(sqrt(6)*x/2), sympy.S(-8)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_187():
    f = sqrt(1 - x**2)/sqrt(x**2 + 1)
    F = -elliptic_e(asin(x), -1) + 2*elliptic_f(asin(x), -1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_188():
    f = sqrt(1 - x**2)/sqrt(3*x**2 + 2)
    F = -sqrt(2)*elliptic_e(asin(x), sympy.S(-3)/2)/3 + 5*sqrt(2)*elliptic_f(asin(x), sympy.S(-3)/2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_189():
    f = sqrt(4 - x**2)/sqrt(3*x**2 + 2)
    F = -sqrt(2)*elliptic_e(asin(x/2), -6)/3 + 7*sqrt(2)*elliptic_f(asin(x/2), -6)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_190():
    f = sqrt(1 - 4*x**2)/sqrt(3*x**2 + 2)
    F = -2*sqrt(2)*elliptic_e(asin(2*x), sympy.S(-3)/8)/3 + 11*sqrt(2)*elliptic_f(asin(2*x), sympy.S(-3)/8)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_191():
    f = sqrt(x**2 + 1)/sqrt(3*x**2 + 2)
    F = x*sqrt(3*x**2 + 2)/(3*sqrt(x**2 + 1)) - sqrt(2)*sqrt(3*x**2 + 2)*elliptic_e(atan(x), sympy.S(-1)/2)/(3*sqrt((3*x**2 + 2)/(x**2 + 1))*sqrt(x**2 + 1)) + sqrt(2)*sqrt(3*x**2 + 2)*elliptic_f(atan(x), sympy.S(-1)/2)/(2*sqrt((3*x**2 + 2)/(x**2 + 1))*sqrt(x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_192():
    f = sqrt(x**2 + 4)/sqrt(3*x**2 + 2)
    F = x*sqrt(3*x**2 + 2)/(3*sqrt(x**2 + 4)) - sqrt(2)*sqrt(3*x**2 + 2)*elliptic_e(atan(x/2), -5)/(3*sqrt((3*x**2 + 2)/(x**2 + 4))*sqrt(x**2 + 4)) + 2*sqrt(2)*sqrt(3*x**2 + 2)*elliptic_f(atan(x/2), -5)/(sqrt((3*x**2 + 2)/(x**2 + 4))*sqrt(x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_193():
    f = sqrt(4*x**2 + 1)/sqrt(3*x**2 + 2)
    F = 4*x*sqrt(3*x**2 + 2)/(3*sqrt(4*x**2 + 1)) - 2*sqrt(2)*sqrt(3*x**2 + 2)*elliptic_e(atan(2*x), sympy.S(5)/8)/(3*sqrt((3*x**2 + 2)/(4*x**2 + 1))*sqrt(4*x**2 + 1)) + sqrt(2)*sqrt(3*x**2 + 2)*elliptic_f(atan(2*x), sympy.S(5)/8)/(4*sqrt((3*x**2 + 2)/(4*x**2 + 1))*sqrt(4*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_194():
    f = sqrt(1 - x**2)/sqrt(2*x**2 - 1)
    F = sqrt(2)*sqrt(1 - 2*x**2)*elliptic_e(asin(sqrt(2)*x), sympy.S.Half)/(2*sqrt(2*x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_195():
    f = (a + b*x**2)**(sympy.S(7)/2)/sqrt(c + d*x**2)
    F = b*x*(a + b*x**2)**(sympy.S(5)/2)*sqrt(c + d*x**2)/(7*d) - 6*b*x*(a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)*(-2*a*d + b*c)/(35*d**2) + b*x*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(71*a**2*d**2 - 71*a*b*c*d + 24*b**2*c**2)/(105*d**3) - sqrt(c)*sqrt(a + b*x**2)*(-7*a*d + 3*b*c)*(15*a**2*d**2 - 11*a*b*c*d + 8*b**2*c**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(105*d**(sympy.S(7)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + 8*sqrt(c)*sqrt(a + b*x**2)*(-2*a*d + b*c)*(11*a**2*d**2 - 11*a*b*c*d + 6*b**2*c**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(105*d**(sympy.S(7)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) - x*sqrt(a + b*x**2)*(-16*a*d + 8*b*c)*(11*a**2*d**2 - 11*a*b*c*d + 6*b**2*c**2)/(105*d**3*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_196():
    f = (a + b*x**2)**(sympy.S(5)/2)/sqrt(c + d*x**2)
    F = b*x*(a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)/(5*d) - 4*b*x*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-2*a*d + b*c)/(15*d**2) + sqrt(c)*sqrt(a + b*x**2)*(15*a**2*d**2 - 11*a*b*c*d + 4*b**2*c**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(15*d**(sympy.S(5)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) - sqrt(c)*sqrt(a + b*x**2)*(23*a**2*d**2 - 23*a*b*c*d + 8*b**2*c**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(15*d**(sympy.S(5)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + x*sqrt(a + b*x**2)*(23*a**2*d**2 - 23*a*b*c*d + 8*b**2*c**2)/(15*d**2*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_197():
    f = (a + b*x**2)**(sympy.S(3)/2)/sqrt(c + d*x**2)
    F = b*x*sqrt(a + b*x**2)*sqrt(c + d*x**2)/(3*d) - sqrt(c)*sqrt(a + b*x**2)*(-3*a*d + b*c)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*d**(sympy.S(3)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + 2*sqrt(c)*sqrt(a + b*x**2)*(-2*a*d + b*c)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*d**(sympy.S(3)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) - x*sqrt(a + b*x**2)*(-4*a*d + 2*b*c)/(3*d*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_198():
    f = sqrt(a + b*x**2)/sqrt(c + d*x**2)
    F = -sqrt(c)*sqrt(a + b*x**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + sqrt(c)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + x*sqrt(a + b*x**2)/sqrt(c + d*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_199():
    f = 1/(sqrt(a + b*x**2)*sqrt(c + d*x**2))
    F = sqrt(c)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(a*sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_200():
    f = 1/((a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2))
    F = b*x*sqrt(c + d*x**2)/(a*sqrt(a + b*x**2)*(-a*d + b*c)) + sqrt(c)*sqrt(d)*sqrt(a + b*x**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(a*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)*(-a*d + b*c)) - sqrt(c)*sqrt(d)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(a*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)*(-a*d + b*c)) - d*x*sqrt(a + b*x**2)/(a*sqrt(c + d*x**2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_201():
    f = 1/((a + b*x**2)**(sympy.S(5)/2)*sqrt(c + d*x**2))
    F = b*x*sqrt(c + d*x**2)/(3*a*(a + b*x**2)**(sympy.S(3)/2)*(-a*d + b*c)) - sqrt(c)*sqrt(d)*sqrt(a + b*x**2)*(-3*a*d + b*c)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*a**2*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)*(-a*d + b*c)**2) + 2*sqrt(b)*sqrt(c + d*x**2)*(-2*a*d + b*c)*elliptic_e(atan(sqrt(b)*x/sqrt(a)), -a*d/(b*c) + 1)/(3*a**(sympy.S(3)/2)*sqrt(a*(c + d*x**2)/(c*(a + b*x**2)))*sqrt(a + b*x**2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_202():
    f = 1/((a + b*x**2)**(sympy.S(7)/2)*sqrt(c + d*x**2))
    F = b*x*sqrt(c + d*x**2)/(5*a*(a + b*x**2)**(sympy.S(5)/2)*(-a*d + b*c)) + 4*b*x*sqrt(c + d*x**2)*(-2*a*d + b*c)/(15*a**2*(a + b*x**2)**(sympy.S(3)/2)*(-a*d + b*c)**2) - sqrt(c)*sqrt(d)*sqrt(a + b*x**2)*(15*a**2*d**2 - 11*a*b*c*d + 4*b**2*c**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(15*a**3*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)*(-a*d + b*c)**3) + sqrt(b)*sqrt(c + d*x**2)*(23*a**2*d**2 - 23*a*b*c*d + 8*b**2*c**2)*elliptic_e(atan(sqrt(b)*x/sqrt(a)), -a*d/(b*c) + 1)/(15*a**(sympy.S(5)/2)*sqrt(a*(c + d*x**2)/(c*(a + b*x**2)))*sqrt(a + b*x**2)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_203():
    f = (a + b*x**2)**(sympy.S(7)/2)/(c + d*x**2)**(sympy.S(3)/2)
    F = b*sqrt(c)*sqrt(a + b*x**2)*(45*a**2*d**2 - 61*a*b*c*d + 24*b**2*c**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(15*d**(sympy.S(7)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + b*x*(a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)*(-5*a*d + 6*b*c)/(5*c*d**2) - b*x*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(15*a**2*d**2 - 43*a*b*c*d + 24*b**2*c**2)/(15*c*d**3) - x*(a + b*x**2)**(sympy.S(5)/2)*(-a*d + b*c)/(c*d*sqrt(c + d*x**2)) + x*sqrt(a + b*x**2)*(-15*a**3*d**3 + 103*a**2*b*c*d**2 - 128*a*b**2*c**2*d + 48*b**3*c**3)/(15*c*d**3*sqrt(c + d*x**2)) - sqrt(a + b*x**2)*(-15*a**3*d**3 + 103*a**2*b*c*d**2 - 128*a*b**2*c**2*d + 48*b**3*c**3)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(15*sqrt(c)*d**(sympy.S(7)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_204():
    f = (a + b*x**2)**(sympy.S(5)/2)/(c + d*x**2)**(sympy.S(3)/2)
    F = -2*b*sqrt(c)*sqrt(a + b*x**2)*(-3*a*d + 2*b*c)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*d**(sympy.S(5)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + b*x*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-3*a*d + 4*b*c)/(3*c*d**2) - x*(a + b*x**2)**(sympy.S(3)/2)*(-a*d + b*c)/(c*d*sqrt(c + d*x**2)) - x*sqrt(a + b*x**2)*(3*a**2*d**2 - 13*a*b*c*d + 8*b**2*c**2)/(3*c*d**2*sqrt(c + d*x**2)) + sqrt(a + b*x**2)*(3*a**2*d**2 - 13*a*b*c*d + 8*b**2*c**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*sqrt(c)*d**(sympy.S(5)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_205():
    f = (a + b*x**2)**(sympy.S(3)/2)/(c + d*x**2)**(sympy.S(3)/2)
    F = b*sqrt(c)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(d**(sympy.S(3)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) - x*sqrt(a + b*x**2)*(-a*d + b*c)/(c*d*sqrt(c + d*x**2)) + x*sqrt(a + b*x**2)*(-a*d + 2*b*c)/(c*d*sqrt(c + d*x**2)) - sqrt(a + b*x**2)*(-a*d + 2*b*c)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(sqrt(c)*d**(sympy.S(3)/2)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_206():
    f = sqrt(a + b*x**2)/(c + d*x**2)**(sympy.S(3)/2)
    F = sqrt(a + b*x**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(sqrt(c)*sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_207():
    f = 1/(sqrt(a + b*x**2)*(c + d*x**2)**(sympy.S(3)/2))
    F = -sqrt(d)*sqrt(a + b*x**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(sqrt(c)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)*(-a*d + b*c)) + b*sqrt(c)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(a*sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_208():
    f = 1/((a + b*x**2)**(sympy.S(3)/2)*(c + d*x**2)**(sympy.S(3)/2))
    F = -2*b*sqrt(c)*sqrt(d)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(a*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)*(-a*d + b*c)**2) + b*x/(a*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-a*d + b*c)) + sqrt(d)*sqrt(a + b*x**2)*(a*d + b*c)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(a*sqrt(c)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_209():
    f = 1/((a + b*x**2)**(sympy.S(5)/2)*(c + d*x**2)**(sympy.S(3)/2))
    F = b*x/(3*a*(a + b*x**2)**(sympy.S(3)/2)*sqrt(c + d*x**2)*(-a*d + b*c)) - b*sqrt(c)*sqrt(d)*sqrt(a + b*x**2)*(-9*a*d + b*c)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*a**2*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)*(-a*d + b*c)**3) + 2*b*x*(-3*a*d + b*c)/(3*a**2*sqrt(a + b*x**2)*sqrt(c + d*x**2)*(-a*d + b*c)**2) + sqrt(d)*sqrt(a + b*x**2)*(-3*a**2*d**2 - 7*a*b*c*d + 2*b**2*c**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(3*a**2*sqrt(c)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)*(-a*d + b*c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_210():
    f = 1/(sqrt(a + b*x**2)*sqrt(c + d*x**2))
    F = sqrt(c)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(a*sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_211():
    f = 1/(sqrt(a - b*x**2)*sqrt(c + d*x**2))
    F = sqrt(a)*sqrt(1 - b*x**2/a)*sqrt(1 + d*x**2/c)*elliptic_f(asin(sqrt(b)*x/sqrt(a)), -a*d/(b*c))/(sqrt(b)*sqrt(a - b*x**2)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_212():
    f = 1/(sqrt(a + b*x**2)*sqrt(c - d*x**2))
    F = sqrt(c)*sqrt(1 + b*x**2/a)*sqrt(1 - d*x**2/c)*elliptic_f(asin(sqrt(d)*x/sqrt(c)), -b*c/(a*d))/(sqrt(d)*sqrt(a + b*x**2)*sqrt(c - d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_213():
    f = 1/(sqrt(a - b*x**2)*sqrt(c - d*x**2))
    F = sqrt(c)*sqrt(1 - b*x**2/a)*sqrt(1 - d*x**2/c)*elliptic_f(asin(sqrt(d)*x/sqrt(c)), b*c/(a*d))/(sqrt(d)*sqrt(a - b*x**2)*sqrt(c - d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_214():
    f = 1/(sqrt(1 - x**2)*sqrt(5*x**2 + 2))
    F = sqrt(2)*elliptic_f(asin(x), sympy.S(-5)/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_215():
    f = 1/(sqrt(1 - x**2)*sqrt(4*x**2 + 2))
    F = sqrt(2)*elliptic_f(asin(x), -2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_216():
    f = 1/(sqrt(1 - x**2)*sqrt(3*x**2 + 2))
    F = sqrt(2)*elliptic_f(asin(x), sympy.S(-3)/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_217():
    f = 1/(sqrt(1 - x**2)*sqrt(2*x**2 + 2))
    F = sqrt(2)*elliptic_f(asin(x), -1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_218():
    f = 1/(sqrt(1 - x**2)*sqrt(x**2 + 2))
    F = sqrt(2)*elliptic_f(asin(x), sympy.S(-1)/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_219():
    f = 1/(sqrt(1 - x**2)*sqrt(2 - x**2))
    F = sqrt(2)*elliptic_f(asin(x), sympy.S.Half)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_220():
    f = 1/(sqrt(1 - x**2)*sqrt(2 - 2*x**2))
    F = sqrt(2)*atanh(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_221():
    f = 1/(sqrt(1 - x**2)*sqrt(2 - 3*x**2))
    F = sqrt(2)*elliptic_f(asin(x), sympy.S(3)/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_222():
    f = 1/(sqrt(1 - x**2)*sqrt(2 - 4*x**2))
    F = sqrt(2)*elliptic_f(asin(x), 2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_223():
    f = 1/(sqrt(1 - x**2)*sqrt(2 - 5*x**2))
    F = sqrt(2)*elliptic_f(asin(x), sympy.S(5)/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_224():
    f = 1/(sqrt(x**2 + 1)*sqrt(5*x**2 + 2))
    F = sqrt(2)*sqrt(5*x**2 + 2)*elliptic_f(atan(x), sympy.S(-3)/2)/(2*sqrt((5*x**2 + 2)/(x**2 + 1))*sqrt(x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_225():
    f = 1/(sqrt(x**2 + 1)*sqrt(4*x**2 + 2))
    F = sqrt(2)*sqrt(2*x**2 + 1)*elliptic_f(atan(x), -1)/(2*sqrt((2*x**2 + 1)/(x**2 + 1))*sqrt(x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_226():
    f = 1/(sqrt(x**2 + 1)*sqrt(3*x**2 + 2))
    F = sqrt(2)*sqrt(3*x**2 + 2)*elliptic_f(atan(x), sympy.S(-1)/2)/(2*sqrt((3*x**2 + 2)/(x**2 + 1))*sqrt(x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_227():
    f = 1/(sqrt(x**2 + 1)*sqrt(2*x**2 + 2))
    F = sqrt(2)*atan(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_228():
    f = 1/(sqrt(x**2 + 1)*sqrt(x**2 + 2))
    F = sqrt(2)*sqrt(x**2 + 2)*elliptic_f(atan(x), sympy.S.Half)/(2*sqrt((x**2 + 2)/(x**2 + 1))*sqrt(x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_229():
    f = 1/(sqrt(2 - x**2)*sqrt(x**2 + 1))
    F = elliptic_f(asin(sqrt(2)*x/2), -2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_230():
    f = 1/(sqrt(2 - 2*x**2)*sqrt(x**2 + 1))
    F = sqrt(2)*elliptic_f(asin(x), -1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_231():
    f = 1/(sqrt(2 - 3*x**2)*sqrt(x**2 + 1))
    F = sqrt(3)*elliptic_f(asin(sqrt(6)*x/2), sympy.S(-2)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_232():
    f = 1/(sqrt(2 - 4*x**2)*sqrt(x**2 + 1))
    F = elliptic_f(asin(sqrt(2)*x), sympy.S(-1)/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_233():
    f = 1/(sqrt(2 - 5*x**2)*sqrt(x**2 + 1))
    F = sqrt(5)*elliptic_f(asin(sqrt(10)*x/2), sympy.S(-2)/5)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_234():
    f = 1/(sqrt(x**2 - 1)*sqrt(5*x**2 + 2))
    F = sqrt(2)*sqrt(1 - x**2)*elliptic_f(asin(x), sympy.S(-5)/2)/(2*sqrt(x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_235():
    f = 1/(sqrt(x**2 - 1)*sqrt(4*x**2 + 2))
    F = sqrt(2)*sqrt(1 - x**2)*elliptic_f(asin(x), -2)/(2*sqrt(x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_236():
    f = 1/(sqrt(x**2 - 1)*sqrt(3*x**2 + 2))
    F = sqrt(2)*sqrt(1 - x**2)*elliptic_f(asin(x), sympy.S(-3)/2)/(2*sqrt(x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_237():
    f = 1/(sqrt(x**2 - 1)*sqrt(2*x**2 + 2))
    F = elliptic_f(asin(sqrt(2)*x/sqrt(x**2 - 1)), sympy.S.Half)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_238():
    f = 1/(sqrt(x**2 - 1)*sqrt(x**2 + 2))
    F = sqrt(2)*sqrt(1 - x**2)*elliptic_f(asin(x), sympy.S(-1)/2)/(2*sqrt(x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_239():
    f = 1/(sqrt(2 - x**2)*sqrt(x**2 - 1))
    F = -elliptic_f(acos(sqrt(2)*x/2), 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_240():
    f = 1/(sqrt(2 - 2*x**2)*sqrt(x**2 - 1))
    F = -sqrt(2)*sqrt(x**2 - 1)*atanh(x)/(2*sqrt(1 - x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_241():
    f = 1/(sqrt(2 - 3*x**2)*sqrt(x**2 - 1))
    F = sqrt(2)*sqrt(1 - x**2)*elliptic_f(asin(x), sympy.S(3)/2)/(2*sqrt(x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_242():
    f = 1/(sqrt(2 - 4*x**2)*sqrt(x**2 - 1))
    F = sqrt(2)*sqrt(1 - x**2)*elliptic_f(asin(x), 2)/(2*sqrt(x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_243():
    f = 1/(sqrt(2 - 5*x**2)*sqrt(x**2 - 1))
    F = sqrt(2)*sqrt(1 - x**2)*elliptic_f(asin(x), sympy.S(5)/2)/(2*sqrt(x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_244():
    f = 1/(sqrt(-x**2 - 1)*sqrt(5*x**2 + 2))
    F = sqrt(2)*sqrt(5*x**2 + 2)*elliptic_f(atan(x), sympy.S(-3)/2)/(2*sqrt((5*x**2 + 2)/(x**2 + 1))*sqrt(-x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_245():
    f = 1/(sqrt(-x**2 - 1)*sqrt(4*x**2 + 2))
    F = sqrt(2)*sqrt(2*x**2 + 1)*elliptic_f(atan(x), -1)/(2*sqrt((2*x**2 + 1)/(x**2 + 1))*sqrt(-x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_246():
    f = 1/(sqrt(-x**2 - 1)*sqrt(3*x**2 + 2))
    F = sqrt(2)*sqrt(3*x**2 + 2)*elliptic_f(atan(x), sympy.S(-1)/2)/(2*sqrt((3*x**2 + 2)/(x**2 + 1))*sqrt(-x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_247():
    f = 1/(sqrt(-x**2 - 1)*sqrt(2*x**2 + 2))
    F = sqrt(2)*sqrt(x**2 + 1)*atan(x)/(2*sqrt(-x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_248():
    f = 1/(sqrt(-x**2 - 1)*sqrt(x**2 + 2))
    F = sqrt(2)*sqrt(x**2 + 2)*elliptic_f(atan(x), sympy.S.Half)/(2*sqrt((x**2 + 2)/(x**2 + 1))*sqrt(-x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_249():
    f = 1/(sqrt(2 - x**2)*sqrt(-x**2 - 1))
    F = sqrt(x**2 + 1)*elliptic_f(asin(sqrt(2)*x/2), -2)/sqrt(-x**2 - 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_250():
    f = 1/(sqrt(2 - 3*x**2)*sqrt(-x**2 - 1))
    F = sqrt(3)*sqrt(x**2 + 1)*elliptic_f(asin(sqrt(6)*x/2), sympy.S(-2)/3)/(3*sqrt(-x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_251():
    f = 1/(sqrt(2 - 4*x**2)*sqrt(-x**2 - 1))
    F = sqrt(x**2 + 1)*elliptic_f(asin(sqrt(2)*x), sympy.S(-1)/2)/(2*sqrt(-x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_252():
    f = 1/(sqrt(2 - 5*x**2)*sqrt(-x**2 - 1))
    F = sqrt(5)*sqrt(x**2 + 1)*elliptic_f(asin(sqrt(10)*x/2), sympy.S(-2)/5)/(5*sqrt(-x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_253():
    f = sqrt(a + b*x**2)/sqrt(c - d*x**2)
    F = sqrt(c)*sqrt(1 - d*x**2/c)*sqrt(a + b*x**2)*elliptic_e(asin(sqrt(d)*x/sqrt(c)), -b*c/(a*d))/(sqrt(d)*sqrt(1 + b*x**2/a)*sqrt(c - d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_254():
    f = sqrt(-a - b*x**2)/sqrt(c - d*x**2)
    F = sqrt(c)*sqrt(1 - d*x**2/c)*sqrt(-a - b*x**2)*elliptic_e(asin(sqrt(d)*x/sqrt(c)), -b*c/(a*d))/(sqrt(d)*sqrt(1 + b*x**2/a)*sqrt(c - d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_255():
    f = sqrt(a + b*x**2)/sqrt(-c + d*x**2)
    F = sqrt(c)*sqrt(1 - d*x**2/c)*sqrt(a + b*x**2)*elliptic_e(asin(sqrt(d)*x/sqrt(c)), -b*c/(a*d))/(sqrt(d)*sqrt(1 + b*x**2/a)*sqrt(-c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_256():
    f = sqrt(-a - b*x**2)/sqrt(-c + d*x**2)
    F = sqrt(c)*sqrt(1 - d*x**2/c)*sqrt(-a - b*x**2)*elliptic_e(asin(sqrt(d)*x/sqrt(c)), -b*c/(a*d))/(sqrt(d)*sqrt(1 + b*x**2/a)*sqrt(-c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_257():
    f = sqrt(a - b*x**2)/sqrt(c - d*x**2)
    F = sqrt(c)*sqrt(1 - d*x**2/c)*sqrt(a - b*x**2)*elliptic_e(asin(sqrt(d)*x/sqrt(c)), b*c/(a*d))/(sqrt(d)*sqrt(1 - b*x**2/a)*sqrt(c - d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_258():
    f = sqrt(-a + b*x**2)/sqrt(c - d*x**2)
    F = sqrt(c)*sqrt(1 - d*x**2/c)*sqrt(-a + b*x**2)*elliptic_e(asin(sqrt(d)*x/sqrt(c)), b*c/(a*d))/(sqrt(d)*sqrt(1 - b*x**2/a)*sqrt(c - d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_259():
    f = sqrt(a - b*x**2)/sqrt(-c + d*x**2)
    F = sqrt(c)*sqrt(1 - d*x**2/c)*sqrt(a - b*x**2)*elliptic_e(asin(sqrt(d)*x/sqrt(c)), b*c/(a*d))/(sqrt(d)*sqrt(1 - b*x**2/a)*sqrt(-c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_260():
    f = sqrt(-a + b*x**2)/sqrt(-c + d*x**2)
    F = sqrt(c)*sqrt(1 - d*x**2/c)*sqrt(-a + b*x**2)*elliptic_e(asin(sqrt(d)*x/sqrt(c)), b*c/(a*d))/(sqrt(d)*sqrt(1 - b*x**2/a)*sqrt(-c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_261():
    f = sqrt(a + b*x**2)/sqrt(c + d*x**2)
    F = -sqrt(c)*sqrt(a + b*x**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + sqrt(c)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + x*sqrt(a + b*x**2)/sqrt(c + d*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_262():
    f = sqrt(-a - b*x**2)/sqrt(c + d*x**2)
    F = -sqrt(c)*sqrt(-a - b*x**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + sqrt(c)*sqrt(-a - b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + x*sqrt(-a - b*x**2)/sqrt(c + d*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_263():
    f = sqrt(a + b*x**2)/sqrt(-c - d*x**2)
    F = -sqrt(c)*sqrt(a + b*x**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(-c - d*x**2)) + sqrt(c)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(-c - d*x**2)) + x*sqrt(a + b*x**2)/sqrt(-c - d*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_264():
    f = sqrt(-a - b*x**2)/sqrt(-c - d*x**2)
    F = -sqrt(c)*sqrt(-a - b*x**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(-c - d*x**2)) + sqrt(c)*sqrt(-a - b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(-c - d*x**2)) + x*sqrt(-a - b*x**2)/sqrt(-c - d*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_265():
    f = sqrt(a - b*x**2)/sqrt(c + d*x**2)
    F = -sqrt(a)*sqrt(b)*sqrt(1 - b*x**2/a)*sqrt(c + d*x**2)*elliptic_e(asin(sqrt(b)*x/sqrt(a)), -a*d/(b*c))/(d*sqrt(1 + d*x**2/c)*sqrt(a - b*x**2)) + sqrt(a)*sqrt(1 - b*x**2/a)*sqrt(1 + d*x**2/c)*(a*d + b*c)*elliptic_f(asin(sqrt(b)*x/sqrt(a)), -a*d/(b*c))/(sqrt(b)*d*sqrt(a - b*x**2)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_266():
    f = sqrt(-a + b*x**2)/sqrt(c + d*x**2)
    F = sqrt(a)*sqrt(b)*sqrt(1 - b*x**2/a)*sqrt(c + d*x**2)*elliptic_e(asin(sqrt(b)*x/sqrt(a)), -a*d/(b*c))/(d*sqrt(1 + d*x**2/c)*sqrt(-a + b*x**2)) - sqrt(a)*sqrt(1 - b*x**2/a)*sqrt(1 + d*x**2/c)*(a*d + b*c)*elliptic_f(asin(sqrt(b)*x/sqrt(a)), -a*d/(b*c))/(sqrt(b)*d*sqrt(-a + b*x**2)*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_267():
    f = sqrt(a - b*x**2)/sqrt(-c - d*x**2)
    F = sqrt(a)*sqrt(b)*sqrt(1 - b*x**2/a)*sqrt(-c - d*x**2)*elliptic_e(asin(sqrt(b)*x/sqrt(a)), -a*d/(b*c))/(d*sqrt(1 + d*x**2/c)*sqrt(a - b*x**2)) + sqrt(a)*sqrt(1 - b*x**2/a)*sqrt(1 + d*x**2/c)*(a*d + b*c)*elliptic_f(asin(sqrt(b)*x/sqrt(a)), -a*d/(b*c))/(sqrt(b)*d*sqrt(a - b*x**2)*sqrt(-c - d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_268():
    f = sqrt(-a + b*x**2)/sqrt(-c - d*x**2)
    F = -sqrt(a)*sqrt(b)*sqrt(1 - b*x**2/a)*sqrt(-c - d*x**2)*elliptic_e(asin(sqrt(b)*x/sqrt(a)), -a*d/(b*c))/(d*sqrt(1 + d*x**2/c)*sqrt(-a + b*x**2)) - sqrt(a)*sqrt(1 - b*x**2/a)*sqrt(1 + d*x**2/c)*(a*d + b*c)*elliptic_f(asin(sqrt(b)*x/sqrt(a)), -a*d/(b*c))/(sqrt(b)*d*sqrt(-a + b*x**2)*sqrt(-c - d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_269():
    f = sqrt(c + d*x**2)/sqrt(a - b*x**2)
    F = sqrt(a)*sqrt(1 - b*x**2/a)*sqrt(c + d*x**2)*elliptic_e(asin(sqrt(b)*x/sqrt(a)), -a*d/(b*c))/(sqrt(b)*sqrt(1 + d*x**2/c)*sqrt(a - b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_270():
    f = sqrt(-c - d*x**2)/sqrt(a - b*x**2)
    F = sqrt(a)*sqrt(1 - b*x**2/a)*sqrt(-c - d*x**2)*elliptic_e(asin(sqrt(b)*x/sqrt(a)), -a*d/(b*c))/(sqrt(b)*sqrt(1 + d*x**2/c)*sqrt(a - b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_271():
    f = sqrt(c + d*x**2)/sqrt(-a + b*x**2)
    F = sqrt(a)*sqrt(1 - b*x**2/a)*sqrt(c + d*x**2)*elliptic_e(asin(sqrt(b)*x/sqrt(a)), -a*d/(b*c))/(sqrt(b)*sqrt(1 + d*x**2/c)*sqrt(-a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_272():
    f = sqrt(-c - d*x**2)/sqrt(-a + b*x**2)
    F = sqrt(a)*sqrt(1 - b*x**2/a)*sqrt(-c - d*x**2)*elliptic_e(asin(sqrt(b)*x/sqrt(a)), -a*d/(b*c))/(sqrt(b)*sqrt(1 + d*x**2/c)*sqrt(-a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_273():
    f = sqrt(c - d*x**2)/sqrt(a - b*x**2)
    F = sqrt(a)*sqrt(1 - b*x**2/a)*sqrt(c - d*x**2)*elliptic_e(asin(sqrt(b)*x/sqrt(a)), a*d/(b*c))/(sqrt(b)*sqrt(1 - d*x**2/c)*sqrt(a - b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_274():
    f = sqrt(-c + d*x**2)/sqrt(a - b*x**2)
    F = sqrt(a)*sqrt(1 - b*x**2/a)*sqrt(-c + d*x**2)*elliptic_e(asin(sqrt(b)*x/sqrt(a)), a*d/(b*c))/(sqrt(b)*sqrt(1 - d*x**2/c)*sqrt(a - b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_275():
    f = sqrt(c - d*x**2)/sqrt(-a + b*x**2)
    F = sqrt(a)*sqrt(1 - b*x**2/a)*sqrt(c - d*x**2)*elliptic_e(asin(sqrt(b)*x/sqrt(a)), a*d/(b*c))/(sqrt(b)*sqrt(1 - d*x**2/c)*sqrt(-a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_276():
    f = sqrt(-c + d*x**2)/sqrt(-a + b*x**2)
    F = sqrt(a)*sqrt(1 - b*x**2/a)*sqrt(-c + d*x**2)*elliptic_e(asin(sqrt(b)*x/sqrt(a)), a*d/(b*c))/(sqrt(b)*sqrt(1 - d*x**2/c)*sqrt(-a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_277():
    f = sqrt(c + d*x**2)/sqrt(a + b*x**2)
    F = -sqrt(c)*sqrt(d)*sqrt(a + b*x**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(b*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) + d*x*sqrt(a + b*x**2)/(b*sqrt(c + d*x**2)) + c**(sympy.S(3)/2)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(a*sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_278():
    f = sqrt(-c - d*x**2)/sqrt(a + b*x**2)
    F = sqrt(c)*sqrt(d)*sqrt(a + b*x**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(b*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(-c - d*x**2)) - d*x*sqrt(a + b*x**2)/(b*sqrt(-c - d*x**2)) - c**(sympy.S(3)/2)*sqrt(a + b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(a*sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(-c - d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_279():
    f = sqrt(c + d*x**2)/sqrt(-a - b*x**2)
    F = sqrt(c)*sqrt(d)*sqrt(-a - b*x**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(b*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2)) - d*x*sqrt(-a - b*x**2)/(b*sqrt(c + d*x**2)) - c**(sympy.S(3)/2)*sqrt(-a - b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(a*sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_280():
    f = sqrt(-c - d*x**2)/sqrt(-a - b*x**2)
    F = -sqrt(c)*sqrt(d)*sqrt(-a - b*x**2)*elliptic_e(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(b*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(-c - d*x**2)) + d*x*sqrt(-a - b*x**2)/(b*sqrt(-c - d*x**2)) + c**(sympy.S(3)/2)*sqrt(-a - b*x**2)*elliptic_f(atan(sqrt(d)*x/sqrt(c)), 1 - b*c/(a*d))/(a*sqrt(d)*sqrt(c*(a + b*x**2)/(a*(c + d*x**2)))*sqrt(-c - d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_281():
    f = sqrt(c - d*x**2)/sqrt(a + b*x**2)
    F = -sqrt(c)*sqrt(d)*sqrt(1 - d*x**2/c)*sqrt(a + b*x**2)*elliptic_e(asin(sqrt(d)*x/sqrt(c)), -b*c/(a*d))/(b*sqrt(1 + b*x**2/a)*sqrt(c - d*x**2)) + sqrt(c)*sqrt(1 + b*x**2/a)*sqrt(1 - d*x**2/c)*(a*d + b*c)*elliptic_f(asin(sqrt(d)*x/sqrt(c)), -b*c/(a*d))/(b*sqrt(d)*sqrt(a + b*x**2)*sqrt(c - d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_282():
    f = sqrt(-c + d*x**2)/sqrt(a + b*x**2)
    F = sqrt(c)*sqrt(d)*sqrt(1 - d*x**2/c)*sqrt(a + b*x**2)*elliptic_e(asin(sqrt(d)*x/sqrt(c)), -b*c/(a*d))/(b*sqrt(1 + b*x**2/a)*sqrt(-c + d*x**2)) - sqrt(c)*sqrt(1 + b*x**2/a)*sqrt(1 - d*x**2/c)*(a*d + b*c)*elliptic_f(asin(sqrt(d)*x/sqrt(c)), -b*c/(a*d))/(b*sqrt(d)*sqrt(a + b*x**2)*sqrt(-c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_283():
    f = sqrt(c - d*x**2)/sqrt(-a - b*x**2)
    F = sqrt(c)*sqrt(d)*sqrt(1 - d*x**2/c)*sqrt(-a - b*x**2)*elliptic_e(asin(sqrt(d)*x/sqrt(c)), -b*c/(a*d))/(b*sqrt(1 + b*x**2/a)*sqrt(c - d*x**2)) + sqrt(c)*sqrt(1 + b*x**2/a)*sqrt(1 - d*x**2/c)*(a*d + b*c)*elliptic_f(asin(sqrt(d)*x/sqrt(c)), -b*c/(a*d))/(b*sqrt(d)*sqrt(-a - b*x**2)*sqrt(c - d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_284():
    f = sqrt(-c + d*x**2)/sqrt(-a - b*x**2)
    F = -sqrt(c)*sqrt(d)*sqrt(1 - d*x**2/c)*sqrt(-a - b*x**2)*elliptic_e(asin(sqrt(d)*x/sqrt(c)), -b*c/(a*d))/(b*sqrt(1 + b*x**2/a)*sqrt(-c + d*x**2)) - sqrt(c)*sqrt(1 + b*x**2/a)*sqrt(1 - d*x**2/c)*(a*d + b*c)*elliptic_f(asin(sqrt(d)*x/sqrt(c)), -b*c/(a*d))/(b*sqrt(d)*sqrt(-a - b*x**2)*sqrt(-c + d*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_285():
    f = 1/(sqrt(b*x**2 + 2)*sqrt(d*x**2 + 3))
    F = sqrt(2)*sqrt(b*x**2 + 2)*elliptic_f(atan(sqrt(3)*sqrt(d)*x/3), -3*b/(2*d) + 1)/(2*sqrt(d)*sqrt((b*x**2 + 2)/(d*x**2 + 3))*sqrt(d*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_286():
    f = 1/(sqrt(4 - x**2)*sqrt(c + d*x**2))
    F = sqrt(1 + d*x**2/c)*elliptic_f(asin(x/2), -4*d/c)/sqrt(c + d*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_287():
    f = 1/(sqrt(c + d*x**2)*sqrt(x**2 + 4))
    F = sqrt(c + d*x**2)*elliptic_f(atan(x/2), 1 - 4*d/c)/(c*sqrt((c + d*x**2)/(c*(x**2 + 4)))*sqrt(x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_288():
    f = 1/(sqrt(1 - x**2)*sqrt(2*x**2 - 1))
    F = -elliptic_f(acos(x), 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_289():
    f = sqrt(-c**2*x**2 + 1)/sqrt(c**2*x**2 + 1)
    F = -elliptic_e(asin(c*x), -1)/c + 2*elliptic_f(asin(c*x), -1)/c
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_290():
    f = sqrt(b*x**2 + 2)/sqrt(d*x**2 + 3)
    F = x*sqrt(b*x**2 + 2)/sqrt(d*x**2 + 3) - sqrt(2)*sqrt(b*x**2 + 2)*elliptic_e(atan(sqrt(3)*sqrt(d)*x/3), -3*b/(2*d) + 1)/(sqrt(d)*sqrt((b*x**2 + 2)/(d*x**2 + 3))*sqrt(d*x**2 + 3)) + sqrt(2)*sqrt(b*x**2 + 2)*elliptic_f(atan(sqrt(3)*sqrt(d)*x/3), -3*b/(2*d) + 1)/(sqrt(d)*sqrt((b*x**2 + 2)/(d*x**2 + 3))*sqrt(d*x**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_291():
    f = sqrt(3*x**2 - 1)/sqrt(2 - 3*x**2)
    F = -sqrt(3)*elliptic_e(acos(sqrt(6)*x/2), 2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_292():
    f = sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)/sqrt(-2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)
    F = sqrt(2)*sqrt(b + sqrt(-4*a*c + b**2))*elliptic_e(asin(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2))), -(b + sqrt(-4*a*c + b**2))/(b - sqrt(-4*a*c + b**2)))/(2*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_293():
    f = sqrt(-2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)/sqrt(-2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)
    F = sqrt(2)*sqrt(b + sqrt(-4*a*c + b**2))*elliptic_e(asin(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2))), (b + sqrt(-4*a*c + b**2))/(b - sqrt(-4*a*c + b**2)))/(2*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_294():
    f = sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)/sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)
    F = x*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)/sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1) - sqrt(2)*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*elliptic_e(atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2))), -2*sqrt(-4*a*c + b**2)/(b - sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt((2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)/(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)) + sqrt(2)*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)*elliptic_f(atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2))), -2*sqrt(-4*a*c + b**2)/(b - sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt((2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)/(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))*sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_295():
    f = sqrt(-2*c*x**2/(b - sqrt(-4*a*c + b**2)) + 1)/sqrt(2*c*x**2/(b + sqrt(-4*a*c + b**2)) + 1)
    F = sqrt(2)*b*elliptic_f(asin(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2))), -(b - sqrt(-4*a*c + b**2))/(b + sqrt(-4*a*c + b**2)))/(sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))) - sqrt(2)*(b + sqrt(-4*a*c + b**2))*elliptic_e(asin(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2))), -(b - sqrt(-4*a*c + b**2))/(b + sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_296():
    f = (1 - 2*x**2)**m/sqrt(1 - x**2)
    F = -2**(-m - 2)*(2 - 4*x**2)**(m + 1)*sqrt(x**2)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), (1 - 2*x**2)**2)/(x*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_297():
    f = 1/(sqrt(x**2 - 1)*sqrt(x**2 - 4*sqrt(3) + 7))
    F = sqrt(1 - x**2)*elliptic_f(asin(x), -7 - 4*sqrt(3))/(sqrt(7 - 4*sqrt(3))*sqrt(x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_298():
    f = 1/(sqrt(x**2*(-3 + sqrt(3)) + 3)*sqrt(2*sqrt(3)*x**2 - 3*sqrt(3) + 3))
    F = -sqrt(sqrt(3) + 3)*elliptic_f(acos(x*sqrt(1 - sqrt(3)/3)), sympy.S.Half + sqrt(3)/2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_299():
    f = 1/((3*x**2 + 2)**(sympy.S(1)/4)*(3*x**2 + 4))
    F = -2**(sympy.S(1)/4)*sqrt(3)*atan(sqrt(3)*(2*2**(sympy.S(1)/4)*sqrt(3*x**2 + 2) + 2*2**(sympy.S(3)/4))/(6*x*(3*x**2 + 2)**(sympy.S(1)/4)))/12 - 2**(sympy.S(1)/4)*sqrt(3)*atanh(sqrt(3)*(-2*2**(sympy.S(1)/4)*sqrt(3*x**2 + 2) + 2*2**(sympy.S(3)/4))/(6*x*(3*x**2 + 2)**(sympy.S(1)/4)))/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_300():
    f = 1/((2 - 3*x**2)**(sympy.S(1)/4)*(4 - 3*x**2))
    F = 2**(sympy.S(1)/4)*sqrt(3)*atan(2**(sympy.S(3)/4)*sqrt(3)*(-sqrt(2)*sqrt(2 - 3*x**2) + 2)/(6*x*(2 - 3*x**2)**(sympy.S(1)/4)))/12 + 2**(sympy.S(1)/4)*sqrt(3)*atanh(2**(sympy.S(3)/4)*sqrt(3)*(sqrt(2)*sqrt(2 - 3*x**2) + 2)/(6*x*(2 - 3*x**2)**(sympy.S(1)/4)))/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_301():
    f = 1/((b*x**2 + 2)**(sympy.S(1)/4)*(b*x**2 + 4))
    F = -2**(sympy.S(1)/4)*atan((2*2**(sympy.S(1)/4)*sqrt(b*x**2 + 2) + 2*2**(sympy.S(3)/4))/(2*sqrt(b)*x*(b*x**2 + 2)**(sympy.S(1)/4)))/(4*sqrt(b)) - 2**(sympy.S(1)/4)*atanh((-2*2**(sympy.S(1)/4)*sqrt(b*x**2 + 2) + 2*2**(sympy.S(3)/4))/(2*sqrt(b)*x*(b*x**2 + 2)**(sympy.S(1)/4)))/(4*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_302():
    f = 1/((-b*x**2 + 2)**(sympy.S(1)/4)*(-b*x**2 + 4))
    F = 2**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*(-sqrt(2)*sqrt(-b*x**2 + 2) + 2)/(2*sqrt(b)*x*(-b*x**2 + 2)**(sympy.S(1)/4)))/(4*sqrt(b)) + 2**(sympy.S(1)/4)*atanh(2**(sympy.S(3)/4)*(sqrt(2)*sqrt(-b*x**2 + 2) + 2)/(2*sqrt(b)*x*(-b*x**2 + 2)**(sympy.S(1)/4)))/(4*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_303():
    f = 1/((a + 3*x**2)**(sympy.S(1)/4)*(2*a + 3*x**2))
    F = -sqrt(3)*atan(sqrt(3)*a**(sympy.S(3)/4)*(1 + sqrt(a + 3*x**2)/sqrt(a))/(3*x*(a + 3*x**2)**(sympy.S(1)/4)))/(6*a**(sympy.S(3)/4)) - sqrt(3)*atanh(sqrt(3)*a**(sympy.S(3)/4)*(1 - sqrt(a + 3*x**2)/sqrt(a))/(3*x*(a + 3*x**2)**(sympy.S(1)/4)))/(6*a**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_304():
    f = 1/((a - 3*x**2)**(sympy.S(1)/4)*(2*a - 3*x**2))
    F = sqrt(3)*atan(sqrt(3)*a**(sympy.S(3)/4)*(1 - sqrt(a - 3*x**2)/sqrt(a))/(3*x*(a - 3*x**2)**(sympy.S(1)/4)))/(6*a**(sympy.S(3)/4)) + sqrt(3)*atanh(sqrt(3)*a**(sympy.S(3)/4)*(1 + sqrt(a - 3*x**2)/sqrt(a))/(3*x*(a - 3*x**2)**(sympy.S(1)/4)))/(6*a**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_305():
    f = 1/((a + b*x**2)**(sympy.S(1)/4)*(2*a + b*x**2))
    F = -atan(a**(sympy.S(3)/4)*(1 + sqrt(a + b*x**2)/sqrt(a))/(sqrt(b)*x*(a + b*x**2)**(sympy.S(1)/4)))/(2*a**(sympy.S(3)/4)*sqrt(b)) - atanh(a**(sympy.S(3)/4)*(1 - sqrt(a + b*x**2)/sqrt(a))/(sqrt(b)*x*(a + b*x**2)**(sympy.S(1)/4)))/(2*a**(sympy.S(3)/4)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_306():
    f = 1/((a - b*x**2)**(sympy.S(1)/4)*(2*a - b*x**2))
    F = atan(a**(sympy.S(3)/4)*(1 - sqrt(a - b*x**2)/sqrt(a))/(sqrt(b)*x*(a - b*x**2)**(sympy.S(1)/4)))/(2*a**(sympy.S(3)/4)*sqrt(b)) + atanh(a**(sympy.S(3)/4)*(1 + sqrt(a - b*x**2)/sqrt(a))/(sqrt(b)*x*(a - b*x**2)**(sympy.S(1)/4)))/(2*a**(sympy.S(3)/4)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_307():
    f = 1/((3*x**2 - 2)*(3*x**2 - 1)**(sympy.S(1)/4))
    F = -sqrt(6)*atan(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/12 - sqrt(6)*atanh(sqrt(6)*x/(2*(3*x**2 - 1)**(sympy.S(1)/4)))/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_308():
    f = 1/((-3*x**2 - 2)*(-3*x**2 - 1)**(sympy.S(1)/4))
    F = -sqrt(6)*atan(sqrt(6)*x/(2*(-3*x**2 - 1)**(sympy.S(1)/4)))/12 - sqrt(6)*atanh(sqrt(6)*x/(2*(-3*x**2 - 1)**(sympy.S(1)/4)))/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_309():
    f = 1/((b*x**2 - 2)*(b*x**2 - 1)**(sympy.S(1)/4))
    F = -sqrt(2)*atan(sqrt(2)*sqrt(b)*x/(2*(b*x**2 - 1)**(sympy.S(1)/4)))/(4*sqrt(b)) - sqrt(2)*atanh(sqrt(2)*sqrt(b)*x/(2*(b*x**2 - 1)**(sympy.S(1)/4)))/(4*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_310():
    f = 1/((-b*x**2 - 2)*(-b*x**2 - 1)**(sympy.S(1)/4))
    F = -sqrt(2)*atan(sqrt(2)*sqrt(b)*x/(2*(-b*x**2 - 1)**(sympy.S(1)/4)))/(4*sqrt(b)) - sqrt(2)*atanh(sqrt(2)*sqrt(b)*x/(2*(-b*x**2 - 1)**(sympy.S(1)/4)))/(4*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_311():
    f = 1/((-2*a + 3*x**2)*(-a + 3*x**2)**(sympy.S(1)/4))
    F = -sqrt(6)*atan(sqrt(6)*x/(2*a**(sympy.S(1)/4)*(-a + 3*x**2)**(sympy.S(1)/4)))/(12*a**(sympy.S(3)/4)) - sqrt(6)*atanh(sqrt(6)*x/(2*a**(sympy.S(1)/4)*(-a + 3*x**2)**(sympy.S(1)/4)))/(12*a**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_312():
    f = 1/((-2*a - 3*x**2)*(-a - 3*x**2)**(sympy.S(1)/4))
    F = -sqrt(6)*atan(sqrt(6)*x/(2*a**(sympy.S(1)/4)*(-a - 3*x**2)**(sympy.S(1)/4)))/(12*a**(sympy.S(3)/4)) - sqrt(6)*atanh(sqrt(6)*x/(2*a**(sympy.S(1)/4)*(-a - 3*x**2)**(sympy.S(1)/4)))/(12*a**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_313():
    f = 1/((-2*a + b*x**2)*(-a + b*x**2)**(sympy.S(1)/4))
    F = -sqrt(2)*atan(sqrt(2)*sqrt(b)*x/(2*a**(sympy.S(1)/4)*(-a + b*x**2)**(sympy.S(1)/4)))/(4*a**(sympy.S(3)/4)*sqrt(b)) - sqrt(2)*atanh(sqrt(2)*sqrt(b)*x/(2*a**(sympy.S(1)/4)*(-a + b*x**2)**(sympy.S(1)/4)))/(4*a**(sympy.S(3)/4)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_314():
    f = 1/((-2*a - b*x**2)*(-a - b*x**2)**(sympy.S(1)/4))
    F = -sqrt(2)*atan(sqrt(2)*sqrt(b)*x/(2*a**(sympy.S(1)/4)*(-a - b*x**2)**(sympy.S(1)/4)))/(4*a**(sympy.S(3)/4)*sqrt(b)) - sqrt(2)*atanh(sqrt(2)*sqrt(b)*x/(2*a**(sympy.S(1)/4)*(-a - b*x**2)**(sympy.S(1)/4)))/(4*a**(sympy.S(3)/4)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_315():
    f = 1/((2 - x**2)*(x**2 - 1)**(sympy.S(1)/4))
    F = sqrt(2)*atan(sqrt(2)*x/(2*(x**2 - 1)**(sympy.S(1)/4)))/4 + sqrt(2)*atanh(sqrt(2)*x/(2*(x**2 - 1)**(sympy.S(1)/4)))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_316():
    f = (a + b*x**2)**(sympy.S(7)/4)/(c + d*x**2)
    F = ((Integer(6) * Symbol('a') * Symbol('b') * x) * ((Integer(5) * Symbol('d') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * x) * (((Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * x * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1))))) * ((Integer(5) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * ((Integer(1) + ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1)))))**((Integer(4))**(Integer(-1))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(5) * Symbol('d') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Integer(1) + ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1)))))**((Integer(4))**(Integer(-1))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * (((Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1)))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * (((Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * x))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * (((Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_317():
    f = (a + b*x**2)**(sympy.S(5)/4)/(c + d*x**2)
    F = ((Integer(2) * Symbol('b') * x * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1)))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * ((Integer(1) + ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(3) * Symbol('d') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Integer(1) + ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * (((Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1)))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * (((Symbol('d'))**(Integer(2)) * x))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * (((Symbol('d'))**(Integer(2)) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_318():
    f = (a + b*x**2)**(sympy.S(3)/4)/(c + d*x**2)
    F = ((Integer(2) * Symbol('b') * x) * ((Symbol('d') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * ((Integer(1) + ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1)))))**((Integer(4))**(Integer(-1))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Symbol('d') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1)))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * (((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * x))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * (((Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_319():
    f = (a + b*x**2)**(sympy.S(1)/4)/(c + d*x**2)
    F = ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * ((Integer(1) + ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Symbol('d') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1)))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Symbol('d') * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Symbol('d') * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_320():
    f = 1/((a + b*x**2)**(sympy.S(1)/4)*(c + d*x**2))
    F = (((Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1)))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * x))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_321():
    f = 1/((a + b*x**2)**(sympy.S(3)/4)*(c + d*x**2))
    F = (((Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1)))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * x))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_322():
    f = 1/((a + b*x**2)**(sympy.S(5)/4)*(c + d*x**2))
    F = ((Integer(2) * sympy.sqrt(Symbol('b')) * ((Integer(1) + ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1)))))**((Integer(4))**(Integer(-1))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((sympy.sqrt(Symbol('a')) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1)))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * x))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_323():
    f = 1/((a + b*x**2)**(sympy.S(7)/4)*(c + d*x**2))
    F = ((Integer(2) * Symbol('b') * x) * ((Integer(3) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt(Symbol('b')) * ((Integer(1) + ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(3) * sympy.sqrt(Symbol('a')) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * Symbol('d') * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1)))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * Symbol('d') * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_324():
    f = 1/((a + b*x**2)**(sympy.S(9)/4)*(c + d*x**2))
    F = ((Integer(2) * Symbol('b') * x) * ((Integer(5) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(5) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt(Symbol('b')) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(8) * Symbol('a') * Symbol('d')))) * ((Integer(1) + ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1)))))**((Integer(4))**(Integer(-1))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(5) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1)))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * x))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_325():
    f = 1/((a + b*x**2)**(sympy.S(11)/4)*(c + d*x**2))
    F = ((Integer(2) * Symbol('b') * x) * ((Integer(7) * Symbol('a') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(7) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(12) * Symbol('a') * Symbol('d')))) * x) * ((Integer(21) * (Symbol('a'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt(Symbol('b')) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(12) * Symbol('a') * Symbol('d')))) * ((Integer(1) + ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(21) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1)))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * x))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_326():
    f = (a + b*x**2)**(sympy.S(7)/4)/(c + d*x**2)**2
    F = ((Symbol('b') * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * x) * ((Integer(2) * Symbol('c') * (Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * x * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1))))) * ((Integer(2) * Symbol('c') * Symbol('d') * (Symbol('c') + (Symbol('d') * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Integer(1) + ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1)))))**((Integer(4))**(Integer(-1))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(2) * Symbol('c') * (Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(2) * Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1)))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * x))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(2) * Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_327():
    f = (a + b*x**2)**(sympy.S(5)/4)/(c + d*x**2)**2
    F = (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * x * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1)))) * ((Integer(2) * Symbol('c') * Symbol('d') * (Symbol('c') + (Symbol('d') * (x)**(Integer(2))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * ((Integer(1) + ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(2) * Symbol('c') * (Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(2) * Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1)))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * (Symbol('d'))**(Integer(2)) * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(2) * Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * (Symbol('d'))**(Integer(2)) * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_328():
    f = (a + b*x**2)**(sympy.S(3)/4)/(c + d*x**2)**2
    F = (Integer(-1) * ((Symbol('b') * x) * ((Integer(2) * Symbol('c') * Symbol('d') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + ((x * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1))))) * ((Integer(2) * Symbol('c') * (Symbol('c') + (Symbol('d') * (x)**(Integer(2))))))**(Integer(-1))) + ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * ((Integer(1) + ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1)))))**((Integer(4))**(Integer(-1))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(2) * Symbol('c') * Symbol('d') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(2) * Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1)))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * x))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(2) * Symbol('a') * Symbol('d'))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))) * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_329():
    f = (a + b*x**2)**(sympy.S(1)/4)/(c + d*x**2)**2
    F = ((x * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1)))) * ((Integer(2) * Symbol('c') * (Symbol('c') + (Symbol('d') * (x)**(Integer(2))))))**(Integer(-1))) + ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * ((Integer(1) + ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(2) * Symbol('c') * Symbol('d') * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1)))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_330():
    f = 1/((a + b*x**2)**(sympy.S(1)/4)*(c + d*x**2)**2)
    F = ((Symbol('b') * x) * ((Integer(2) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * x * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1))))) * ((Integer(2) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('c') + (Symbol('d') * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * ((Integer(1) + ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1)))))**((Integer(4))**(Integer(-1))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(2) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1)))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * sympy.sqrt(Symbol('d')) * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * x))**(Integer(-1)))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * sympy.sqrt(Symbol('d')) * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_331():
    f = 1/((a + b*x**2)**(sympy.S(3)/4)*(c + d*x**2)**2)
    F = (Integer(-1) * ((Symbol('d') * x * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1)))) * ((Integer(2) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('c') + (Symbol('d') * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * ((Integer(1) + ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(2) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1)))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * x))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_332():
    f = 1/((a + b*x**2)**(sympy.S(5)/4)*(c + d*x**2)**2)
    F = (Integer(-1) * ((Symbol('d') * x) * ((Integer(2) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * (Symbol('c') + (Symbol('d') * (x)**(Integer(2))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('b')) * ((Integer(4) * Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * ((Integer(1) + ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1)))))**((Integer(4))**(Integer(-1))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('d')) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1)))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * x))**(Integer(-1)))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('d')) * ((Integer(7) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_333():
    f = 1/((a + b*x**2)**(sympy.S(7)/4)*(c + d*x**2)**2)
    F = ((Symbol('b') * ((Integer(4) * Symbol('b') * Symbol('c')) + (Integer(3) * Symbol('a') * Symbol('d'))) * x) * ((Integer(6) * Symbol('a') * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * x) * ((Integer(2) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * (x)**(Integer(2))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('b')) * ((Integer(4) * Symbol('b') * Symbol('c')) + (Integer(3) * Symbol('a') * Symbol('d'))) * ((Integer(1) + ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(6) * sympy.sqrt(Symbol('a')) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * Symbol('d') * ((Integer(9) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1)))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * Symbol('d') * ((Integer(9) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_334():
    f = 1/((a + b*x**2)**(sympy.S(9)/4)*(c + d*x**2)**2)
    F = ((Symbol('b') * ((Integer(4) * Symbol('b') * Symbol('c')) + (Integer(5) * Symbol('a') * Symbol('d'))) * x) * ((Integer(10) * Symbol('a') * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(5) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * x) * ((Integer(2) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * (x)**(Integer(2))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('b')) * ((Integer(12) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(52) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))) * ((Integer(1) + ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1)))))**((Integer(4))**(Integer(-1))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(10) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(11) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1)))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(7) * (Integer(2))**(Integer(-1)))) * x))**(Integer(-1)))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(11) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * ((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))))**((Integer(7) * (Integer(2))**(Integer(-1)))) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_335():
    f = 1/((a + b*x**2)**(sympy.S(11)/4)*(c + d*x**2)**2)
    F = ((Symbol('b') * ((Integer(4) * Symbol('b') * Symbol('c')) + (Integer(7) * Symbol('a') * Symbol('d'))) * x) * ((Integer(14) * Symbol('a') * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(7) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(20) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(76) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(-1) * (Integer(21) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))) * x) * ((Integer(42) * (Symbol('a'))**(Integer(2)) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * x) * ((Integer(2) * Symbol('c') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * (x)**(Integer(2))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('b')) * ((Integer(20) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(76) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(-1) * (Integer(21) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))) * ((Integer(1) + ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))), Integer(2))) * ((Integer(42) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**(Integer(2)) * ((Integer(13) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1)))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * x))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**(Integer(2)) * ((Integer(13) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Symbol('a'))**(Integer(-1))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('d'))) * (sympy.sqrt((((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d')))))**(Integer(-1))), sympy.asin((((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**((Integer(4))**(Integer(-1))) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(4) * Symbol('c') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_336():
    f = (a + b*x**2)**p*(c + d*x**2)**q
    F = x*(a + b*x**2)**p*(c + d*x**2)**q*appellf1(sympy.S.Half, -p, -q, sympy.S(3)/2, -b*x**2/a, -d*x**2/c)/((1 + b*x**2/a)**p*(1 + d*x**2/c)**q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_337():
    f = (a + b*x**2)**p*(c + d*x**2)**3
    F = d*x*(a + b*x**2)**(p + 1)*(c + d*x**2)**2/(b*(2*p + 7)) - d*x*(a + b*x**2)**(p + 1)*(c + d*x**2)*(5*a*d - b*c*(2*p + 11))/(b**2*(2*p + 5)*(2*p + 7)) + d*x*(a + b*x**2)**(p + 1)*(15*a**2*d**2 - 8*a*b*c*d*(p + 6) + b**2*c**2*(4*p**2 + 28*p + 57))/(b**3*(2*p + 3)*(2*p + 5)*(2*p + 7)) - x*(a + b*x**2)**p*(15*a**3*d**3 - 9*a**2*b*c*d**2*(2*p + 7) + 3*a*b**2*c**2*d*(4*p**2 + 24*p + 35) - b**3*c**3*(8*p**3 + 60*p**2 + 142*p + 105))*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*x**2/a)/(b**3*(1 + b*x**2/a)**p*(2*p + 3)*(2*p + 5)*(2*p + 7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_338():
    f = (a + b*x**2)**p*(c + d*x**2)**2
    F = d*x*(a + b*x**2)**(p + 1)*(c + d*x**2)/(b*(2*p + 5)) - d*x*(a + b*x**2)**(p + 1)*(3*a*d - b*c*(2*p + 7))/(b**2*(2*p + 3)*(2*p + 5)) + x*(a + b*x**2)**p*(3*a**2*d**2 - 2*a*b*c*d*(2*p + 5) + b**2*c**2*(4*p**2 + 16*p + 15))*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*x**2/a)/(b**2*(1 + b*x**2/a)**p*(2*p + 3)*(2*p + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_339():
    f = (a + b*x**2)**p*(c + d*x**2)
    F = d*x*(a + b*x**2)**(p + 1)/(b*(2*p + 3)) - x*(a + b*x**2)**p*(a*d - b*c*(2*p + 3))*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*x**2/a)/(b*(1 + b*x**2/a)**p*(2*p + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_340():
    f = (a + b*x**2)**p
    F = x*(a + b*x**2)**p*hyper((sympy.S.Half, -p), (sympy.S(3)/2,), -b*x**2/a)/(1 + b*x**2/a)**p
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_341():
    f = (a + b*x**2)**p/(c + d*x**2)
    F = x*(a + b*x**2)**p*appellf1(sympy.S.Half, 1, -p, sympy.S(3)/2, -d*x**2/c, -b*x**2/a)/(c*(1 + b*x**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_342():
    f = (a + b*x**2)**p/(c + d*x**2)**2
    F = x*(a + b*x**2)**p*appellf1(sympy.S.Half, 2, -p, sympy.S(3)/2, -d*x**2/c, -b*x**2/a)/(c**2*(1 + b*x**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_343():
    f = (a + b*x**2)**p/(c + d*x**2)**3
    F = x*(a + b*x**2)**p*appellf1(sympy.S.Half, 3, -p, sympy.S(3)/2, -d*x**2/c, -b*x**2/a)/(c**3*(1 + b*x**2/a)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_1_Binomial_products_1_1_2_Quadratic_1_1_2_3_a_plus_b_x_pow_2_pow_p_c_plus_d_x_pow_2_pow_q_344():
    f = (a + b*x**2)**(-b*c/(-2*a*d + 2*b*c) - 1)*(c + d*x**2)**(a*d/(-2*a*d + 2*b*c) - 1)
    F = x*(c + d*x**2)**(a*d/(-2*a*d + 2*b*c))/(a*c*(a + b*x**2)**(b*c/(-2*a*d + 2*b*c)))
    assert integrate(f, x) == F

