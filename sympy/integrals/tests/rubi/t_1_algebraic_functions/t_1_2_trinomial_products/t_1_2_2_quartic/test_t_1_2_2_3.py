"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.2 Trinomial products/1.2.2 Quartic/1.2.2.3 (d+e x^2)^m (a+b x^2+c x^4)^p.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, a, b, c, d, e, f, g, p, q = symbols('A B a b c d e f g p q')

def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_1():
    f = (c + d*x**2)/(a + b*x**4)
    F = -sqrt(2)*(-sqrt(a)*d + sqrt(b)*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(-sqrt(a)*d + sqrt(b)*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)) - sqrt(2)*(sqrt(a)*d + sqrt(b)*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(sqrt(a)*d + sqrt(b)*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_2():
    f = (c - d*x**2)/(a + b*x**4)
    F = -sqrt(2)*(-sqrt(a)*d + sqrt(b)*c)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(-sqrt(a)*d + sqrt(b)*c)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)) - sqrt(2)*(sqrt(a)*d + sqrt(b)*c)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)) + sqrt(2)*(sqrt(a)*d + sqrt(b)*c)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_3():
    f = (c + d*x**2)/(a - b*x**4)
    F = (-sqrt(a)*d + sqrt(b)*c)*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)) + (sqrt(a)*d + sqrt(b)*c)*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_4():
    f = (c - d*x**2)/(a - b*x**4)
    F = (-sqrt(a)*d + sqrt(b)*c)*atanh(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)) + (sqrt(a)*d + sqrt(b)*c)*atan(b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*b**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_5():
    f = (3*x**2 + 2)/(9*x**4 + 4)
    F = sqrt(3)*atan(sqrt(3)*x - 1)/6 + sqrt(3)*atan(sqrt(3)*x + 1)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_6():
    f = (2 - 3*x**2)/(9*x**4 + 4)
    F = -sqrt(3)*log(3*x**2 - 2*sqrt(3)*x + 2)/12 + sqrt(3)*log(3*x**2 + 2*sqrt(3)*x + 2)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_7():
    f = (3*x**2 + 2)/(4 - 9*x**4)
    F = sqrt(6)*atanh(sqrt(6)*x/2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_8():
    f = (2 - 3*x**2)/(4 - 9*x**4)
    F = sqrt(6)*atan(sqrt(6)*x/2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_9():
    f = (sqrt(a)*sqrt(b) + b*x**2)/(a + b*x**4)
    F = -sqrt(2)*b**(sympy.S(1)/4)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)) + sqrt(2)*b**(sympy.S(1)/4)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_10():
    f = (sqrt(a)*sqrt(b) - b*x**2)/(a + b*x**4)
    F = -sqrt(2)*b**(sympy.S(1)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(4*a**(sympy.S(1)/4)) + sqrt(2)*b**(sympy.S(1)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(4*a**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_11():
    f = (d + e*x**2)/(d**2 + e**2*x**4)
    F = -sqrt(2)*atan(1 - sqrt(2)*sqrt(e)*x/sqrt(d))/(2*sqrt(d)*sqrt(e)) + sqrt(2)*atan(1 + sqrt(2)*sqrt(e)*x/sqrt(d))/(2*sqrt(d)*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_12():
    f = (d - e*x**2)/(d**2 + e**2*x**4)
    F = -sqrt(2)*log(-sqrt(2)*sqrt(d)*sqrt(e)*x + d + e*x**2)/(4*sqrt(d)*sqrt(e)) + sqrt(2)*log(sqrt(2)*sqrt(d)*sqrt(e)*x + d + e*x**2)/(4*sqrt(d)*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_13():
    f = (2*x**2 + 5)/(x**4 - 1)
    F = -3*atan(x)/2 - 7*atanh(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_14():
    f = (b*x**2 + 1)/sqrt(-b**2*x**4 + 1)
    F = elliptic_e(asin(sqrt(b)*x), -1)/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_15():
    f = (-b*x**2 + 1)/sqrt(-b**2*x**4 + 1)
    F = -elliptic_e(asin(sqrt(b)*x), -1)/sqrt(b) + 2*elliptic_f(asin(sqrt(b)*x), -1)/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_16():
    f = (b*x**2 + 1)/sqrt(b**2*x**4 - 1)
    F = sqrt(-b**2*x**4 + 1)*elliptic_e(asin(sqrt(b)*x), -1)/(sqrt(b)*sqrt(b**2*x**4 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_17():
    f = (-b*x**2 + 1)/sqrt(b**2*x**4 - 1)
    F = -sqrt(-b**2*x**4 + 1)*elliptic_e(asin(sqrt(b)*x), -1)/(sqrt(b)*sqrt(b**2*x**4 - 1)) + 2*sqrt(-b**2*x**4 + 1)*elliptic_f(asin(sqrt(b)*x), -1)/(sqrt(b)*sqrt(b**2*x**4 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_18():
    f = (-b*x**2 + 1)/sqrt(b**2*x**4 + 1)
    F = -x*sqrt(b**2*x**4 + 1)/(b*x**2 + 1) + sqrt((b**2*x**4 + 1)/(b*x**2 + 1)**2)*(b*x**2 + 1)*elliptic_e(2*atan(sqrt(b)*x), sympy.S.Half)/(sqrt(b)*sqrt(b**2*x**4 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_19():
    f = (b*x**2 + 1)/sqrt(b**2*x**4 + 1)
    F = x*sqrt(b**2*x**4 + 1)/(b*x**2 + 1) - sqrt((b**2*x**4 + 1)/(b*x**2 + 1)**2)*(b*x**2 + 1)*elliptic_e(2*atan(sqrt(b)*x), sympy.S.Half)/(sqrt(b)*sqrt(b**2*x**4 + 1)) + sqrt((b**2*x**4 + 1)/(b*x**2 + 1)**2)*(b*x**2 + 1)*elliptic_f(2*atan(sqrt(b)*x), sympy.S.Half)/(sqrt(b)*sqrt(b**2*x**4 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_20():
    f = (-b*x**2 + 1)/sqrt(-b**2*x**4 - 1)
    F = x*sqrt(-b**2*x**4 - 1)/(b*x**2 + 1) + sqrt((b**2*x**4 + 1)/(b*x**2 + 1)**2)*(b*x**2 + 1)*elliptic_e(2*atan(sqrt(b)*x), sympy.S.Half)/(sqrt(b)*sqrt(-b**2*x**4 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_21():
    f = (b*x**2 + 1)/sqrt(-b**2*x**4 - 1)
    F = -x*sqrt(-b**2*x**4 - 1)/(b*x**2 + 1) - sqrt((b**2*x**4 + 1)/(b*x**2 + 1)**2)*(b*x**2 + 1)*elliptic_e(2*atan(sqrt(b)*x), sympy.S.Half)/(sqrt(b)*sqrt(-b**2*x**4 - 1)) + sqrt((b**2*x**4 + 1)/(b*x**2 + 1)**2)*(b*x**2 + 1)*elliptic_f(2*atan(sqrt(b)*x), sympy.S.Half)/(sqrt(b)*sqrt(-b**2*x**4 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_22():
    f = sqrt(c**2*x**2 + 1)/sqrt(-c**2*x**2 + 1)
    F = elliptic_e(asin(c*x), -1)/c
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_23():
    f = (c**2*x**2 + 1)/sqrt(-c**4*x**4 + 1)
    F = elliptic_e(asin(c*x), -1)/c
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_24():
    f = sqrt(-c**2*x**2 + 1)/sqrt(c**2*x**2 + 1)
    F = -elliptic_e(asin(c*x), -1)/c + 2*elliptic_f(asin(c*x), -1)/c
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_25():
    f = (-c**2*x**2 + 1)/sqrt(-c**4*x**4 + 1)
    F = -elliptic_e(asin(c*x), -1)/c + 2*elliptic_f(asin(c*x), -1)/c
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_26():
    f = (d + e*x**2)/(b*x**2 + d**2 + e**2*x**4)
    F = -atan((-2*e*x + sqrt(-b + 2*d*e))/sqrt(b + 2*d*e))/sqrt(b + 2*d*e) + atan((2*e*x + sqrt(-b + 2*d*e))/sqrt(b + 2*d*e))/sqrt(b + 2*d*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_27():
    f = (d + e*x**2)/(d**2 + e**2*x**4 + f*x**2)
    F = -atan((-2*e*x + sqrt(2*d*e - f))/sqrt(2*d*e + f))/sqrt(2*d*e + f) + atan((2*e*x + sqrt(2*d*e - f))/sqrt(2*d*e + f))/sqrt(2*d*e + f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_28():
    f = (d + e*x**2)/(-b*x**2 + d**2 + e**2*x**4)
    F = atanh((-2*e*x + sqrt(b + 2*d*e))/sqrt(b - 2*d*e))/sqrt(b - 2*d*e) - atanh((2*e*x + sqrt(b + 2*d*e))/sqrt(b - 2*d*e))/sqrt(b - 2*d*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_29():
    f = (d + e*x**2)/(d**2 + e**2*x**4 - f*x**2)
    F = -atan((-2*e*x + sqrt(2*d*e + f))/sqrt(2*d*e - f))/sqrt(2*d*e - f) + atan((2*e*x + sqrt(2*d*e + f))/sqrt(2*d*e - f))/sqrt(2*d*e - f)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_30():
    f = (d - e*x**2)/(b*x**2 + d**2 + e**2*x**4)
    F = -log(d + e*x**2 - x*sqrt(-b + 2*d*e))/(2*sqrt(-b + 2*d*e)) + log(d + e*x**2 + x*sqrt(-b + 2*d*e))/(2*sqrt(-b + 2*d*e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_31():
    f = (d - e*x**2)/(d**2 + e**2*x**4 + f*x**2)
    F = -log(d + e*x**2 - x*sqrt(2*d*e - f))/(2*sqrt(2*d*e - f)) + log(d + e*x**2 + x*sqrt(2*d*e - f))/(2*sqrt(2*d*e - f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_32():
    f = (d - e*x**2)/(-b*x**2 + d**2 + e**2*x**4)
    F = -log(d + e*x**2 - x*sqrt(b + 2*d*e))/(2*sqrt(b + 2*d*e)) + log(d + e*x**2 + x*sqrt(b + 2*d*e))/(2*sqrt(b + 2*d*e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_33():
    f = (d - e*x**2)/(d**2 + e**2*x**4 - f*x**2)
    F = -log(d + e*x**2 - x*sqrt(2*d*e + f))/(2*sqrt(2*d*e + f)) + log(d + e*x**2 + x*sqrt(2*d*e + f))/(2*sqrt(2*d*e + f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_34():
    f = (d - e*x**2)/(b*x**2 + c*d**2/e**2 + c*x**4)
    F = -e**(sympy.S(3)/2)*log(sqrt(c)*d + sqrt(c)*e*x**2 - sqrt(e)*x*sqrt(-b*e + 2*c*d))/(2*sqrt(c)*sqrt(-b*e + 2*c*d)) + e**(sympy.S(3)/2)*log(sqrt(c)*d + sqrt(c)*e*x**2 + sqrt(e)*x*sqrt(-b*e + 2*c*d))/(2*sqrt(c)*sqrt(-b*e + 2*c*d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_35():
    f = (d + e*x**2)/(b*x**2 + c*d**2/e**2 + c*x**4)
    F = -e**(sympy.S(3)/2)*atan((-2*sqrt(c)*sqrt(e)*x + sqrt(-b*e + 2*c*d))/sqrt(b*e + 2*c*d))/(sqrt(c)*sqrt(b*e + 2*c*d)) + e**(sympy.S(3)/2)*atan((2*sqrt(c)*sqrt(e)*x + sqrt(-b*e + 2*c*d))/sqrt(b*e + 2*c*d))/(sqrt(c)*sqrt(b*e + 2*c*d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_36():
    f = (d + e*x**2)/(b*x**2 + c*(d**2/e**2 + x**4))
    F = -e**(sympy.S(3)/2)*atan((-2*sqrt(c)*sqrt(e)*x + sqrt(-b*e + 2*c*d))/sqrt(b*e + 2*c*d))/(sqrt(c)*sqrt(b*e + 2*c*d)) + e**(sympy.S(3)/2)*atan((2*sqrt(c)*sqrt(e)*x + sqrt(-b*e + 2*c*d))/sqrt(b*e + 2*c*d))/(sqrt(c)*sqrt(b*e + 2*c*d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_37():
    f = (a - b*x**2)/(a**2 + b**2*x**4 + x**2*(2*a*b - 1))
    F = -log(a + b*x**2 - x)/2 + log(a + b*x**2 + x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_38():
    f = (a + b*x**2)/(a**2 + b**2*x**4 + x**2*(2*a*b - 1))
    F = atanh((-2*b*x + 1)/sqrt(-4*a*b + 1))/sqrt(-4*a*b + 1) - atanh((2*b*x + 1)/sqrt(-4*a*b + 1))/sqrt(-4*a*b + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_39():
    f = (2*x**2 + 1)/(b*x**2 + 4*x**4 + 1)
    F = -atan((-4*x + sqrt(4 - b))/sqrt(b + 4))/sqrt(b + 4) + atan((4*x + sqrt(4 - b))/sqrt(b + 4))/sqrt(b + 4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_40():
    f = (2*x**2 + 1)/(-b*x**2 + 4*x**4 + 1)
    F = -atan((-4*x + sqrt(b + 4))/sqrt(4 - b))/sqrt(4 - b) + atan((4*x + sqrt(b + 4))/sqrt(4 - b))/sqrt(4 - b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_41():
    f = (2*x**2 + 1)/(4*x**4 + 6*x**2 + 1)
    F = sqrt(10)*atan(2*x/sqrt(3 - sqrt(5)))/10 + sqrt(10)*atan(2*x/sqrt(sqrt(5) + 3))/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_42():
    f = (2*x**2 + 1)/(4*x**4 + 5*x**2 + 1)
    F = atan(x)/3 + atan(2*x)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_43():
    f = (2*x**2 + 1)/(4*x**4 + 4*x**2 + 1)
    F = sqrt(2)*atan(sqrt(2)*x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_44():
    f = (2*x**2 + 1)/(4*x**4 + 3*x**2 + 1)
    F = -sqrt(7)*atan(sqrt(7)*(1 - 4*x)/7)/7 + sqrt(7)*atan(sqrt(7)*(4*x + 1)/7)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_45():
    f = (2*x**2 + 1)/(4*x**4 + 2*x**2 + 1)
    F = -sqrt(6)*atan(sqrt(3)*(-2*sqrt(2)*x + 1)/3)/6 + sqrt(6)*atan(sqrt(3)*(2*sqrt(2)*x + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_46():
    f = (2*x**2 + 1)/(4*x**4 + x**2 + 1)
    F = -sqrt(5)*atan(sqrt(5)*(-4*x + sqrt(3))/5)/5 + sqrt(5)*atan(sqrt(5)*(4*x + sqrt(3))/5)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_47():
    f = (2*x**2 + 1)/(4*x**4 + 1)
    F = atan(2*x - 1)/2 + atan(2*x + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_48():
    f = (2*x**2 + 1)/(4*x**4 - x**2 + 1)
    F = -sqrt(3)*atan(sqrt(3)*(-4*x + sqrt(5))/3)/3 + sqrt(3)*atan(sqrt(3)*(4*x + sqrt(5))/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_49():
    f = (2*x**2 + 1)/(4*x**4 - 2*x**2 + 1)
    F = sqrt(2)*atan(2*sqrt(2)*x - sqrt(3))/2 + sqrt(2)*atan(2*sqrt(2)*x + sqrt(3))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_50():
    f = (2*x**2 + 1)/(4*x**4 - 3*x**2 + 1)
    F = atan(4*x - sqrt(7)) + atan(4*x + sqrt(7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_51():
    f = (2*x**2 + 1)/(4*x**4 - 4*x**2 + 1)
    F = x/(1 - 2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_52():
    f = (2*x**2 + 1)/(4*x**4 - 5*x**2 + 1)
    F = -log(1 - 2*x)/2 + log(1 - x)/2 - log(x + 1)/2 + log(2*x + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_53():
    f = (2*x**2 + 1)/(4*x**4 - 6*x**2 + 1)
    F = -sqrt(2)*atanh(2*sqrt(2)*x - sqrt(5))/2 - sqrt(2)*atanh(2*sqrt(2)*x + sqrt(5))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_54():
    f = (1 - 2*x**2)/(b*x**2 + 4*x**4 + 1)
    F = -log(2*x**2 - x*sqrt(4 - b) + 1)/(2*sqrt(4 - b)) + log(2*x**2 + x*sqrt(4 - b) + 1)/(2*sqrt(4 - b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_55():
    f = (1 - 2*x**2)/(4*x**4 + 6*x**2 + 1)
    F = sqrt(2)*atan(2*x/sqrt(3 - sqrt(5)))/2 - sqrt(2)*atan(2*x/sqrt(sqrt(5) + 3))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_56():
    f = (1 - 2*x**2)/(4*x**4 + 5*x**2 + 1)
    F = -atan(x) + atan(2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_57():
    f = (1 - 2*x**2)/(4*x**4 + 4*x**2 + 1)
    F = x/(2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_58():
    f = (1 - 2*x**2)/(4*x**4 + 3*x**2 + 1)
    F = -log(2*x**2 - x + 1)/2 + log(2*x**2 + x + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_59():
    f = (1 - 2*x**2)/(4*x**4 + 2*x**2 + 1)
    F = -sqrt(2)*log(2*x**2 - sqrt(2)*x + 1)/4 + sqrt(2)*log(2*x**2 + sqrt(2)*x + 1)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_60():
    f = (1 - 2*x**2)/(4*x**4 + x**2 + 1)
    F = -sqrt(3)*log(2*x**2 - sqrt(3)*x + 1)/6 + sqrt(3)*log(2*x**2 + sqrt(3)*x + 1)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_61():
    f = (1 - 2*x**2)/(4*x**4 + 1)
    F = -log(2*x**2 - 2*x + 1)/4 + log(2*x**2 + 2*x + 1)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_62():
    f = (1 - 2*x**2)/(4*x**4 - x**2 + 1)
    F = -sqrt(5)*log(2*x**2 - sqrt(5)*x + 1)/10 + sqrt(5)*log(2*x**2 + sqrt(5)*x + 1)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_63():
    f = (1 - 2*x**2)/(4*x**4 - 2*x**2 + 1)
    F = -sqrt(6)*log(2*x**2 - sqrt(6)*x + 1)/12 + sqrt(6)*log(2*x**2 + sqrt(6)*x + 1)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_64():
    f = (1 - 2*x**2)/(4*x**4 - 3*x**2 + 1)
    F = -sqrt(7)*log(2*x**2 - sqrt(7)*x + 1)/14 + sqrt(7)*log(2*x**2 + sqrt(7)*x + 1)/14
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_65():
    f = (1 - 2*x**2)/(4*x**4 - 4*x**2 + 1)
    F = sqrt(2)*atanh(sqrt(2)*x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_66():
    f = (1 - 2*x**2)/(4*x**4 - 5*x**2 + 1)
    F = -log(1 - 2*x)/6 - log(1 - x)/6 + log(x + 1)/6 + log(2*x + 1)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_67():
    f = (1 - 2*x**2)/(4*x**4 - 6*x**2 + 1)
    F = -sqrt(10)*atanh(sqrt(5)*(-2*sqrt(2)*x + 1)/5)/10 + sqrt(10)*atanh(sqrt(5)*(2*sqrt(2)*x + 1)/5)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_68():
    f = (x**2 + 1)/(b*x**2 + x**4 + 1)
    F = -atan((-2*x + sqrt(2 - b))/sqrt(b + 2))/sqrt(b + 2) + atan((2*x + sqrt(2 - b))/sqrt(b + 2))/sqrt(b + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_69():
    f = (x**2 + 1)/(x**4 + 5*x**2 + 1)
    F = sqrt(7)*atan(x*sqrt(sqrt(21)/2 + sympy.S(5)/2))/7 + sqrt(7)*atan(sqrt(2)*x/sqrt(sqrt(21) + 5))/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_70():
    f = (x**2 + 1)/(x**4 + 4*x**2 + 1)
    F = sqrt(6)*atan(x/sqrt(2 - sqrt(3)))/6 + sqrt(6)*atan(x/sqrt(sqrt(3) + 2))/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_71():
    f = (x**2 + 1)/(x**4 + 3*x**2 + 1)
    F = sqrt(5)*atan(x*sqrt(sqrt(5)/2 + sympy.S(3)/2))/5 + sqrt(5)*atan(sqrt(2)*x/sqrt(sqrt(5) + 3))/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_72():
    f = (x**2 + 1)/(x**4 + 2*x**2 + 1)
    F = atan(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_73():
    f = (x**2 + 1)/(x**4 + x**2 + 1)
    F = -sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/3 + sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_74():
    f = (x**2 + 1)/(x**4 + 1)
    F = sqrt(2)*atan(sqrt(2)*x - 1)/2 + sqrt(2)*atan(sqrt(2)*x + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_75():
    f = (x**2 + 1)/(x**4 - x**2 + 1)
    F = atan(2*x - sqrt(3)) + atan(2*x + sqrt(3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_76():
    f = (x**2 + 1)/(x**4 - 2*x**2 + 1)
    F = x/(1 - x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_77():
    f = (x**2 + 1)/(x**4 - 3*x**2 + 1)
    F = log(-2*x + 1 + sqrt(5))/2 + log(-2*x - sqrt(5) + 1)/2 - log(2*x + 1 + sqrt(5))/2 - log(2*x - sqrt(5) + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_78():
    f = (x**2 + 1)/(x**4 - 4*x**2 + 1)
    F = -sqrt(2)*atanh(sqrt(2)*x - sqrt(3))/2 - sqrt(2)*atanh(sqrt(2)*x + sqrt(3))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_79():
    f = (x**2 + 1)/(x**4 - 5*x**2 + 1)
    F = sqrt(3)*atanh(sqrt(3)*(-2*x + sqrt(7))/3)/3 - sqrt(3)*atanh(sqrt(3)*(2*x + sqrt(7))/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_80():
    f = (1 - x**2)/(b*x**2 + x**4 + 1)
    F = -log(x**2 - x*sqrt(2 - b) + 1)/(2*sqrt(2 - b)) + log(x**2 + x*sqrt(2 - b) + 1)/(2*sqrt(2 - b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_81():
    f = (1 - x**2)/(x**4 + 5*x**2 + 1)
    F = sqrt(3)*atan(x*sqrt(sqrt(21)/2 + sympy.S(5)/2))/3 - sqrt(3)*atan(sqrt(2)*x/sqrt(sqrt(21) + 5))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_82():
    f = (1 - x**2)/(x**4 + 4*x**2 + 1)
    F = sqrt(2)*atan(x/sqrt(2 - sqrt(3)))/2 - sqrt(2)*atan(x/sqrt(sqrt(3) + 2))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_83():
    f = (1 - x**2)/(x**4 + 3*x**2 + 1)
    F = atan(x*sqrt(sqrt(5)/2 + sympy.S(3)/2)) - atan(sqrt(2)*x/sqrt(sqrt(5) + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_84():
    f = (1 - x**2)/(x**4 + 2*x**2 + 1)
    F = x/(x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_85():
    f = (1 - x**2)/(x**4 + x**2 + 1)
    F = -log(x**2 - x + 1)/2 + log(x**2 + x + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_86():
    f = (1 - x**2)/(x**4 + 1)
    F = -sqrt(2)*log(x**2 - sqrt(2)*x + 1)/4 + sqrt(2)*log(x**2 + sqrt(2)*x + 1)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_87():
    f = (1 - x**2)/(x**4 - x**2 + 1)
    F = -sqrt(3)*log(x**2 - sqrt(3)*x + 1)/6 + sqrt(3)*log(x**2 + sqrt(3)*x + 1)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_88():
    f = (1 - x**2)/(x**4 - 2*x**2 + 1)
    F = atanh(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_89():
    f = (1 - x**2)/(x**4 - 3*x**2 + 1)
    F = -sqrt(5)*atanh(sqrt(5)*(1 - 2*x)/5)/5 + sqrt(5)*atanh(sqrt(5)*(2*x + 1)/5)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_90():
    f = (1 - x**2)/(x**4 - 4*x**2 + 1)
    F = -sqrt(6)*atanh(sqrt(3)*(-sqrt(2)*x + 1)/3)/6 + sqrt(6)*atanh(sqrt(3)*(sqrt(2)*x + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_91():
    f = (1 - x**2)/(x**4 - 5*x**2 + 1)
    F = -sqrt(7)*atanh(sqrt(7)*(-2*x + sqrt(3))/7)/7 + sqrt(7)*atanh(sqrt(7)*(2*x + sqrt(3))/7)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_92():
    f = (-3*x**2 - 1)/(9*x**4 + 2*x**2 + 1)
    F = sqrt(2)*atan(sqrt(2)*(1 - 3*x)/2)/4 - sqrt(2)*atan(sqrt(2)*(3*x + 1)/2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_93():
    f = (3*x**2 + 1)/(-9*x**4 - 2*x**2 - 1)
    F = sqrt(2)*atan(sqrt(2)*(1 - 3*x)/2)/4 - sqrt(2)*atan(sqrt(2)*(3*x + 1)/2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_94():
    f = (2*x**2 + 3)/(x**4 - 2*x**2 + 1)
    F = 5*x/(2 - 2*x**2) + atanh(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_95():
    f = (3*x**2 + 2)/(3*x**4 - 8*x**2 + 5)
    F = 5*atanh(x)/2 - 7*sqrt(15)*atanh(sqrt(15)*x/5)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_96():
    f = (d + e*x**2)/(3*x**4 - 8*x**2 + 5)
    F = (d/2 + e/2)*atanh(x) - sqrt(15)*(3*d + 5*e)*atanh(sqrt(15)*x/5)/30
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_97():
    f = (x**2 + 3)/(x**4 + 3*x**2 + 1)
    F = sqrt(10)*(sqrt(5) + 3)**(sympy.S(3)/2)*atan(x*sqrt(sqrt(5)/2 + sympy.S(3)/2))/20 - sqrt(180 - 80*sqrt(5))*atan(sqrt(2)*x/sqrt(sqrt(5) + 3))/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_98():
    f = (a + b*x**2)/(x**4 + x**2 + 1)
    F = -(a/4 - b/4)*log(x**2 - x + 1) + (a/4 - b/4)*log(x**2 + x + 1) - sqrt(3)*(a + b)*atan(sqrt(3)*(1 - 2*x)/3)/6 + sqrt(3)*(a + b)*atan(sqrt(3)*(2*x + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_99():
    f = (a + b*x**2)/(x**4 + x**2 + 1)**2
    F = x*(a + b - x**2*(a - 2*b))/(6*x**4 + 6*x**2 + 6) - (a/4 - b/8)*log(x**2 - x + 1) + (a/4 - b/8)*log(x**2 + x + 1) - sqrt(3)*(4*a + b)*atan(sqrt(3)*(1 - 2*x)/3)/36 + sqrt(3)*(4*a + b)*atan(sqrt(3)*(2*x + 1)/3)/36
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_100():
    f = (a + b*x**2)/(x**4 + x**2 + 2)
    F = -(a - sqrt(2)*b)*log(x**2 - x*sqrt(-1 + 2*sqrt(2)) + sqrt(2))/(4*sqrt(-2 + 4*sqrt(2))) + (a - sqrt(2)*b)*log(x**2 + x*sqrt(-1 + 2*sqrt(2)) + sqrt(2))/(4*sqrt(-2 + 4*sqrt(2))) - sqrt(sympy.S(-1)/14 + sqrt(2)/7)*(a + sqrt(2)*b)*atan((-2*x + sqrt(-1 + 2*sqrt(2)))/sqrt(1 + 2*sqrt(2)))/2 + sqrt(sympy.S(-1)/14 + sqrt(2)/7)*(a + sqrt(2)*b)*atan((2*x + sqrt(-1 + 2*sqrt(2)))/sqrt(1 + 2*sqrt(2)))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_101():
    f = (a + b*x**2)/(x**4 + x**2 + 2)**2
    F = x*(3*a + 2*b - x**2*(a - 4*b))/(28*x**4 + 28*x**2 + 56) - sqrt(sympy.S(-1)/14 + sqrt(2)/7)*(a*(11 - sqrt(2)) - b*(2 - 4*sqrt(2)))*atan((-2*x + sqrt(-1 + 2*sqrt(2)))/sqrt(1 + 2*sqrt(2)))/56 + sqrt(sympy.S(-1)/14 + sqrt(2)/7)*(a*(11 - sqrt(2)) - b*(2 - 4*sqrt(2)))*atan((2*x + sqrt(-1 + 2*sqrt(2)))/sqrt(1 + 2*sqrt(2)))/56 - (11*a - 2*b + sqrt(2)*(a - 4*b))*log(x**2 - x*sqrt(-1 + 2*sqrt(2)) + sqrt(2))/(112*sqrt(-2 + 4*sqrt(2))) + (a*(sqrt(2) + 11) - 4*sqrt(2)*b - 2*b)*log(x**2 + x*sqrt(-1 + 2*sqrt(2)) + sqrt(2))/(112*sqrt(-2 + 4*sqrt(2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_102():
    f = (-x**2 + sqrt(2))/(x**4 - sqrt(2)*x**2 + 1)
    F = -sqrt(sqrt(2)/2 + 1)*log(x**2 - x*sqrt(sqrt(2) + 2) + 1)/4 + sqrt(sqrt(2)/2 + 1)*log(x**2 + x*sqrt(sqrt(2) + 2) + 1)/4 - atan((-2*x + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(2*sqrt(sqrt(2) + 2)) + atan((2*x + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(2*sqrt(sqrt(2) + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_103():
    f = (x**2 + sqrt(2))/(x**4 + sqrt(2)*x**2 + 1)
    F = -sqrt(1 - sqrt(2)/2)*log(x**2 - x*sqrt(2 - sqrt(2)) + 1)/4 + sqrt(1 - sqrt(2)/2)*log(x**2 + x*sqrt(2 - sqrt(2)) + 1)/4 - atan((-2*x + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(2*sqrt(2 - sqrt(2))) + atan((2*x + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(2*sqrt(2 - sqrt(2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_104():
    f = (-x**2 + sqrt(2))/(b*x**2 + x**4 + 1)
    F = (1 - sqrt(2))*atan((-2*x + sqrt(2 - b))/sqrt(b + 2))/(2*sqrt(b + 2)) - (1 - sqrt(2))*atan((2*x + sqrt(2 - b))/sqrt(b + 2))/(2*sqrt(b + 2)) - (1 + sqrt(2))*log(x**2 - x*sqrt(2 - b) + 1)/(4*sqrt(2 - b)) + (1 + sqrt(2))*log(x**2 + x*sqrt(2 - b) + 1)/(4*sqrt(2 - b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_105():
    f = (x**2 + sqrt(2))/(b*x**2 + x**4 + 1)
    F = -(1 + sqrt(2))*atan((-2*x + sqrt(2 - b))/sqrt(b + 2))/(2*sqrt(b + 2)) + (1 + sqrt(2))*atan((2*x + sqrt(2 - b))/sqrt(b + 2))/(2*sqrt(b + 2)) + (1 - sqrt(2))*log(x**2 - x*sqrt(2 - b) + 1)/(4*sqrt(2 - b)) - (1 - sqrt(2))*log(x**2 + x*sqrt(2 - b) + 1)/(4*sqrt(2 - b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_106():
    f = (2*a - x**2)/(a**2 - a*x**2 + x**4)
    F = -sqrt(3)*log(-sqrt(3)*sqrt(a)*x + a + x**2)/(4*sqrt(a)) + sqrt(3)*log(sqrt(3)*sqrt(a)*x + a + x**2)/(4*sqrt(a)) - atan(sqrt(3) - 2*x/sqrt(a))/(2*sqrt(a)) + atan(sqrt(3) + 2*x/sqrt(a))/(2*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_107():
    f = (2*sqrt(a) - x**2)/(-sqrt(a)*x**2 + a + x**4)
    F = -sqrt(3)*log(-sqrt(3)*a**(sympy.S(1)/4)*x + sqrt(a) + x**2)/(4*a**(sympy.S(1)/4)) + sqrt(3)*log(sqrt(3)*a**(sympy.S(1)/4)*x + sqrt(a) + x**2)/(4*a**(sympy.S(1)/4)) - atan(sqrt(3) - 2*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4)) + atan(sqrt(3) + 2*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_108():
    f = (2*b**(sympy.S(2)/3) + x**2)/(b**(sympy.S(4)/3) + b**(sympy.S(2)/3)*x**2 + x**4)
    F = -log(b**(sympy.S(2)/3) - b**(sympy.S(1)/3)*x + x**2)/(4*b**(sympy.S(1)/3)) + log(b**(sympy.S(2)/3) + b**(sympy.S(1)/3)*x + x**2)/(4*b**(sympy.S(1)/3)) - sqrt(3)*atan(sqrt(3)*(b**(sympy.S(1)/3) - 2*x)/(3*b**(sympy.S(1)/3)))/(2*b**(sympy.S(1)/3)) + sqrt(3)*atan(sqrt(3)*(b**(sympy.S(1)/3) + 2*x)/(3*b**(sympy.S(1)/3)))/(2*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_109():
    f = (A + B*x**2)/(a**2 - a*x**2 + x**4)
    F = -sqrt(3)*(A - B*a)*log(-sqrt(3)*sqrt(a)*x + a + x**2)/(12*a**(sympy.S(3)/2)) + sqrt(3)*(A - B*a)*log(sqrt(3)*sqrt(a)*x + a + x**2)/(12*a**(sympy.S(3)/2)) - (A + B*a)*atan(sqrt(3) - 2*x/sqrt(a))/(2*a**(sympy.S(3)/2)) + (A + B*a)*atan(sqrt(3) + 2*x/sqrt(a))/(2*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_110():
    f = (A + B*x**2)/(-sqrt(a)*x**2 + a + x**4)
    F = -sqrt(3)*(A - B*sqrt(a))*log(-sqrt(3)*a**(sympy.S(1)/4)*x + sqrt(a) + x**2)/(12*a**(sympy.S(3)/4)) + sqrt(3)*(A - B*sqrt(a))*log(sqrt(3)*a**(sympy.S(1)/4)*x + sqrt(a) + x**2)/(12*a**(sympy.S(3)/4)) - (A + B*sqrt(a))*atan(sqrt(3) - 2*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)) + (A + B*sqrt(a))*atan(sqrt(3) + 2*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_111():
    f = (A + B*x**2)/(a + c*x**4 - x**2*sqrt(a*c))
    F = -(A - B*sqrt(a)/sqrt(c))*log(sqrt(a) + sqrt(c)*x**2 - x*sqrt(2*sqrt(a)*sqrt(c) + sqrt(a*c)))/(4*sqrt(a)*sqrt(2*sqrt(a)*sqrt(c) + sqrt(a*c))) + (A - B*sqrt(a)/sqrt(c))*log(sqrt(a) + sqrt(c)*x**2 + x*sqrt(2*sqrt(a)*sqrt(c) + sqrt(a*c)))/(4*sqrt(a)*sqrt(2*sqrt(a)*sqrt(c) + sqrt(a*c))) - (A*sqrt(c) + B*sqrt(a))*atan((-2*sqrt(c)*x + sqrt(2*sqrt(a)*sqrt(c) + sqrt(a*c)))/sqrt(2*sqrt(a)*sqrt(c) - sqrt(a*c)))/(2*sqrt(a)*sqrt(c)*sqrt(2*sqrt(a)*sqrt(c) - sqrt(a*c))) + (A*sqrt(c) + B*sqrt(a))*atan((2*sqrt(c)*x + sqrt(2*sqrt(a)*sqrt(c) + sqrt(a*c)))/sqrt(2*sqrt(a)*sqrt(c) - sqrt(a*c)))/(2*sqrt(a)*sqrt(c)*sqrt(2*sqrt(a)*sqrt(c) - sqrt(a*c)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_112():
    f = (A + B*x**2)/(-sqrt(a)*sqrt(c)*x**2 + a + c*x**4)
    F = -sqrt(3)*(A - B*sqrt(a)/sqrt(c))*log(-sqrt(3)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(12*a**(sympy.S(3)/4)*c**(sympy.S(1)/4)) + sqrt(3)*(A - B*sqrt(a)/sqrt(c))*log(sqrt(3)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(12*a**(sympy.S(3)/4)*c**(sympy.S(1)/4)) - (A*sqrt(c) + B*sqrt(a))*atan(sqrt(3) - 2*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*c**(sympy.S(3)/4)) + (A*sqrt(c) + B*sqrt(a))*atan(sqrt(3) + 2*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(2*a**(sympy.S(3)/4)*c**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_113():
    f = (3 - x**2)/sqrt(-x**4 + x**2 + 3)
    F = -sqrt(sympy.S(-1)/2 + sqrt(13)/2)*elliptic_e(asin(sqrt(2)*x/sqrt(1 + sqrt(13))), sympy.S(-7)/6 - sqrt(13)/6) + sqrt(7 + 2*sqrt(13))*elliptic_f(asin(sqrt(2)*x/sqrt(1 + sqrt(13))), sympy.S(-7)/6 - sqrt(13)/6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_114():
    f = (3 - x**2)/sqrt(-x**4 + 2*x**2 + 3)
    F = -elliptic_e(asin(sqrt(3)*x/3), -3) + 4*elliptic_f(asin(sqrt(3)*x/3), -3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_115():
    f = (3 - x**2)/sqrt(-x**4 + 3*x**2 + 3)
    F = -sqrt(sympy.S(-3)/2 + sqrt(21)/2)*elliptic_e(asin(sqrt(2)*x/sqrt(3 + sqrt(21))), sympy.S(-5)/2 - sqrt(21)/2) + sqrt(9 + 2*sqrt(21))*elliptic_f(asin(sqrt(2)*x/sqrt(3 + sqrt(21))), sympy.S(-5)/2 - sqrt(21)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_116():
    f = (3 - x**2)/sqrt(-x**4 - x**2 + 3)
    F = -sqrt(sympy.S.Half + sqrt(13)/2)*elliptic_e(asin(sqrt(2)*x/sqrt(-1 + sqrt(13))), sympy.S(-7)/6 + sqrt(13)/6) + sqrt(5 + 2*sqrt(13))*elliptic_f(asin(sqrt(2)*x/sqrt(-1 + sqrt(13))), sympy.S(-7)/6 + sqrt(13)/6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_117():
    f = (3 - x**2)/sqrt(-x**4 - 2*x**2 + 3)
    F = -sqrt(3)*elliptic_e(asin(x), sympy.S(-1)/3) + 2*sqrt(3)*elliptic_f(asin(x), sympy.S(-1)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_118():
    f = (3 - x**2)/sqrt(-x**4 - 3*x**2 + 3)
    F = -sqrt(sympy.S(3)/2 + sqrt(21)/2)*elliptic_e(asin(sqrt(2)*x/sqrt(-3 + sqrt(21))), sympy.S(-5)/2 + sqrt(21)/2) + sqrt(3 + 2*sqrt(21))*elliptic_f(asin(sqrt(2)*x/sqrt(-3 + sqrt(21))), sympy.S(-5)/2 + sqrt(21)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_119():
    f = (b + 2*c*x**2 - sqrt(-4*a*c + b**2))/sqrt(a + b*x**2 + c*x**4)
    F = -2*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/sqrt(a + b*x**2 + c*x**4) + 2*sqrt(c)*x*sqrt(a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2) + sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(2*sqrt(a)*sqrt(c) + b - sqrt(-4*a*c + b**2))*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(2*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_120():
    f = (a + c*x**4)*(d + e*x**2)**4
    F = a*d**4*x + 4*a*d**3*e*x**3/3 + 4*c*d*e**3*x**11/11 + c*e**4*x**13/13 + d**2*x**5*(6*a*e**2 + c*d**2)/5 + 4*d*e*x**7*(a*e**2 + c*d**2)/7 + e**2*x**9*(a*e**2 + 6*c*d**2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_121():
    f = (a + c*x**4)*(d + e*x**2)**3
    F = a*d**3*x + a*d**2*e*x**3 + c*d*e**2*x**9/3 + c*e**3*x**11/11 + d*x**5*(3*a*e**2 + c*d**2)/5 + e*x**7*(a*e**2 + 3*c*d**2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_122():
    f = (a + c*x**4)*(d + e*x**2)**2
    F = a*d**2*x + 2*a*d*e*x**3/3 + 2*c*d*e*x**7/7 + c*e**2*x**9/9 + x**5*(a*e**2 + c*d**2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_123():
    f = (a + c*x**4)*(d + e*x**2)
    F = a*d*x + a*e*x**3/3 + c*d*x**5/5 + c*e*x**7/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_124():
    f = (a + c*x**4)/(d + e*x**2)
    F = -c*d*x/e**2 + c*x**3/(3*e) + (a*e**2 + c*d**2)*atan(sqrt(e)*x/sqrt(d))/(sqrt(d)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_125():
    f = (a + c*x**4)/(d + e*x**2)**2
    F = c*x/e**2 + x*(a + c*d**2/e**2)/(2*d*(d + e*x**2)) - (-a*e**2 + 3*c*d**2)*atan(sqrt(e)*x/sqrt(d))/(2*d**(sympy.S(3)/2)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_126():
    f = (a + c*x**4)/(d + e*x**2)**3
    F = x*(3*a/d**2 - 5*c/e**2)/(8*d + 8*e*x**2) + x*(a + c*d**2/e**2)/(4*d*(d + e*x**2)**2) + (3*a*e**2 + 3*c*d**2)*atan(sqrt(e)*x/sqrt(d))/(8*d**(sympy.S(5)/2)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_127():
    f = (a + c*x**4)/(d + e*x**2)**4
    F = x*(5*a/d**2 - 7*c/e**2)/(24*(d + e*x**2)**2) + x*(a + c*d**2/e**2)/(6*d*(d + e*x**2)**3) + x*(5*a/d**2 + c/e**2)/(16*d*(d + e*x**2)) + (5*a*e**2 + c*d**2)*atan(sqrt(e)*x/sqrt(d))/(16*d**(sympy.S(7)/2)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_128():
    f = (a + c*x**4)**2*(d + e*x**2)**3
    F = a**2*d**3*x + a**2*d**2*e*x**3 + a*d*x**5*(3*a*e**2 + 2*c*d**2)/5 + a*e*x**7*(a*e**2 + 6*c*d**2)/7 + 3*c**2*d*e**2*x**13/13 + c**2*e**3*x**15/15 + c*d*x**9*(6*a*e**2 + c*d**2)/9 + c*e*x**11*(2*a*e**2 + 3*c*d**2)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_129():
    f = (a + c*x**4)**2*(d + e*x**2)**2
    F = a**2*d**2*x + 2*a**2*d*e*x**3/3 + 4*a*c*d*e*x**7/7 + a*x**5*(a*e**2 + 2*c*d**2)/5 + 2*c**2*d*e*x**11/11 + c**2*e**2*x**13/13 + c*x**9*(2*a*e**2 + c*d**2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_130():
    f = (a + c*x**4)**2*(d + e*x**2)
    F = a**2*d*x + a**2*e*x**3/3 + 2*a*c*d*x**5/5 + 2*a*c*e*x**7/7 + c**2*d*x**9/9 + c**2*e*x**11/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_131():
    f = (a + c*x**4)**2
    F = a**2*x + 2*a*c*x**5/5 + c**2*x**9/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_132():
    f = (a + c*x**4)**2/(d + e*x**2)
    F = -c**2*d*x**5/(5*e**2) + c**2*x**7/(7*e) - c*d*x*(2*a*e**2 + c*d**2)/e**4 + c*x**3*(2*a*e**2 + c*d**2)/(3*e**3) + (a*e**2 + c*d**2)**2*atan(sqrt(e)*x/sqrt(d))/(sqrt(d)*e**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_133():
    f = (a + c*x**4)**2/(d + e*x**2)**2
    F = -2*c**2*d*x**3/(3*e**3) + c**2*x**5/(5*e**2) + c*x*(2*a*e**2 + 3*c*d**2)/e**4 + x*(a*e**2 + c*d**2)**2/(2*d*e**4*(d + e*x**2)) - (-a*e**2 + 7*c*d**2)*(a*e**2 + c*d**2)*atan(sqrt(e)*x/sqrt(d))/(2*d**(sympy.S(3)/2)*e**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_134():
    f = (a + c*x**4)**2/(d + e*x**2)**3
    F = -3*c**2*d*x/e**4 + c**2*x**3/(3*e**3) + x*(a*e**2 + c*d**2)**2/(4*d*e**4*(d + e*x**2)**2) + x*(3*a**2 - 10*a*c*d**2/e**2 - 13*c**2*d**4/e**4)/(8*d**2*(d + e*x**2)) + (3*a**2*e**4 + 6*a*c*d**2*e**2 + 35*c**2*d**4)*atan(sqrt(e)*x/sqrt(d))/(8*d**(sympy.S(5)/2)*e**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_135():
    f = (a + c*x**4)**2/(d + e*x**2)**4
    F = c**2*x/e**4 + x*(a*e**2 + c*d**2)**2/(6*d*e**4*(d + e*x**2)**3) + x*(5*a**2 - 14*a*c*d**2/e**2 - 19*c**2*d**4/e**4)/(24*d**2*(d + e*x**2)**2) + x*(5*a**2 + 2*a*c*d**2/e**2 + 29*c**2*d**4/e**4)/(16*d**3*(d + e*x**2)) - (-5*a**2*e**4 - 2*a*c*d**2*e**2 + 35*c**2*d**4)*atan(sqrt(e)*x/sqrt(d))/(16*d**(sympy.S(7)/2)*e**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_136():
    f = (a + c*x**4)**2/(d + e*x**2)**5
    F = x*(a*e**2 + c*d**2)**2/(8*d*e**4*(d + e*x**2)**4) + x*(7*a**2 - 18*a*c*d**2/e**2 - 25*c**2*d**4/e**4)/(48*d**2*(d + e*x**2)**3) + x*(35*a**2 + 6*a*c*d**2/e**2 + 163*c**2*d**4/e**4)/(192*d**3*(d + e*x**2)**2) - x*(-35*a**2*e**4 - 6*a*c*d**2*e**2 + 93*c**2*d**4)/(128*d**4*e**4*(d + e*x**2)) + (35*a**2*e**4 + 6*a*c*d**2*e**2 + 35*c**2*d**4)*atan(sqrt(e)*x/sqrt(d))/(128*d**(sympy.S(9)/2)*e**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_137():
    f = (d + e*x**2)**4/(a + c*x**4)
    F = 4*d*e**3*x**3/(3*c) + e**4*x**5/(5*c) + e**2*x*(-a*e**2 + 6*c*d**2)/c**2 - sqrt(2)*(-4*sqrt(a)*sqrt(c)*d*e*(-a*e**2 + c*d**2) + a**2*e**4 - 6*a*c*d**2*e**2 + c**2*d**4)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*c**(sympy.S(9)/4)) + sqrt(2)*(-4*sqrt(a)*sqrt(c)*d*e*(-a*e**2 + c*d**2) + a**2*e**4 - 6*a*c*d**2*e**2 + c**2*d**4)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*c**(sympy.S(9)/4)) - sqrt(2)*(4*sqrt(a)*sqrt(c)*d*e*(-a*e**2 + c*d**2) + a**2*e**4 - 6*a*c*d**2*e**2 + c**2*d**4)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*c**(sympy.S(9)/4)) + sqrt(2)*(4*sqrt(a)*sqrt(c)*d*e*(-a*e**2 + c*d**2) + a**2*e**4 - 6*a*c*d**2*e**2 + c**2*d**4)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*c**(sympy.S(9)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_138():
    f = (d + e*x**2)**3/(a + c*x**4)
    F = 3*d*e**2*x/c + e**3*x**3/(3*c) - sqrt(2)*(-sqrt(a)*e*(-a*e**2 + 3*c*d**2) + sqrt(c)*d*(-3*a*e**2 + c*d**2))*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*c**(sympy.S(7)/4)) + sqrt(2)*(-sqrt(a)*e*(-a*e**2 + 3*c*d**2) + sqrt(c)*d*(-3*a*e**2 + c*d**2))*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*c**(sympy.S(7)/4)) - sqrt(2)*(sqrt(a)*e*(-a*e**2 + 3*c*d**2) + sqrt(c)*d*(-3*a*e**2 + c*d**2))*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*c**(sympy.S(7)/4)) + sqrt(2)*(sqrt(a)*e*(-a*e**2 + 3*c*d**2) + sqrt(c)*d*(-3*a*e**2 + c*d**2))*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*c**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_139():
    f = (d + e*x**2)**2/(a + c*x**4)
    F = e**2*x/c - sqrt(2)*(-2*sqrt(a)*sqrt(c)*d*e - a*e**2 + c*d**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*c**(sympy.S(5)/4)) + sqrt(2)*(-2*sqrt(a)*sqrt(c)*d*e - a*e**2 + c*d**2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*c**(sympy.S(5)/4)) - sqrt(2)*(2*sqrt(a)*sqrt(c)*d*e - a*e**2 + c*d**2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*c**(sympy.S(5)/4)) + sqrt(2)*(2*sqrt(a)*sqrt(c)*d*e - a*e**2 + c*d**2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*c**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_140():
    f = (d + e*x**2)/(a + c*x**4)
    F = -sqrt(2)*(-sqrt(a)*e + sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*c**(sympy.S(3)/4)) + sqrt(2)*(-sqrt(a)*e + sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*c**(sympy.S(3)/4)) - sqrt(2)*(sqrt(a)*e + sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*c**(sympy.S(3)/4)) + sqrt(2)*(sqrt(a)*e + sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*c**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_141():
    f = 1/(a + c*x**4)
    F = -sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*c**(sympy.S(1)/4)) + sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*c**(sympy.S(1)/4)) - sqrt(2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*c**(sympy.S(1)/4)) + sqrt(2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*c**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_142():
    f = 1/((a + c*x**4)*(d + e*x**2))
    F = e**(sympy.S(3)/2)*atan(sqrt(e)*x/sqrt(d))/(sqrt(d)*(a*e**2 + c*d**2)) - sqrt(2)*c**(sympy.S(1)/4)*(-sqrt(a)*e + sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)) + sqrt(2)*c**(sympy.S(1)/4)*(-sqrt(a)*e + sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)) - sqrt(2)*c**(sympy.S(1)/4)*(sqrt(a)*e + sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)) + sqrt(2)*c**(sympy.S(1)/4)*(sqrt(a)*e + sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_143():
    f = 1/((a + c*x**4)*(d + e*x**2)**2)
    F = 2*c*sqrt(d)*e**(sympy.S(3)/2)*atan(sqrt(e)*x/sqrt(d))/(a*e**2 + c*d**2)**2 + e**2*x/(2*d*(d + e*x**2)*(a*e**2 + c*d**2)) + e**(sympy.S(3)/2)*atan(sqrt(e)*x/sqrt(d))/(2*d**(sympy.S(3)/2)*(a*e**2 + c*d**2)) - sqrt(2)*c**(sympy.S(3)/4)*(-2*sqrt(a)*sqrt(c)*d*e - a*e**2 + c*d**2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) + sqrt(2)*c**(sympy.S(3)/4)*(-2*sqrt(a)*sqrt(c)*d*e - a*e**2 + c*d**2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) - sqrt(2)*c**(sympy.S(3)/4)*(2*sqrt(a)*sqrt(c)*d*e - a*e**2 + c*d**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) + sqrt(2)*c**(sympy.S(3)/4)*(2*sqrt(a)*sqrt(c)*d*e - a*e**2 + c*d**2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_144():
    f = (d + e*x**2)**3/(a + c*x**4)**2
    F = -e**3*x**3/(c*(a + c*x**4)) + x*(d*(-3*a*e**2 + c*d**2) + 3*e*x**2*(a*e**2 + c*d**2))/(4*a*c*(a + c*x**4)) - sqrt(2)*(-3*sqrt(a)*e + 3*sqrt(c)*d)*(a*e**2 + c*d**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*c**(sympy.S(7)/4)) + sqrt(2)*(-3*sqrt(a)*e + 3*sqrt(c)*d)*(a*e**2 + c*d**2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*c**(sympy.S(7)/4)) - sqrt(2)*(3*sqrt(a)*e + 3*sqrt(c)*d)*(a*e**2 + c*d**2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*c**(sympy.S(7)/4)) + sqrt(2)*(3*sqrt(a)*e + 3*sqrt(c)*d)*(a*e**2 + c*d**2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*c**(sympy.S(7)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_145():
    f = (d + e*x**2)**2/(a + c*x**4)**2
    F = -e**2*x/(3*c*(a + c*x**4)) + x*(a*e**2 + 3*c*d**2 + 6*c*d*e*x**2)/(12*a*c*(a + c*x**4)) - sqrt(2)*(-2*sqrt(a)*sqrt(c)*d*e + a*e**2 + 3*c*d**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*c**(sympy.S(5)/4)) + sqrt(2)*(-2*sqrt(a)*sqrt(c)*d*e + a*e**2 + 3*c*d**2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*c**(sympy.S(5)/4)) - sqrt(2)*(2*sqrt(a)*sqrt(c)*d*e + a*e**2 + 3*c*d**2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*c**(sympy.S(5)/4)) + sqrt(2)*(2*sqrt(a)*sqrt(c)*d*e + a*e**2 + 3*c*d**2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*c**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_146():
    f = (d + e*x**2)/(a + c*x**4)**2
    F = x*(d + e*x**2)/(4*a*(a + c*x**4)) - sqrt(2)*(-sqrt(a)*e + 3*sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*c**(sympy.S(3)/4)) + sqrt(2)*(-sqrt(a)*e + 3*sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*c**(sympy.S(3)/4)) - sqrt(2)*(sqrt(a)*e + 3*sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*c**(sympy.S(3)/4)) + sqrt(2)*(sqrt(a)*e + 3*sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*c**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_147():
    f = (a + c*x**4)**(-2)
    F = x/(4*a*(a + c*x**4)) - 3*sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*c**(sympy.S(1)/4)) + 3*sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*c**(sympy.S(1)/4)) - 3*sqrt(2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*c**(sympy.S(1)/4)) + 3*sqrt(2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*c**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_148():
    f = 1/((a + c*x**4)**2*(d + e*x**2))
    F = e**(sympy.S(7)/2)*atan(sqrt(e)*x/sqrt(d))/(sqrt(d)*(a*e**2 + c*d**2)**2) + c*x*(d - e*x**2)/(4*a*(a + c*x**4)*(a*e**2 + c*d**2)) - sqrt(2)*c**(sympy.S(1)/4)*e**2*(-sqrt(a)*e + sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) + sqrt(2)*c**(sympy.S(1)/4)*e**2*(-sqrt(a)*e + sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) - sqrt(2)*c**(sympy.S(1)/4)*e**2*(sqrt(a)*e + sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) + sqrt(2)*c**(sympy.S(1)/4)*e**2*(sqrt(a)*e + sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**2) - sqrt(2)*c**(sympy.S(1)/4)*(-sqrt(a)*e + 3*sqrt(c)*d)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(a*e**2 + c*d**2)) + sqrt(2)*c**(sympy.S(1)/4)*(-sqrt(a)*e + 3*sqrt(c)*d)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(a*e**2 + c*d**2)) - sqrt(2)*c**(sympy.S(1)/4)*(sqrt(a)*e + 3*sqrt(c)*d)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*(a*e**2 + c*d**2)) + sqrt(2)*c**(sympy.S(1)/4)*(sqrt(a)*e + 3*sqrt(c)*d)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*(a*e**2 + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_149():
    f = 1/((a + c*x**4)**2*(d + e*x**2)**2)
    F = 4*c*sqrt(d)*e**(sympy.S(7)/2)*atan(sqrt(e)*x/sqrt(d))/(a*e**2 + c*d**2)**3 + e**4*x/(2*d*(d + e*x**2)*(a*e**2 + c*d**2)**2) + e**(sympy.S(7)/2)*atan(sqrt(e)*x/sqrt(d))/(2*d**(sympy.S(3)/2)*(a*e**2 + c*d**2)**2) + c*x*(-a*e**2 + c*d**2 - 2*c*d*e*x**2)/(4*a*(a + c*x**4)*(a*e**2 + c*d**2)**2) - sqrt(2)*c**(sympy.S(3)/4)*e**2*(-4*sqrt(a)*sqrt(c)*d*e - a*e**2 + 3*c*d**2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**3) + sqrt(2)*c**(sympy.S(3)/4)*e**2*(-4*sqrt(a)*sqrt(c)*d*e - a*e**2 + 3*c*d**2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**3) - sqrt(2)*c**(sympy.S(3)/4)*e**2*(4*sqrt(a)*sqrt(c)*d*e - a*e**2 + 3*c*d**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**3) + sqrt(2)*c**(sympy.S(3)/4)*e**2*(4*sqrt(a)*sqrt(c)*d*e - a*e**2 + 3*c*d**2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**2 + c*d**2)**3) - sqrt(2)*c**(sympy.S(3)/4)*(-2*sqrt(a)*sqrt(c)*d*e - 3*a*e**2 + 3*c*d**2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(a*e**2 + c*d**2)**2) + sqrt(2)*c**(sympy.S(3)/4)*(-2*sqrt(a)*sqrt(c)*d*e - 3*a*e**2 + 3*c*d**2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(a*e**2 + c*d**2)**2) - sqrt(2)*c**(sympy.S(3)/4)*(2*sqrt(a)*sqrt(c)*d*e - 3*a*e**2 + 3*c*d**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*(a*e**2 + c*d**2)**2) + sqrt(2)*c**(sympy.S(3)/4)*(2*sqrt(a)*sqrt(c)*d*e - 3*a*e**2 + 3*c*d**2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*(a*e**2 + c*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_150():
    f = (d + e*x**2)**4/sqrt(a + c*x**4)
    F = -4*a**(sympy.S(1)/4)*d*e*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-3*a*e**2 + 5*c*d**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*c**(sympy.S(7)/4)*sqrt(a + c*x**4)) + 4*d*e**3*x**3*sqrt(a + c*x**4)/(5*c) + e**4*x**5*sqrt(a + c*x**4)/(7*c) + e**2*x*sqrt(a + c*x**4)*(-5*a*e**2 + 42*c*d**2)/(21*c**2) + 4*d*e*x*sqrt(a + c*x**4)*(-3*a*e**2 + 5*c*d**2)/(5*c**(sympy.S(3)/2)*(sqrt(a) + sqrt(c)*x**2)) + sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-252*a**(sympy.S(3)/2)*sqrt(c)*d*e**3 + 420*sqrt(a)*c**(sympy.S(3)/2)*d**3*e + 25*a**2*e**4 - 210*a*c*d**2*e**2 + 105*c**2*d**4)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(210*a**(sympy.S(1)/4)*c**(sympy.S(9)/4)*sqrt(a + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_151():
    f = (d + e*x**2)**3/sqrt(a + c*x**4)
    F = -3*a**(sympy.S(1)/4)*e*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-a*e**2 + 5*c*d**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*c**(sympy.S(7)/4)*sqrt(a + c*x**4)) + a**(sympy.S(1)/4)*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-3*a*e**3 + 15*c*d**2*e + 5*sqrt(c)*d*(-a*e**2 + c*d**2)/sqrt(a))*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(10*c**(sympy.S(7)/4)*sqrt(a + c*x**4)) + d*e**2*x*sqrt(a + c*x**4)/c + e**3*x**3*sqrt(a + c*x**4)/(5*c) + 3*e*x*sqrt(a + c*x**4)*(-a*e**2 + 5*c*d**2)/(5*c**(sympy.S(3)/2)*(sqrt(a) + sqrt(c)*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_152():
    f = (d + e*x**2)**2/sqrt(a + c*x**4)
    F = -2*a**(sympy.S(1)/4)*d*e*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(c**(sympy.S(3)/4)*sqrt(a + c*x**4)) + e**2*x*sqrt(a + c*x**4)/(3*c) + 2*d*e*x*sqrt(a + c*x**4)/(sqrt(c)*(sqrt(a) + sqrt(c)*x**2)) + sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(6*sqrt(a)*sqrt(c)*d*e - a*e**2 + 3*c*d**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(6*a**(sympy.S(1)/4)*c**(sympy.S(5)/4)*sqrt(a + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_153():
    f = (d + e*x**2)/sqrt(a + c*x**4)
    F = -a**(sympy.S(1)/4)*e*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(c**(sympy.S(3)/4)*sqrt(a + c*x**4)) + a**(sympy.S(1)/4)*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(e + sqrt(c)*d/sqrt(a))*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*c**(sympy.S(3)/4)*sqrt(a + c*x**4)) + e*x*sqrt(a + c*x**4)/(sqrt(c)*(sqrt(a) + sqrt(c)*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_154():
    f = 1/(sqrt(a + c*x**4)*(d + e*x**2))
    F = ((sympy.sqrt(Symbol('e')) * sympy.atan(((sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2))))) * x) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('d')) * sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('c')) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * Symbol('e')))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt(Symbol('a')))**(Integer(-1))) + Symbol('e')))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * Symbol('e')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * Symbol('d') * Symbol('e')))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('d') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_155():
    f = 1/(sqrt(a + c*x**4)*(d + e*x**2)**2)
    F = (Integer(-1) * ((sympy.sqrt(Symbol('c')) * Symbol('e') * x * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))) * ((Integer(2) * Symbol('d') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2))))))**(Integer(-1)))) + (((Symbol('e'))**(Integer(2)) * x * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))) * ((Integer(2) * Symbol('d') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))) + ((sympy.sqrt(Symbol('e')) * ((Integer(3) * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * sympy.atan(((sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2))))) * x) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('e') * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * Symbol('d') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * Symbol('d') * ((sympy.sqrt(Symbol('c')) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * Symbol('e')))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) + (sympy.sqrt(Symbol('a')) * Symbol('e'))) * ((Integer(3) * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * Symbol('e')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * Symbol('d') * Symbol('e')))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**(Integer(2)) * ((sympy.sqrt(Symbol('c')) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * Symbol('e')))) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_156():
    f = 1/(sqrt(a + c*x**4)*(d + e*x**2)**3)
    F = (Integer(-1) * ((Integer(3) * sympy.sqrt(Symbol('c')) * Symbol('e') * ((Integer(3) * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * x * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))) * ((Integer(8) * (Symbol('d'))**(Integer(2)) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2))))))**(Integer(-1)))) + (((Symbol('e'))**(Integer(2)) * x * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))) * ((Integer(4) * Symbol('d') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * (Symbol('e'))**(Integer(2)) * ((Integer(3) * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * x * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))) * ((Integer(8) * (Symbol('d'))**(Integer(2)) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(3) * sympy.sqrt(Symbol('e')) * ((Integer(5) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(4))) + (Integer(2) * Symbol('a') * Symbol('c') * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + ((Symbol('a'))**(Integer(2)) * (Symbol('e'))**(Integer(4)))) * sympy.atan(((sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2))))) * x) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(16) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(3) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('e') * ((Integer(3) * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('d'))**(Integer(2)) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(4) * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * Symbol('d') * Symbol('e'))) + (Integer(3) * Symbol('a') * (Symbol('e'))**(Integer(2)))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(8) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**(Integer(2)) * ((sympy.sqrt(Symbol('c')) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * Symbol('e')))) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * ((sympy.sqrt(Symbol('c')) * Symbol('d')) + (sympy.sqrt(Symbol('a')) * Symbol('e'))) * ((Integer(5) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(4))) + (Integer(2) * Symbol('a') * Symbol('c') * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + ((Symbol('a'))**(Integer(2)) * (Symbol('e'))**(Integer(4)))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * Symbol('e')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * Symbol('d') * Symbol('e')))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(32) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**(Integer(3)) * ((sympy.sqrt(Symbol('c')) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * Symbol('e')))) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_157():
    f = (d + e*x**2)**3/sqrt(a - c*x**4)
    F = 3*a**(sympy.S(3)/4)*e*sqrt(1 - c*x**4/a)*(a*e**2 + 5*c*d**2)*elliptic_e(asin(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), -1)/(5*c**(sympy.S(7)/4)*sqrt(a - c*x**4)) + a**(sympy.S(3)/4)*sqrt(1 - c*x**4/a)*(-3*e*(a*e**2 + 5*c*d**2) + 5*sqrt(c)*d*(a*e**2 + c*d**2)/sqrt(a))*elliptic_f(asin(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), -1)/(5*c**(sympy.S(7)/4)*sqrt(a - c*x**4)) - d*e**2*x*sqrt(a - c*x**4)/c - e**3*x**3*sqrt(a - c*x**4)/(5*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_158():
    f = (d + e*x**2)**2/sqrt(a - c*x**4)
    F = 2*a**(sympy.S(3)/4)*d*e*sqrt(1 - c*x**4/a)*elliptic_e(asin(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), -1)/(c**(sympy.S(3)/4)*sqrt(a - c*x**4)) + a**(sympy.S(1)/4)*sqrt(1 - c*x**4/a)*(-6*sqrt(a)*sqrt(c)*d*e + a*e**2 + 3*c*d**2)*elliptic_f(asin(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), -1)/(3*c**(sympy.S(5)/4)*sqrt(a - c*x**4)) - e**2*x*sqrt(a - c*x**4)/(3*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_159():
    f = (d + e*x**2)/sqrt(a - c*x**4)
    F = a**(sympy.S(3)/4)*e*sqrt(1 - c*x**4/a)*elliptic_e(asin(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), -1)/(c**(sympy.S(3)/4)*sqrt(a - c*x**4)) + a**(sympy.S(3)/4)*sqrt(1 - c*x**4/a)*(-e + sqrt(c)*d/sqrt(a))*elliptic_f(asin(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), -1)/(c**(sympy.S(3)/4)*sqrt(a - c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_160():
    f = 1/(sqrt(a - c*x**4)*(d + e*x**2))
    F = ((Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * Symbol('e')) * ((sympy.sqrt(Symbol('c')) * Symbol('d')))**(Integer(-1)))), sympy.asin((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_161():
    f = 1/(sqrt(a - c*x**4)*(d + e*x**2)**2)
    F = (Integer(-1) * (((Symbol('e'))**(Integer(2)) * x * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))) * ((Integer(2) * Symbol('d') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))) * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('e') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('d') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * Symbol('d') * ((sympy.sqrt(Symbol('c')) * Symbol('d')) + (sympy.sqrt(Symbol('a')) * Symbol('e'))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))))**(Integer(-1)))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Integer(3) * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * Symbol('e')) * ((sympy.sqrt(Symbol('c')) * Symbol('d')))**(Integer(-1)))), sympy.asin((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(2) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**(Integer(2)) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_162():
    f = 1/(sqrt(a - c*x**4)*(d + e*x**2)**3)
    F = (Integer(-1) * (((Symbol('e'))**(Integer(2)) * x * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))) * ((Integer(4) * Symbol('d') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('e'))**(Integer(2)) * ((Integer(3) * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))) * x * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))) * ((Integer(8) * (Symbol('d'))**(Integer(2)) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))))**(Integer(2)) * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('e') * ((Integer(3) * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('d'))**(Integer(2)) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(7) * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * Symbol('d') * Symbol('e'))) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('e'))**(Integer(2))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('d'))**(Integer(2)) * ((sympy.sqrt(Symbol('c')) * Symbol('d')) + (sympy.sqrt(Symbol('a')) * Symbol('e'))) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Integer(5) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('c') * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2)))) + ((Symbol('a'))**(Integer(2)) * (Symbol('e'))**(Integer(4)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * Symbol('e')) * ((sympy.sqrt(Symbol('c')) * Symbol('d')))**(Integer(-1)))), sympy.asin((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(8) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**(Integer(3)) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_163():
    f = 1/(sqrt(a - c*x**4)*(d + e*x**2)**4)
    F = (Integer(-1) * (((Symbol('e'))**(Integer(2)) * x * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))) * ((Integer(6) * Symbol('d') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Symbol('e'))**(Integer(2)) * ((Integer(3) * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))) * x * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))) * ((Integer(24) * (Symbol('d'))**(Integer(2)) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))))**(Integer(2)) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * ((Integer(29) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Integer(14) * Symbol('a') * Symbol('c') * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2)))) + (Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('e'))**(Integer(4)))) * x * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))) * ((Integer(16) * (Symbol('d'))**(Integer(3)) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))))**(Integer(3)) * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('e') * ((Integer(29) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Integer(14) * Symbol('a') * Symbol('c') * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2)))) + (Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('e'))**(Integer(4)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(16) * (Symbol('d'))**(Integer(3)) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Integer(57) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Integer(30) * sympy.sqrt(Symbol('a')) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**(Integer(3)) * Symbol('e'))) + (Integer(-1) * (Integer(32) * Symbol('a') * Symbol('c') * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2)))) + (Integer(10) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('c')) * Symbol('d') * (Symbol('e'))**(Integer(3))) + (Integer(15) * (Symbol('a'))**(Integer(2)) * (Symbol('e'))**(Integer(4)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(48) * (Symbol('d'))**(Integer(3)) * (((sympy.sqrt(Symbol('c')) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * Symbol('e')))))**(Integer(2)) * (((sympy.sqrt(Symbol('c')) * Symbol('d')) + (sympy.sqrt(Symbol('a')) * Symbol('e'))))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))))**(Integer(-1)))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Integer(35) * (Symbol('c'))**(Integer(3)) * (Symbol('d'))**(Integer(6))) + (Integer(-1) * (Integer(7) * Symbol('a') * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(4)) * (Symbol('e'))**(Integer(2)))) + (Integer(17) * (Symbol('a'))**(Integer(2)) * Symbol('c') * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(4))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * (Symbol('e'))**(Integer(6))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * Symbol('e')) * ((sympy.sqrt(Symbol('c')) * Symbol('d')))**(Integer(-1)))), sympy.asin((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * ((Integer(16) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**(Integer(4)) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_164():
    f = (d + e*x**2)/sqrt(-a + c*x**4)
    F = a**(sympy.S(3)/4)*e*sqrt(1 - c*x**4/a)*elliptic_e(asin(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), -1)/(c**(sympy.S(3)/4)*sqrt(-a + c*x**4)) + a**(sympy.S(3)/4)*sqrt(1 - c*x**4/a)*(-e + sqrt(c)*d/sqrt(a))*elliptic_f(asin(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), -1)/(c**(sympy.S(3)/4)*sqrt(-a + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_165():
    f = 1/(sqrt(-a + c*x**4)*(d + e*x**2))
    F = ((Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * ((sympy.sqrt(Symbol('a')) * Symbol('e')) * ((sympy.sqrt(Symbol('c')) * Symbol('d')))**(Integer(-1)))), sympy.asin((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('d') * sympy.sqrt(((Integer(-1) * Symbol('a')) + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_166():
    f = (sqrt(a) + sqrt(c)*x**2)/sqrt(-a + c*x**4)
    F = a**(sympy.S(3)/4)*sqrt(1 - c*x**4/a)*elliptic_e(asin(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), -1)/(c**(sympy.S(1)/4)*sqrt(-a + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_167():
    f = (x**2*sqrt(c/a) + 1)/sqrt(-a + c*x**4)
    F = sqrt(1 - c*x**4/a)*elliptic_e(asin(x*(c/a)**(sympy.S(1)/4)), -1)/((c/a)**(sympy.S(1)/4)*sqrt(-a + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_168():
    f = (d + e*x**2)/sqrt(-a - c*x**4)
    F = -a**(sympy.S(1)/4)*e*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(c**(sympy.S(3)/4)*sqrt(-a - c*x**4)) + a**(sympy.S(1)/4)*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(e + sqrt(c)*d/sqrt(a))*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*c**(sympy.S(3)/4)*sqrt(-a - c*x**4)) - e*x*sqrt(-a - c*x**4)/(sqrt(c)*(sqrt(a) + sqrt(c)*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_169():
    f = 1/(sqrt(-a - c*x**4)*(d + e*x**2))
    F = ((sympy.sqrt(Symbol('e')) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2)))))) * x) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(Symbol('e')) * sympy.sqrt(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('d')) * sympy.sqrt((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('c')) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * Symbol('e')))) * sympy.sqrt(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt(Symbol('a')))**(Integer(-1))) + Symbol('e')))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * Symbol('e')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * Symbol('d') * Symbol('e')))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('d') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))) * sympy.sqrt(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_170():
    f = 1/(sqrt(4 - 5*x**4)*(a + b*x**2))
    F = sympy.Function('EllipticPi')((Integer(-1) * ((Integer(2) * Symbol('b')) * ((sympy.sqrt(Integer(5)) * Symbol('a')))**(Integer(-1)))), sympy.asin((((Integer(5))**((Integer(4))**(Integer(-1))) * x) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-1)) * ((sympy.sqrt(Integer(2)) * (Integer(5))**((Integer(4))**(Integer(-1))) * Symbol('a')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_171():
    f = 1/((a + b*x**2)*sqrt(5*x**4 + 4))
    F = ((sympy.sqrt(Symbol('b')) * sympy.atan(((sympy.sqrt(((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(4) * (Symbol('b'))**(Integer(2))))) * x) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(4) + (Integer(5) * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt(((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(4) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))) + (((Integer(5))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Integer(5)) * Symbol('a')) + (Integer(2) * Symbol('b'))) * (Integer(2) + (sympy.sqrt(Integer(5)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(4) + (Integer(5) * (x)**(Integer(4)))) * (((Integer(2) + (sympy.sqrt(Integer(5)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Integer(5))**((Integer(4))**(Integer(-1))) * x) * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * sympy.sqrt(Integer(2)) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt((Integer(4) + (Integer(5) * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((((sympy.sqrt(Integer(5)) * Symbol('a')) + (Integer(2) * Symbol('b'))))**(Integer(2)) * (Integer(2) + (sympy.sqrt(Integer(5)) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(4) + (Integer(5) * (x)**(Integer(4)))) * (((Integer(2) + (sympy.sqrt(Integer(5)) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Integer(5)) * Symbol('a')) + (Integer(-1) * (Integer(2) * Symbol('b')))))**(Integer(2)) * ((Integer(8) * sympy.sqrt(Integer(5)) * Symbol('a') * Symbol('b')))**(Integer(-1)))), (Integer(2) * sympy.atan((((Integer(5))**((Integer(4))**(Integer(-1))) * x) * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * sympy.sqrt(Integer(2)) * (Integer(5))**((Integer(4))**(Integer(-1))) * Symbol('a') * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt((Integer(4) + (Integer(5) * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_172():
    f = 1/((a + b*x**2)*sqrt(-d*x**4 + 4))
    F = sympy.Function('EllipticPi')((Integer(-1) * ((Integer(2) * Symbol('b')) * ((Symbol('a') * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), sympy.asin((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-1)) * ((sympy.sqrt(Integer(2)) * Symbol('a') * (Symbol('d'))**((Integer(4))**(Integer(-1)))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_173():
    f = 1/((a + b*x**2)*sqrt(d*x**4 + 4))
    F = ((sympy.sqrt(Symbol('b')) * sympy.atan(((sympy.sqrt(((Integer(4) * (Symbol('b'))**(Integer(2))) + ((Symbol('a'))**(Integer(2)) * Symbol('d')))) * x) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * sympy.sqrt((Integer(4) + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt(((Integer(4) * (Symbol('b'))**(Integer(2))) + ((Symbol('a'))**(Integer(2)) * Symbol('d'))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**((Integer(4))**(Integer(-1))) * (Integer(2) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(4) + (Symbol('d') * (x)**(Integer(4)))) * (((Integer(2) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * sympy.sqrt(Integer(2)) * ((Integer(2) * Symbol('b')) + (Integer(-1) * (Symbol('a') * sympy.sqrt(Symbol('d'))))) * sympy.sqrt((Integer(4) + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))) + ((((Integer(2) * Symbol('b')) + (Symbol('a') * sympy.sqrt(Symbol('d')))) * (Integer(2) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))) * sympy.sqrt(((Integer(4) + (Symbol('d') * (x)**(Integer(4)))) * (((Integer(2) + (sympy.sqrt(Symbol('d')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((Integer(2) * Symbol('b')) + (Integer(-1) * (Symbol('a') * sympy.sqrt(Symbol('d'))))))**(Integer(2)) * ((Integer(8) * Symbol('a') * Symbol('b') * sympy.sqrt(Symbol('d'))))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('d'))**((Integer(4))**(Integer(-1))) * x) * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * sympy.sqrt(Integer(2)) * Symbol('a') * ((Integer(2) * Symbol('b')) + (Integer(-1) * (Symbol('a') * sympy.sqrt(Symbol('d'))))) * (Symbol('d'))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(4) + (Symbol('d') * (x)**(Integer(4)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_174():
    f = sqrt(a + b*x**2)/sqrt(1 - x**4)
    F = (Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(2))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (x)**(Integer(2)))) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**(Integer(-1)))) * sympy.Function('EllipticPi')((Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * x) * (sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Integer(1) + (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * (x)**(Integer(2))))) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_175():
    f = (a + b*x**4)**p*(c + e*x**2)**q
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('e') * (x)**(Integer(2)))))**(Symbol('q')) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_176():
    f = (a + b*x**4)**p*(c + e*x**2)**3
    F = c**3*x*(a + b*x**4)**p*hyper((sympy.S(1)/4, -p), (sympy.S(5)/4,), -b*x**4/a)/(1 + b*x**4/a)**p + 3*c*e**2*x**5*(a + b*x**4)**p*hyper((sympy.S(5)/4, -p), (sympy.S(9)/4,), -b*x**4/a)/(5*(1 + b*x**4/a)**p) + e**3*x**3*(a + b*x**4)**(p + 1)/(b*(4*p + 7)) - e*x**3*(a + b*x**4)**p*(a*e**2 - b*c**2*(4*p + 7))*hyper((sympy.S(3)/4, -p), (sympy.S(7)/4,), -b*x**4/a)/(b*(1 + b*x**4/a)**p*(4*p + 7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_177():
    f = (a + b*x**4)**p*(c + e*x**2)**2
    F = 2*c*e*x**3*(a + b*x**4)**p*hyper((sympy.S(3)/4, -p), (sympy.S(7)/4,), -b*x**4/a)/(3*(1 + b*x**4/a)**p) + e**2*x*(a + b*x**4)**(p + 1)/(b*(4*p + 5)) - x*(a + b*x**4)**p*(a*e**2 - b*c**2*(4*p + 5))*hyper((sympy.S(1)/4, -p), (sympy.S(5)/4,), -b*x**4/a)/(b*(1 + b*x**4/a)**p*(4*p + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_178():
    f = (a + b*x**4)**p*(c + e*x**2)
    F = c*x*(a + b*x**4)**p*hyper((sympy.S(1)/4, -p), (sympy.S(5)/4,), -b*x**4/a)/(1 + b*x**4/a)**p + e*x**3*(a + b*x**4)**p*hyper((sympy.S(3)/4, -p), (sympy.S(7)/4,), -b*x**4/a)/(3*(1 + b*x**4/a)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_179():
    f = (a + b*x**4)**p
    F = x*(a + b*x**4)**p*hyper((sympy.S(1)/4, -p), (sympy.S(5)/4,), -b*x**4/a)/(1 + b*x**4/a)**p
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_180():
    f = (a + b*x**4)**p/(c + e*x**2)
    F = x*(a + b*x**4)**p*appellf1(sympy.S(1)/4, 1, -p, sympy.S(5)/4, e**2*x**4/c**2, -b*x**4/a)/(c*(1 + b*x**4/a)**p) - e*x**3*(a + b*x**4)**p*appellf1(sympy.S(3)/4, 1, -p, sympy.S(7)/4, e**2*x**4/c**2, -b*x**4/a)/(3*c**2*(1 + b*x**4/a)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_181():
    f = (a + b*x**4)**p/(c + e*x**2)**2
    F = x*(a + b*x**4)**p*appellf1(sympy.S(1)/4, 2, -p, sympy.S(5)/4, e**2*x**4/c**2, -b*x**4/a)/(c**2*(1 + b*x**4/a)**p) - 2*e*x**3*(a + b*x**4)**p*appellf1(sympy.S(3)/4, 2, -p, sympy.S(7)/4, e**2*x**4/c**2, -b*x**4/a)/(3*c**3*(1 + b*x**4/a)**p) + e**2*x**5*(a + b*x**4)**p*appellf1(sympy.S(5)/4, 2, -p, sympy.S(9)/4, e**2*x**4/c**2, -b*x**4/a)/(5*c**4*(1 + b*x**4/a)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_182():
    f = (1 - x**2)**3*(b*x**4 + 1)**p
    F = 3*x**5*hyper((sympy.S(5)/4, -p), (sympy.S(9)/4,), -b*x**4)/5 + x*hyper((sympy.S(1)/4, -p), (sympy.S(5)/4,), -b*x**4) - x**3*(b*x**4 + 1)**(p + 1)/(b*(4*p + 7)) + x**3*(-b*(4*p + 7) + 1)*hyper((sympy.S(3)/4, -p), (sympy.S(7)/4,), -b*x**4)/(b*(4*p + 7))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_183():
    f = (1 - x**2)**2*(b*x**4 + 1)**p
    F = -2*x**3*hyper((sympy.S(3)/4, -p), (sympy.S(7)/4,), -b*x**4)/3 + x*(b*x**4 + 1)**(p + 1)/(b*(4*p + 5)) - x*(-b*(4*p + 5) + 1)*hyper((sympy.S(1)/4, -p), (sympy.S(5)/4,), -b*x**4)/(b*(4*p + 5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_184():
    f = (1 - x**2)*(b*x**4 + 1)**p
    F = -x**3*hyper((sympy.S(3)/4, -p), (sympy.S(7)/4,), -b*x**4)/3 + x*hyper((sympy.S(1)/4, -p), (sympy.S(5)/4,), -b*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_185():
    f = (b*x**4 + 1)**p
    F = x*hyper((sympy.S(1)/4, -p), (sympy.S(5)/4,), -b*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_186():
    f = (b*x**4 + 1)**p/(1 - x**2)
    F = x**3*appellf1(sympy.S(3)/4, 1, -p, sympy.S(7)/4, x**4, -b*x**4)/3 + x*appellf1(sympy.S(1)/4, 1, -p, sympy.S(5)/4, x**4, -b*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_187():
    f = (b*x**4 + 1)**p/(1 - x**2)**2
    F = x**5*appellf1(sympy.S(5)/4, 2, -p, sympy.S(9)/4, x**4, -b*x**4)/5 + 2*x**3*appellf1(sympy.S(3)/4, 2, -p, sympy.S(7)/4, x**4, -b*x**4)/3 + x*appellf1(sympy.S(1)/4, 2, -p, sympy.S(5)/4, x**4, -b*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_188():
    f = (b*x**4 + 1)**p/(1 - x**2)**3
    F = x**7*appellf1(sympy.S(7)/4, 3, -p, sympy.S(11)/4, x**4, -b*x**4)/7 + 3*x**5*appellf1(sympy.S(5)/4, 3, -p, sympy.S(9)/4, x**4, -b*x**4)/5 + x**3*appellf1(sympy.S(3)/4, 3, -p, sympy.S(7)/4, x**4, -b*x**4) + x*appellf1(sympy.S(1)/4, 3, -p, sympy.S(5)/4, x**4, -b*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_189():
    f = (d + e*x**2)**4/(d**2 - e**2*x**4)
    F = 8*d**(sympy.S(5)/2)*atanh(sqrt(e)*x/sqrt(d))/sqrt(e) - 7*d**2*x - 4*d*e*x**3/3 - e**2*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_190():
    f = (d + e*x**2)**3/(d**2 - e**2*x**4)
    F = 4*d**(sympy.S(3)/2)*atanh(sqrt(e)*x/sqrt(d))/sqrt(e) - 3*d*x - e*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_191():
    f = (d + e*x**2)**2/(d**2 - e**2*x**4)
    F = 2*sqrt(d)*atanh(sqrt(e)*x/sqrt(d))/sqrt(e) - x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_192():
    f = (d + e*x**2)/(d**2 - e**2*x**4)
    F = atanh(sqrt(e)*x/sqrt(d))/(sqrt(d)*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_193():
    f = 1/((d + e*x**2)*(d**2 - e**2*x**4))
    F = x/(4*d**2*(d + e*x**2)) + atan(sqrt(e)*x/sqrt(d))/(2*d**(sympy.S(5)/2)*sqrt(e)) + atanh(sqrt(e)*x/sqrt(d))/(4*d**(sympy.S(5)/2)*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_194():
    f = 1/((d + e*x**2)**2*(d**2 - e**2*x**4))
    F = x/(8*d**2*(d + e*x**2)**2) + 5*x/(16*d**3*(d + e*x**2)) + 7*atan(sqrt(e)*x/sqrt(d))/(16*d**(sympy.S(7)/2)*sqrt(e)) + atanh(sqrt(e)*x/sqrt(d))/(8*d**(sympy.S(7)/2)*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_195():
    f = (d + e*x**2)**(sympy.S(3)/2)/(d**2 - e**2*x**4)
    F = -atanh(sqrt(e)*x/sqrt(d + e*x**2))/sqrt(e) + sqrt(2)*atanh(sqrt(2)*sqrt(e)*x/sqrt(d + e*x**2))/sqrt(e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_196():
    f = sqrt(d + e*x**2)/(d**2 - e**2*x**4)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(e)*x/sqrt(d + e*x**2))/(2*d*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_197():
    f = 1/(sqrt(d + e*x**2)*(d**2 - e**2*x**4))
    F = x/(2*d**2*sqrt(d + e*x**2)) + sqrt(2)*atanh(sqrt(2)*sqrt(e)*x/sqrt(d + e*x**2))/(4*d**2*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_198():
    f = 1/((d + e*x**2)**(sympy.S(3)/2)*(d**2 - e**2*x**4))
    F = x/(6*d**2*(d + e*x**2)**(sympy.S(3)/2)) + 7*x/(12*d**3*sqrt(d + e*x**2)) + sqrt(2)*atanh(sqrt(2)*sqrt(e)*x/sqrt(d + e*x**2))/(8*d**3*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_199():
    f = (a + b*x**2)**(sympy.S(5)/2)/sqrt(a**2 - b**2*x**4)
    F = 19*a**2*sqrt(a - b*x**2)*sqrt(a + b*x**2)*atan(sqrt(b)*x/sqrt(a - b*x**2))/(8*sqrt(b)*sqrt(a**2 - b**2*x**4)) - 9*a*x*(a - b*x**2)*sqrt(a + b*x**2)/(8*sqrt(a**2 - b**2*x**4)) - x*(a - b*x**2)*(a + b*x**2)**(sympy.S(3)/2)/(4*sqrt(a**2 - b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_200():
    f = (a + b*x**2)**(sympy.S(3)/2)/sqrt(a**2 - b**2*x**4)
    F = 3*a*sqrt(a - b*x**2)*sqrt(a + b*x**2)*atan(sqrt(b)*x/sqrt(a - b*x**2))/(2*sqrt(b)*sqrt(a**2 - b**2*x**4)) - x*(a - b*x**2)*sqrt(a + b*x**2)/(2*sqrt(a**2 - b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_201():
    f = sqrt(a + b*x**2)/sqrt(a**2 - b**2*x**4)
    F = sqrt(a - b*x**2)*sqrt(a + b*x**2)*atan(sqrt(b)*x/sqrt(a - b*x**2))/(sqrt(b)*sqrt(a**2 - b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_202():
    f = 1/(sqrt(a + b*x**2)*sqrt(a**2 - b**2*x**4))
    F = sqrt(2)*sqrt(a - b*x**2)*sqrt(a + b*x**2)*atan(sqrt(2)*sqrt(b)*x/sqrt(a - b*x**2))/(2*a*sqrt(b)*sqrt(a**2 - b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_203():
    f = 1/((a + b*x**2)**(sympy.S(3)/2)*sqrt(a**2 - b**2*x**4))
    F = x*(a - b*x**2)/(4*a**2*sqrt(a + b*x**2)*sqrt(a**2 - b**2*x**4)) + 3*sqrt(2)*sqrt(a - b*x**2)*sqrt(a + b*x**2)*atan(sqrt(2)*sqrt(b)*x/sqrt(a - b*x**2))/(8*a**2*sqrt(b)*sqrt(a**2 - b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_204():
    f = 1/((a + b*x**2)**(sympy.S(5)/2)*sqrt(a**2 - b**2*x**4))
    F = x*(a - b*x**2)/(8*a**2*(a + b*x**2)**(sympy.S(3)/2)*sqrt(a**2 - b**2*x**4)) + 9*x*(a - b*x**2)/(32*a**3*sqrt(a + b*x**2)*sqrt(a**2 - b**2*x**4)) + 19*sqrt(2)*sqrt(a - b*x**2)*sqrt(a + b*x**2)*atan(sqrt(2)*sqrt(b)*x/sqrt(a - b*x**2))/(64*a**3*sqrt(b)*sqrt(a**2 - b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_205():
    f = (a - b*x**2)**(sympy.S(5)/2)/sqrt(a**2 - b**2*x**4)
    F = 19*a**2*sqrt(a - b*x**2)*sqrt(a + b*x**2)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(8*sqrt(b)*sqrt(a**2 - b**2*x**4)) - 9*a*x*sqrt(a - b*x**2)*(a + b*x**2)/(8*sqrt(a**2 - b**2*x**4)) - x*(a - b*x**2)**(sympy.S(3)/2)*(a + b*x**2)/(4*sqrt(a**2 - b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_206():
    f = (a - b*x**2)**(sympy.S(3)/2)/sqrt(a**2 - b**2*x**4)
    F = 3*a*sqrt(a - b*x**2)*sqrt(a + b*x**2)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(2*sqrt(b)*sqrt(a**2 - b**2*x**4)) - x*sqrt(a - b*x**2)*(a + b*x**2)/(2*sqrt(a**2 - b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_207():
    f = sqrt(a - b*x**2)/sqrt(a**2 - b**2*x**4)
    F = sqrt(a - b*x**2)*sqrt(a + b*x**2)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(sqrt(b)*sqrt(a**2 - b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_208():
    f = 1/(sqrt(a - b*x**2)*sqrt(a**2 - b**2*x**4))
    F = sqrt(2)*sqrt(a - b*x**2)*sqrt(a + b*x**2)*atanh(sqrt(2)*sqrt(b)*x/sqrt(a + b*x**2))/(2*a*sqrt(b)*sqrt(a**2 - b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_209():
    f = 1/((a - b*x**2)**(sympy.S(3)/2)*sqrt(a**2 - b**2*x**4))
    F = x*(a + b*x**2)/(4*a**2*sqrt(a - b*x**2)*sqrt(a**2 - b**2*x**4)) + 3*sqrt(2)*sqrt(a - b*x**2)*sqrt(a + b*x**2)*atanh(sqrt(2)*sqrt(b)*x/sqrt(a + b*x**2))/(8*a**2*sqrt(b)*sqrt(a**2 - b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_210():
    f = 1/((a - b*x**2)**(sympy.S(5)/2)*sqrt(a**2 - b**2*x**4))
    F = x*(a + b*x**2)/(8*a**2*(a - b*x**2)**(sympy.S(3)/2)*sqrt(a**2 - b**2*x**4)) + 9*x*(a + b*x**2)/(32*a**3*sqrt(a - b*x**2)*sqrt(a**2 - b**2*x**4)) + 19*sqrt(2)*sqrt(a - b*x**2)*sqrt(a + b*x**2)*atanh(sqrt(2)*sqrt(b)*x/sqrt(a + b*x**2))/(64*a**3*sqrt(b)*sqrt(a**2 - b**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_211():
    f = sqrt(x**2 - 1)/sqrt(x**4 - 1)
    F = sqrt(x**2 - 1)*sqrt(x**2 + 1)*asinh(x)/sqrt(x**4 - 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_212():
    f = (d + e*x**2)**4/(b*d*e + b*e**2*x**2 - c*d**2 + c*e**2*x**4)
    F = e**2*x**5/(5*c) + e*x**3*(-b*e + 4*c*d)/(3*c**2) + x*(b**2*e**2 - 5*b*c*d*e + 7*c**2*d**2)/c**3 - (-b*e + 2*c*d)**3*atanh(sqrt(c)*sqrt(e)*x/sqrt(-b*e + c*d))/(c**(sympy.S(7)/2)*sqrt(e)*sqrt(-b*e + c*d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_213():
    f = (d + e*x**2)**3/(b*d*e + b*e**2*x**2 - c*d**2 + c*e**2*x**4)
    F = e*x**3/(3*c) + x*(-b*e + 3*c*d)/c**2 - (-b*e + 2*c*d)**2*atanh(sqrt(c)*sqrt(e)*x/sqrt(-b*e + c*d))/(c**(sympy.S(5)/2)*sqrt(e)*sqrt(-b*e + c*d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_214():
    f = (d + e*x**2)**2/(b*d*e + b*e**2*x**2 - c*d**2 + c*e**2*x**4)
    F = x/c - (-b*e + 2*c*d)*atanh(sqrt(c)*sqrt(e)*x/sqrt(-b*e + c*d))/(c**(sympy.S(3)/2)*sqrt(e)*sqrt(-b*e + c*d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_215():
    f = (d + e*x**2)/(b*d*e + b*e**2*x**2 - c*d**2 + c*e**2*x**4)
    F = -atanh(sqrt(c)*sqrt(e)*x/sqrt(-b*e + c*d))/(sqrt(c)*sqrt(e)*sqrt(-b*e + c*d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_216():
    f = 1/((d + e*x**2)*(b*d*e + b*e**2*x**2 - c*d**2 + c*e**2*x**4))
    F = -c**(sympy.S(3)/2)*atanh(sqrt(c)*sqrt(e)*x/sqrt(-b*e + c*d))/(sqrt(e)*sqrt(-b*e + c*d)*(-b*e + 2*c*d)**2) - x/(2*d*(d + e*x**2)*(-b*e + 2*c*d)) - (-b*e + 4*c*d)*atan(sqrt(e)*x/sqrt(d))/(2*d**(sympy.S(3)/2)*sqrt(e)*(-b*e + 2*c*d)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_217():
    f = 1/((d + e*x**2)**2*(b*d*e + b*e**2*x**2 - c*d**2 + c*e**2*x**4))
    F = -c**(sympy.S(5)/2)*atanh(sqrt(c)*sqrt(e)*x/sqrt(-b*e + c*d))/(sqrt(e)*sqrt(-b*e + c*d)*(-b*e + 2*c*d)**3) - x/(4*d*(d + e*x**2)**2*(-b*e + 2*c*d)) - x*(-3*b*e + 10*c*d)/(8*d**2*(d + e*x**2)*(-b*e + 2*c*d)**2) - (3*b**2*e**2 - 16*b*c*d*e + 28*c**2*d**2)*atan(sqrt(e)*x/sqrt(d))/(8*d**(sympy.S(5)/2)*sqrt(e)*(-b*e + 2*c*d)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_218():
    f = (d + e*x**2)**(sympy.S(5)/2)/(b*d*e + b*e**2*x**2 - c*d**2 + c*e**2*x**4)
    F = x*sqrt(d + e*x**2)/(2*c) + (-2*b*e + 5*c*d)*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(2*c**2*sqrt(e)) - (-b*e + 2*c*d)**(sympy.S(3)/2)*atanh(sqrt(e)*x*sqrt(-b*e + 2*c*d)/(sqrt(d + e*x**2)*sqrt(-b*e + c*d)))/(c**2*sqrt(e)*sqrt(-b*e + c*d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_219():
    f = (d + e*x**2)**(sympy.S(3)/2)/(b*d*e + b*e**2*x**2 - c*d**2 + c*e**2*x**4)
    F = atanh(sqrt(e)*x/sqrt(d + e*x**2))/(c*sqrt(e)) - sqrt(-b*e + 2*c*d)*atanh(sqrt(e)*x*sqrt(-b*e + 2*c*d)/(sqrt(d + e*x**2)*sqrt(-b*e + c*d)))/(c*sqrt(e)*sqrt(-b*e + c*d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_220():
    f = sqrt(d + e*x**2)/(b*d*e + b*e**2*x**2 - c*d**2 + c*e**2*x**4)
    F = -atanh(sqrt(e)*x*sqrt(-b*e + 2*c*d)/(sqrt(d + e*x**2)*sqrt(-b*e + c*d)))/(sqrt(e)*sqrt(-b*e + c*d)*sqrt(-b*e + 2*c*d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_221():
    f = 1/(sqrt(d + e*x**2)*(b*d*e + b*e**2*x**2 - c*d**2 + c*e**2*x**4))
    F = -c*atanh(sqrt(e)*x*sqrt(-b*e + 2*c*d)/(sqrt(d + e*x**2)*sqrt(-b*e + c*d)))/(sqrt(e)*sqrt(-b*e + c*d)*(-b*e + 2*c*d)**(sympy.S(3)/2)) - x/(d*sqrt(d + e*x**2)*(-b*e + 2*c*d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_222():
    f = 1/((d + e*x**2)**(sympy.S(3)/2)*(b*d*e + b*e**2*x**2 - c*d**2 + c*e**2*x**4))
    F = -c**2*atanh(sqrt(e)*x*sqrt(-b*e + 2*c*d)/(sqrt(d + e*x**2)*sqrt(-b*e + c*d)))/(sqrt(e)*sqrt(-b*e + c*d)*(-b*e + 2*c*d)**(sympy.S(5)/2)) - x/(3*d*(d + e*x**2)**(sympy.S(3)/2)*(-b*e + 2*c*d)) - x*(-2*b*e + 7*c*d)/(3*d**2*sqrt(d + e*x**2)*(-b*e + 2*c*d)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_223():
    f = (x**2 + 1)**3*sqrt(x**4 + x**2 + 1)
    F = x**3*(x**4 + x**2 + 1)**(sympy.S(3)/2)/9 + 2*x*(6*x**2 + 7)*sqrt(x**4 + x**2 + 1)/45 + x*(x**4 + x**2 + 1)**(sympy.S(3)/2)/3 + 26*x*sqrt(x**4 + x**2 + 1)/(45*x**2 + 45) + sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(7*x**2 + 7)*elliptic_f(2*atan(x), sympy.S(1)/4)/(15*sqrt(x**4 + x**2 + 1)) - sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(26*x**2 + 26)*elliptic_e(2*atan(x), sympy.S(1)/4)/(45*sqrt(x**4 + x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_224():
    f = (x**2 + 1)**2*sqrt(x**4 + x**2 + 1)
    F = 2*x*(3*x**2 + 4)*sqrt(x**4 + x**2 + 1)/21 + x*(x**4 + x**2 + 1)**(sympy.S(3)/2)/7 + 2*x*sqrt(x**4 + x**2 + 1)/(3*x**2 + 3) - sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(2*x**2 + 2)*elliptic_e(2*atan(x), sympy.S(1)/4)/(3*sqrt(x**4 + x**2 + 1)) + sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(4*x**2 + 4)*elliptic_f(2*atan(x), sympy.S(1)/4)/(7*sqrt(x**4 + x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_225():
    f = (x**2 + 1)*sqrt(x**4 + x**2 + 1)
    F = x*(x**2 + 2)*sqrt(x**4 + x**2 + 1)/5 + 3*x*sqrt(x**4 + x**2 + 1)/(5*x**2 + 5) - sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(3*x**2 + 3)*elliptic_e(2*atan(x), sympy.S(1)/4)/(5*sqrt(x**4 + x**2 + 1)) + sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(3*x**2 + 3)*elliptic_f(2*atan(x), sympy.S(1)/4)/(5*sqrt(x**4 + x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_226():
    f = sqrt(x**4 + x**2 + 1)/(x**2 + 1)
    F = x*sqrt(x**4 + x**2 + 1)/(x**2 + 1) - sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(x**2 + 1)*elliptic_e(2*atan(x), sympy.S(1)/4)/sqrt(x**4 + x**2 + 1) + sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(3*x**2 + 3)*elliptic_f(2*atan(x), sympy.S(1)/4)/(4*sqrt(x**4 + x**2 + 1)) + atan(x/sqrt(x**4 + x**2 + 1))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_227():
    f = sqrt(x**4 + x**2 + 1)/(x**2 + 1)**2
    F = sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(x**2 + 1)*elliptic_e(2*atan(x), sympy.S(1)/4)/(2*sqrt(x**4 + x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_228():
    f = sqrt(x**4 + x**2 + 1)/(x**2 + 1)**3
    F = x*sqrt(x**4 + x**2 + 1)/(4*(x**2 + 1)**2) + sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(x**2 + 1)*elliptic_e(2*atan(x), sympy.S(1)/4)/(4*sqrt(x**4 + x**2 + 1)) + atan(x/sqrt(x**4 + x**2 + 1))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_229():
    f = sqrt(x**4 + x**2 + 1)/(x**2 + 1)**4
    F = x*sqrt(x**4 + x**2 + 1)/(6*(x**2 + 1)**2) + x*sqrt(x**4 + x**2 + 1)/(6*(x**2 + 1)**3) + sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(x**2 + 1)*elliptic_e(2*atan(x), sympy.S(1)/4)/(3*sqrt(x**4 + x**2 + 1)) - sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(x**2 + 1)*elliptic_f(2*atan(x), sympy.S(1)/4)/(8*sqrt(x**4 + x**2 + 1)) + atan(x/sqrt(x**4 + x**2 + 1))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_230():
    f = (x**2 + 1)**3/sqrt(x**4 + x**2 + 1)
    F = x**3*sqrt(x**4 + x**2 + 1)/5 + 11*x*sqrt(x**4 + x**2 + 1)/15 + 14*x*sqrt(x**4 + x**2 + 1)/(15*x**2 + 15) + sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(3*x**2 + 3)*elliptic_f(2*atan(x), sympy.S(1)/4)/(5*sqrt(x**4 + x**2 + 1)) - sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(14*x**2 + 14)*elliptic_e(2*atan(x), sympy.S(1)/4)/(15*sqrt(x**4 + x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_231():
    f = (x**2 + 1)**2/sqrt(x**4 + x**2 + 1)
    F = x*sqrt(x**4 + x**2 + 1)/3 + 4*x*sqrt(x**4 + x**2 + 1)/(3*x**2 + 3) + sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(x**2 + 1)*elliptic_f(2*atan(x), sympy.S(1)/4)/sqrt(x**4 + x**2 + 1) - sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(4*x**2 + 4)*elliptic_e(2*atan(x), sympy.S(1)/4)/(3*sqrt(x**4 + x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_232():
    f = (x**2 + 1)/sqrt(x**4 + x**2 + 1)
    F = x*sqrt(x**4 + x**2 + 1)/(x**2 + 1) - sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(x**2 + 1)*elliptic_e(2*atan(x), sympy.S(1)/4)/sqrt(x**4 + x**2 + 1) + sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(x**2 + 1)*elliptic_f(2*atan(x), sympy.S(1)/4)/sqrt(x**4 + x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_233():
    f = 1/((x**2 + 1)*sqrt(x**4 + x**2 + 1))
    F = sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(x**2 + 1)*elliptic_f(2*atan(x), sympy.S(1)/4)/(4*sqrt(x**4 + x**2 + 1)) + atan(x/sqrt(x**4 + x**2 + 1))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_234():
    f = 1/((x**2 + 1)**2*sqrt(x**4 + x**2 + 1))
    F = sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(x**2 + 1)*elliptic_e(2*atan(x), sympy.S(1)/4)/(2*sqrt(x**4 + x**2 + 1)) - sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(x**2 + 1)*elliptic_f(2*atan(x), sympy.S(1)/4)/(4*sqrt(x**4 + x**2 + 1)) + atan(x/sqrt(x**4 + x**2 + 1))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_235():
    f = 1/((x**2 + 1)**3*sqrt(x**4 + x**2 + 1))
    F = x*sqrt(x**4 + x**2 + 1)/(4*(x**2 + 1)**2) - sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(x**2 + 1)*elliptic_f(2*atan(x), sympy.S(1)/4)/(2*sqrt(x**4 + x**2 + 1)) + sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(3*x**2 + 3)*elliptic_e(2*atan(x), sympy.S(1)/4)/(4*sqrt(x**4 + x**2 + 1)) + atan(x/sqrt(x**4 + x**2 + 1))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_236():
    f = (x**2 + 1)**3/(x**4 + x**2 + 1)**(sympy.S(3)/2)
    F = -x*(1 - x**2)/(3*sqrt(x**4 + x**2 + 1)) + 2*x*sqrt(x**4 + x**2 + 1)/(3*x**2 + 3) + sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(x**2 + 1)*elliptic_f(2*atan(x), sympy.S(1)/4)/sqrt(x**4 + x**2 + 1) - sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(2*x**2 + 2)*elliptic_e(2*atan(x), sympy.S(1)/4)/(3*sqrt(x**4 + x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_237():
    f = (x**2 + 1)**2/(x**4 + x**2 + 1)**(sympy.S(3)/2)
    F = x*(2*x**2 + 1)/(3*sqrt(x**4 + x**2 + 1)) - 2*x*sqrt(x**4 + x**2 + 1)/(3*x**2 + 3) + sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(2*x**2 + 2)*elliptic_e(2*atan(x), sympy.S(1)/4)/(3*sqrt(x**4 + x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_238():
    f = (x**2 + 1)/(x**4 + x**2 + 1)**(sympy.S(3)/2)
    F = x*(x**2 + 2)/(3*sqrt(x**4 + x**2 + 1)) - x*sqrt(x**4 + x**2 + 1)/(3*x**2 + 3) + sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(x**2 + 1)*elliptic_e(2*atan(x), sympy.S(1)/4)/(3*sqrt(x**4 + x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_239():
    f = 1/((x**2 + 1)*(x**4 + x**2 + 1)**(sympy.S(3)/2))
    F = -x*(2*x**2 + 1)/(3*sqrt(x**4 + x**2 + 1)) + 2*x*sqrt(x**4 + x**2 + 1)/(3*x**2 + 3) - sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(2*x**2 + 2)*elliptic_e(2*atan(x), sympy.S(1)/4)/(3*sqrt(x**4 + x**2 + 1)) + sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(3*x**2 + 3)*elliptic_f(2*atan(x), sympy.S(1)/4)/(4*sqrt(x**4 + x**2 + 1)) + atan(x/sqrt(x**4 + x**2 + 1))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_240():
    f = 1/((x**2 + 1)**2*(x**4 + x**2 + 1)**(sympy.S(3)/2))
    F = -x*(x**2 + 2)/(3*sqrt(x**4 + x**2 + 1)) + x*sqrt(x**4 + x**2 + 1)/(3*x**2 + 3) + sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(x**2 + 1)*elliptic_e(2*atan(x), sympy.S(1)/4)/(6*sqrt(x**4 + x**2 + 1)) + atan(x/sqrt(x**4 + x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_241():
    f = 1/((x**2 + 1)**3*(x**4 + x**2 + 1)**(sympy.S(3)/2))
    F = -x*(1 - x**2)/(3*sqrt(x**4 + x**2 + 1)) - x*sqrt(x**4 + x**2 + 1)/(3*x**2 + 3) + x*sqrt(x**4 + x**2 + 1)/(4*(x**2 + 1)**2) - sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(5*x**2 + 5)*elliptic_f(2*atan(x), sympy.S(1)/4)/(4*sqrt(x**4 + x**2 + 1)) + sqrt((x**4 + x**2 + 1)/(x**2 + 1)**2)*(19*x**2 + 19)*elliptic_e(2*atan(x), sympy.S(1)/4)/(12*sqrt(x**4 + x**2 + 1)) + 3*atan(x/sqrt(x**4 + x**2 + 1))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_242():
    f = (d + e*x**2)**4*(a + b*x**2 + c*x**4)
    F = a*d**4*x + c*e**4*x**13/13 + d**3*x**3*(4*a*e + b*d)/3 + d**2*x**5*(6*a*e**2 + 4*b*d*e + c*d**2)/5 + 2*d*e*x**7*(2*c*d**2 + e*(2*a*e + 3*b*d))/7 + e**3*x**11*(b*e + 4*c*d)/11 + e**2*x**9*(6*c*d**2 + e*(a*e + 4*b*d))/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_243():
    f = (d + e*x**2)**3*(a + b*x**2 + c*x**4)
    F = a*d**3*x + c*e**3*x**11/11 + d**2*x**3*(3*a*e + b*d)/3 + d*x**5*(c*d**2 + 3*e*(a*e + b*d))/5 + e**2*x**9*(b*e + 3*c*d)/9 + e*x**7*(3*c*d**2 + e*(a*e + 3*b*d))/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_244():
    f = (d + e*x**2)**2*(a + b*x**2 + c*x**4)
    F = a*d**2*x + c*e**2*x**9/9 + d*x**3*(2*a*e + b*d)/3 + e*x**7*(b*e + 2*c*d)/7 + x**5*(c*d**2 + e*(a*e + 2*b*d))/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_245():
    f = (d + e*x**2)*(a + b*x**2 + c*x**4)
    F = a*d*x + c*e*x**7/7 + x**5*(b*e + c*d)/5 + x**3*(a*e + b*d)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_246():
    f = (a + b*x**2 + c*x**4)/(d + e*x**2)
    F = c*x**3/(3*e) - x*(-b*e + c*d)/e**2 + (a*e**2 - b*d*e + c*d**2)*atan(sqrt(e)*x/sqrt(d))/(sqrt(d)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_247():
    f = (a + b*x**2 + c*x**4)/(d + e*x**2)**2
    F = c*x/e**2 + x*(a + d*(-b*e + c*d)/e**2)/(2*d*(d + e*x**2)) - (3*c*d**2 - e*(a*e + b*d))*atan(sqrt(e)*x/sqrt(d))/(2*d**(sympy.S(3)/2)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_248():
    f = (a + b*x**2 + c*x**4)/(d + e*x**2)**3
    F = x*(a + d*(-b*e + c*d)/e**2)/(4*d*(d + e*x**2)**2) - x*(5*c*d**2 - e*(3*a*e + b*d))/(8*d**2*e**2*(d + e*x**2)) + (3*c*d**2 + e*(3*a*e + b*d))*atan(sqrt(e)*x/sqrt(d))/(8*d**(sympy.S(5)/2)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_249():
    f = (a + b*x**2 + c*x**4)/(d + e*x**2)**4
    F = x*(a + d*(-b*e + c*d)/e**2)/(6*d*(d + e*x**2)**3) - x*(7*c*d**2 - e*(5*a*e + b*d))/(24*d**2*e**2*(d + e*x**2)**2) + x*(c*d**2 + e*(5*a*e + b*d))/(16*d**3*e**2*(d + e*x**2)) + (c*d**2 + e*(5*a*e + b*d))*atan(sqrt(e)*x/sqrt(d))/(16*d**(sympy.S(7)/2)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_250():
    f = (d + e*x**2)**3*(a + b*x**2 + c*x**4)**2
    F = a**2*d**3*x + a*d**2*x**3*(3*a*e + 2*b*d)/3 + c**2*e**3*x**15/15 + c*e**2*x**13*(2*b*e + 3*c*d)/13 + d*x**5*(6*a*b*d*e + a*(3*a*e**2 + 2*c*d**2) + b**2*d**2)/5 + e*x**11*(b**2*e**2 + 3*c**2*d**2 + 2*c*e*(a*e + 3*b*d))/11 + x**9*(b*e**2*(2*a*e + 3*b*d)/9 + c**2*d**3/9 + 2*c*d*e*(a*e + b*d)/3) + x**7*(a**2*e**3/7 + 6*a*b*d*e**2/7 + 6*a*c*d**2*e/7 + 3*b**2*d**2*e/7 + 2*b*c*d**3/7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_251():
    f = (d + e*x**2)**2*(a + b*x**2 + c*x**4)**2
    F = a**2*d**2*x + 2*a*d*x**3*(a*e + b*d)/3 + c**2*e**2*x**13/13 + 2*c*e*x**11*(b*e + c*d)/11 + x**9*(b**2*e**2/9 + c**2*d**2/9 + 2*c*e*(a*e + 2*b*d)/9) + x**7*(2*a*b*e**2/7 + 4*a*c*d*e/7 + 2*b**2*d*e/7 + 2*b*c*d**2/7) + x**5*(4*a*b*d*e/5 + a*(a*e**2 + 2*c*d**2)/5 + b**2*d**2/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_252():
    f = (d + e*x**2)*(a + b*x**2 + c*x**4)**2
    F = a**2*d*x + a*x**3*(a*e + 2*b*d)/3 + c**2*e*x**11/11 + c*x**9*(2*b*e + c*d)/9 + x**7*(2*a*c*e/7 + b**2*e/7 + 2*b*c*d/7) + x**5*(2*a*b*e/5 + 2*a*c*d/5 + b**2*d/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_253():
    f = (a + b*x**2 + c*x**4)**2
    F = a**2*x + 2*a*b*x**3/3 + 2*b*c*x**7/7 + c**2*x**9/9 + x**5*(2*a*c/5 + b**2/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_254():
    f = (a + b*x**2 + c*x**4)**2/(d + e*x**2)
    F = c**2*x**7/(7*e) - c*x**5*(-2*b*e + c*d)/(5*e**2) + x**3*(b**2*e**2 + c**2*d**2 - 2*c*e*(-a*e + b*d))/(3*e**3) - x*(-b*e + c*d)*(c*d**2 - e*(-2*a*e + b*d))/e**4 + (a*e**2 - b*d*e + c*d**2)**2*atan(sqrt(e)*x/sqrt(d))/(sqrt(d)*e**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_255():
    f = (a + b*x**2 + c*x**4)**2/(d + e*x**2)**2
    F = c**2*x**5/(5*e**2) - 2*c*x**3*(-b*e + c*d)/(3*e**3) + x*(b**2*e**2 + 3*c**2*d**2 - 2*c*e*(-a*e + 2*b*d))/e**4 + x*(a*e**2 - b*d*e + c*d**2)**2/(2*d*e**4*(d + e*x**2)) - (7*c*d**2 - e*(a*e + 3*b*d))*(a*e**2 - b*d*e + c*d**2)*atan(sqrt(e)*x/sqrt(d))/(2*d**(sympy.S(3)/2)*e**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_256():
    f = (a + b*x**2 + c*x**4)**2/(d + e*x**2)**3
    F = c**2*x**3/(3*e**3) - c*x*(-2*b*e + 3*c*d)/e**4 + x*(a*e**2 - b*d*e + c*d**2)**2/(4*d*e**4*(d + e*x**2)**2) - x*(-3*a*e**2 - 5*b*d*e + 13*c*d**2)*(a*e**2 - b*d*e + c*d**2)/(8*d**2*e**4*(d + e*x**2)) + (35*c**2*d**4 - 6*c*d**2*e*(-a*e + 5*b*d) + e**2*(3*a**2*e**2 + 2*a*b*d*e + 3*b**2*d**2))*atan(sqrt(e)*x/sqrt(d))/(8*d**(sympy.S(5)/2)*e**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_257():
    f = (a + b*x**2 + c*x**4)**2/(d + e*x**2)**4
    F = c**2*x/e**4 + x*(a*e**2 - b*d*e + c*d**2)**2/(6*d*e**4*(d + e*x**2)**3) - x*(-5*a*e**2 - 7*b*d*e + 19*c*d**2)*(a*e**2 - b*d*e + c*d**2)/(24*d**2*e**4*(d + e*x**2)**2) + x*(29*c**2*d**4 - 2*c*d**2*e*(-a*e + 11*b*d) + e**2*(5*a**2*e**2 + 2*a*b*d*e + b**2*d**2))/(16*d**3*e**4*(d + e*x**2)) - (35*c**2*d**4 - 2*c*d**2*e*(a*e + 5*b*d) - e**2*(5*a**2*e**2 + 2*a*b*d*e + b**2*d**2))*atan(sqrt(e)*x/sqrt(d))/(16*d**(sympy.S(7)/2)*e**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_258():
    f = (a + b*x**2 + c*x**4)**2/(d + e*x**2)**5
    F = x*(a*e**2 - b*d*e + c*d**2)**2/(8*d*e**4*(d + e*x**2)**4) - x*(-7*a*e**2 - 9*b*d*e + 25*c*d**2)*(a*e**2 - b*d*e + c*d**2)/(48*d**2*e**4*(d + e*x**2)**3) + x*(163*c**2*d**4 - 2*c*d**2*e*(-3*a*e + 59*b*d) + e**2*(35*a**2*e**2 + 10*a*b*d*e + 3*b**2*d**2))/(192*d**3*e**4*(d + e*x**2)**2) - x*(93*c**2*d**4 - 2*c*d**2*e*(3*a*e + 5*b*d) - e**2*(35*a**2*e**2 + 10*a*b*d*e + 3*b**2*d**2))/(128*d**4*e**4*(d + e*x**2)) + (35*c**2*d**4 + 2*c*d**2*e*(3*a*e + 5*b*d) + e**2*(35*a**2*e**2 + 10*a*b*d*e + 3*b**2*d**2))*atan(sqrt(e)*x/sqrt(d))/(128*d**(sympy.S(9)/2)*e**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_259():
    f = (a + b*x**2 + c*x**4)/(d + e*x**2)**2
    F = c*x/e**2 + x*(a + d*(-b*e + c*d)/e**2)/(2*d*(d + e*x**2)) - (3*c*d**2 - e*(a*e + b*d))*atan(sqrt(e)*x/sqrt(d))/(2*d**(sympy.S(3)/2)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_260():
    f = (a + x**2*(b + c*x**2))/(d + e*x**2)**2
    F = c*x/e**2 + x*(a + d*(-b*e + c*d)/e**2)/(2*d*(d + e*x**2)) - (3*c*d**2 - e*(a*e + b*d))*atan(sqrt(e)*x/sqrt(d))/(2*d**(sympy.S(3)/2)*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_261():
    f = (d + e*x**2)**4/(a + b*x**2 + c*x**4)
    F = e**4*x**5/(5*c) + e**3*x**3*(-b*e + 4*c*d)/(3*c**2) + e**2*x*(b**2*e**2 + 6*c**2*d**2 - c*e*(a*e + 4*b*d))/c**3 + sqrt(2)*(e*(-b*e + 2*c*d)*(b**2*e**2 + 2*c**2*d**2 - 2*c*e*(a*e + b*d)) - (b**4*e**4 - 4*b**2*c*e**3*(a*e + b*d) + 2*c**4*d**4 - 4*c**3*d**2*e*(3*a*e + b*d) + 2*c**2*e**2*(a**2*e**2 + 6*a*b*d*e + 3*b**2*d**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(7)/2)*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*(e*(-b*e + 2*c*d)*(b**2*e**2 + 2*c**2*d**2 - 2*c*e*(a*e + b*d)) + (b**4*e**4 - 4*b**2*c*e**3*(a*e + b*d) + 2*c**4*d**4 - 4*c**3*d**2*e*(3*a*e + b*d) + 2*c**2*e**2*(a**2*e**2 + 6*a*b*d*e + 3*b**2*d**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(7)/2)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_262():
    f = (d + e*x**2)**3/(a + b*x**2 + c*x**4)
    F = e**3*x**3/(3*c) + e**2*x*(-b*e + 3*c*d)/c**2 + sqrt(2)*(e*(b**2*e**2 + 3*c**2*d**2 - c*e*(a*e + 3*b*d)) - (-b*e + 2*c*d)*(b**2*e**2 + c**2*d**2 - c*e*(3*a*e + b*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(5)/2)*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*(e*(b**2*e**2 + 3*c**2*d**2 - c*e*(a*e + 3*b*d)) + (-b*e + 2*c*d)*(b**2*e**2 + c**2*d**2 - c*e*(3*a*e + b*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(5)/2)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_263():
    f = (d + e*x**2)**2/(a + b*x**2 + c*x**4)
    F = e**2*x/c + sqrt(2)*(e*(-b*e + 2*c*d) - (b**2*e**2 + 2*c**2*d**2 - 2*c*e*(a*e + b*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*(e*(-b*e + 2*c*d) + (b**2*e**2 + 2*c**2*d**2 - 2*c*e*(a*e + b*d))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_264():
    f = (d + e*x**2)/(a + b*x**2 + c*x**4)
    F = sqrt(2)*(e - (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*(e + (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_265():
    f = 1/(a + b*x**2 + c*x**4)
    F = -sqrt(2)*sqrt(c)*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)) + sqrt(2)*sqrt(c)*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_266():
    f = 1/((d + e*x**2)*(a + b*x**2 + c*x**4))
    F = -sqrt(2)*sqrt(c)*(e + (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(b + sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)) - sqrt(2)*sqrt(c)*(e - (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(b - sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)) + e**(sympy.S(3)/2)*atan(sqrt(e)*x/sqrt(d))/(sqrt(d)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_267():
    f = 1/((d + e*x**2)**2*(a + b*x**2 + c*x**4))
    F = -sqrt(2)*sqrt(c)*(b*e**2*(b - sqrt(-4*a*c + b**2)) + 2*c**2*d**2 - 2*c*e*(a*e + b*d - d*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(b + sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2)**2) + sqrt(2)*sqrt(c)*(b*e**2*(b + sqrt(-4*a*c + b**2)) + 2*c**2*d**2 - 2*c*e*(a*e + b*d + d*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(b - sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2)**2) + e**2*x/(2*d*(d + e*x**2)*(a*e**2 - b*d*e + c*d**2)) + e**(sympy.S(3)/2)*(-b*e + 2*c*d)*atan(sqrt(e)*x/sqrt(d))/(sqrt(d)*(a*e**2 - b*d*e + c*d**2)**2) + e**(sympy.S(3)/2)*atan(sqrt(e)*x/sqrt(d))/(2*d**(sympy.S(3)/2)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_268():
    f = (d + e*x**2)**3/(a + b*x**2 + c*x**4)**2
    F = x*(c*(-a*b*e*(a*e**2 + 3*c*d**2)/c - 2*a*d*(-3*a*e**2 + c*d**2) + b**2*d**3) - x**2*(a*b**2*e**3 + 2*a*c*e*(-a*e**2 + 3*c*d**2) - b*c*d*(3*a*e**2 + c*d**2)))/(2*a*c*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*(a*b**3*e**3 + 6*a*c*(a*e**2 + c*d**2)*(2*c*d - e*sqrt(-4*a*c + b**2)) - b**2*(-3*a*c*d*e**2 - a*e**3*sqrt(-4*a*c + b**2) + c**2*d**3) + b*c*(a*e**2*(-8*a*e + 3*d*sqrt(-4*a*c + b**2)) + c*d**2*(-12*a*e + d*sqrt(-4*a*c + b**2))))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*c**(sympy.S(3)/2)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) - sqrt(2)*(a*b**3*e**3 + 6*a*c*(a*e**2 + c*d**2)*(2*c*d + e*sqrt(-4*a*c + b**2)) - b**2*(-3*a*c*d*e**2 + a*e**3*sqrt(-4*a*c + b**2) + c**2*d**3) - b*c*(a*e**2*(8*a*e + 3*d*sqrt(-4*a*c + b**2)) + c*d**2*(12*a*e + d*sqrt(-4*a*c + b**2))))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*c**(sympy.S(3)/2)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_269():
    f = (d + e*x**2)**2/(a + b*x**2 + c*x**4)**2
    F = x*(-2*a*b*d*e - 2*a*(-a*e**2 + c*d**2) + b**2*d**2 + x**2*(a*b*e**2 - 4*a*c*d*e + b*c*d**2))/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)) + sqrt(2)*(a*b*e**2 - 4*a*c*d*e + b*c*d**2 - (8*a*b*c*d*e - 4*a*c*(a*e**2 + 3*c*d**2) + b**2*(-a*e**2 + c*d**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*(a*b*e**2 - 4*a*c*d*e + b*c*d**2 + (8*a*b*c*d*e - 4*a*c*(a*e**2 + 3*c*d**2) + b**2*(-a*e**2 + c*d**2))/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_270():
    f = (d + e*x**2)/(a + b*x**2 + c*x**4)**2
    F = sqrt(2)*sqrt(c)*(-2*a*e + b*d - (4*a*b*e - 12*a*c*d + b**2*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + sqrt(2)*sqrt(c)*(-2*a*e + b*d + (4*a*b*e - 12*a*c*d + b**2*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)) + x*(-a*b*e - 2*a*c*d + b**2*d + c*x**2*(-2*a*e + b*d))/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_271():
    f = (a + b*x**2 + c*x**4)**(-2)
    F = -sqrt(2)*sqrt(c)*(-12*a*c + b**2 - b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + sqrt(2)*sqrt(c)*(-12*a*c + b**2 + b*sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)) + x*(-2*a*c + b**2 + b*c*x**2)/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_272():
    f = 1/((d + e*x**2)*(a + b*x**2 + c*x**4)**2)
    F = -sqrt(2)*sqrt(c)*e**2*(e + (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(b + sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)**2) - sqrt(2)*sqrt(c)*e**2*(e - (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(b - sqrt(-4*a*c + b**2))*(a*e**2 - b*d*e + c*d**2)**2) + e**(sympy.S(7)/2)*atan(sqrt(e)*x/sqrt(d))/(sqrt(d)*(a*e**2 - b*d*e + c*d**2)**2) + sqrt(2)*sqrt(c)*(2*a*c*e - b**2*e + b*c*d - (8*a*b*c*e - 12*a*c**2*d - b**3*e + b**2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2)) + sqrt(2)*sqrt(c)*(2*a*c*e - b**2*e + b*c*d + (8*a*b*c*e - 12*a*c**2*d - b**3*e + b**2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2)) + x*(3*a*b*c*e - 2*a*c**2*d - b**3*e + b**2*c*d + c*x**2*(2*a*c*e - b**2*e + b*c*d))/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)*(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_273():
    f = 1/((d + e*x**2)**2*(a + b*x**2 + c*x**4)**2)
    F = -sqrt(2)*sqrt(c)*e**2*(b*e**2*(b - sqrt(-4*a*c + b**2)) + 3*c**2*d**2 - c*e*(a*e + 3*b*d - 2*d*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(sqrt(b + sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2)**3) + sqrt(2)*sqrt(c)*e**2*(b*e**2*(b + sqrt(-4*a*c + b**2)) + 3*c**2*d**2 - c*e*(a*e + 3*b*d + 2*d*sqrt(-4*a*c + b**2)))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(sqrt(b - sqrt(-4*a*c + b**2))*sqrt(-4*a*c + b**2)*(a*e**2 - b*d*e + c*d**2)**3) + e**4*x/(2*d*(d + e*x**2)*(a*e**2 - b*d*e + c*d**2)**2) + 2*e**(sympy.S(7)/2)*(-b*e + 2*c*d)*atan(sqrt(e)*x/sqrt(d))/(sqrt(d)*(a*e**2 - b*d*e + c*d**2)**3) + e**(sympy.S(7)/2)*atan(sqrt(e)*x/sqrt(d))/(2*d**(sympy.S(3)/2)*(a*e**2 - b*d*e + c*d**2)**2) - sqrt(2)*sqrt(c)*(-4*a*c**2*(3*c*d**2 + e*(-3*a*e + d*sqrt(-4*a*c + b**2))) + b**4*e**2 - b**3*e*(2*c*d + e*sqrt(-4*a*c + b**2)) + b**2*c*(c*d**2 + e*(-9*a*e + 2*d*sqrt(-4*a*c + b**2))) + b*c*(3*a*e**2*sqrt(-4*a*c + b**2) - c*d*(-16*a*e + d*sqrt(-4*a*c + b**2))))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(4*a*sqrt(b + sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)*(a*e**2 - b*d*e + c*d**2)**2) + sqrt(2)*sqrt(c)*(-4*a*c**2*(3*c*d**2 - e*(3*a*e + d*sqrt(-4*a*c + b**2))) + b**4*e**2 - b**3*e*(2*c*d - e*sqrt(-4*a*c + b**2)) + b**2*c*(c*d**2 - e*(9*a*e + 2*d*sqrt(-4*a*c + b**2))) - b*c*(3*a*e**2*sqrt(-4*a*c + b**2) - c*d*(16*a*e + d*sqrt(-4*a*c + b**2))))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(4*a*sqrt(b - sqrt(-4*a*c + b**2))*(-4*a*c + b**2)**(sympy.S(3)/2)*(a*e**2 - b*d*e + c*d**2)**2) + x*(a*b*c*e*(-b*e + 2*c*d) - c*x**2*(-4*a*c**2*d*e - b**3*e**2 + 2*b**2*c*d*e - b*c*(-3*a*e**2 + c*d**2)) + (-2*a*c + b**2)*(b**2*e**2 + c**2*d**2 - c*e*(a*e + 2*b*d)))/(2*a*(-4*a*c + b**2)*(a + b*x**2 + c*x**4)*(a*e**2 - b*d*e + c*d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_274():
    f = (d + e*x**2)**(sympy.S(5)/2)*(a + b*x**2 + c*x**4)
    F = c*x**3*(d + e*x**2)**(sympy.S(7)/2)/(10*e) + d**3*(80*a*e**2 - 10*b*d*e + 3*c*d**2)*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(256*e**(sympy.S(5)/2)) + d**2*x*sqrt(d + e*x**2)*(80*a*e**2 - 10*b*d*e + 3*c*d**2)/(256*e**2) + d*x*(d + e*x**2)**(sympy.S(3)/2)*(80*a*e**2 - 10*b*d*e + 3*c*d**2)/(384*e**2) - x*(d + e*x**2)**(sympy.S(7)/2)*(-10*b*e + 3*c*d)/(80*e**2) + x*(d + e*x**2)**(sympy.S(5)/2)*(80*a*e**2 - 10*b*d*e + 3*c*d**2)/(480*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_275():
    f = (d + e*x**2)**(sympy.S(3)/2)*(a + b*x**2 + c*x**4)
    F = c*x**3*(d + e*x**2)**(sympy.S(5)/2)/(8*e) + d**2*(48*a*e**2 - 8*b*d*e + 3*c*d**2)*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(128*e**(sympy.S(5)/2)) + d*x*sqrt(d + e*x**2)*(48*a*e**2 - 8*b*d*e + 3*c*d**2)/(128*e**2) - x*(d + e*x**2)**(sympy.S(5)/2)*(-8*b*e + 3*c*d)/(48*e**2) + x*(d + e*x**2)**(sympy.S(3)/2)*(48*a*e**2 - 8*b*d*e + 3*c*d**2)/(192*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_276():
    f = sqrt(d + e*x**2)*(a + b*x**2 + c*x**4)
    F = c*x**3*(d + e*x**2)**(sympy.S(3)/2)/(6*e) + d*(8*a*e**2 - 2*b*d*e + c*d**2)*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(16*e**(sympy.S(5)/2)) - x*(d + e*x**2)**(sympy.S(3)/2)*(-2*b*e + c*d)/(8*e**2) + x*sqrt(d + e*x**2)*(8*a*e**2 - 2*b*d*e + c*d**2)/(16*e**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_277():
    f = (a + b*x**2 + c*x**4)/sqrt(d + e*x**2)
    F = c*x**3*sqrt(d + e*x**2)/(4*e) - x*sqrt(d + e*x**2)*(-4*b*e + 3*c*d)/(8*e**2) + (8*a*e**2 - 4*b*d*e + 3*c*d**2)*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(8*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_278():
    f = (a + b*x**2 + c*x**4)/(d + e*x**2)**(sympy.S(3)/2)
    F = c*x*sqrt(d + e*x**2)/(2*e**2) - (-2*b*e + 3*c*d)*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(2*e**(sympy.S(5)/2)) + x*(a + d*(-b*e + c*d)/e**2)/(d*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_279():
    f = (a + b*x**2 + c*x**4)/(d + e*x**2)**(sympy.S(5)/2)
    F = c*atanh(sqrt(e)*x/sqrt(d + e*x**2))/e**(sympy.S(5)/2) + x*(a + d*(-b*e + c*d)/e**2)/(3*d*(d + e*x**2)**(sympy.S(3)/2)) - x*(4*c*d**2 - e*(2*a*e + b*d))/(3*d**2*e**2*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_280():
    f = (a + b*x**2 + c*x**4)/(d + e*x**2)**(sympy.S(7)/2)
    F = a*x/(d*(d + e*x**2)**(sympy.S(5)/2)) + x**3*(4*a*e + b*d)/(3*d**2*(d + e*x**2)**(sympy.S(5)/2)) + x**5*(3*c*d**2 + 2*e*(4*a*e + b*d))/(15*d**3*(d + e*x**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_281():
    f = (a + b*x**2 + c*x**4)/(d + e*x**2)**(sympy.S(9)/2)
    F = a*x/(d*(d + e*x**2)**(sympy.S(7)/2)) + x**3*(6*a*e + b*d)/(3*d**2*(d + e*x**2)**(sympy.S(7)/2)) + x**5*(3*c*d**2 + 4*e*(6*a*e + b*d))/(15*d**3*(d + e*x**2)**(sympy.S(7)/2)) + 2*e*x**7*(3*c*d**2 + 4*e*(6*a*e + b*d))/(105*d**4*(d + e*x**2)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_282():
    f = (a + b*x**2 + c*x**4)/(d + e*x**2)**(sympy.S(11)/2)
    F = a*x/(d*(d + e*x**2)**(sympy.S(9)/2)) + x**3*(8*a*e + b*d)/(3*d**2*(d + e*x**2)**(sympy.S(9)/2)) + x**5*(c*d**2 + 2*e*(8*a*e + b*d))/(5*d**3*(d + e*x**2)**(sympy.S(9)/2)) + 4*e*x**7*(c*d**2 + 2*e*(8*a*e + b*d))/(35*d**4*(d + e*x**2)**(sympy.S(9)/2)) + 8*e**2*x**9*(c*d**2 + 2*e*(8*a*e + b*d))/(315*d**5*(d + e*x**2)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_283():
    f = (a + b*x**2 + c*x**4)/(d + e*x**2)**(sympy.S(13)/2)
    F = a*x/(d*(d + e*x**2)**(sympy.S(11)/2)) + x**3*(10*a*e + b*d)/(3*d**2*(d + e*x**2)**(sympy.S(11)/2)) + x**5*(3*c*d**2 + 8*e*(10*a*e + b*d))/(15*d**3*(d + e*x**2)**(sympy.S(11)/2)) + 2*e*x**7*(3*c*d**2 + 8*e*(10*a*e + b*d))/(35*d**4*(d + e*x**2)**(sympy.S(11)/2)) + 8*e**2*x**9*(3*c*d**2 + 8*e*(10*a*e + b*d))/(315*d**5*(d + e*x**2)**(sympy.S(11)/2)) + 16*e**3*x**11*(3*c*d**2 + 8*e*(10*a*e + b*d))/(3465*d**6*(d + e*x**2)**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_284():
    f = (5*x**2 + 7)**3*sqrt(x**4 + 3*x**2 + 2)
    F = 125*x**3*(x**4 + 3*x**2 + 2)**(sympy.S(3)/2)/9 + 577*x*(x**2 + 2)/(3*sqrt(x**4 + 3*x**2 + 2)) + x*(757*x**2 + 2608)*sqrt(x**4 + 3*x**2 + 2)/21 + 275*x*(x**4 + 3*x**2 + 2)**(sympy.S(3)/2)/7 - 577*sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_e(atan(x), sympy.S.Half)/(3*sqrt(x**4 + 3*x**2 + 2)) + 2945*sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_f(atan(x), sympy.S.Half)/(21*sqrt(x**4 + 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_285():
    f = (5*x**2 + 7)**2*sqrt(x**4 + 3*x**2 + 2)
    F = 31*x*(x**2 + 2)/sqrt(x**4 + 3*x**2 + 2) + x*(114*x**2 + 407)*sqrt(x**4 + 3*x**2 + 2)/21 + 25*x*(x**4 + 3*x**2 + 2)**(sympy.S(3)/2)/7 - 31*sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_e(atan(x), sympy.S.Half)/sqrt(x**4 + 3*x**2 + 2) + 472*sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_f(atan(x), sympy.S.Half)/(21*sqrt(x**4 + 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_286():
    f = (5*x**2 + 7)*sqrt(x**4 + 3*x**2 + 2)
    F = 5*x*(x**2 + 2)/sqrt(x**4 + 3*x**2 + 2) + x*(3*x**2 + 10)*sqrt(x**4 + 3*x**2 + 2)/3 - 5*sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_e(atan(x), sympy.S.Half)/sqrt(x**4 + 3*x**2 + 2) + 11*sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_f(atan(x), sympy.S.Half)/(3*sqrt(x**4 + 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_287():
    f = sqrt(x**4 + 3*x**2 + 2)
    F = x*(x**2 + 2)/sqrt(x**4 + 3*x**2 + 2) + x*sqrt(x**4 + 3*x**2 + 2)/3 - sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_e(atan(x), sympy.S.Half)/sqrt(x**4 + 3*x**2 + 2) + 2*sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_f(atan(x), sympy.S.Half)/(3*sqrt(x**4 + 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_288():
    f = sqrt(x**4 + 3*x**2 + 2)/(5*x**2 + 7)**2
    F = (Integer(-1) * ((x * (Integer(2) + (x)**(Integer(2)))) * ((Integer(70) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + ((x * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(14) * (Integer(7) + (Integer(5) * (x)**(Integer(2))))))**(Integer(-1))) + (((Integer(1) + (x)**(Integer(2))) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.elliptic_e(sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(35) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + ((Integer(3) * (Integer(1) + (x)**(Integer(2))) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.elliptic_f(sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(140) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * (((Integer(2) + (x)**(Integer(2))) * sympy.Function('EllipticPi')((Integer(2) * (Integer(7))**(Integer(-1))), sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(980) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_289():
    f = sqrt(x**4 + 3*x**2 + 2)/(5*x**2 + 7)**3
    F = (Integer(-1) * ((Integer(11) * x * (Integer(2) + (x)**(Integer(2)))) * ((Integer(11760) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + ((x * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(28) * ((Integer(7) + (Integer(5) * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))) + ((Integer(11) * x * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(2352) * (Integer(7) + (Integer(5) * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(11) * (Integer(1) + (x)**(Integer(2))) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.elliptic_e(sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(5880) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + ((Integer(81) * (Integer(1) + (x)**(Integer(2))) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.elliptic_f(sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(7840) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * ((Integer(1201) * (Integer(2) + (x)**(Integer(2))) * sympy.Function('EllipticPi')((Integer(2) * (Integer(7))**(Integer(-1))), sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(164640) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_290():
    f = (5*x**2 + 7)**3*(x**4 + 3*x**2 + 2)**(sympy.S(3)/2)
    F = 125*x**3*(x**4 + 3*x**2 + 2)**(sympy.S(5)/2)/13 + 20884*x*(x**2 + 2)/(65*sqrt(x**4 + 3*x**2 + 2)) + x*(65345*x**2 + 208212)*(x**4 + 3*x**2 + 2)**(sympy.S(3)/2)/3003 + x*(297911*x**2 + 1032541)*sqrt(x**4 + 3*x**2 + 2)/5005 + 3825*x*(x**4 + 3*x**2 + 2)**(sympy.S(5)/2)/143 - 20884*sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_e(atan(x), sympy.S.Half)/(65*sqrt(x**4 + 3*x**2 + 2)) + 1171349*sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_f(atan(x), sympy.S.Half)/(5005*sqrt(x**4 + 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_291():
    f = (5*x**2 + 7)**2*(x**4 + 3*x**2 + 2)**(sympy.S(3)/2)
    F = 742*x*(x**2 + 2)/(15*sqrt(x**4 + 3*x**2 + 2)) + x*(2240*x**2 + 7281)*(x**4 + 3*x**2 + 2)**(sympy.S(3)/2)/693 + x*(10643*x**2 + 36783)*sqrt(x**4 + 3*x**2 + 2)/1155 + 25*x*(x**4 + 3*x**2 + 2)**(sympy.S(5)/2)/11 - 742*sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_e(atan(x), sympy.S.Half)/(15*sqrt(x**4 + 3*x**2 + 2)) + 13879*sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_f(atan(x), sympy.S.Half)/(385*sqrt(x**4 + 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_292():
    f = (5*x**2 + 7)*(x**4 + 3*x**2 + 2)**(sympy.S(3)/2)
    F = 116*x*(x**2 + 2)/(15*sqrt(x**4 + 3*x**2 + 2)) + x*(35*x**2 + 108)*(x**4 + 3*x**2 + 2)**(sympy.S(3)/2)/63 + x*(149*x**2 + 519)*sqrt(x**4 + 3*x**2 + 2)/105 - 116*sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_e(atan(x), sympy.S.Half)/(15*sqrt(x**4 + 3*x**2 + 2)) + 197*sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_f(atan(x), sympy.S.Half)/(35*sqrt(x**4 + 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_293():
    f = (x**4 + 3*x**2 + 2)**(sympy.S(3)/2)
    F = 6*x*(x**2 + 2)/(5*sqrt(x**4 + 3*x**2 + 2)) + x*(9*x**2 + 29)*sqrt(x**4 + 3*x**2 + 2)/35 + x*(x**4 + 3*x**2 + 2)**(sympy.S(3)/2)/7 - 6*sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_e(atan(x), sympy.S.Half)/(5*sqrt(x**4 + 3*x**2 + 2)) + 31*sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_f(atan(x), sympy.S.Half)/(35*sqrt(x**4 + 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_294():
    f = (x**4 + 3*x**2 + 2)**(sympy.S(3)/2)/(5*x**2 + 7)
    F = ((Integer(24) * x * (Integer(2) + (x)**(Integer(2)))) * ((Integer(125) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + ((Integer(75))**(Integer(-1)) * x * (Integer(11) + (Integer(3) * (x)**(Integer(2)))) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) + (Integer(-1) * ((Integer(24) * sympy.sqrt(Integer(2)) * (Integer(1) + (x)**(Integer(2))) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.elliptic_e(sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(125) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + ((Integer(56) * sympy.sqrt(Integer(2)) * (Integer(1) + (x)**(Integer(2))) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.elliptic_f(sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(375) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * sympy.sqrt(Integer(2)) * (Integer(2) + (x)**(Integer(2))) * sympy.Function('EllipticPi')((Integer(2) * (Integer(7))**(Integer(-1))), sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(875) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_295():
    f = (5*x**2 + 7)**3/sqrt(x**4 + 3*x**2 + 2)
    F = 25*x**3*sqrt(x**4 + 3*x**2 + 2) + 135*x*(x**2 + 2)/sqrt(x**4 + 3*x**2 + 2) + 75*x*sqrt(x**4 + 3*x**2 + 2) - 135*sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_e(atan(x), sympy.S.Half)/sqrt(x**4 + 3*x**2 + 2) + sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(193*x**2 + 193)*elliptic_f(atan(x), sympy.S.Half)/(2*sqrt(x**4 + 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_296():
    f = (5*x**2 + 7)**2/sqrt(x**4 + 3*x**2 + 2)
    F = 20*x*(x**2 + 2)/sqrt(x**4 + 3*x**2 + 2) + 25*x*sqrt(x**4 + 3*x**2 + 2)/3 - 20*sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_e(atan(x), sympy.S.Half)/sqrt(x**4 + 3*x**2 + 2) + sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(97*x**2 + 97)*elliptic_f(atan(x), sympy.S.Half)/(6*sqrt(x**4 + 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_297():
    f = (5*x**2 + 7)/sqrt(x**4 + 3*x**2 + 2)
    F = 5*x*(x**2 + 2)/sqrt(x**4 + 3*x**2 + 2) - 5*sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_e(atan(x), sympy.S.Half)/sqrt(x**4 + 3*x**2 + 2) + sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(7*x**2 + 7)*elliptic_f(atan(x), sympy.S.Half)/(2*sqrt(x**4 + 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_298():
    f = 1/sqrt(x**4 + 3*x**2 + 2)
    F = sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_f(atan(x), sympy.S.Half)/(2*sqrt(x**4 + 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_299():
    f = 1/((5*x**2 + 7)*sqrt(x**4 + 3*x**2 + 2))
    F = (((Integer(1) + (x)**(Integer(2))) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.elliptic_f(sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Integer(2) + (x)**(Integer(2))) * sympy.Function('EllipticPi')((Integer(2) * (Integer(7))**(Integer(-1))), sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(14) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_300():
    f = 1/((5*x**2 + 7)**2*sqrt(x**4 + 3*x**2 + 2))
    F = ((Integer(5) * x * (Integer(2) + (x)**(Integer(2)))) * ((Integer(84) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * ((Integer(25) * x * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(84) * (Integer(7) + (Integer(5) * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Integer(1) + (x)**(Integer(2))) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.elliptic_e(sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(42) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + ((Integer(9) * (Integer(1) + (x)**(Integer(2))) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.elliptic_f(sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(56) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * ((Integer(65) * (Integer(2) + (x)**(Integer(2))) * sympy.Function('EllipticPi')((Integer(2) * (Integer(7))**(Integer(-1))), sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(1176) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_301():
    f = 1/((5*x**2 + 7)**3*sqrt(x**4 + 3*x**2 + 2))
    F = ((Integer(65) * x * (Integer(2) + (x)**(Integer(2)))) * ((Integer(4704) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * ((Integer(25) * x * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(168) * ((Integer(7) + (Integer(5) * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(325) * x * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(4704) * (Integer(7) + (Integer(5) * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(65) * (Integer(1) + (x)**(Integer(2))) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.elliptic_e(sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(2352) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + ((Integer(631) * (Integer(1) + (x)**(Integer(2))) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.elliptic_f(sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(9408) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2525) * (Integer(2) + (x)**(Integer(2))) * sympy.Function('EllipticPi')((Integer(2) * (Integer(7))**(Integer(-1))), sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(65856) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_302():
    f = (5*x**2 + 7)**5/(x**4 + 3*x**2 + 2)**(sympy.S(3)/2)
    F = 625*x**3*sqrt(x**4 + 3*x**2 + 2) + 7679*x*(x**2 + 2)/(2*sqrt(x**4 + 3*x**2 + 2)) - x*(179*x**2 + 115)/(2*sqrt(x**4 + 3*x**2 + 2)) + 5000*x*sqrt(x**4 + 3*x**2 + 2)/3 - sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(7679*x**2 + 7679)*elliptic_e(atan(x), sympy.S.Half)/(2*sqrt(x**4 + 3*x**2 + 2)) + sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(15383*x**2 + 15383)*elliptic_f(atan(x), sympy.S.Half)/(6*sqrt(x**4 + 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_303():
    f = (5*x**2 + 7)**4/(x**4 + 3*x**2 + 2)**(sympy.S(3)/2)
    F = 637*x*(x**2 + 2)/(2*sqrt(x**4 + 3*x**2 + 2)) + x*(113*x**2 + 145)/(2*sqrt(x**4 + 3*x**2 + 2)) + 625*x*sqrt(x**4 + 3*x**2 + 2)/3 + 1067*sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_f(atan(x), sympy.S.Half)/(3*sqrt(x**4 + 3*x**2 + 2)) - sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(637*x**2 + 637)*elliptic_e(atan(x), sympy.S.Half)/(2*sqrt(x**4 + 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_304():
    f = (5*x**2 + 7)**3/(x**4 + 3*x**2 + 2)**(sympy.S(3)/2)
    F = x*(5 - 11*x**2)/(2*sqrt(x**4 + 3*x**2 + 2)) + 261*x*(x**2 + 2)/(2*sqrt(x**4 + 3*x**2 + 2)) + sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(169*x**2 + 169)*elliptic_f(atan(x), sympy.S.Half)/(2*sqrt(x**4 + 3*x**2 + 2)) - sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(261*x**2 + 261)*elliptic_e(atan(x), sympy.S.Half)/(2*sqrt(x**4 + 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_305():
    f = (5*x**2 + 7)**2/(x**4 + 3*x**2 + 2)**(sympy.S(3)/2)
    F = -17*x*(x**2 + 2)/(2*sqrt(x**4 + 3*x**2 + 2)) + x*(17*x**2 + 25)/(2*sqrt(x**4 + 3*x**2 + 2)) + 6*sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_f(atan(x), sympy.S.Half)/sqrt(x**4 + 3*x**2 + 2) + sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(17*x**2 + 17)*elliptic_e(atan(x), sympy.S.Half)/(2*sqrt(x**4 + 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_306():
    f = (5*x**2 + 7)/(x**4 + 3*x**2 + 2)**(sympy.S(3)/2)
    F = -x*(x**2 + 2)/(2*sqrt(x**4 + 3*x**2 + 2)) + x*(x**2 + 5)/(2*sqrt(x**4 + 3*x**2 + 2)) + sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_e(atan(x), sympy.S.Half)/(2*sqrt(x**4 + 3*x**2 + 2)) + sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_f(atan(x), sympy.S.Half)/(2*sqrt(x**4 + 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_307():
    f = (x**4 + 3*x**2 + 2)**(sympy.S(-3)/2)
    F = -3*x*(x**2 + 2)/(2*sqrt(x**4 + 3*x**2 + 2)) + x*(3*x**2 + 5)/(2*sqrt(x**4 + 3*x**2 + 2)) - sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_f(atan(x), sympy.S.Half)/sqrt(x**4 + 3*x**2 + 2) + sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(3*x**2 + 3)*elliptic_e(atan(x), sympy.S.Half)/(2*sqrt(x**4 + 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_308():
    f = 1/((5*x**2 + 7)**2*(x**4 + 3*x**2 + 2)**(sympy.S(3)/2))
    F = (Integer(-1) * ((Integer(31) * x * (Integer(2) + (x)**(Integer(2)))) * ((Integer(56) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + ((x * (Integer(20) + (Integer(11) * (x)**(Integer(2))))) * ((Integer(36) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + ((Integer(625) * x * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(504) * (Integer(7) + (Integer(5) * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(31) * (Integer(1) + (x)**(Integer(2))) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.elliptic_e(sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(28) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * ((Integer(463) * (Integer(1) + (x)**(Integer(2))) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.elliptic_f(sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(336) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + ((Integer(375) * (Integer(2) + (x)**(Integer(2))) * sympy.Function('EllipticPi')((Integer(2) * (Integer(7))**(Integer(-1))), sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(784) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_309():
    f = 1/((5*x**2 + 7)**3*(x**4 + 3*x**2 + 2)**(sympy.S(3)/2))
    F = (Integer(-1) * ((Integer(5797) * x * (Integer(2) + (x)**(Integer(2)))) * ((Integer(28224) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + ((x * (Integer(50) + (Integer(23) * (x)**(Integer(2))))) * ((Integer(216) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + ((Integer(625) * x * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(1008) * ((Integer(7) + (Integer(5) * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))) + ((Integer(41875) * x * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(84672) * (Integer(7) + (Integer(5) * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(5797) * (Integer(1) + (x)**(Integer(2))) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.elliptic_e(sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(14112) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * ((Integer(49907) * (Integer(1) + (x)**(Integer(2))) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.elliptic_f(sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(56448) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + ((Integer(192625) * (Integer(2) + (x)**(Integer(2))) * sympy.Function('EllipticPi')((Integer(2) * (Integer(7))**(Integer(-1))), sympy.atan(x), (Integer(2))**(Integer(-1)))) * ((Integer(395136) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Integer(2) + (x)**(Integer(2))) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1)))) * sympy.sqrt((Integer(2) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_310():
    f = (5*x**2 + 7)**4*sqrt(-x**4 + x**2 + 2)
    F = -625*x**5*(-x**4 + x**2 + 2)**(sympy.S(3)/2)/11 - 14500*x**3*(-x**4 + x**2 + 2)**(sympy.S(3)/2)/33 + x*(717372*x**2 + 177953)*sqrt(-x**4 + x**2 + 2)/231 - 116100*x*(-x**4 + x**2 + 2)**(sympy.S(3)/2)/77 + 3764813*elliptic_e(asin(sqrt(2)*x/2), -2)/231 - 539419*elliptic_f(asin(sqrt(2)*x/2), -2)/77
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_311():
    f = (5*x**2 + 7)**3*sqrt(-x**4 + x**2 + 2)
    F = -125*x**3*(-x**4 + x**2 + 2)**(sympy.S(3)/2)/9 + x*(14691*x**2 + 5956)*sqrt(-x**4 + x**2 + 2)/63 - 1825*x*(-x**4 + x**2 + 2)**(sympy.S(3)/2)/21 + 79411*elliptic_e(asin(sqrt(2)*x/2), -2)/63 - 8735*elliptic_f(asin(sqrt(2)*x/2), -2)/21
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_312():
    f = (5*x**2 + 7)**2*sqrt(-x**4 + x**2 + 2)
    F = x*(354*x**2 + 275)*sqrt(-x**4 + x**2 + 2)/21 - 25*x*(-x**4 + x**2 + 2)**(sympy.S(3)/2)/7 + 2045*elliptic_e(asin(sqrt(2)*x/2), -2)/21 - 79*elliptic_f(asin(sqrt(2)*x/2), -2)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_313():
    f = (5*x**2 + 7)*sqrt(-x**4 + x**2 + 2)
    F = x*(x**2 + 2)*sqrt(-x**4 + x**2 + 2) + 7*elliptic_e(asin(sqrt(2)*x/2), -2) + 3*elliptic_f(asin(sqrt(2)*x/2), -2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_314():
    f = sqrt(-x**4 + x**2 + 2)
    F = x*sqrt(-x**4 + x**2 + 2)/3 + elliptic_e(asin(sqrt(2)*x/2), -2)/3 + elliptic_f(asin(sqrt(2)*x/2), -2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_315():
    f = sqrt(-x**4 + x**2 + 2)/(5*x**2 + 7)
    F = ((Integer(-1) * (Integer(5))**(Integer(-1))) * sympy.elliptic_e(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) + ((Integer(17) * (Integer(25))**(Integer(-1))) * sympy.elliptic_f(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) + (Integer(-1) * ((Integer(34) * (Integer(175))**(Integer(-1))) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(10) * (Integer(7))**(Integer(-1)))), sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_316():
    f = sqrt(-x**4 + x**2 + 2)/(5*x**2 + 7)**2
    F = ((x * sympy.sqrt((Integer(2) + (x)**(Integer(2)) + (Integer(-1) * (x)**(Integer(4)))))) * ((Integer(14) * (Integer(7) + (Integer(5) * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(70))**(Integer(-1)) * sympy.elliptic_e(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) + (Integer(-1) * ((Integer(6) * (Integer(175))**(Integer(-1))) * sympy.elliptic_f(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2)))) + ((Integer(99) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(10) * (Integer(7))**(Integer(-1)))), sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) * (Integer(2450))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_317():
    f = sqrt(-x**4 + x**2 + 2)/(5*x**2 + 7)**3
    F = ((x * sympy.sqrt((Integer(2) + (x)**(Integer(2)) + (Integer(-1) * (x)**(Integer(4)))))) * ((Integer(28) * ((Integer(7) + (Integer(5) * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(31) * x * sympy.sqrt((Integer(2) + (x)**(Integer(2)) + (Integer(-1) * (x)**(Integer(4)))))) * ((Integer(13328) * (Integer(7) + (Integer(5) * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(31) * sympy.elliptic_e(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) * (Integer(66640))**(Integer(-1)))) + (Integer(-1) * ((Integer(269) * sympy.elliptic_f(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) * (Integer(166600))**(Integer(-1)))) + ((Integer(16601) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(10) * (Integer(7))**(Integer(-1)))), sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) * (Integer(2332400))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_318():
    f = (5*x**2 + 7)**4*(-x**4 + x**2 + 2)**(sympy.S(3)/2)
    F = -125*x**5*(-x**4 + x**2 + 2)**(sympy.S(5)/2)/3 - 11750*x**3*(-x**4 + x**2 + 2)**(sympy.S(5)/2)/39 - x*(69817 - 1581440*x**2)*(-x**4 + x**2 + 2)**(sympy.S(3)/2)/1001 + 3*x*(7837383*x**2 + 2193559)*sqrt(-x**4 + x**2 + 2)/5005 - 132300*x*(-x**4 + x**2 + 2)**(sympy.S(5)/2)/143 + 124141422*elliptic_e(asin(sqrt(2)*x/2), -2)/5005 - 50794416*elliptic_f(asin(sqrt(2)*x/2), -2)/5005
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_319():
    f = (5*x**2 + 7)**3*(-x**4 + x**2 + 2)**(sympy.S(3)/2)
    F = -125*x**3*(-x**4 + x**2 + 2)**(sympy.S(5)/2)/13 + x*(374045*x**2 + 33792)*(-x**4 + x**2 + 2)**(sympy.S(3)/2)/3003 + x*(5712051*x**2 + 2512273)*sqrt(-x**4 + x**2 + 2)/15015 - 7825*x*(-x**4 + x**2 + 2)**(sympy.S(5)/2)/143 + 31072528*elliptic_e(asin(sqrt(2)*x/2), -2)/15015 - 3199778*elliptic_f(asin(sqrt(2)*x/2), -2)/5005
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_320():
    f = (5*x**2 + 7)**2*(-x**4 + x**2 + 2)**(sympy.S(3)/2)
    F = x*(920*x**2 + 363)*(-x**4 + x**2 + 2)**(sympy.S(3)/2)/99 + x*(14889*x**2 + 11497)*sqrt(-x**4 + x**2 + 2)/495 - 25*x*(-x**4 + x**2 + 2)**(sympy.S(5)/2)/11 + 85942*elliptic_e(asin(sqrt(2)*x/2), -2)/495 - 3392*elliptic_f(asin(sqrt(2)*x/2), -2)/165
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_321():
    f = (5*x**2 + 7)*(-x**4 + x**2 + 2)**(sympy.S(3)/2)
    F = x*(35*x**2 + 48)*(-x**4 + x**2 + 2)**(sympy.S(3)/2)/63 + x*(669*x**2 + 1087)*sqrt(-x**4 + x**2 + 2)/315 + 4432*elliptic_e(asin(sqrt(2)*x/2), -2)/315 + 418*elliptic_f(asin(sqrt(2)*x/2), -2)/105
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_322():
    f = (-x**4 + x**2 + 2)**(sympy.S(3)/2)
    F = x*(3*x**2 + 19)*sqrt(-x**4 + x**2 + 2)/35 + x*(-x**4 + x**2 + 2)**(sympy.S(3)/2)/7 + 34*elliptic_e(asin(sqrt(2)*x/2), -2)/35 + 48*elliptic_f(asin(sqrt(2)*x/2), -2)/35
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_323():
    f = (-x**4 + x**2 + 2)**(sympy.S(3)/2)/(5*x**2 + 7)
    F = ((Integer(75))**(Integer(-1)) * x * (Integer(13) + (Integer(-1) * (Integer(3) * (x)**(Integer(2))))) * sympy.sqrt((Integer(2) + (x)**(Integer(2)) + (Integer(-1) * (x)**(Integer(4)))))) + ((Integer(92) * (Integer(375))**(Integer(-1))) * sympy.elliptic_e(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) + (Integer(-1) * ((Integer(178) * (Integer(625))**(Integer(-1))) * sympy.elliptic_f(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2)))) + ((Integer(1156) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(10) * (Integer(7))**(Integer(-1)))), sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) * (Integer(4375))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_324():
    f = (-x**4 + x**2 + 2)**(sympy.S(3)/2)/(5*x**2 + 7)**2
    F = ((Integer(-1) * (Integer(75))**(Integer(-1))) * x * sympy.sqrt((Integer(2) + (x)**(Integer(2)) + (Integer(-1) * (x)**(Integer(4)))))) + (Integer(-1) * ((Integer(17) * x * sympy.sqrt((Integer(2) + (x)**(Integer(2)) + (Integer(-1) * (x)**(Integer(4)))))) * ((Integer(175) * (Integer(7) + (Integer(5) * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(97) * (Integer(525))**(Integer(-1))) * sympy.elliptic_e(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2)))) + ((Integer(458) * (Integer(875))**(Integer(-1))) * sympy.elliptic_f(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) + (Integer(-1) * ((Integer(1241) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(10) * (Integer(7))**(Integer(-1)))), sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) * (Integer(6125))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_325():
    f = (-x**4 + x**2 + 2)**(sympy.S(3)/2)/(5*x**2 + 7)**3
    F = (Integer(-1) * ((Integer(17) * x * sympy.sqrt((Integer(2) + (x)**(Integer(2)) + (Integer(-1) * (x)**(Integer(4)))))) * ((Integer(350) * ((Integer(7) + (Integer(5) * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1)))) + ((Integer(563) * x * sympy.sqrt((Integer(2) + (x)**(Integer(2)) + (Integer(-1) * (x)**(Integer(4)))))) * ((Integer(9800) * (Integer(7) + (Integer(5) * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(191) * sympy.elliptic_e(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) * (Integer(9800))**(Integer(-1))) + (Integer(-1) * ((Integer(1251) * sympy.elliptic_f(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) * (Integer(24500))**(Integer(-1)))) + ((Integer(9879) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(10) * (Integer(7))**(Integer(-1)))), sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) * (Integer(343000))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_326():
    f = (5*x**2 + 7)**3/sqrt(-x**4 + x**2 + 2)
    F = -25*x**3*sqrt(-x**4 + x**2 + 2) - 625*x*sqrt(-x**4 + x**2 + 2)/3 + 3905*elliptic_e(asin(sqrt(2)*x/2), -2)/3 - 542*elliptic_f(asin(sqrt(2)*x/2), -2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_327():
    f = (5*x**2 + 7)**2/sqrt(-x**4 + x**2 + 2)
    F = -25*x*sqrt(-x**4 + x**2 + 2)/3 + 260*elliptic_e(asin(sqrt(2)*x/2), -2)/3 - 21*elliptic_f(asin(sqrt(2)*x/2), -2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_328():
    f = (5*x**2 + 7)/sqrt(-x**4 + x**2 + 2)
    F = 5*elliptic_e(asin(sqrt(2)*x/2), -2) + 2*elliptic_f(asin(sqrt(2)*x/2), -2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_329():
    f = 1/sqrt(-x**4 + x**2 + 2)
    F = elliptic_f(asin(sqrt(2)*x/2), -2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_330():
    f = 1/((5*x**2 + 7)*sqrt(-x**4 + x**2 + 2))
    F = (Integer(7))**(Integer(-1)) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(10) * (Integer(7))**(Integer(-1)))), sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_331():
    f = 1/((5*x**2 + 7)**2*sqrt(-x**4 + x**2 + 2))
    F = (Integer(-1) * ((Integer(25) * x * sympy.sqrt((Integer(2) + (x)**(Integer(2)) + (Integer(-1) * (x)**(Integer(4)))))) * ((Integer(476) * (Integer(7) + (Integer(5) * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Integer(476))**(Integer(-1))) * sympy.elliptic_e(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2)))) + (Integer(-1) * ((Integer(238))**(Integer(-1)) * sympy.elliptic_f(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2)))) + ((Integer(167) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(10) * (Integer(7))**(Integer(-1)))), sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) * (Integer(3332))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_332():
    f = 1/((5*x**2 + 7)**3*sqrt(-x**4 + x**2 + 2))
    F = (Integer(-1) * ((Integer(25) * x * sympy.sqrt((Integer(2) + (x)**(Integer(2)) + (Integer(-1) * (x)**(Integer(4)))))) * ((Integer(952) * ((Integer(7) + (Integer(5) * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(12525) * x * sympy.sqrt((Integer(2) + (x)**(Integer(2)) + (Integer(-1) * (x)**(Integer(4)))))) * ((Integer(453152) * (Integer(7) + (Integer(5) * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2505) * sympy.elliptic_e(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) * (Integer(453152))**(Integer(-1)))) + (Integer(-1) * ((Integer(263) * sympy.elliptic_f(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) * (Integer(226576))**(Integer(-1)))) + ((Integer(58915) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(10) * (Integer(7))**(Integer(-1)))), sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) * (Integer(3172064))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_333():
    f = (5*x**2 + 7)**5/(-x**4 + x**2 + 2)**(sympy.S(3)/2)
    F = 625*x**3*sqrt(-x**4 + x**2 + 2) + x*(1419793*x**2 + 1419985)/(18*sqrt(-x**4 + x**2 + 2)) + 27500*x*sqrt(-x**4 + x**2 + 2)/3 - 3482293*elliptic_e(asin(sqrt(2)*x/2), -2)/18 + 627857*elliptic_f(asin(sqrt(2)*x/2), -2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_334():
    f = (5*x**2 + 7)**4/(-x**4 + x**2 + 2)**(sympy.S(3)/2)
    F = x*(83489*x**2 + 83585)/(18*sqrt(-x**4 + x**2 + 2)) + 625*x*sqrt(-x**4 + x**2 + 2)/3 - 165239*elliptic_e(asin(sqrt(2)*x/2), -2)/18 + 31921*elliptic_f(asin(sqrt(2)*x/2), -2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_335():
    f = (5*x**2 + 7)**3/(-x**4 + x**2 + 2)**(sympy.S(3)/2)
    F = x*(4897*x**2 + 4945)/(18*sqrt(-x**4 + x**2 + 2)) - 7147*elliptic_e(asin(sqrt(2)*x/2), -2)/18 + 1763*elliptic_f(asin(sqrt(2)*x/2), -2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_336():
    f = (5*x**2 + 7)**2/(-x**4 + x**2 + 2)**(sympy.S(3)/2)
    F = x*(281*x**2 + 305)/(18*sqrt(-x**4 + x**2 + 2)) - 281*elliptic_e(asin(sqrt(2)*x/2), -2)/18 + 139*elliptic_f(asin(sqrt(2)*x/2), -2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_337():
    f = (5*x**2 + 7)/(-x**4 + x**2 + 2)**(sympy.S(3)/2)
    F = x*(13*x**2 + 25)/(18*sqrt(-x**4 + x**2 + 2)) - 13*elliptic_e(asin(sqrt(2)*x/2), -2)/18 + 17*elliptic_f(asin(sqrt(2)*x/2), -2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_338():
    f = (-x**4 + x**2 + 2)**(sympy.S(-3)/2)
    F = x*(5 - x**2)/(18*sqrt(-x**4 + x**2 + 2)) + elliptic_e(asin(sqrt(2)*x/2), -2)/18 + elliptic_f(asin(sqrt(2)*x/2), -2)/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_339():
    f = 1/((5*x**2 + 7)*(-x**4 + x**2 + 2)**(sympy.S(3)/2))
    F = ((x * (Integer(35) + (Integer(-1) * (Integer(16) * (x)**(Integer(2)))))) * ((Integer(306) * sympy.sqrt((Integer(2) + (x)**(Integer(2)) + (Integer(-1) * (x)**(Integer(4)))))))**(Integer(-1))) + ((Integer(8) * (Integer(153))**(Integer(-1))) * sympy.elliptic_e(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) + ((Integer(102))**(Integer(-1)) * sympy.elliptic_f(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) + (Integer(-1) * ((Integer(25) * (Integer(238))**(Integer(-1))) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(10) * (Integer(7))**(Integer(-1)))), sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_340():
    f = 1/((5*x**2 + 7)**2*(-x**4 + x**2 + 2)**(sympy.S(3)/2))
    F = ((x * (Integer(580) + (Integer(-1) * (Integer(287) * (x)**(Integer(2)))))) * ((Integer(10404) * sympy.sqrt((Integer(2) + (x)**(Integer(2)) + (Integer(-1) * (x)**(Integer(4)))))))**(Integer(-1))) + ((Integer(625) * x * sympy.sqrt((Integer(2) + (x)**(Integer(2)) + (Integer(-1) * (x)**(Integer(4)))))) * ((Integer(16184) * (Integer(7) + (Integer(5) * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(5143) * sympy.elliptic_e(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) * (Integer(145656))**(Integer(-1))) + ((Integer(89) * sympy.elliptic_f(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) * (Integer(24276))**(Integer(-1))) + (Integer(-1) * ((Integer(10825) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(10) * (Integer(7))**(Integer(-1)))), sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) * (Integer(113288))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_341():
    f = 1/((5*x**2 + 7)**3*(-x**4 + x**2 + 2)**(sympy.S(3)/2))
    F = ((x * (Integer(9830) + (Integer(-1) * (Integer(4909) * (x)**(Integer(2)))))) * ((Integer(353736) * sympy.sqrt((Integer(2) + (x)**(Integer(2)) + (Integer(-1) * (x)**(Integer(4)))))))**(Integer(-1))) + ((Integer(625) * x * sympy.sqrt((Integer(2) + (x)**(Integer(2)) + (Integer(-1) * (x)**(Integer(4)))))) * ((Integer(32368) * ((Integer(7) + (Integer(5) * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))) + ((Integer(645625) * x * sympy.sqrt((Integer(2) + (x)**(Integer(2)) + (Integer(-1) * (x)**(Integer(4)))))) * ((Integer(15407168) * (Integer(7) + (Integer(5) * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(3086453) * sympy.elliptic_e(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) * (Integer(138664512))**(Integer(-1))) + ((Integer(60409) * sympy.elliptic_f(sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) * (Integer(23110752))**(Integer(-1))) + (Integer(-1) * ((Integer(6898575) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(10) * (Integer(7))**(Integer(-1)))), sympy.asin((x * (sympy.sqrt(Integer(2)))**(Integer(-1)))), Integer(-2))) * (Integer(107850176))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_342():
    f = (5*x**2 + 7)**4*sqrt(x**4 + 3*x**2 + 4)
    F = 625*x**5*(x**4 + 3*x**2 + 4)**(sympy.S(3)/2)/11 + 23500*x**3*(x**4 + 3*x**2 + 4)**(sympy.S(3)/2)/99 + x*(4516*x**2 + 18727)*sqrt(x**4 + 3*x**2 + 4)/33 + 3050*x*(x**4 + 3*x**2 + 4)**(sympy.S(3)/2)/11 + 51665*x*sqrt(x**4 + 3*x**2 + 4)/(33*x**2 + 66) - 51665*sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(x**2 + 2)*elliptic_e(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(33*sqrt(x**4 + 3*x**2 + 4)) + sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(33159*x**2 + 66318)*elliptic_f(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(22*sqrt(x**4 + 3*x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_343():
    f = (5*x**2 + 7)**3*sqrt(x**4 + 3*x**2 + 4)
    F = 125*x**3*(x**4 + 3*x**2 + 4)**(sympy.S(3)/2)/9 + x*(407*x**2 + 1708)*sqrt(x**4 + 3*x**2 + 4)/21 + 275*x*(x**4 + 3*x**2 + 4)**(sympy.S(3)/2)/7 + 4717*x*sqrt(x**4 + 3*x**2 + 4)/(21*x**2 + 42) - 4717*sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(x**2 + 2)*elliptic_e(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(21*sqrt(x**4 + 3*x**2 + 4)) + sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(1301*x**2 + 2602)*elliptic_f(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(6*sqrt(x**4 + 3*x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_344():
    f = (5*x**2 + 7)**2*sqrt(x**4 + 3*x**2 + 4)
    F = x*(38*x**2 + 119)*sqrt(x**4 + 3*x**2 + 4)/7 + 25*x*(x**4 + 3*x**2 + 4)**(sympy.S(3)/2)/7 + 319*x*sqrt(x**4 + 3*x**2 + 4)/(7*x**2 + 14) - 319*sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(x**2 + 2)*elliptic_e(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(7*sqrt(x**4 + 3*x**2 + 4)) + sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(81*x**2 + 162)*elliptic_f(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(2*sqrt(x**4 + 3*x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_345():
    f = (5*x**2 + 7)*sqrt(x**4 + 3*x**2 + 4)
    F = x*(3*x**2 + 10)*sqrt(x**4 + 3*x**2 + 4)/3 + 9*x*sqrt(x**4 + 3*x**2 + 4)/(x**2 + 2) - 9*sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(x**2 + 2)*elliptic_e(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/sqrt(x**4 + 3*x**2 + 4) + sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(49*x**2 + 98)*elliptic_f(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(6*sqrt(x**4 + 3*x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_346():
    f = sqrt(x**4 + 3*x**2 + 4)
    F = x*sqrt(x**4 + 3*x**2 + 4)/3 + x*sqrt(x**4 + 3*x**2 + 4)/(x**2 + 2) - sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(x**2 + 2)*elliptic_e(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/sqrt(x**4 + 3*x**2 + 4) + sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(7*x**2 + 14)*elliptic_f(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(6*sqrt(x**4 + 3*x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_347():
    f = sqrt(x**4 + 3*x**2 + 4)/(5*x**2 + 7)
    F = ((x * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(5) * (Integer(2) + (x)**(Integer(2)))))**(Integer(-1))) + ((Integer(5))**(Integer(-1)) * sympy.sqrt((Integer(11) * (Integer(35))**(Integer(-1)))) * sympy.atan(((Integer(2) * sympy.sqrt((Integer(11) * (Integer(35))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4)))))**(Integer(-1))))) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(5) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + ((Integer(9) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(25) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * ((Integer(11) * sympy.sqrt(Integer(2)) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(75) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + ((Integer(187) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(9) * (Integer(280))**(Integer(-1)))), (Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(525) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_348():
    f = sqrt(x**4 + 3*x**2 + 4)/(5*x**2 + 7)**2
    F = (Integer(-1) * ((x * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(70) * (Integer(2) + (x)**(Integer(2)))))**(Integer(-1)))) + ((x * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(14) * (Integer(7) + (Integer(5) * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(51) * sympy.atan(((Integer(2) * sympy.sqrt((Integer(11) * (Integer(35))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4)))))**(Integer(-1))))) * ((Integer(280) * sympy.sqrt(Integer(385))))**(Integer(-1))) + (((Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(35) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * (((Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(35) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + ((Integer(289) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(9) * (Integer(280))**(Integer(-1)))), (Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(9800) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_349():
    f = sqrt(x**4 + 3*x**2 + 4)/(5*x**2 + 7)**3
    F = (Integer(-1) * ((Integer(139) * x * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(86240) * (Integer(2) + (x)**(Integer(2)))))**(Integer(-1)))) + ((x * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(28) * ((Integer(7) + (Integer(5) * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))) + ((Integer(139) * x * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(17248) * (Integer(7) + (Integer(5) * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(14999) * sympy.atan(((Integer(2) * sympy.sqrt((Integer(11) * (Integer(35))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4)))))**(Integer(-1))))) * ((Integer(344960) * sympy.sqrt(Integer(385))))**(Integer(-1))) + ((Integer(139) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(43120) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * ((Integer(23) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(2940) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + ((Integer(254983) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(9) * (Integer(280))**(Integer(-1)))), (Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(36220800) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_350():
    f = (5*x**2 + 7)**4*(x**4 + 3*x**2 + 4)**(sympy.S(3)/2)
    F = 125*x**5*(x**4 + 3*x**2 + 4)**(sympy.S(5)/2)/3 + 2250*x**3*(x**4 + 3*x**2 + 4)**(sympy.S(5)/2)/13 + x*(131080*x**2 + 452001)*(x**4 + 3*x**2 + 4)**(sympy.S(3)/2)/1287 + 7*x*(174989*x**2 + 661429)*sqrt(x**4 + 3*x**2 + 4)/2145 + 92150*x*(x**4 + 3*x**2 + 4)**(sympy.S(5)/2)/429 + 12665086*x*sqrt(x**4 + 3*x**2 + 4)/(2145*x**2 + 4290) - 12665086*sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(x**2 + 2)*elliptic_e(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(2145*sqrt(x**4 + 3*x**2 + 4)) + 2383556*sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(x**2 + 2)*elliptic_f(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(429*sqrt(x**4 + 3*x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_351():
    f = (5*x**2 + 7)**3*(x**4 + 3*x**2 + 4)**(sympy.S(3)/2)
    F = 125*x**3*(x**4 + 3*x**2 + 4)**(sympy.S(5)/2)/13 + x*(15365*x**2 + 53504)*(x**4 + 3*x**2 + 4)**(sympy.S(3)/2)/1001 + x*(435441*x**2 + 1653701)*sqrt(x**4 + 3*x**2 + 4)/5005 + 3825*x*(x**4 + 3*x**2 + 4)**(sympy.S(5)/2)/143 + 4525662*x*sqrt(x**4 + 3*x**2 + 4)/(5005*x**2 + 10010) - 4525662*sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(x**2 + 2)*elliptic_e(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(5005*sqrt(x**4 + 3*x**2 + 4)) + 121826*sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(x**2 + 2)*elliptic_f(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(143*sqrt(x**4 + 3*x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_352():
    f = (5*x**2 + 7)**2*(x**4 + 3*x**2 + 4)**(sympy.S(3)/2)
    F = x*(2240*x**2 + 6831)*(x**4 + 3*x**2 + 4)**(sympy.S(3)/2)/693 + x*(18253*x**2 + 64533)*sqrt(x**4 + 3*x**2 + 4)/1155 + 25*x*(x**4 + 3*x**2 + 4)**(sympy.S(5)/2)/11 + 175346*x*sqrt(x**4 + 3*x**2 + 4)/(1155*x**2 + 2310) - 175346*sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(x**2 + 2)*elliptic_e(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(1155*sqrt(x**4 + 3*x**2 + 4)) + 4628*sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(x**2 + 2)*elliptic_f(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(33*sqrt(x**4 + 3*x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_353():
    f = (5*x**2 + 7)*(x**4 + 3*x**2 + 4)**(sympy.S(3)/2)
    F = x*(35*x**2 + 108)*(x**4 + 3*x**2 + 4)**(sympy.S(3)/2)/63 + x*(289*x**2 + 1029)*sqrt(x**4 + 3*x**2 + 4)/105 + 2798*x*sqrt(x**4 + 3*x**2 + 4)/(105*x**2 + 210) - 2798*sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(x**2 + 2)*elliptic_e(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(105*sqrt(x**4 + 3*x**2 + 4)) + 74*sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(x**2 + 2)*elliptic_f(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(3*sqrt(x**4 + 3*x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_354():
    f = (x**4 + 3*x**2 + 4)**(sympy.S(3)/2)
    F = x*(9*x**2 + 49)*sqrt(x**4 + 3*x**2 + 4)/35 + x*(x**4 + 3*x**2 + 4)**(sympy.S(3)/2)/7 + 138*x*sqrt(x**4 + 3*x**2 + 4)/(35*x**2 + 70) - 138*sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(x**2 + 2)*elliptic_e(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(35*sqrt(x**4 + 3*x**2 + 4)) + 4*sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(x**2 + 2)*elliptic_f(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/sqrt(x**4 + 3*x**2 + 4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_355():
    f = (x**4 + 3*x**2 + 4)**(sympy.S(3)/2)/(5*x**2 + 7)
    F = ((Integer(94) * x * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(125) * (Integer(2) + (x)**(Integer(2)))))**(Integer(-1))) + ((Integer(75))**(Integer(-1)) * x * (Integer(11) + (Integer(3) * (x)**(Integer(2)))) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) + ((Integer(44) * (Integer(125))**(Integer(-1))) * sympy.sqrt((Integer(11) * (Integer(35))**(Integer(-1)))) * sympy.atan(((Integer(2) * sympy.sqrt((Integer(11) * (Integer(35))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4)))))**(Integer(-1))))) + (Integer(-1) * ((Integer(94) * sympy.sqrt(Integer(2)) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(125) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + ((Integer(54) * sympy.sqrt(Integer(2)) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(125) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + ((Integer(4114) * sympy.sqrt(Integer(2)) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(9) * (Integer(280))**(Integer(-1)))), (Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(13125) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_356():
    f = (x**4 + 3*x**2 + 4)**(sympy.S(3)/2)/(5*x**2 + 7)**3
    F = ((Integer(9) * x * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(1960) * (Integer(2) + (x)**(Integer(2)))))**(Integer(-1))) + ((Integer(11) * x * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(175) * ((Integer(7) + (Integer(5) * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))) + ((Integer(167) * x * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(9800) * (Integer(7) + (Integer(5) * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(1347) * sympy.atan(((Integer(2) * sympy.sqrt((Integer(11) * (Integer(35))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4)))))**(Integer(-1))))) * ((Integer(7840) * sympy.sqrt(Integer(385))))**(Integer(-1))) + ((Integer(111) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(24500) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.sqrt(Integer(2)) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(875) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(817) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(91875) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(22) * sympy.sqrt(Integer(2)) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(13125) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + ((Integer(7633) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(9) * (Integer(280))**(Integer(-1)))), (Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(274400) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_357():
    f = (5*x**2 + 7)**3/sqrt(x**4 + 3*x**2 + 4)
    F = 25*x**3*sqrt(x**4 + 3*x**2 + 4) + 75*x*sqrt(x**4 + 3*x**2 + 4) - 15*x*sqrt(x**4 + 3*x**2 + 4)/(x**2 + 2) + 15*sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(x**2 + 2)*elliptic_e(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/sqrt(x**4 + 3*x**2 + 4) + sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(13*x**2 + 26)*elliptic_f(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(4*sqrt(x**4 + 3*x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_358():
    f = (5*x**2 + 7)**2/sqrt(x**4 + 3*x**2 + 4)
    F = 25*x*sqrt(x**4 + 3*x**2 + 4)/3 + 20*x*sqrt(x**4 + 3*x**2 + 4)/(x**2 + 2) - 20*sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(x**2 + 2)*elliptic_e(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/sqrt(x**4 + 3*x**2 + 4) + sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(167*x**2 + 334)*elliptic_f(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(12*sqrt(x**4 + 3*x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_359():
    f = (5*x**2 + 7)/sqrt(x**4 + 3*x**2 + 4)
    F = 5*x*sqrt(x**4 + 3*x**2 + 4)/(x**2 + 2) - 5*sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(x**2 + 2)*elliptic_e(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/sqrt(x**4 + 3*x**2 + 4) + sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(17*x**2 + 34)*elliptic_f(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(4*sqrt(x**4 + 3*x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_360():
    f = 1/sqrt(x**4 + 3*x**2 + 4)
    F = sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(x**2 + 2)*elliptic_f(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(4*sqrt(x**4 + 3*x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_361():
    f = 1/((5*x**2 + 7)*sqrt(x**4 + 3*x**2 + 4))
    F = ((Integer(4))**(Integer(-1)) * sympy.sqrt((Integer(5) * (Integer(77))**(Integer(-1)))) * sympy.atan(((Integer(2) * sympy.sqrt((Integer(11) * (Integer(35))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4)))))**(Integer(-1))))) + (Integer(-1) * (((Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(6) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + ((Integer(17) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(9) * (Integer(280))**(Integer(-1)))), (Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(84) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_362():
    f = 1/((5*x**2 + 7)**2*sqrt(x**4 + 3*x**2 + 4))
    F = (Integer(-1) * ((Integer(5) * x * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(616) * (Integer(2) + (x)**(Integer(2)))))**(Integer(-1)))) + ((Integer(25) * x * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(616) * (Integer(7) + (Integer(5) * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(37) * sympy.sqrt((Integer(5) * (Integer(77))**(Integer(-1)))) * sympy.atan(((Integer(2) * sympy.sqrt((Integer(11) * (Integer(35))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4)))))**(Integer(-1))))) * (Integer(2464))**(Integer(-1))) + ((Integer(5) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(308) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * (((Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(42) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + ((Integer(629) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(9) * (Integer(280))**(Integer(-1)))), (Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(51744) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_363():
    f = 1/((5*x**2 + 7)**3*sqrt(x**4 + 3*x**2 + 4))
    F = (Integer(-1) * ((Integer(555) * x * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(758912) * (Integer(2) + (x)**(Integer(2)))))**(Integer(-1)))) + ((Integer(25) * x * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(1232) * ((Integer(7) + (Integer(5) * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))) + ((Integer(2775) * x * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(758912) * (Integer(7) + (Integer(5) * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3285) * sympy.sqrt((Integer(5) * (Integer(77))**(Integer(-1)))) * sympy.atan(((Integer(2) * sympy.sqrt((Integer(11) * (Integer(35))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4)))))**(Integer(-1))))) * (Integer(3035648))**(Integer(-1)))) + ((Integer(555) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(379456) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * (((Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(8624) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(18615) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(9) * (Integer(280))**(Integer(-1)))), (Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(21249536) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_364():
    f = (5*x**2 + 7)**5/(x**4 + 3*x**2 + 4)**(sympy.S(3)/2)
    F = 625*x**3*sqrt(x**4 + 3*x**2 + 4) + x*(45779*x**2 + 99493)/(28*sqrt(x**4 + 3*x**2 + 4)) + 5000*x*sqrt(x**4 + 3*x**2 + 4)/3 - 220779*x*sqrt(x**4 + 3*x**2 + 4)/(28*x**2 + 56) - sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(130729*x**2 + 261458)*elliptic_f(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(24*sqrt(x**4 + 3*x**2 + 4)) + sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(220779*x**2 + 441558)*elliptic_e(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(28*sqrt(x**4 + 3*x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_365():
    f = (5*x**2 + 7)**4/(x**4 + 3*x**2 + 4)**(sympy.S(3)/2)
    F = x*(2719 - 4023*x**2)/(28*sqrt(x**4 + 3*x**2 + 4)) + 625*x*sqrt(x**4 + 3*x**2 + 4)/3 + 14523*x*sqrt(x**4 + 3*x**2 + 4)/(28*x**2 + 56) + sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(4243*x**2 + 8486)*elliptic_f(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(24*sqrt(x**4 + 3*x**2 + 4)) - sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(14523*x**2 + 29046)*elliptic_e(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(28*sqrt(x**4 + 3*x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_366():
    f = (5*x**2 + 7)**3/(x**4 + 3*x**2 + 4)**(sympy.S(3)/2)
    F = -x*(949*x**2 + 2323)/(28*sqrt(x**4 + 3*x**2 + 4)) + 4449*x*sqrt(x**4 + 3*x**2 + 4)/(28*x**2 + 56) + sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(973*x**2 + 1946)*elliptic_f(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(8*sqrt(x**4 + 3*x**2 + 4)) - sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(4449*x**2 + 8898)*elliptic_e(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(28*sqrt(x**4 + 3*x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_367():
    f = (5*x**2 + 7)**2/(x**4 + 3*x**2 + 4)**(sympy.S(3)/2)
    F = -x*(9 - 113*x**2)/(28*sqrt(x**4 + 3*x**2 + 4)) - 113*x*sqrt(x**4 + 3*x**2 + 4)/(28*x**2 + 56) + sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(9*x**2 + 18)*elliptic_f(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(8*sqrt(x**4 + 3*x**2 + 4)) + sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(113*x**2 + 226)*elliptic_e(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(28*sqrt(x**4 + 3*x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_368():
    f = (5*x**2 + 7)/(x**4 + 3*x**2 + 4)**(sympy.S(3)/2)
    F = x*(19*x**2 + 53)/(28*sqrt(x**4 + 3*x**2 + 4)) - 19*x*sqrt(x**4 + 3*x**2 + 4)/(28*x**2 + 56) - sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(3*x**2 + 6)*elliptic_f(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(8*sqrt(x**4 + 3*x**2 + 4)) + sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(19*x**2 + 38)*elliptic_e(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(28*sqrt(x**4 + 3*x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_369():
    f = (x**4 + 3*x**2 + 4)**(sympy.S(-3)/2)
    F = -x*(3*x**2 + 1)/(28*sqrt(x**4 + 3*x**2 + 4)) + 3*x*sqrt(x**4 + 3*x**2 + 4)/(28*x**2 + 56) + sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(x**2 + 2)*elliptic_f(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(8*sqrt(x**4 + 3*x**2 + 4)) - sqrt(2)*sqrt((x**4 + 3*x**2 + 4)/(x**2 + 2)**2)*(3*x**2 + 6)*elliptic_e(2*atan(sqrt(2)*x/2), sympy.S(1)/8)/(28*sqrt(x**4 + 3*x**2 + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_370():
    f = 1/((5*x**2 + 7)*(x**4 + 3*x**2 + 4)**(sympy.S(3)/2))
    F = (Integer(-1) * ((x * (Integer(13) + (Integer(4) * (x)**(Integer(2))))) * ((Integer(308) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + ((x * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(77) * (Integer(2) + (x)**(Integer(2)))))**(Integer(-1))) + ((Integer(25) * (Integer(176))**(Integer(-1))) * sympy.sqrt((Integer(5) * (Integer(77))**(Integer(-1)))) * sympy.atan(((Integer(2) * sympy.sqrt((Integer(11) * (Integer(35))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4)))))**(Integer(-1))))) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(77) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + (Integer(-1) * (((Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(12) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + ((Integer(425) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(9) * (Integer(280))**(Integer(-1)))), (Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(3696) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_371():
    f = 1/((5*x**2 + 7)**2*(x**4 + 3*x**2 + 4)**(sympy.S(3)/2))
    F = ((x * (Integer(24) + (Integer(37) * (x)**(Integer(2))))) * ((Integer(13552) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * ((Integer(199) * x * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(27104) * (Integer(2) + (x)**(Integer(2)))))**(Integer(-1)))) + ((Integer(625) * x * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(27104) * (Integer(7) + (Integer(5) * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(575) * sympy.sqrt((Integer(5) * (Integer(77))**(Integer(-1)))) * sympy.atan(((Integer(2) * sympy.sqrt((Integer(11) * (Integer(35))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4)))))**(Integer(-1))))) * (Integer(108416))**(Integer(-1))) + ((Integer(199) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(13552) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Integer(2)) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(231) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))) + ((Integer(9775) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(9) * (Integer(280))**(Integer(-1)))), (Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(2276736) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_372():
    f = 1/((5*x**2 + 7)**3*(x**4 + 3*x**2 + 4)**(sympy.S(3)/2))
    F = ((x * (Integer(548) + (Integer(139) * (x)**(Integer(2))))) * ((Integer(596288) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * ((Integer(18159) * x * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(33392128) * (Integer(2) + (x)**(Integer(2)))))**(Integer(-1)))) + ((Integer(625) * x * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(54208) * ((Integer(7) + (Integer(5) * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))) + ((Integer(51875) * x * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))) * ((Integer(33392128) * (Integer(7) + (Integer(5) * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((Integer(529425) * sympy.sqrt((Integer(5) * (Integer(77))**(Integer(-1)))) * sympy.atan(((Integer(2) * sympy.sqrt((Integer(11) * (Integer(35))**(Integer(-1)))) * x) * (sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4)))))**(Integer(-1))))) * (Integer(133568512))**(Integer(-1)))) + ((Integer(18159) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(16696064) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + ((Integer(843) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(379456) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3000075) * (Integer(2) + (x)**(Integer(2))) * sympy.sqrt(((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))) * (((Integer(2) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(9) * (Integer(280))**(Integer(-1)))), (Integer(2) * sympy.atan((x * (sympy.sqrt(Integer(2)))**(Integer(-1))))), (Integer(8))**(Integer(-1)))) * ((Integer(934979584) * sympy.sqrt(Integer(2)) * sympy.sqrt((Integer(4) + (Integer(3) * (x)**(Integer(2))) + (x)**(Integer(4))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_373():
    f = (d + e*x**2)**3/sqrt(a + b*x**2 + c*x**4)
    F = -a**(sympy.S(1)/4)*e*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(8*b**2*e**2 + 45*c**2*d**2 - 3*c*e*(3*a*e + 10*b*d))*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(15*c**(sympy.S(11)/4)*sqrt(a + b*x**2 + c*x**4)) + a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(e*(8*b**2*e**2 + 45*c**2*d**2 - 3*c*e*(3*a*e + 10*b*d)) + sqrt(c)*(4*a*b*e**3 - 15*a*c*d*e**2 + 15*c**2*d**3)/sqrt(a))*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(30*c**(sympy.S(11)/4)*sqrt(a + b*x**2 + c*x**4)) + e**3*x**3*sqrt(a + b*x**2 + c*x**4)/(5*c) + e**2*x*(-4*b*e + 15*c*d)*sqrt(a + b*x**2 + c*x**4)/(15*c**2) + e*x*sqrt(a + b*x**2 + c*x**4)*(8*b**2*e**2 + 45*c**2*d**2 - 3*c*e*(3*a*e + 10*b*d))/(15*c**(sympy.S(5)/2)*(sqrt(a) + sqrt(c)*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_374():
    f = (d + e*x**2)**2/sqrt(a + b*x**2 + c*x**4)
    F = -2*a**(sympy.S(1)/4)*e*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-b*e + 3*c*d)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(3*c**(sympy.S(7)/4)*sqrt(a + b*x**2 + c*x**4)) + a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(2*e*(-b*e + 3*c*d) + sqrt(c)*(-a*e**2 + 3*c*d**2)/sqrt(a))*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(6*c**(sympy.S(7)/4)*sqrt(a + b*x**2 + c*x**4)) + e**2*x*sqrt(a + b*x**2 + c*x**4)/(3*c) + 2*e*x*(-b*e + 3*c*d)*sqrt(a + b*x**2 + c*x**4)/(3*c**(sympy.S(3)/2)*(sqrt(a) + sqrt(c)*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_375():
    f = (d + e*x**2)/sqrt(a + b*x**2 + c*x**4)
    F = -a**(sympy.S(1)/4)*e*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(c**(sympy.S(3)/4)*sqrt(a + b*x**2 + c*x**4)) + a**(sympy.S(1)/4)*sqrt((a + b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(e + sqrt(c)*d/sqrt(a))*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half - b/(4*sqrt(a)*sqrt(c)))/(2*c**(sympy.S(3)/4)*sqrt(a + b*x**2 + c*x**4)) + e*x*sqrt(a + b*x**2 + c*x**4)/(sqrt(c)*(sqrt(a) + sqrt(c)*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_376():
    f = 1/((d + e*x**2)*sqrt(a + b*x**2 + c*x**4))
    F = ((sympy.sqrt(Symbol('e')) * sympy.atan(((sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2))))) * x) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('d')) * sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * (Symbol('b') * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))))))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('c')) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * Symbol('e')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt(Symbol('a')))**(Integer(-1))) + Symbol('e')))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * Symbol('e')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * Symbol('d') * Symbol('e')))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * (Symbol('b') * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))))))) * ((Integer(4) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('d') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_377():
    f = 1/((d + e*x**2)**2*sqrt(a + b*x**2 + c*x**4))
    F = (Integer(-1) * ((sympy.sqrt(Symbol('c')) * Symbol('e') * x * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))))) * ((Integer(2) * Symbol('d') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2))))))**(Integer(-1)))) + (((Symbol('e'))**(Integer(2)) * x * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))))) * ((Integer(2) * Symbol('d') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))) + ((sympy.sqrt(Symbol('e')) * ((Integer(3) * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('e') * ((Integer(2) * Symbol('b') * Symbol('d')) + (Integer(-1) * (Symbol('a') * Symbol('e'))))))) * sympy.atan(((sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2))))) * x) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('e') * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * (Symbol('b') * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))))))) * ((Integer(2) * Symbol('d') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * (Symbol('b') * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))))))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * Symbol('d') * ((sympy.sqrt(Symbol('c')) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * Symbol('e')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) + (sympy.sqrt(Symbol('a')) * Symbol('e'))) * ((Integer(3) * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('e') * ((Integer(2) * Symbol('b') * Symbol('d')) + (Integer(-1) * (Symbol('a') * Symbol('e'))))))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * Symbol('e')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * Symbol('d') * Symbol('e')))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * (Symbol('b') * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))))))) * ((Integer(8) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**(Integer(2)) * ((sympy.sqrt(Symbol('c')) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * Symbol('e')))) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('e'))) + (Symbol('a') * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_378():
    f = (d + e*x**2)**3/sqrt(a + b*x**2 - c*x**4)
    F = -e**3*x**3*sqrt(a + b*x**2 - c*x**4)/(5*c) - e**2*x*(4*b*e + 15*c*d)*sqrt(a + b*x**2 - c*x**4)/(15*c**2) - sqrt(2)*e*(b - sqrt(4*a*c + b**2))*sqrt(b + sqrt(4*a*c + b**2))*sqrt(-2*c*x**2/(b - sqrt(4*a*c + b**2)) + 1)*sqrt(-2*c*x**2/(b + sqrt(4*a*c + b**2)) + 1)*(8*b**2*e**2 + 45*c**2*d**2 + 3*c*e*(3*a*e + 10*b*d))*elliptic_e(asin(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(4*a*c + b**2))), (b + sqrt(4*a*c + b**2))/(b - sqrt(4*a*c + b**2)))/(60*c**(sympy.S(7)/2)*sqrt(a + b*x**2 - c*x**4)) + sqrt(2)*(b - sqrt(4*a*c + b**2))*sqrt(b + sqrt(4*a*c + b**2))*sqrt(-2*c*x**2/(b - sqrt(4*a*c + b**2)) + 1)*sqrt(-2*c*x**2/(b + sqrt(4*a*c + b**2)) + 1)*(2*c*(4*a*b*e**3 + 15*a*c*d*e**2 + 15*c**2*d**3)/(b - sqrt(4*a*c + b**2)) + e*(8*b**2*e**2 + 45*c**2*d**2 + 3*c*e*(3*a*e + 10*b*d)))*elliptic_f(asin(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(4*a*c + b**2))), (b + sqrt(4*a*c + b**2))/(b - sqrt(4*a*c + b**2)))/(60*c**(sympy.S(7)/2)*sqrt(a + b*x**2 - c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_379():
    f = (d + e*x**2)**2/sqrt(a + b*x**2 - c*x**4)
    F = -e**2*x*sqrt(a + b*x**2 - c*x**4)/(3*c) - sqrt(2)*e*(b - sqrt(4*a*c + b**2))*sqrt(b + sqrt(4*a*c + b**2))*(b*e + 3*c*d)*sqrt(-2*c*x**2/(b - sqrt(4*a*c + b**2)) + 1)*sqrt(-2*c*x**2/(b + sqrt(4*a*c + b**2)) + 1)*elliptic_e(asin(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(4*a*c + b**2))), (b + sqrt(4*a*c + b**2))/(b - sqrt(4*a*c + b**2)))/(6*c**(sympy.S(5)/2)*sqrt(a + b*x**2 - c*x**4)) + sqrt(2)*sqrt(b + sqrt(4*a*c + b**2))*sqrt(-2*c*x**2/(b - sqrt(4*a*c + b**2)) + 1)*sqrt(-2*c*x**2/(b + sqrt(4*a*c + b**2)) + 1)*(b*e**2*(b - sqrt(4*a*c + b**2)) + 3*c**2*d**2 + c*e*(a*e + 3*b*d - 3*d*sqrt(4*a*c + b**2)))*elliptic_f(asin(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(4*a*c + b**2))), (b + sqrt(4*a*c + b**2))/(b - sqrt(4*a*c + b**2)))/(6*c**(sympy.S(5)/2)*sqrt(a + b*x**2 - c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_380():
    f = (d + e*x**2)/sqrt(a + b*x**2 - c*x**4)
    F = -sqrt(2)*e*(b - sqrt(4*a*c + b**2))*sqrt(b + sqrt(4*a*c + b**2))*sqrt(-2*c*x**2/(b - sqrt(4*a*c + b**2)) + 1)*sqrt(-2*c*x**2/(b + sqrt(4*a*c + b**2)) + 1)*elliptic_e(asin(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(4*a*c + b**2))), (b + sqrt(4*a*c + b**2))/(b - sqrt(4*a*c + b**2)))/(4*c**(sympy.S(3)/2)*sqrt(a + b*x**2 - c*x**4)) + sqrt(2)*sqrt(b + sqrt(4*a*c + b**2))*(2*c*d + e*(b - sqrt(4*a*c + b**2)))*sqrt(-2*c*x**2/(b - sqrt(4*a*c + b**2)) + 1)*sqrt(-2*c*x**2/(b + sqrt(4*a*c + b**2)) + 1)*elliptic_f(asin(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(4*a*c + b**2))), (b + sqrt(4*a*c + b**2))/(b - sqrt(4*a*c + b**2)))/(4*c**(sympy.S(3)/2)*sqrt(a + b*x**2 - c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_381():
    f = 1/((d + e*x**2)*sqrt(a + b*x**2 - c*x**4))
    F = (sympy.sqrt((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (x)**(Integer(2))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (x)**(Integer(2))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * (((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e')) * ((Integer(2) * Symbol('c') * Symbol('d')))**(Integer(-1)))), sympy.asin(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))), ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c'))))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1))))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_382():
    f = 1/((d + e*x**2)**2*sqrt(a + b*x**2 - c*x**4))
    F = (Integer(-1) * (((Symbol('e'))**(Integer(2)) * x * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))) * ((Integer(2) * Symbol('d') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('b') * Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))) * (Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1)))) + (((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))) * sympy.sqrt((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (x)**(Integer(2))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (x)**(Integer(2))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))), ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c'))))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * Symbol('d') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('e') * ((Symbol('b') * Symbol('d')) + (Integer(-1) * (Symbol('a') * Symbol('e')))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))) * ((Integer(2) * Symbol('c') * Symbol('d')) + ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (x)**(Integer(2))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (x)**(Integer(2))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))), ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c'))))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * Symbol('d') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('e') * ((Symbol('b') * Symbol('d')) + (Integer(-1) * (Symbol('a') * Symbol('e')))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))) * ((Integer(3) * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('e') * ((Integer(2) * Symbol('b') * Symbol('d')) + (Integer(-1) * (Symbol('a') * Symbol('e')))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (x)**(Integer(2))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (x)**(Integer(2))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))))) * sympy.Function('EllipticPi')((Integer(-1) * (((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e')) * ((Integer(2) * Symbol('c') * Symbol('d')))**(Integer(-1)))), sympy.asin(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))), ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c'))))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2)) * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('e') * ((Symbol('b') * Symbol('d')) + (Integer(-1) * (Symbol('a') * Symbol('e')))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_383():
    f = (d + e*x**2)/sqrt(-a + b*x**2 + c*x**4)
    F = e*x*(b - sqrt(4*a*c + b**2))*(2*c*x**2/(b - sqrt(4*a*c + b**2)) + 1)/(2*c*sqrt(-a + b*x**2 + c*x**4)) + sqrt(2)*d*sqrt(b + sqrt(4*a*c + b**2))*(2*c*x**2/(b - sqrt(4*a*c + b**2)) + 1)*elliptic_f(atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(4*a*c + b**2))), -2*sqrt(4*a*c + b**2)/(b - sqrt(4*a*c + b**2)))/(2*sqrt(c)*sqrt((2*c*x**2/(b - sqrt(4*a*c + b**2)) + 1)/(2*c*x**2/(b + sqrt(4*a*c + b**2)) + 1))*sqrt(-a + b*x**2 + c*x**4)) - sqrt(2)*e*(b - sqrt(4*a*c + b**2))*sqrt(b + sqrt(4*a*c + b**2))*(2*c*x**2/(b - sqrt(4*a*c + b**2)) + 1)*elliptic_e(atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(4*a*c + b**2))), -2*sqrt(4*a*c + b**2)/(b - sqrt(4*a*c + b**2)))/(4*c**(sympy.S(3)/2)*sqrt((2*c*x**2/(b - sqrt(4*a*c + b**2)) + 1)/(2*c*x**2/(b + sqrt(4*a*c + b**2)) + 1))*sqrt(-a + b*x**2 + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_384():
    f = 1/((d + e*x**2)*sqrt(-a + b*x**2 + c*x**4))
    F = (sympy.sqrt(((Integer(-1) * Symbol('b')) + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))) * sympy.sqrt((Integer(1) + ((Integer(2) * Symbol('c') * (x)**(Integer(2))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1))))) * sympy.sqrt((Integer(1) + ((Integer(2) * Symbol('c') * (x)**(Integer(2))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * sympy.Function('EllipticPi')((((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')) * ((Integer(2) * Symbol('c') * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt(((Integer(-1) * Symbol('b')) + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))), ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * Symbol('d') * sympy.sqrt(((Integer(-1) * Symbol('a')) + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_385():
    f = (d + e*x**2)/sqrt(-a + b*x**2 - c*x**4)
    F = -a**(sympy.S(1)/4)*e*sqrt((a - b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half + b/(4*sqrt(a)*sqrt(c)))/(c**(sympy.S(3)/4)*sqrt(-a + b*x**2 - c*x**4)) + a**(sympy.S(1)/4)*sqrt((a - b*x**2 + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(e + sqrt(c)*d/sqrt(a))*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half + b/(4*sqrt(a)*sqrt(c)))/(2*c**(sympy.S(3)/4)*sqrt(-a + b*x**2 - c*x**4)) - e*x*sqrt(-a + b*x**2 - c*x**4)/(sqrt(c)*(sqrt(a) + sqrt(c)*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_386():
    f = 1/((d + e*x**2)*sqrt(-a + b*x**2 - c*x**4))
    F = ((sympy.sqrt(Symbol('e')) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('e') * ((Symbol('b') * Symbol('d')) + (Symbol('a') * Symbol('e'))))))) * x) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt(Symbol('e')) * sympy.sqrt(((Integer(-1) * Symbol('a')) + (Symbol('b') * (x)**(Integer(2))) + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('d')) * sympy.sqrt((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('e') * ((Symbol('b') * Symbol('d')) + (Symbol('a') * Symbol('e')))))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))) + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Symbol('b') * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('c')) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * Symbol('e')))) * sympy.sqrt(((Integer(-1) * Symbol('a')) + (Symbol('b') * (x)**(Integer(2))) + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) * (sympy.sqrt(Symbol('a')))**(Integer(-1))) + Symbol('e')))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))) + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(-1) * ((((sympy.sqrt(Symbol('c')) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * Symbol('e')))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * Symbol('d') * Symbol('e')))**(Integer(-1)))), (Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Symbol('b') * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))))) * ((Integer(4) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('d') * ((Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(2))))) * sympy.sqrt(((Integer(-1) * Symbol('a')) + (Symbol('b') * (x)**(Integer(2))) + (Integer(-1) * (Symbol('c') * (x)**(Integer(4))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_387():
    f = (d + e*x**2)**3/sqrt(x**4 + 3*x**2 + 2)
    F = e**3*x**3*sqrt(x**4 + 3*x**2 + 2)/5 + e**2*x*(d - 4*e/5)*sqrt(x**4 + 3*x**2 + 2) + 3*e*x*(x**2 + 2)*(5*d**2 - 10*d*e + 6*e**2)/(5*sqrt(x**4 + 3*x**2 + 2)) - 3*sqrt(2)*e*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*(5*d**2 - 10*d*e + 6*e**2)*elliptic_e(atan(x), sympy.S.Half)/(5*sqrt(x**4 + 3*x**2 + 2)) + sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*(5*d**3 - 10*d*e**2 + 8*e**3)*elliptic_f(atan(x), sympy.S.Half)/(10*sqrt(x**4 + 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_388():
    f = (d + e*x**2)**2/sqrt(x**4 + 3*x**2 + 2)
    F = e**2*x*sqrt(x**4 + 3*x**2 + 2)/3 + e*x*(2*d - 2*e)*(x**2 + 2)/sqrt(x**4 + 3*x**2 + 2) - 2*sqrt(2)*e*sqrt((x**2 + 2)/(x**2 + 1))*(d - e)*(x**2 + 1)*elliptic_e(atan(x), sympy.S.Half)/sqrt(x**4 + 3*x**2 + 2) + sqrt(2)*sqrt((x**2 + 2)/(x**2 + 1))*(3*d**2 - 2*e**2)*(x**2 + 1)*elliptic_f(atan(x), sympy.S.Half)/(6*sqrt(x**4 + 3*x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_389():
    f = (d + e*x**2)/sqrt(x**4 + 3*x**2 + 2)
    F = sqrt(2)*d*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_f(atan(x), sympy.S.Half)/(2*sqrt(x**4 + 3*x**2 + 2)) + e*x*(x**2 + 2)/sqrt(x**4 + 3*x**2 + 2) - sqrt(2)*e*sqrt((x**2 + 2)/(x**2 + 1))*(x**2 + 1)*elliptic_e(atan(x), sympy.S.Half)/sqrt(x**4 + 3*x**2 + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_390():
    f = (c + e*x**2)**q*(a + b*x**4 + c*x**2)**p
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('e') * (x)**(Integer(2)))))**(Symbol('q')) * ((Symbol('a') + (Symbol('c') * (x)**(Integer(2))) + (Symbol('b') * (x)**(Integer(4)))))**(Symbol('p'))), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_391():
    f = (c + e*x**2)**3*(a + b*x**4 + c*x**2)**p
    F = e**3*x**3*(a + b*x**4 + c*x**2)**(p + 1)/(b*(4*p + 7)) + c*e**2*x*(a + b*x**4 + c*x**2)**(p + 1)*(12*b*p + 21*b - 2*e*p - 5*e)/(b**2*(4*p + 5)*(4*p + 7)) + c*x*(a + b*x**4 + c*x**2)**p*(-3*a*b*e**2*(4*p + 7) + a*e**3*(2*p + 5) + b**2*c**2*(16*p**2 + 48*p + 35))*appellf1(sympy.S.Half, -p, -p, sympy.S(3)/2, -2*b*x**2/(c - sqrt(-4*a*b + c**2)), -2*b*x**2/(c + sqrt(-4*a*b + c**2)))/(b**2*(4*p + 5)*(4*p + 7)*(2*b*x**2/(c - sqrt(-4*a*b + c**2)) + 1)**p*(2*b*x**2/(c + sqrt(-4*a*b + c**2)) + 1)**p) + e*x**3*(a + b*x**4 + c*x**2)**p*(3*b**2*c**2*(16*p**2 + 48*p + 35) - 3*b*e*(a*e*(4*p + 5) + c**2*(8*p**2 + 26*p + 21)) + c**2*e**2*(4*p**2 + 16*p + 15))*appellf1(sympy.S(3)/2, -p, -p, sympy.S(5)/2, -2*b*x**2/(c - sqrt(-4*a*b + c**2)), -2*b*x**2/(c + sqrt(-4*a*b + c**2)))/(3*b**2*(4*p + 5)*(4*p + 7)*(2*b*x**2/(c - sqrt(-4*a*b + c**2)) + 1)**p*(2*b*x**2/(c + sqrt(-4*a*b + c**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_392():
    f = (c + e*x**2)**2*(a + b*x**4 + c*x**2)**p
    F = c*e*x**3*(a + b*x**4 + c*x**2)**p*(8*b*p + 10*b - 2*e*p - 3*e)*appellf1(sympy.S(3)/2, -p, -p, sympy.S(5)/2, -2*b*x**2/(c - sqrt(-4*a*b + c**2)), -2*b*x**2/(c + sqrt(-4*a*b + c**2)))/(3*b*(4*p + 5)*(2*b*x**2/(c - sqrt(-4*a*b + c**2)) + 1)**p*(2*b*x**2/(c + sqrt(-4*a*b + c**2)) + 1)**p) + e**2*x*(a + b*x**4 + c*x**2)**(p + 1)/(b*(4*p + 5)) - x*(a*e**2 - b*c**2*(4*p + 5))*(a + b*x**4 + c*x**2)**p*appellf1(sympy.S.Half, -p, -p, sympy.S(3)/2, -2*b*x**2/(c - sqrt(-4*a*b + c**2)), -2*b*x**2/(c + sqrt(-4*a*b + c**2)))/(b*(4*p + 5)*(2*b*x**2/(c - sqrt(-4*a*b + c**2)) + 1)**p*(2*b*x**2/(c + sqrt(-4*a*b + c**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_393():
    f = (c + e*x**2)*(a + b*x**4 + c*x**2)**p
    F = c*x*(a + b*x**4 + c*x**2)**p*appellf1(sympy.S.Half, -p, -p, sympy.S(3)/2, -2*b*x**2/(c - sqrt(-4*a*b + c**2)), -2*b*x**2/(c + sqrt(-4*a*b + c**2)))/((2*b*x**2/(c - sqrt(-4*a*b + c**2)) + 1)**p*(2*b*x**2/(c + sqrt(-4*a*b + c**2)) + 1)**p) + e*x**3*(a + b*x**4 + c*x**2)**p*appellf1(sympy.S(3)/2, -p, -p, sympy.S(5)/2, -2*b*x**2/(c - sqrt(-4*a*b + c**2)), -2*b*x**2/(c + sqrt(-4*a*b + c**2)))/(3*(2*b*x**2/(c - sqrt(-4*a*b + c**2)) + 1)**p*(2*b*x**2/(c + sqrt(-4*a*b + c**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_394():
    f = (a + b*x**4 + c*x**2)**p
    F = x*(a + b*x**4 + c*x**2)**p*appellf1(sympy.S.Half, -p, -p, sympy.S(3)/2, -2*b*x**2/(c - sqrt(-4*a*b + c**2)), -2*b*x**2/(c + sqrt(-4*a*b + c**2)))/((2*b*x**2/(c - sqrt(-4*a*b + c**2)) + 1)**p*(2*b*x**2/(c + sqrt(-4*a*b + c**2)) + 1)**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_395():
    f = (a + b*x**4 + c*x**2)**p/(c + e*x**2)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('c') * (x)**(Integer(2))) + (Symbol('b') * (x)**(Integer(4)))))**(Symbol('p')) * ((Symbol('c') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_396():
    f = (a + b*x**4 + c*x**2)**p/(c + e*x**2)**2
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('c') * (x)**(Integer(2))) + (Symbol('b') * (x)**(Integer(4)))))**(Symbol('p')) * (((Symbol('c') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_397():
    f = (f + g*x)/(sqrt(a + c*x**4)*(d + e*x))
    F = ((((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4)))))) * x) * ((Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4))))))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.atanh((((Symbol('a') * (Symbol('e'))**(Integer(2))) + (Symbol('c') * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2)))) * ((sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4))))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))))))**(Integer(-1)))) + ((((sympy.sqrt(Symbol('c')) * Symbol('d') * Symbol('f')) + (sympy.sqrt(Symbol('a')) * Symbol('e') * Symbol('g'))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2))))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('d') * Symbol('e') * ((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_398():
    f = (f + g*x)/(sqrt(-a + c*x**4)*(d + e*x))
    F = ((((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.atanh((((Symbol('a') * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Symbol('c') * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2))))) * ((sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4)))))) * sympy.sqrt(((Integer(-1) * Symbol('a')) + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4))))))))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * Symbol('g') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.elliptic_f(sympy.asin((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('e') * sympy.sqrt(((Integer(-1) * Symbol('a')) + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(4))) * (Symbol('a'))**(Integer(-1)))))) * sympy.Function('EllipticPi')(((sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2))) * ((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))))**(Integer(-1))), sympy.asin((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1)))), Integer(-1))) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('d') * Symbol('e') * sympy.sqrt(((Integer(-1) * Symbol('a')) + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_399():
    f = (x - sqrt(3) + 1)/((x + 1 + sqrt(3))*sqrt(x**4 + 4*sqrt(3)*x**2 - 4))
    F = sqrt(-3 + 2*sqrt(3))*atanh((x - sqrt(3) + 1)**2/(sqrt(-9 + 6*sqrt(3))*sqrt(x**4 + 4*sqrt(3)*x**2 - 4)))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_400():
    f = (x + 1 + sqrt(3))/((x - sqrt(3) + 1)*sqrt(x**4 - 4*sqrt(3)*x**2 - 4))
    F = -sqrt(3 + 2*sqrt(3))*atan((x + 1 + sqrt(3))**2/(sqrt(9 + 6*sqrt(3))*sqrt(x**4 - 4*sqrt(3)*x**2 - 4)))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_401():
    f = (2*x - sqrt(3) + 1)/((2*x + 1 + sqrt(3))*sqrt(4*x**4 + 4*sqrt(3)*x**2 - 1))
    F = sqrt(-3 + 2*sqrt(3))*atanh((2*x - sqrt(3) + 1)**2/(2*sqrt(-9 + 6*sqrt(3))*sqrt(4*x**4 + 4*sqrt(3)*x**2 - 1)))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_402():
    f = (2*x + 1 + sqrt(3))/((2*x - sqrt(3) + 1)*sqrt(4*x**4 - 4*sqrt(3)*x**2 - 1))
    F = -sqrt(3 + 2*sqrt(3))*atan((2*x + 1 + sqrt(3))**2/(2*sqrt(9 + 6*sqrt(3))*sqrt(4*x**4 - 4*sqrt(3)*x**2 - 1)))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_403():
    f = (f + g*x)/((d + e*x)*sqrt(a + b*x**2 + c*x**4))
    F = ((((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('b') * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2)))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4)))))) * x) * ((Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * ((Symbol('b') * (Symbol('d'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(2))))))))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.atanh((((Symbol('b') * (Symbol('d'))**(Integer(2))) + (Integer(2) * Symbol('a') * (Symbol('e'))**(Integer(2))) + (((Integer(2) * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('b') * (Symbol('e'))**(Integer(2)))) * (x)**(Integer(2)))) * ((Integer(2) * sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('b') * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(4))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('b') * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))))))**(Integer(-1)))) + ((((sympy.sqrt(Symbol('c')) * Symbol('d') * Symbol('f')) + (sympy.sqrt(Symbol('a')) * Symbol('e') * Symbol('g'))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * (Symbol('b') * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))))))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2))))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), ((Integer(4))**(Integer(-1)) * (Integer(2) + (Integer(-1) * (Symbol('b') * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))))))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('d') * Symbol('e') * ((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_2_Trinomial_products_1_2_2_Quartic_1_2_2_3_d_plus_e_x_pow_2_pow_m_a_plus_b_x_pow_2_plus_c_x_pow_4_pow_p_404():
    f = (f + g*x)/((d + e*x)*sqrt(-a + b*x**2 + c*x**4))
    F = (Integer(-1) * ((((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.atanh((((Symbol('b') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('a') * (Symbol('e'))**(Integer(2)))) + (((Integer(2) * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Symbol('b') * (Symbol('e'))**(Integer(2)))) * (x)**(Integer(2)))) * ((Integer(2) * sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('b') * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4)))))) * sympy.sqrt(((Integer(-1) * Symbol('a')) + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('b') * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4))))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g') * (Integer(1) + ((Integer(2) * Symbol('c') * (x)**(Integer(2))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))) * sympy.elliptic_f(sympy.atan(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))), (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c'))))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * Symbol('e') * sympy.sqrt(((Integer(1) + ((Integer(2) * Symbol('c') * (x)**(Integer(2))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))) * ((Integer(1) + ((Integer(2) * Symbol('c') * (x)**(Integer(2))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))))**(Integer(-1)))) * sympy.sqrt(((Integer(-1) * Symbol('a')) + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + ((sympy.sqrt(((Integer(-1) * Symbol('b')) + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))) * ((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))) * sympy.sqrt((Integer(1) + ((Integer(2) * Symbol('c') * (x)**(Integer(2))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1))))) * sympy.sqrt((Integer(1) + ((Integer(2) * Symbol('c') * (x)**(Integer(2))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * (((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))) * (Symbol('e'))**(Integer(2))) * ((Integer(2) * Symbol('c') * (Symbol('d'))**(Integer(2))))**(Integer(-1)))), sympy.asin(((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt(((Integer(-1) * Symbol('b')) + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))), ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c')))))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('c')) * Symbol('d') * Symbol('e') * sympy.sqrt(((Integer(-1) * Symbol('a')) + (Symbol('b') * (x)**(Integer(2))) + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1)))
    assert integrate(f, x) == F

